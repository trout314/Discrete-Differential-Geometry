#!/usr/bin/env python3
"""Replica exchange (parallel tempering) over VDV coefficient.

Runs K replicas at geometrically-spaced codim3_degree_variance_coef values.
Between sampling rounds, attempts swaps of adjacent replicas' VDV coefficients.
The highest-VDV replica benefits from mixing at lower VDV, producing
low-degree-variance samples that would otherwise be inaccessible.

Usage:
    # 8 replicas from VDV=0.5 to VDV=20000, 100K facet seed
    python scripts/replica_exchange_vdv.py \
        --seed-file seeds/S3_N1e5_1e-1_ED5p0043_1e-1_VDV_5e-1_s000.mfd \
        --vdv-min 0.5 --vdv-max 20000 --num-replicas 8

    # Custom ladder
    python scripts/replica_exchange_vdv.py \
        --seed-file seeds/S3_N1e4_1e-1_ED5p0043_1e-1_VDV_5e-1_s000.mfd \
        --vdv-ladder 0.5,2,8,32,128,512,2048,8192
"""

import argparse
import csv
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'python'))
import discrete_differential_geometry as ddg


# --- Memory safety ---

MIN_FREE_MEMORY_GB = 4.0


def get_free_memory_gb():
    """Return available system memory in GB from /proc/meminfo."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / (1024 * 1024)
    except (FileNotFoundError, OSError):
        return float("inf")
    return float("inf")


def check_memory(context=""):
    """Abort if free memory is below the safety threshold."""
    free = get_free_memory_gb()
    if free < MIN_FREE_MEMORY_GB:
        print(f"ABORTING: free memory {free:.1f} GB < {MIN_FREE_MEMORY_GB:.0f} GB threshold"
              f"{' (' + context + ')' if context else ''}")
        sys.exit(42)


# --- VDV penalty computation ---

def codim3_penalty(mfd, dim):
    """Compute the intensive codim-3 degree variance penalty.

    This matches the localSolidAngleCurvPenalty in sampler.d:
      (degTarget^2 * n - 2*degTarget*c3pf*F + totSqDeg - minPenalty) / n

    Which simplifies to:
      degree_variance - (x - x^2)
    where x = frac(degTarget) and degTarget = c3pf * F / n.
    """
    if dim <= 2:
        return 0.0

    codim3_dim = dim - 3  # simplex dimension of codim-3 faces
    fv = mfd.f_vector
    n_facets = fv[dim]
    n_codim3 = fv[codim3_dim]

    if n_codim3 == 0:
        return 0.0

    # codim3 simplices per facet = binomial(dim+1, dim-2)
    c3pf = math.comb(dim + 1, dim - 2)
    deg_target = c3pf * n_facets / n_codim3
    x = deg_target - math.floor(deg_target)
    min_penalty = x * (1 - x)

    dv = mfd.degree_variance(codim3_dim)
    return dv - min_penalty


def swap_accept(beta_i, beta_j, penalty_i, penalty_j):
    """Compute replica exchange swap acceptance probability.

    Returns True if the swap should be accepted.
    """
    delta = (beta_j - beta_i) * (penalty_i - penalty_j)
    if delta <= 0:
        return True
    return np.random.random() < math.exp(-delta)


def main():
    parser = argparse.ArgumentParser(
        description='Replica exchange over VDV coefficient')
    parser.add_argument('--seed-file', type=str, required=True,
                        help='Path to seed .mfd file')
    parser.add_argument('--dim', type=int, default=3,
                        help='Manifold dimension (default: 3)')
    parser.add_argument('--hinge-target', type=float, default=5.0043,
                        help='Hinge degree target (default: 5.0043)')
    parser.add_argument('--num-facets-coef', type=float, default=0.1,
                        help='Volume penalty coefficient (default: 0.1)')
    parser.add_argument('--num-hinges-coef', type=float, default=0.1,
                        help='Hinge count penalty coefficient (default: 0.1)')

    # Ladder specification (two ways)
    parser.add_argument('--vdv-ladder', type=str, default=None,
                        help='Explicit comma-separated VDV coefficient ladder')
    parser.add_argument('--vdv-min', type=float, default=0.5,
                        help='Minimum VDV coefficient (default: 0.5)')
    parser.add_argument('--vdv-max', type=float, default=20000,
                        help='Maximum VDV coefficient (default: 20000)')
    parser.add_argument('--num-replicas', type=int, default=8,
                        help='Number of replicas (default: 8)')

    # Sampling parameters
    parser.add_argument('--sweeps-per-exchange', type=float, default=5,
                        help='Sweeps between swap attempts (default: 5)')
    parser.add_argument('--num-exchanges', type=int, default=1000,
                        help='Total number of exchange rounds (default: 1000)')
    parser.add_argument('--log-interval', type=int, default=50,
                        help='Log stats every N exchanges (default: 50)')

    # Output
    parser.add_argument('--output-dir', default='data/replica_exchange',
                        help='Directory for output files')
    parser.add_argument('--save-interval', type=int, default=0,
                        help='Save top replica every N exchanges (0=only at end)')
    parser.add_argument('--min-free-memory-gb', type=float, default=4.0,
                        help='Abort if free memory drops below this (default: 4 GB)')

    args = parser.parse_args()

    global MIN_FREE_MEMORY_GB
    MIN_FREE_MEMORY_GB = args.min_free_memory_gb

    dim = args.dim
    check_memory("pre-flight")

    # --- Build VDV ladder ---
    if args.vdv_ladder:
        ladder = sorted(float(x) for x in args.vdv_ladder.split(','))
    else:
        # Geometric spacing
        K = args.num_replicas
        ratio = (args.vdv_max / args.vdv_min) ** (1.0 / (K - 1))
        ladder = [args.vdv_min * ratio ** i for i in range(K)]

    K = len(ladder)
    print(f"Replica exchange with {K} replicas")
    print(f"VDV ladder: {', '.join(f'{v:.1f}' for v in ladder)}")

    # --- Memory estimate ---
    # Rough: 36 MB base + 2.15 MB per 1K facets, times K replicas
    m_probe = ddg.Manifold.load(args.seed_file, dim)
    size = m_probe.num_facets
    del m_probe
    ddg.gc_collect()

    per_replica_mb = 36 + 2.15 * size / 1000
    total_mb = per_replica_mb * K
    free_gb = get_free_memory_gb()
    print(f"Estimated memory: {total_mb:.0f} MB ({per_replica_mb:.0f} MB x {K}), "
          f"free: {free_gb:.1f} GB")
    if total_mb / 1024 > free_gb - MIN_FREE_MEMORY_GB:
        print(f"ERROR: estimated {total_mb/1024:.1f} GB needed but only "
              f"{free_gb:.1f} GB free ({MIN_FREE_MEMORY_GB:.0f} GB reserved)")
        sys.exit(1)

    # --- Create replicas ---
    # Each replica loads its own copy of the seed
    samplers = []
    for i, vdv in enumerate(ladder):
        check_memory(f"creating replica {i}")
        print(f"  Replica {i}: loading seed, VDV={vdv:.1f}...")
        m = ddg.Manifold.load(args.seed_file, dim)
        params = ddg.SamplerParams(
            num_facets_target=m.num_facets,
            hinge_degree_target=args.hinge_target,
            num_facets_coef=args.num_facets_coef,
            num_hinges_coef=args.num_hinges_coef,
            hinge_degree_variance_coef=0.0,
            codim3_degree_variance_coef=vdv,
        )
        sampler = ddg.ManifoldSampler(m, params)
        samplers.append(sampler)
        ddg.gc_collect()

    # beta[i] is the current VDV coefficient assigned to sampler i
    beta = list(ladder)

    # --- Prepare output ---
    os.makedirs(args.output_dir, exist_ok=True)
    tag = f'{dim}d_N{size}'
    csv_path = os.path.join(args.output_dir, f'{tag}_replica_exchange.csv')
    csv_file = open(csv_path, 'w', newline='')
    writer = csv.writer(csv_file)
    writer.writerow([
        'exchange', 'elapsed_s',
        *[f'beta_{i}' for i in range(K)],
        *[f'accept_{i}' for i in range(K)],
        *[f'degvar0_{i}' for i in range(K)],
        *[f'degvar_hinge_{i}' for i in range(K)],
        *[f'swap_accept_{i}_{i+1}' for i in range(K - 1)],
    ])

    # Swap statistics
    swap_tried = np.zeros(K - 1, dtype=np.int64)
    swap_accepted = np.zeros(K - 1, dtype=np.int64)

    # --- Main loop ---
    t_start = time.monotonic()

    # Print header
    print(f"\n{'Exch':>6} {'Time':>7} ", end='')
    for i in range(K):
        print(f"{'R'+str(i)+' β':>10}", end='')
    print(f"  {'TopDV0':>8} {'TopDVh':>8} {'TopAcc%':>8} {'Swaps':>20}")
    print("-" * (30 + 10 * K + 50))

    for exch in range(args.num_exchanges):
        check_memory(f"exchange {exch}")

        # 1. Run each replica for some sweeps
        for i, sampler in enumerate(samplers):
            sampler.reset_stats()
            sampler.run(sweeps=args.sweeps_per_exchange)

        # 2. Attempt swaps between adjacent replicas (even/odd alternation)
        parity = exch % 2
        for i in range(parity, K - 1, 2):
            j = i + 1
            p_i = codim3_penalty(samplers[i].manifold, dim)
            p_j = codim3_penalty(samplers[j].manifold, dim)

            swap_tried[i] += 1
            if swap_accept(beta[i], beta[j], p_i, p_j):
                swap_accepted[i] += 1
                # Swap VDV coefficients between samplers i and j
                beta[i], beta[j] = beta[j], beta[i]
                samplers[i].set_codim3_degree_variance_coef(beta[i])
                samplers[j].set_codim3_degree_variance_coef(beta[j])

        # 3. Logging
        if (exch + 1) % args.log_interval == 0 or exch == 0:
            elapsed = time.monotonic() - t_start

            # Find which sampler currently holds the highest VDV coefficient
            top_idx = max(range(K), key=lambda i: beta[i])
            top_mfd = samplers[top_idx].manifold
            top_dv0 = top_mfd.degree_variance(0)
            top_dvh = top_mfd.degree_variance(max(0, dim - 2))
            top_stats = samplers[top_idx].get_stats()
            top_acc = (top_stats.total_accepted / top_stats.total_tried * 100
                       if top_stats.total_tried > 0 else 0)

            # Per-replica stats
            accept_rates = []
            degvars_0 = []
            degvars_h = []
            for i, sampler in enumerate(samplers):
                st = sampler.get_stats()
                acc = st.total_accepted / st.total_tried if st.total_tried > 0 else 0
                accept_rates.append(acc)
                degvars_0.append(sampler.manifold.degree_variance(0))
                degvars_h.append(sampler.manifold.degree_variance(max(0, dim - 2)))

            # Swap rates for log
            swap_strs = []
            for i in range(K - 1):
                if swap_tried[i] > 0:
                    r = swap_accepted[i] / swap_tried[i]
                    swap_strs.append(f"{100*r:.0f}%")
                else:
                    swap_strs.append("--")

            # CSV row
            writer.writerow([
                exch + 1, f'{elapsed:.1f}',
                *[f'{beta[i]:.2f}' for i in range(K)],
                *[f'{accept_rates[i]:.4f}' for i in range(K)],
                *[f'{degvars_0[i]:.2f}' for i in range(K)],
                *[f'{degvars_h[i]:.2f}' for i in range(K)],
                *swap_strs,
            ])
            csv_file.flush()

            # Console
            print(f"{exch+1:6d} {elapsed:6.0f}s ", end='')
            for i in range(K):
                print(f"{beta[i]:10.1f}", end='')
            print(f"  {top_dv0:8.2f} {top_dvh:8.2f} {top_acc:7.1f}%", end='')
            print(f"  swaps: {'/'.join(swap_strs)}")

        # 4. Periodic save
        if args.save_interval > 0 and (exch + 1) % args.save_interval == 0:
            top_idx = max(range(K), key=lambda i: beta[i])
            mfd_path = os.path.join(args.output_dir,
                                    f'{tag}_re_exch{exch+1}.mfd')
            samplers[top_idx].manifold.dup().save(mfd_path)
            print(f"  -> Saved checkpoint: {mfd_path}")

    csv_file.close()
    print(f"\nResults saved to {csv_path}")

    # --- Final save ---
    top_idx = max(range(K), key=lambda i: beta[i])
    top_mfd = samplers[top_idx].manifold
    top_beta = beta[top_idx]
    mfd_path = os.path.join(args.output_dir, f'{tag}_re_final.mfd')
    top_mfd.dup().save(mfd_path)
    print(f"Final manifold saved to {mfd_path}")
    print(f"  VDV coef: {top_beta:.1f}")
    print(f"  Vertex degree variance: {top_mfd.degree_variance(0):.4f}")
    print(f"  Hinge degree variance: {top_mfd.degree_variance(max(0, dim - 2)):.4f}")
    print(f"  Mean edge degree: {top_mfd.mean_degree(1):.4f}")
    print(f"  Facets: {top_mfd.num_facets}")

    # Swap summary
    print(f"\nSwap acceptance rates:")
    for i in range(K - 1):
        if swap_tried[i] > 0:
            r = 100 * swap_accepted[i] / swap_tried[i]
            print(f"  {ladder[i]:.1f} <-> {ladder[i+1]:.1f}: "
                  f"{swap_accepted[i]}/{swap_tried[i]} ({r:.1f}%)")

    ddg.gc_collect()
    ddg.gc_minimize()


if __name__ == '__main__':
    main()
