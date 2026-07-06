#!/usr/bin/env python3
"""Ramp codim-3 degree variance coefficient until acceptance drops to zero.

Grows a manifold at fixed hinge_degree_target, then incrementally increases
codim3_degree_variance_coef, equilibrating at each step and logging the
acceptance rate. Stops when acceptance falls below a threshold.

The goal is to find the maximum VDV coefficient that still allows moves,
producing samples with minimal degree variance.

Usage:
    python scripts/ramp_vdv.py --size 1000
    python scripts/ramp_vdv.py --size 10000 --vdv-start 0.0 --vdv-step 0.05 --vdv-max 5.0
    python scripts/ramp_vdv.py --size 1000 --eq-sweeps 200 --probe-sweeps 50

Output:
    Prints a table of (vdv_coef, acceptance_rate, degree_variance, objective)
    and saves the final manifold when acceptance drops below --min-accept.
"""

import argparse
import csv
import os
import sys
import time

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


def main():
    parser = argparse.ArgumentParser(
        description='Ramp VDV coefficient until acceptance drops to zero')
    parser.add_argument('--seed-file', type=str, default=None,
                        help='Path to an existing .mfd seed file to load')
    parser.add_argument('--size', type=int, default=1000,
                        help='Target number of facets when growing from scratch (default: 1000)')
    parser.add_argument('--dim', type=int, default=3,
                        help='Manifold dimension (default: 3)')
    parser.add_argument('--hinge-target', type=float, default=5.0043,
                        help='Hinge degree target (default: 5.0043)')
    parser.add_argument('--num-facets-coef', type=float, default=0.1,
                        help='Volume penalty coefficient (default: 0.1)')
    parser.add_argument('--num-hinges-coef', type=float, default=0.05,
                        help='Hinge count penalty coefficient (default: 0.05)')
    parser.add_argument('--vdv-start', type=float, default=0.0,
                        help='Starting VDV coefficient (default: 0.0)')
    parser.add_argument('--vdv-step', type=float, default=0.05,
                        help='VDV coefficient additive increment (default: 0.05)')
    parser.add_argument('--vdv-factor', type=float, default=None,
                        help='VDV coefficient multiplicative factor (overrides --vdv-step)')
    parser.add_argument('--vdv-max', type=float, default=10.0,
                        help='Maximum VDV coefficient to try (default: 10.0)')
    parser.add_argument('--eq-sweeps', type=int, default=100,
                        help='Equilibration sweeps at each VDV step (default: 100)')
    parser.add_argument('--probe-sweeps', type=int, default=20,
                        help='Sweeps to measure acceptance rate after equilibrating (default: 20)')
    parser.add_argument('--min-accept', type=float, default=0.01,
                        help='Stop when acceptance rate falls below this (default: 0.01 = 1%%)')
    parser.add_argument('--output-dir', default='data/vdv_ramp',
                        help='Directory for output files (default: data/vdv_ramp)')
    parser.add_argument('--save-manifold', action='store_true',
                        help='Save the manifold at the last good VDV step')
    parser.add_argument('--min-free-memory-gb', type=float, default=4.0,
                        help='Abort if free memory drops below this (default: 4 GB)')
    args = parser.parse_args()

    global MIN_FREE_MEMORY_GB
    MIN_FREE_MEMORY_GB = args.min_free_memory_gb

    dim = args.dim
    check_memory("pre-flight")

    # --- Load or grow manifold ---
    if args.seed_file:
        print(f"Loading seed from {args.seed_file}...")
        m = ddg.Manifold.load(args.seed_file, dim)
        size = m.num_facets
        print(f"Loaded {size} facets")
    else:
        size = args.size
        print(f"Growing {dim}d manifold to {size} facets (hinge_target={args.hinge_target})...")
        m = ddg.Manifold.standard_sphere(dim)

    params = ddg.SamplerParams(
        num_facets_target=size,
        hinge_degree_target=args.hinge_target,
        num_facets_coef=args.num_facets_coef,
        num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=args.vdv_start,
    )
    sampler = ddg.ManifoldSampler(m, params)

    if not args.seed_file:
        step_size = max(100, size // 20)

        def growth_callback(cur, tgt):
            check_memory(f"growth {cur}/{tgt}")

        sampler.ramped_grow(size, step_size=step_size, eq_sweeps_per_step=3,
                            callback=growth_callback)
        print(f"Grown to {sampler.manifold.num_facets} facets")

    # --- Initial equilibration at starting VDV ---
    print(f"Equilibrating at VDV={args.vdv_start:.3f} for {args.eq_sweeps} sweeps...")
    sampler.run(sweeps=args.eq_sweeps)

    # --- Prepare output ---
    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, f'{dim}d_N{size}_vdv_ramp.csv')
    csv_file = open(csv_path, 'w', newline='')
    writer = csv.writer(csv_file)
    writer.writerow([
        'vdv_coef', 'acceptance_rate', 'total_tried', 'total_accepted',
        'num_facets', 'degree_variance_0', 'degree_variance_hinge',
        'mean_degree_0', 'mean_degree_hinge', 'objective',
    ])

    print(f"\n{'VDV coef':>10} {'Accept%':>9} {'DegVar(vtx)':>12} "
          f"{'DegVar(hinge)':>14} {'Objective':>10} {'Facets':>8}")
    print("-" * 72)

    last_good_vdv = args.vdv_start
    vdv = args.vdv_start

    while vdv <= args.vdv_max + 1e-9:
        check_memory(f"VDV={vdv:.3f}")

        # Set the new VDV coefficient
        sampler.set_codim3_degree_variance_coef(vdv)

        # Equilibrate
        sampler.run(sweeps=args.eq_sweeps)

        # Probe: measure acceptance over a clean window
        sampler.reset_stats()
        sampler.run(sweeps=args.probe_sweeps)
        stats = sampler.get_stats()

        tried = stats.total_tried
        accepted = stats.total_accepted
        accept_rate = accepted / tried if tried > 0 else 0.0

        mfd = sampler.manifold
        nf = mfd.num_facets
        dv0 = mfd.degree_variance(0)
        dv_hinge = mfd.degree_variance(max(0, dim - 2))
        md0 = mfd.mean_degree(0)
        md_hinge = mfd.mean_degree(max(0, dim - 2))
        obj = sampler.current_objective

        # Log
        writer.writerow([
            f'{vdv:.4f}', f'{accept_rate:.6f}', tried, accepted,
            nf, f'{dv0:.4f}', f'{dv_hinge:.4f}',
            f'{md0:.4f}', f'{md_hinge:.4f}', f'{obj:.4f}',
        ])
        csv_file.flush()

        print(f"{vdv:10.4f} {100*accept_rate:8.2f}% {dv0:12.2f} "
              f"{dv_hinge:14.2f} {obj:10.1f} {nf:8d}")

        if accept_rate < args.min_accept:
            print(f"\nAcceptance rate {100*accept_rate:.2f}% < {100*args.min_accept:.1f}% "
                  f"at VDV={vdv:.4f}. Stopping.")
            print(f"Last good VDV: {last_good_vdv:.4f}")
            break

        last_good_vdv = vdv
        if args.vdv_factor is not None:
            vdv *= args.vdv_factor
        else:
            vdv += args.vdv_step

    csv_file.close()
    print(f"\nResults saved to {csv_path}")

    # Optionally save the manifold at the last good VDV
    if args.save_manifold:
        # Re-equilibrate at the last good VDV
        sampler.set_codim3_degree_variance_coef(last_good_vdv)
        sampler.run(sweeps=args.eq_sweeps)
        mfd_path = os.path.join(args.output_dir,
                                f'{dim}d_N{size}_vdv{last_good_vdv:.4f}.mfd')
        sampler.manifold.dup().save(mfd_path)
        print(f"Manifold saved to {mfd_path}")

    ddg.gc_collect()
    ddg.gc_minimize()


if __name__ == '__main__':
    main()
