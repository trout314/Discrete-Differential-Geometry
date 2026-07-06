#!/usr/bin/env python3
"""Replica exchange (parallel tempering) over the VDV coefficient.

Runs K replicas at different codim3_degree_variance_coef (VDV) values and
attempts swaps between adjacent replicas so that the high-VDV replica mixes
via excursions to lower VDV. This lets us sample low-vertex-degree-variance
triangulations that a single high-VDV chain cannot reach (acceptance there
collapses to ~0).

Two things make this version actually mix, unlike a naive geometric ladder:

1. Thermodynamic-length ladder. A geometric ladder in VDV starves swaps at
   BOTH ends: at the bottom Delta(beta) is tiny but the vertex-degree-variance
   penalty E fluctuates hugely; at the top E barely moves but Delta(beta) is
   enormous. Healthy swaps need Delta(beta) * sigma_E ~ const. A short pilot
   measures sigma_E(beta) and rungs are placed at equal thermodynamic length
   L(beta) = integral sigma_E dbeta, which equalizes swap acceptance.

2. Correct adjacency. Swaps are attempted between adjacent *rungs* (temperature
   values), tracked explicitly, not adjacent array slots (which drift out of
   temperature order as betas shuffle between replicas).

Sampling / uniformity. The default sampler is pure Metropolis with no Hastings
correction, so its stationary distribution is exp(-objective) * V(x), where
V(x) is the number of valid Pachner moves. To recover uniform-within-ensemble
statistics, each recorded sample carries weight 1/V(x) (importance_weight());
reweight observables afterward. The diagnostics report both the autocorrelation
ESS and the weight-variance (Kish) ESS so the true independent-sample count is
visible.

Usage:
    # Auto-sized thermodynamic ladder targeting ~23% swaps
    python scripts/replica_exchange_vdv.py \
        --seed-file seeds/S3_N1e5_..._VDV_5e-1_s000.mfd \
        --vdv-min 0.5 --vdv-max 20000 --target-swap 0.23

    # Fix the replica count instead of auto-sizing
    python scripts/replica_exchange_vdv.py \
        --seed-file seeds/S3_N1e4_..._s000.mfd \
        --vdv-min 0.5 --vdv-max 2000 --num-replicas 12

    # Skip the pilot and use an explicit ladder
    python scripts/replica_exchange_vdv.py \
        --seed-file seeds/... --vdv-ladder 0.5,2,8,32,128,512,2048,8192
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
from discrete_differential_geometry.convergence import (
    integrated_autocorrelation_time,
    weighted_ess,
)


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
    """Intensive codim-3 degree variance penalty E, matching sampler.d.

    The sampler's objective contribution is coef * E where (from
    localSolidAngleCurvPenalty, divided by nCoDim3):
        E = degree_variance(codim3_dim) - min_penalty
    with min_penalty = x*(1-x), x = frac(coDim3DegTarget). For dim=3 the
    codim-3 faces are vertices, so E tracks vertex degree variance.
    """
    if dim <= 2:
        return 0.0

    codim3_dim = dim - 3
    fv = mfd.f_vector
    n_facets = fv[dim]
    n_codim3 = fv[codim3_dim]
    if n_codim3 == 0:
        return 0.0

    c3pf = math.comb(dim + 1, dim - 2)
    deg_target = c3pf * n_facets / n_codim3
    x = deg_target - math.floor(deg_target)
    min_penalty = x * (1 - x)

    return mfd.degree_variance(codim3_dim) - min_penalty


def swap_accept(beta_i, beta_j, penalty_i, penalty_j):
    """Replica-exchange swap acceptance for swapping VDV coefficients.

    Only the codim-3 term carries beta and differs between replicas; all other
    objective terms share identical coefficients and cancel. So
        Delta = (beta_j - beta_i) * (penalty_i - penalty_j)
    and the swap is accepted with probability min(1, exp(-Delta)).
    """
    delta = (beta_j - beta_i) * (penalty_i - penalty_j)
    if delta <= 0:
        return True
    return np.random.random() < math.exp(-delta)


# --- Thermodynamic-length ladder construction ---

def _erfcinv(y):
    """Inverse complementary error function for y in (0, 2), via bisection."""
    if y <= 0:
        return float('inf')
    if y >= 2:
        return float('-inf')
    lo, hi = -6.0, 6.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        # erfc is decreasing in its argument
        if math.erfc(mid) > y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def target_delta_for_swap(target_swap):
    """Per-rung thermodynamic length delta = Delta(beta)*sigma_E that yields
    roughly the requested mean swap acceptance.

    Rule of thumb (Gaussian energy, adjacent replicas): acceptance ~
    erfc(delta / 2), so delta = 2 * erfcinv(target_swap). At 23% -> ~1.7.
    The achieved acceptance is measured and reported; this only sizes the
    rung count.
    """
    return 2.0 * _erfcinv(target_swap)


def pilot_sigma(seed_file, dim, base_params, beta_grid,
                eq_sweeps, sample_sweeps, n_samples):
    """Warm-start ramp over beta_grid measuring sigma_E(beta).

    Reuses a single manifold, ramping beta upward and equilibrating at each
    grid point (a warm start is far cheaper than an independent chain per
    point). Returns (betas, sigma_E, mean_E) as numpy arrays.
    """
    m = ddg.Manifold.load(seed_file, dim)
    params = ddg.SamplerParams(
        num_facets_target=m.num_facets,
        hinge_degree_target=base_params['hinge_target'],
        num_facets_coef=base_params['num_facets_coef'],
        num_hinges_coef=base_params['num_hinges_coef'],
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=beta_grid[0],
    )
    sampler = ddg.ManifoldSampler(m, params)

    betas, sig, mean = [], [], []
    print(f"  {'beta':>10} {'mean_E':>10} {'sigma_E':>10} {'accept%':>8}")
    for b in beta_grid:
        check_memory(f"pilot beta={b:.2f}")
        sampler.set_codim3_degree_variance_coef(b)
        sampler.run(sweeps=eq_sweeps)

        sampler.reset_stats()
        samples = []
        for _ in range(n_samples):
            sampler.run(sweeps=sample_sweeps)
            samples.append(codim3_penalty(sampler.manifold, dim))
        st = sampler.get_stats()
        acc = st.total_accepted / st.total_tried if st.total_tried > 0 else 0.0

        samples = np.asarray(samples)
        betas.append(b)
        mean.append(float(samples.mean()))
        sig.append(float(samples.std()))
        print(f"  {b:10.2f} {samples.mean():10.3f} {samples.std():10.4f} {100*acc:7.1f}%")

    del sampler, m
    ddg.gc_collect()
    return np.asarray(betas), np.asarray(sig), np.asarray(mean)


def build_ladder(betas, sigma, vdv_min, vdv_max, target_delta, num_replicas=None):
    """Place rungs at equal thermodynamic length L(beta) = integral sigma dbeta.

    Equal spacing in L means equal Delta(beta)*sigma_E per rung, hence roughly
    constant swap acceptance. If num_replicas is given it is used directly;
    otherwise the count is chosen so each rung spans ~target_delta of length.

    If the penalty freezes (sigma_E -> 0) below vdv_max -- the triangulation
    has crystallized into a rigid low-variance state -- the ladder is capped at
    that freeze point: rungs above it would all sample the same frozen state.

    Returns (ladder, total_L, K, beta_top) where beta_top is the effective top
    rung (== vdv_max unless a freeze was detected).
    """
    betas = np.asarray(betas, dtype=float)
    sigma = np.asarray(sigma, dtype=float)

    # Detect the freeze point: the last grid point with non-negligible
    # fluctuation. Everything above it is rigid and not worth a rung.
    thresh = max(1e-9, 1e-3 * sigma.max())
    live = np.where(sigma > thresh)[0]
    if len(live) < 1:
        # No fluctuation anywhere -- degenerate; just span the endpoints.
        return [vdv_min, vdv_max], 0.0, 2, vdv_max

    # Include one grid point past the last live one so the top rung sits at the
    # onset of freezing rather than short of it, but never above vdv_max.
    hi = min(live[-1] + 1, len(betas) - 1)
    b = betas[:hi + 1].copy()
    s = np.maximum(sigma[:hi + 1], 1e-9)
    beta_top = min(vdv_max, float(b[-1]))
    b[0] = vdv_min
    b[-1] = beta_top

    # Cumulative thermodynamic length via the trapezoid rule.
    L = np.zeros(len(b))
    for k in range(1, len(b)):
        L[k] = L[k - 1] + 0.5 * (s[k] + s[k - 1]) * (b[k] - b[k - 1])
    total_L = L[-1]

    if num_replicas is not None:
        K = max(2, num_replicas)
    else:
        K = max(2, int(round(total_L / target_delta)) + 1)

    targets = np.linspace(0.0, total_L, K)
    ladder = np.interp(targets, L, b)
    ladder[0], ladder[-1] = vdv_min, beta_top

    # Drop any duplicate rungs (can happen if L saturates near the top).
    ladder = np.unique(np.round(ladder, 6))
    return [float(v) for v in ladder], float(total_L), len(ladder), beta_top


def main():
    parser = argparse.ArgumentParser(
        description='Replica exchange over VDV coefficient (thermodynamic ladder)')
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

    # Ladder specification
    parser.add_argument('--vdv-ladder', type=str, default=None,
                        help='Explicit comma-separated VDV ladder (skips the pilot)')
    parser.add_argument('--vdv-min', type=float, default=0.5,
                        help='Minimum VDV coefficient (default: 0.5)')
    parser.add_argument('--vdv-max', type=float, default=2000.0,
                        help='Maximum VDV coefficient (default: 2000)')
    parser.add_argument('--num-replicas', type=int, default=None,
                        help='Fix replica count (default: auto-size from --target-swap)')
    parser.add_argument('--target-swap', type=float, default=0.23,
                        help='Target mean swap acceptance for ladder sizing (default: 0.23)')

    # Pilot
    parser.add_argument('--pilot-grid', type=int, default=16,
                        help='Number of geometric pilot beta points (default: 16)')
    parser.add_argument('--pilot-eq-sweeps', type=int, default=20,
                        help='Equilibration sweeps per pilot point (default: 20)')
    parser.add_argument('--pilot-sample-sweeps', type=int, default=2,
                        help='Sweeps between pilot E samples (default: 2)')
    parser.add_argument('--pilot-samples', type=int, default=30,
                        help='E samples per pilot point (default: 30)')

    # Sampling parameters
    parser.add_argument('--sweeps-per-exchange', type=float, default=5,
                        help='Sweeps between swap attempts (default: 5)')
    parser.add_argument('--num-exchanges', type=int, default=1000,
                        help='Total number of exchange rounds (default: 1000)')
    parser.add_argument('--warmup-exchanges', type=int, default=None,
                        help='Exchanges to discard before recording (default: 20%% of total)')
    parser.add_argument('--target-rung', type=int, default=-1,
                        help='Rung whose ensemble is recorded for ESS (default: -1 = top)')
    parser.add_argument('--log-interval', type=int, default=50,
                        help='Log stats every N exchanges (default: 50)')
    parser.add_argument('--gc-interval', type=int, default=10,
                        help='Reclaim D-side GC garbage every N exchanges (default: 10). '
                             'Essential for long runs: each exchange generates '
                             'move temporaries that otherwise pile up unbounded.')

    # Output
    parser.add_argument('--output-dir', default='data/replica_exchange',
                        help='Directory for output files')
    parser.add_argument('--save-interval', type=int, default=0,
                        help='Save target replica every N exchanges (0=only at end)')
    parser.add_argument('--min-free-memory-gb', type=float, default=4.0,
                        help='Abort if free memory drops below this (default: 4 GB)')

    args = parser.parse_args()

    global MIN_FREE_MEMORY_GB
    MIN_FREE_MEMORY_GB = args.min_free_memory_gb

    dim = args.dim
    check_memory("pre-flight")

    base_params = dict(
        hinge_target=args.hinge_target,
        num_facets_coef=args.num_facets_coef,
        num_hinges_coef=args.num_hinges_coef,
    )

    # --- Build the VDV ladder ---
    if args.vdv_ladder:
        ladder = sorted(float(x) for x in args.vdv_ladder.split(','))
        print(f"Using explicit ladder ({len(ladder)} rungs), skipping pilot.")
    else:
        print(f"Pilot: measuring sigma_E over {args.pilot_grid} beta points "
              f"[{args.vdv_min:g}, {args.vdv_max:g}]...")
        grid = np.geomspace(args.vdv_min, args.vdv_max, args.pilot_grid)
        betas, sigma, mean_E = pilot_sigma(
            args.seed_file, dim, base_params, grid,
            args.pilot_eq_sweeps, args.pilot_sample_sweeps, args.pilot_samples)
        target_delta = target_delta_for_swap(args.target_swap)
        ladder, total_L, K, beta_top = build_ladder(
            betas, sigma, args.vdv_min, args.vdv_max,
            target_delta, args.num_replicas)
        print(f"\n  Total thermodynamic length L = {total_L:.2f}, "
              f"target delta/rung = {target_delta:.2f} "
              f"(~{100*args.target_swap:.0f}% swaps) -> {K} rungs")
        if beta_top < args.vdv_max * (1 - 1e-6):
            print(f"  NOTE: penalty froze at VDV~{beta_top:.2f} "
                  f"(< requested max {args.vdv_max:g}); ladder capped there. "
                  f"The triangulation crystallizes above this coefficient.")

    K = len(ladder)
    print(f"\nReplica exchange with {K} replicas")
    print(f"VDV ladder: {', '.join(f'{v:.2f}' for v in ladder)}")

    # --- Memory estimate ---
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
    samplers = []
    for i, vdv in enumerate(ladder):
        check_memory(f"creating replica {i}")
        print(f"  Replica {i}: loading seed, VDV={vdv:.2f}...")
        m = ddg.Manifold.load(args.seed_file, dim)
        params = ddg.SamplerParams(
            num_facets_target=m.num_facets,
            hinge_degree_target=args.hinge_target,
            num_facets_coef=args.num_facets_coef,
            num_hinges_coef=args.num_hinges_coef,
            hinge_degree_variance_coef=0.0,
            codim3_degree_variance_coef=vdv,
        )
        samplers.append(ddg.ManifoldSampler(m, params))
        ddg.gc_collect()

    # Rung bookkeeping. slot i holds temperature beta[i]; rung_of_slot/slot_of_rung
    # track the mapping so swaps are always between adjacent *rungs*, not slots.
    beta = list(ladder)
    rung_of_slot = list(range(K))
    slot_of_rung = list(range(K))

    target_rung = args.target_rung % K
    warmup = args.warmup_exchanges
    if warmup is None:
        warmup = max(1, args.num_exchanges // 5)

    # Round-trip tracking: each slot is a walker labelled by the last ladder
    # extreme it touched. A walker only becomes 'down' via a genuine bottom->top
    # traversal, so a bottom->top->bottom loop counts exactly one round trip.
    # Only the bottom walker is seeded 'up'; others start unlabelled so that a
    # replica beginning at the top does not fabricate a trip it never made.
    walker_label = [None] * K   # 'up' (came from bottom), 'down' (reached top)
    walker_label[slot_of_rung[0]] = 'up'
    up_trips = 0                 # completed bottom->top traversals
    round_trips = 0             # completed bottom->top->bottom loops

    # Recorded time series at the target rung (post-warmup). The importance
    # weight 1/V(x) is computed exactly once per kept sample here (it needs
    # V(x) at that configuration and cannot be deferred without storing every
    # manifold); the reweighting itself happens post-run in the diagnostics.
    # last_weight caches the most recent kept weight so the periodic CSV log
    # can report it without a second (needless) countValidBistellarMoves pass.
    rec_E, rec_weight, rec_meandeg = [], [], []
    last_weight = float('nan')

    swap_tried = np.zeros(K - 1, dtype=np.int64)
    swap_accepted = np.zeros(K - 1, dtype=np.int64)

    # --- Prepare output ---
    os.makedirs(args.output_dir, exist_ok=True)
    tag = f'{dim}d_N{size}'
    csv_path = os.path.join(args.output_dir, f'{tag}_replica_exchange.csv')
    csv_file = open(csv_path, 'w', newline='')
    writer = csv.writer(csv_file)
    writer.writerow([
        'exchange', 'elapsed_s', 'target_dv0', 'target_meandeg',
        'target_weight', 'up_trips', 'round_trips',
        *[f'swap_accept_{r}_{r+1}' for r in range(K - 1)],
    ])

    t_start = time.monotonic()
    print(f"\n{'Exch':>6} {'Time':>7} {'tgtDV0':>8} {'tgtMdeg':>8} "
          f"{'swap(min/mean/max)':>20} {'rtrips':>7}")
    print("-" * 68)

    for exch in range(args.num_exchanges):
        check_memory(f"exchange {exch}")

        # 1. Advance each replica.
        for sampler in samplers:
            sampler.reset_stats()
            sampler.run(sweeps=args.sweeps_per_exchange)

        # 2. Attempt swaps between adjacent rungs (even/odd alternation).
        parity = exch % 2
        for r in range(parity, K - 1, 2):
            a = slot_of_rung[r]      # slot currently at lower rung r
            b = slot_of_rung[r + 1]  # slot currently at higher rung r+1
            p_a = codim3_penalty(samplers[a].manifold, dim)
            p_b = codim3_penalty(samplers[b].manifold, dim)

            swap_tried[r] += 1
            if swap_accept(beta[a], beta[b], p_a, p_b):
                swap_accepted[r] += 1
                beta[a], beta[b] = beta[b], beta[a]
                rung_of_slot[a], rung_of_slot[b] = r + 1, r
                slot_of_rung[r], slot_of_rung[r + 1] = b, a
                samplers[a].set_codim3_degree_variance_coef(beta[a])
                samplers[b].set_codim3_degree_variance_coef(beta[b])

        # 3. Round-trip accounting on the two extreme rungs. A walker becomes
        #    'down' only by reaching the top while labelled 'up' (came from the
        #    bottom); reaching the bottom while 'down' closes a round trip.
        bottom_slot = slot_of_rung[0]
        top_slot = slot_of_rung[K - 1]
        if walker_label[top_slot] == 'up':
            up_trips += 1
            walker_label[top_slot] = 'down'
        if walker_label[bottom_slot] == 'down':
            round_trips += 1
        if walker_label[bottom_slot] is None or walker_label[bottom_slot] == 'down':
            # First visit to the bottom (or return after a completed trip)
            # (re)starts the cycle.
            walker_label[bottom_slot] = 'up'

        # 4. Record the target-rung ensemble (post-warmup).
        if exch >= warmup:
            tslot = slot_of_rung[target_rung]
            tmfd = samplers[tslot].manifold
            last_weight = tmfd.importance_weight()
            rec_E.append(codim3_penalty(tmfd, dim))
            rec_weight.append(last_weight)
            rec_meandeg.append(tmfd.mean_degree(0))

        # 5. Logging.
        if (exch + 1) % args.log_interval == 0 or exch == 0:
            elapsed = time.monotonic() - t_start
            tslot = slot_of_rung[target_rung]
            tmfd = samplers[tslot].manifold
            tdv0 = tmfd.degree_variance(0)
            tmdeg = tmfd.mean_degree(0)
            tweight = last_weight  # most recent kept-sample weight (nan pre-warmup)

            rates = np.where(swap_tried > 0, swap_accepted / np.maximum(swap_tried, 1), 0.0)
            writer.writerow([
                exch + 1, f'{elapsed:.1f}', f'{tdv0:.4f}', f'{tmdeg:.4f}',
                f'{tweight:.6e}', up_trips, round_trips,
                *[f'{r:.4f}' for r in rates],
            ])
            csv_file.flush()

            print(f"{exch+1:6d} {elapsed:6.0f}s {tdv0:8.2f} {tmdeg:8.3f} "
                  f"  {100*rates.min():4.0f}/{100*rates.mean():4.0f}/{100*rates.max():4.0f}%   "
                  f"{round_trips:7d}")

        # 6. Reclaim accumulated move temporaries so memory stays bounded.
        if args.gc_interval > 0 and (exch + 1) % args.gc_interval == 0:
            ddg.gc_collect()

        # 7. Periodic checkpoint save.
        if args.save_interval > 0 and (exch + 1) % args.save_interval == 0:
            tslot = slot_of_rung[target_rung]
            mfd_path = os.path.join(args.output_dir, f'{tag}_re_exch{exch+1}.mfd')
            samplers[tslot].manifold.dup().save(mfd_path)
            print(f"  -> Saved checkpoint: {mfd_path}")

    csv_file.close()
    elapsed = time.monotonic() - t_start

    # --- Diagnostics ---
    print(f"\n{'='*68}")
    print("Diagnostics")
    print(f"{'='*68}")

    # Per-rung swap acceptance.
    rates = np.where(swap_tried > 0, swap_accepted / np.maximum(swap_tried, 1), 0.0)
    print(f"\nSwap acceptance per rung (target ~{100*args.target_swap:.0f}%):")
    for r in range(K - 1):
        bar = '#' * int(40 * rates[r])
        print(f"  {ladder[r]:9.2f} <-> {ladder[r+1]:9.2f}: "
              f"{100*rates[r]:5.1f}%  {bar}")
    print(f"  min/mean/max: {100*rates.min():.1f}% / "
          f"{100*rates.mean():.1f}% / {100*rates.max():.1f}%")

    # Round trips / mixing.
    print(f"\nMixing:")
    print(f"  up-traversals (bottom->top):     {up_trips}")
    print(f"  round trips (bottom->top->bot):  {round_trips}")
    if round_trips > 0:
        rt_time = (args.num_exchanges - warmup) / round_trips
        print(f"  mean round-trip time:            {rt_time:.0f} exchanges")
    else:
        print(f"  mean round-trip time:            n/a (no full round trips)")

    # Target-rung ESS (autocorrelation x weight variance).
    print(f"\nTarget rung {target_rung} (VDV={ladder[target_rung]:.2f}), "
          f"{len(rec_E)} recorded samples:")
    if len(rec_E) >= 4:
        E = np.asarray(rec_E)
        w = np.asarray(rec_weight)
        md = np.asarray(rec_meandeg)

        tau = integrated_autocorrelation_time(E)
        ess_ac = len(E) / tau if tau and tau == tau else float('nan')
        ess_w = weighted_ess(w)
        # Combined: autocorrelation thins the chain, weights further deflate it.
        ess_combined = ess_ac * (ess_w / len(w)) if len(w) else float('nan')

        # Reweighted mean vertex degree variance and mean vertex degree.
        wsum = w.sum()
        E_rw = float((w * E).sum() / wsum) if wsum > 0 else float('nan')
        md_rw = float((w * md).sum() / wsum) if wsum > 0 else float('nan')
        cv_w = float(w.std() / w.mean()) if w.mean() > 0 else float('nan')

        print(f"  vertex deg variance:  raw mean {E.mean():.4f}, "
              f"reweighted {E_rw:.4f}")
        print(f"  mean vertex degree:   raw mean {md.mean():.4f}, "
              f"reweighted {md_rw:.4f}")
        print(f"  autocorr time tau:    {tau:.1f}  -> ESS_ac  = {ess_ac:.0f}")
        print(f"  weight CV:            {cv_w:.3f}  -> ESS_wt  = {ess_w:.0f} "
              f"/ {len(w)} (Kish)")
        print(f"  combined ESS:         {ess_combined:.0f} "
              f"independent samples from {len(E)} recorded")
    else:
        print("  (too few samples for ESS; increase --num-exchanges "
              "or lower --warmup-exchanges)")

    print(f"\nTotal time: {elapsed:.0f}s")
    print(f"Results saved to {csv_path}")

    # --- Final save ---
    tslot = slot_of_rung[target_rung]
    tmfd = samplers[tslot].manifold
    mfd_path = os.path.join(args.output_dir, f'{tag}_re_final.mfd')
    tmfd.dup().save(mfd_path)
    print(f"\nFinal target-rung manifold saved to {mfd_path}")
    print(f"  VDV coef: {ladder[target_rung]:.2f}")
    print(f"  Vertex degree variance: {tmfd.degree_variance(0):.4f}")
    print(f"  Mean vertex degree:     {tmfd.mean_degree(0):.4f}")
    print(f"  Mean edge degree:       {tmfd.mean_degree(1):.4f}")
    print(f"  Facets: {tmfd.num_facets}")

    ddg.gc_collect()
    ddg.gc_minimize()


if __name__ == '__main__':
    main()
