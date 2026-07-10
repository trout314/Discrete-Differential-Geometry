#!/usr/bin/env python3
"""Two-stage VDV annealing to produce minimal degree variance samples.

Stage 1 (coarse): Ramp VDV coefficient up quickly (1.5x geometric) to find
the approximate ceiling where acceptance drops below a threshold.

Stage 2 (fine): Back off to the last good VDV, then ramp slowly with many
sweeps per step, squeezing out remaining variance. Stops on a VDV *plateau*
(best degree variance improves less than --fine-vdv-tol over --fine-patience
steps), not on an acceptance threshold, so it keeps annealing while VDV is
still dropping even at sub-1% acceptance.

Usage:
    # From a saved seed
    python scripts/anneal_vdv.py \
        --seed-file seeds/S3_N1e3_1e-1_ED5p0043_1e-1_VDV_5e-1_s000.mfd

    # Grow from scratch
    python scripts/anneal_vdv.py --size 1000

    # Custom parameters
    python scripts/anneal_vdv.py \
        --seed-file seeds/S3_N1e4_1e-1_ED5p0043_1e-1_VDV_5e-1_s000.mfd \
        --coarse-factor 1.3 --fine-factor 1.05 --fine-sweeps 200
"""

import argparse
import csv
import glob
import json
import os
import statistics
import sys
import time
from collections import deque

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))   # discrete_differential_geometry
sys.path.insert(0, os.path.join(_ROOT, "tools"))    # seed_utils
import discrete_differential_geometry as ddg
from seed_utils import get_free_memory_gb


# --- Memory safety ---

MIN_FREE_MEMORY_GB = 4.0


def check_memory(context=""):
    """Abort if free memory is below the safety threshold."""
    free = get_free_memory_gb()
    if free < MIN_FREE_MEMORY_GB:
        print(f"ABORTING: free memory {free:.1f} GB < {MIN_FREE_MEMORY_GB:.0f} GB threshold"
              f"{' (' + context + ')' if context else ''}")
        sys.exit(42)


def even_floor(dbar):
    """Minimum vertex-degree variance achievable at mean degree ``dbar``.

    Codim-3 (vertex) degrees are always even, so the minimum-variance degree
    distribution splits across the two nearest even integers, giving
    4 y (1-y) with y = frac(dbar/2).
    """
    y = (dbar / 2.0) % 1.0
    return 4.0 * y * (1.0 - y)


def save_checkpoint(outdir, tag, sampler, state):
    """Atomically write a resumable checkpoint: manifold + fine-stage state.

    Written to <name>.tmp then os.replace'd so an interrupted write never leaves
    a half-file. The JSON is written last and references the manifold it pairs
    with, so resume trusts the JSON.
    """
    mfd_path = os.path.join(outdir, f'{tag}_checkpoint.mfd')
    json_path = os.path.join(outdir, f'{tag}_checkpoint.json')
    sampler.manifold.dup().save(mfd_path + '.tmp')
    os.replace(mfd_path + '.tmp', mfd_path)
    state = dict(state, manifold=os.path.basename(mfd_path))
    with open(json_path + '.tmp', 'w') as f:
        json.dump(state, f, indent=2)
    os.replace(json_path + '.tmp', json_path)


def load_checkpoint(resume_arg, dim):
    """Load a checkpoint written by save_checkpoint.

    Returns (state, manifold, outdir). resume_arg may be the output directory or
    the *_checkpoint.json path itself.
    """
    if os.path.isdir(resume_arg):
        found = sorted(glob.glob(os.path.join(resume_arg, '*_checkpoint.json')))
        if not found:
            raise SystemExit(f"No *_checkpoint.json found in {resume_arg}")
        json_path = found[0]
    else:
        json_path = resume_arg
    with open(json_path) as f:
        state = json.load(f)
    outdir = os.path.dirname(json_path) or '.'
    m = ddg.Manifold.load(os.path.join(outdir, state['manifold']), dim)
    return state, m, outdir


def probe_acceptance(sampler, sweeps):
    """Run sweeps and return acceptance rate."""
    sampler.reset_stats()
    sampler.run(sweeps=sweeps)
    stats = sampler.get_stats()
    if stats.total_tried == 0:
        return 0.0
    return stats.total_accepted / stats.total_tried


def log_state(writer, vdv, sampler, dim, stage, accept_rate):
    """Log one row to CSV and print to console."""
    mfd = sampler.manifold
    nf = mfd.num_facets
    dv0 = mfd.degree_variance(0)
    dv_hinge = mfd.degree_variance(max(0, dim - 2))
    md0 = mfd.mean_degree(0)
    md_hinge = mfd.mean_degree(max(0, dim - 2))
    obj = sampler.current_objective

    writer.writerow([
        stage, f'{vdv:.4f}', f'{accept_rate:.6f}',
        nf, f'{dv0:.4f}', f'{dv_hinge:.4f}',
        f'{md0:.4f}', f'{md_hinge:.4f}', f'{obj:.4f}',
    ])

    print(f"  {stage:>5s}  {vdv:12.2f} {100*accept_rate:8.2f}% {dv0:12.2f} "
          f"{dv_hinge:14.2f} {obj:10.1f} {nf:8d}")

    return dv0, dv_hinge


def main():
    parser = argparse.ArgumentParser(
        description='Two-stage VDV annealing for minimal degree variance')
    parser.add_argument('--seed-file', type=str, default=None,
                        help='Path to an existing .mfd seed file to load')
    parser.add_argument('--size', type=int, default=1000,
                        help='Target number of facets when growing from scratch')
    parser.add_argument('--dim', type=int, default=3,
                        help='Manifold dimension (default: 3)')
    parser.add_argument('--hinge-target', type=float, default=5.0043,
                        help='Hinge degree target (default: 5.0043)')
    parser.add_argument('--num-facets-coef', type=float, default=0.1,
                        help='Volume penalty coefficient (default: 0.1)')
    parser.add_argument('--num-hinges-coef', type=float, default=0.1,
                        help='Hinge count penalty coefficient (default: 0.1)')
    parser.add_argument('--vdv-start', type=float, default=0.5,
                        help='Starting VDV coefficient (default: 0.5)')

    # Stage 1: coarse
    parser.add_argument('--coarse-factor', type=float, default=1.5,
                        help='Coarse stage multiplicative factor (default: 1.5)')
    parser.add_argument('--coarse-sweeps', type=int, default=5,
                        help='Sweeps per coarse step (default: 5)')
    parser.add_argument('--coarse-probe', type=int, default=5,
                        help='Probe sweeps for coarse acceptance / logging (default: 5)')
    parser.add_argument('--coarse-threshold', type=float, default=0.05,
                        help='Acceptance floor that also ends the coarse stage '
                             '(default: 0.05 = 5%%). With the floor-subtracted (excess) '
                             'objective, acceptance stays high at the floor, so the '
                             'VDV-plateau test below is what normally stops it.')
    parser.add_argument('--coarse-vdv-tol', type=float, default=0.02,
                        help='End the coarse stage when the best VDV improves by less '
                             'than this fraction over --coarse-patience steps (default: '
                             '0.02). This is the primary coarse stop: the excess '
                             'objective flattens acceptance near the floor, so ramping '
                             'on acceptance alone runs the coefficient away.')
    parser.add_argument('--coarse-patience', type=int, default=4,
                        help='Window (coarse steps) for the VDV-plateau test (default: 4)')
    parser.add_argument('--coarse-max-steps', type=int, default=80,
                        help='Hard cap on coarse steps (default: 80)')

    # Stage 2: fine
    parser.add_argument('--fine-factor', type=float, default=1.05,
                        help='Fine stage multiplicative factor (default: 1.05)')
    parser.add_argument('--fine-sweeps', type=int, default=100,
                        help='Sweeps per fine step (default: 100)')
    parser.add_argument('--fine-probe', type=int, default=20,
                        help='Probe sweeps for fine acceptance / logging (default: 20)')
    parser.add_argument('--fine-vdv-tol', type=float, default=0.01,
                        help='Stop when the best VDV improves by less than this fraction '
                             'over --fine-patience steps (default: 0.01 = 1%%)')
    parser.add_argument('--fine-patience', type=int, default=20,
                        help='Window (steps) over which --fine-vdv-tol improvement is '
                             'required to keep annealing (default: 20)')
    parser.add_argument('--fine-min-accept', type=float, default=0.0,
                        help='Optional hard acceptance floor to stop the fine stage '
                             '(default: 0 = disabled; stopping is governed by VDV plateau)')
    parser.add_argument('--fine-max-steps', type=int, default=500,
                        help='Hard cap on fine steps (default: 500)')

    # Final equilibrium measurement
    parser.add_argument('--measure-samples', type=int, default=30,
                        help='VDV samples to average for the final equilibrium '
                             'measurement at the best coefficient (default: 30)')
    parser.add_argument('--measure-sweeps', type=int, default=20,
                        help='Sweeps between equilibrium-measurement samples (default: 20)')

    # Checkpoint / resume
    parser.add_argument('--checkpoint-interval', type=int, default=25,
                        help='Save a resumable checkpoint every N fine steps '
                             '(default: 25; 0 disables). Robustness for long runs.')
    parser.add_argument('--resume', type=str, default=None,
                        help='Resume the fine stage from a checkpoint: pass the '
                             'output dir (or the *_checkpoint.json). Skips the coarse '
                             'stage and continues from the saved coefficient. Schedule '
                             'flags (--fine-sweeps, --fine-factor, --fine-patience, '
                             '--fine-max-steps) may be changed to continue with heavier '
                             'equilibration; the physical params (pin coefs, hinge '
                             'target) are taken from the checkpoint.')

    # Output
    parser.add_argument('--output-dir', default='data/anneal_vdv',
                        help='Directory for output files')
    parser.add_argument('--min-free-memory-gb', type=float, default=4.0,
                        help='Abort if free memory drops below this (default: 4 GB)')
    args = parser.parse_args()

    global MIN_FREE_MEMORY_GB
    MIN_FREE_MEMORY_GB = args.min_free_memory_gb

    dim = args.dim
    check_memory("pre-flight")

    # --- Load manifold: resume a checkpoint, load a seed, or grow ---
    resume_state = None
    if args.resume:
        resume_state, m, args.output_dir = load_checkpoint(args.resume, dim)
        size = resume_state['size']
        params_dict = dict(resume_state['params'])
        params_dict['num_facets_target'] = size
        print(f"Resuming from checkpoint in {args.output_dir}: {size} facets, "
              f"fine step {resume_state['fine_steps']}, "
              f"VDV coef {resume_state['vdv']:.1f}, best VDV {resume_state['best_dv0']:.4f}")
    elif args.seed_file:
        print(f"Loading seed from {args.seed_file}...")
        m = ddg.Manifold.load(args.seed_file, dim)
        size = m.num_facets
        print(f"Loaded {size} facets")
        params_dict = None
    else:
        size = args.size
        print(f"Growing {dim}d manifold to {size} facets...")
        m = ddg.Manifold.standard_sphere(dim)
        params_dict = None

    if params_dict is None:
        params_dict = {
            'num_facets_target': size,
            'hinge_degree_target': args.hinge_target,
            'num_facets_coef': args.num_facets_coef,
            'num_hinges_coef': args.num_hinges_coef,
        }

    params = ddg.SamplerParams(
        num_facets_target=params_dict['num_facets_target'],
        hinge_degree_target=params_dict['hinge_degree_target'],
        num_facets_coef=params_dict['num_facets_coef'],
        num_hinges_coef=params_dict['num_hinges_coef'],
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=(resume_state['vdv'] if resume_state else args.vdv_start),
    )
    sampler = ddg.ManifoldSampler(m, params)

    if resume_state is None and not args.seed_file:
        step_size = max(100, size // 20)
        def growth_callback(cur, tgt):
            check_memory(f"growth {cur}/{tgt}")
        sampler.ramped_grow(size, step_size=step_size, eq_sweeps_per_step=3,
                            callback=growth_callback)
        print(f"Grown to {sampler.manifold.num_facets} facets")

    # --- Prepare output ---
    os.makedirs(args.output_dir, exist_ok=True)
    tag = f'{dim}d_N{size}'
    csv_path = os.path.join(args.output_dir, f'{tag}_anneal.csv')
    csv_file = open(csv_path, 'a' if resume_state else 'w', newline='')
    writer = csv.writer(csv_file)
    if resume_state is None:
        writer.writerow([
            'stage', 'vdv_coef', 'acceptance_rate',
            'num_facets', 'degree_variance_0', 'degree_variance_hinge',
            'mean_degree_0', 'mean_degree_hinge', 'objective',
        ])

    print(f"\n{'Stage':>7} {'VDV coef':>12} {'Accept%':>9} {'DegVar(vtx)':>12} "
          f"{'DegVar(hinge)':>14} {'Objective':>10} {'Facets':>8}")
    print("-" * 82)

    t_start = time.monotonic()

    if resume_state is not None:
        # Resumed: skip the coarse stage and restore the fine-stage state.
        coarse_steps = resume_state.get('coarse_steps', 0)
        vdv = resume_state['vdv']
        best_dv0 = resume_state['best_dv0']
        best_vdv = resume_state['best_vdv']
        fine_steps = resume_state['fine_steps']
        best_history = deque(resume_state['best_history'], maxlen=args.fine_patience + 1)
        print(f"\n  Fine stage (resumed): VDV={vdv:.1f}, "
              f"factor={args.fine_factor}, sweeps/step={args.fine_sweeps}, "
              f"from fine step {fine_steps}")
    else:
        # =================================================================
        # Stage 1: Coarse ramp — fast-forward to where VDV reaches its floor
        # =================================================================
        # Stop on a VDV plateau, not an acceptance threshold. With the
        # floor-subtracted (excess) objective, once the triangulation sits on
        # its degree-variance floor the excess is ~0 and floor-preserving moves
        # stay accepted, so acceptance never falls and an acceptance-thresholded
        # ramp runs the coefficient away to absurd values. Ramping until VDV
        # stops dropping lands us near the floor at a sane coefficient.
        vdv = args.vdv_start
        last_good_vdv = vdv
        coarse_steps = 0
        coarse_best = float('inf')
        coarse_hist = deque(maxlen=args.coarse_patience + 1)
        coarse_reason = f"reached max {args.coarse_max_steps} coarse steps"

        while coarse_steps < args.coarse_max_steps:
            check_memory(f"coarse VDV={vdv:.1f}")
            sampler.set_codim3_degree_variance_coef(vdv)
            sampler.run(sweeps=args.coarse_sweeps)
            accept = probe_acceptance(sampler, args.coarse_probe)
            dv0, _ = log_state(writer, vdv, sampler, dim, 'coarse', accept)
            csv_file.flush()
            coarse_steps += 1

            if dv0 < coarse_best:
                coarse_best = dv0
                last_good_vdv = vdv
            coarse_hist.append(coarse_best)

            # Primary stop: VDV plateaued over the patience window.
            if len(coarse_hist) > args.coarse_patience:
                past = coarse_hist[0]
                if past > 0 and (past - coarse_best) / past < args.coarse_vdv_tol:
                    coarse_reason = (f"VDV plateaued at {coarse_best:.4f} "
                                     f"(<{100*args.coarse_vdv_tol:.1f}% gain over "
                                     f"{args.coarse_patience} steps)")
                    break

            # Secondary stop: acceptance genuinely collapsed (rarely fires here).
            if accept < args.coarse_threshold:
                coarse_reason = f"acceptance {100*accept:.1f}% below floor"
                break

            vdv *= args.coarse_factor

        print(f"\n  Coarse stage done: {coarse_reason}")
        print(f"  Backing off to last improving VDV: {last_good_vdv:.1f}")

        # =================================================================
        # Stage 2: Fine ramp — squeeze variance near the ceiling
        # =================================================================
        # Back off to last good VDV and re-equilibrate
        vdv = last_good_vdv
        sampler.set_codim3_degree_variance_coef(vdv)
        print(f"\n  Fine stage: starting at VDV={vdv:.1f}, "
              f"factor={args.fine_factor}, sweeps/step={args.fine_sweeps}")
        sampler.run(sweeps=args.fine_sweeps)

        best_dv0 = float('inf')
        best_vdv = vdv
        fine_steps = 0
        best_history = deque(maxlen=args.fine_patience + 1)

    # We stop when the best VDV has improved by less than --fine-vdv-tol over the
    # last --fine-patience steps (a VDV *plateau*), rather than when acceptance
    # crosses a threshold. This keeps annealing while VDV is still dropping even
    # at sub-1% acceptance, where the old acceptance-thresholded stop quit early.
    stop_reason = (f"reached max {args.fine_max_steps} steps "
                   f"(VDV may still be improving — raise --fine-max-steps)")

    while fine_steps < args.fine_max_steps:
        check_memory(f"fine VDV={vdv:.1f}")
        sampler.set_codim3_degree_variance_coef(vdv)
        sampler.run(sweeps=args.fine_sweeps)
        accept = probe_acceptance(sampler, args.fine_probe)
        dv0, _ = log_state(writer, vdv, sampler, dim, 'fine', accept)
        csv_file.flush()
        fine_steps += 1

        if dv0 < best_dv0:
            best_dv0 = dv0
            best_vdv = vdv
        best_history.append(best_dv0)

        # Plateau test: best VDV gained < tol (relative) over the patience window.
        if len(best_history) > args.fine_patience:
            past = best_history[0]
            if past > 0 and (past - best_dv0) / past < args.fine_vdv_tol:
                stop_reason = (f"VDV plateaued at {best_dv0:.4f} "
                               f"(<{100*args.fine_vdv_tol:.1f}% gain over "
                               f"{args.fine_patience} steps; acceptance {100*accept:.2f}%)")
                break

        # Optional hard acceptance floor (disabled by default).
        if args.fine_min_accept > 0 and accept < args.fine_min_accept:
            stop_reason = (f"acceptance {100*accept:.3f}% below floor "
                           f"at VDV={vdv:.1f}")
            break

        vdv *= args.fine_factor

        # Resumable checkpoint (before advancing to the next coefficient).
        if args.checkpoint_interval > 0 and fine_steps % args.checkpoint_interval == 0:
            save_checkpoint(args.output_dir, tag, sampler, {
                'stage': 'fine', 'vdv': vdv, 'fine_steps': fine_steps,
                'best_vdv': best_vdv, 'best_dv0': best_dv0,
                'best_history': list(best_history), 'coarse_steps': coarse_steps,
                'dim': dim, 'size': size, 'params': params_dict,
            })

    print(f"\n  Fine stage done: {stop_reason}")

    elapsed = time.monotonic() - t_start

    csv_file.close()
    print(f"\nResults saved to {csv_path}")
    print(f"Total time: {elapsed:.0f}s "
          f"({coarse_steps} coarse + {fine_steps} fine steps)")

    # --- Final equilibrium measurement at the best coefficient ---
    # Record VDV, mean vertex degree, mean edge degree and the even-granularity
    # excess *per sample*, so the excess is the internally consistent
    # VDV - floor(dbar) of the *same* configuration. (Averaging VDV over the
    # measure phase but taking dbar/floor from the final config mismatches the
    # two and biases the excess -- badly at small N, where dbar fluctuates
    # during measurement; this produced spurious negative excesses earlier.)
    # The full per-sample series is written to <tag>_equilibrium.csv for
    # downstream mean-offset / density-of-states analysis.
    sampler.set_codim3_degree_variance_coef(best_vdv)
    sampler.run(sweeps=args.fine_sweeps)
    vdv_s, dbar_s, edeg_s, excess_s = [], [], [], []
    samples_path = os.path.join(args.output_dir, f'{tag}_equilibrium.csv')
    with open(samples_path, 'w', newline='') as sf:
        sw = csv.writer(sf)
        sw.writerow(['sample', 'vdv', 'mean_vertex_degree', 'mean_edge_degree',
                     'even_floor', 'excess', 'hinge_degree_variance'])
        for i in range(args.measure_samples):
            sampler.run(sweeps=args.measure_sweeps)
            m = sampler.manifold
            v = m.degree_variance(0)
            d = m.mean_degree(0)
            e = m.mean_degree(1)
            fl = even_floor(d)
            hv = m.degree_variance(max(0, dim - 2))
            vdv_s.append(v); dbar_s.append(d); edeg_s.append(e); excess_s.append(v - fl)
            sw.writerow([i, f'{v:.6f}', f'{d:.6f}', f'{e:.6f}',
                         f'{fl:.6f}', f'{v - fl:.6f}', f'{hv:.6f}'])

    def _ms(x):
        return (statistics.fmean(x),
                statistics.pstdev(x) if len(x) > 1 else 0.0)
    eq_vdv, eq_std = _ms(vdv_s)
    eq_dbar, eq_dbar_std = _ms(dbar_s)
    eq_edeg, eq_edeg_std = _ms(edeg_s)
    eq_exc, eq_exc_std = _ms(excess_s)

    mfd = sampler.manifold
    mfd_path = os.path.join(args.output_dir, f'{tag}_annealed.mfd')
    mfd.dup().save(mfd_path)
    print(f"\nSaved: {mfd_path}")
    print(f"  Best VDV coef:           {best_vdv:.2f}")
    print(f"  Equilibrium VDV:         {eq_vdv:.4f} +/- {eq_std:.4f} "
          f"(mean of {args.measure_samples} samples)")
    print(f"  Mean vertex degree:      {eq_dbar:.4f} +/- {eq_dbar_std:.4f}")
    print(f"  Mean edge degree:        {eq_edeg:.5f} +/- {eq_edeg_std:.5f}")
    print(f"  Excess (per-sample):     {eq_exc:.4f} +/- {eq_exc_std:.4f}   "
          f"[floor(<dbar>) = {even_floor(eq_dbar):.4f}]")
    print(f"  Hinge degree variance:   {mfd.degree_variance(max(0, dim-2)):.4f}")
    print(f"  Facets:                  {mfd.num_facets}")
    print(f"  Equilibrium samples:     {samples_path}")

    ddg.gc_collect()
    ddg.gc_minimize()


if __name__ == '__main__':
    main()
