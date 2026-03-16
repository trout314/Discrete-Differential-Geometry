#!/usr/bin/env python3
"""Generate a single equilibrated seed triangulation."""

import argparse
import os
import sys
import time

# Allow importing the library from the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from seed_utils import build_metadata_comments, build_seed_filename


def main():
    parser = argparse.ArgumentParser(
        description="Equilibrate a seed triangulation under a given energy function."
    )
    parser.add_argument("--topology", default="S3", help="Topology label (default: S3)")
    parser.add_argument("--dimension", type=int, default=3, help="Manifold dimension (default: 3)")
    parser.add_argument("--num-facets-target", type=int, default=1000)
    parser.add_argument("--num-facets-coef", type=float, default=0.1)
    parser.add_argument("--hinge-degree-target", type=float, default=5.1)
    parser.add_argument("--num-hinges-coef", type=float, default=0.1)
    parser.add_argument("--hinge-degree-variance-coef", type=float, default=0.0)
    parser.add_argument("--codim3-degree-variance-coef", type=float, default=0.1)
    parser.add_argument("--equilibration-sweeps", type=int, default=500)
    parser.add_argument("--growth-step-size", type=int, default=500)
    parser.add_argument("--eq-sweeps-per-step", type=int, default=5)
    parser.add_argument("--output-dir", default="seeds")
    parser.add_argument(
        "--seed-file", default=None,
        help="Path to existing .mfd file to use as starting triangulation"
    )
    parser.add_argument(
        "--seed-index", type=int, default=None,
        help="Replica index. Appends _s{N:03d} to the output filename."
    )
    parser.add_argument(
        "--batch-mode", action="store_true",
        help="Print one-line status updates instead of interactive display."
    )
    args = parser.parse_args()

    # Build the output filename early so we can use it as a prefix in batch mode
    params = SamplerParams(
        num_facets_target=args.num_facets_target,
        hinge_degree_target=args.hinge_degree_target,
        num_facets_coef=args.num_facets_coef,
        num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hinge_degree_variance_coef,
        codim3_degree_variance_coef=args.codim3_degree_variance_coef,
    )
    filename = build_seed_filename(args.topology, params, seed_index=args.seed_index)
    tag = filename.removesuffix(".mfd")

    def log(msg):
        if args.batch_mode:
            print(f"[{tag}] {msg}", flush=True)
        else:
            print(msg)

    # 1. Load or create initial manifold
    if args.seed_file:
        log(f"Loading seed from {args.seed_file}")
        mfd = Manifold.load(args.seed_file, args.dimension)
        initial_triangulation = args.seed_file
    else:
        log(f"Creating standard sphere of dimension {args.dimension}")
        mfd = Manifold.standard_sphere(args.dimension)
        initial_triangulation = f"standard_sphere({args.dimension})"

    # 2. Create sampler
    sampler = ManifoldSampler(mfd, params)

    # 3. Ramped growth
    log(f"Growing to {args.num_facets_target} facets...")

    def growth_callback(cur, tgt):
        log(f"growing: {cur}/{tgt} facets")

    sampler.ramped_grow(
        target_facets=args.num_facets_target,
        step_size=args.growth_step_size,
        eq_sweeps_per_step=args.eq_sweeps_per_step,
        callback=growth_callback,
    )
    log(f"Growth complete: {sampler.manifold.num_facets} facets")

    # 4. Equilibration
    log(f"Equilibrating for {args.equilibration_sweeps} sweeps...")
    sampler.reset_stats()

    if args.batch_mode:
        t_start = time.monotonic()
        last_print = [0.0]
        stats_before = sampler.get_stats()

        def batch_callback(done, total):
            now = time.monotonic()
            if now - last_print[0] < 1.0:
                return False
            last_print[0] = now
            elapsed = now - t_start
            stats = sampler.get_stats()
            tried = stats.total_tried - stats_before.total_tried
            accepted = stats.total_accepted - stats_before.total_accepted
            accept_pct = 100.0 * accepted / tried if tried > 0 else 0.0
            pct = 100.0 * done / total if total > 0 else 0.0
            nf = sampler.manifold.num_facets
            obj = sampler.current_objective
            log(f"{pct:.1f}% facets={nf} obj={obj:.1f} accept={accept_pct:.1f}% {elapsed:.0f}s")
            return False

        sampler.run(sweeps=args.equilibration_sweeps, callback=batch_callback)
    else:
        sampler.run_with_display(sweeps=args.equilibration_sweeps)

    # 5. Build metadata and save
    comments = build_metadata_comments(
        topology=args.topology,
        dimension=args.dimension,
        initial_triangulation=initial_triangulation,
        num_facets_target=args.num_facets_target,
        num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_degree_target,
        num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hinge_degree_variance_coef,
        codim3_degree_variance_coef=args.codim3_degree_variance_coef,
        growth_step_size=args.growth_step_size,
        eq_sweeps_per_step=args.eq_sweeps_per_step,
        equilibration_sweeps=args.equilibration_sweeps,
        manifold_view=sampler.manifold,
        objective=sampler.current_objective,
    )

    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, filename)

    sampler.manifold.dup().save(output_path, comments=comments)
    log(f"Saved: {output_path}")

    # 6. Summary
    mfd_view = sampler.manifold
    dim = mfd_view.dimension
    log(
        f"facets={mfd_view.num_facets} "
        f"hinge_deg_avg={mfd_view.mean_degree(dim - 2):.3f} "
        f"hinge_deg_var={mfd_view.degree_variance(dim - 2):.2f} "
        f"vtx_deg_var={mfd_view.degree_variance(0):.1f} "
        f"obj={sampler.current_objective:.1f}"
    )


if __name__ == "__main__":
    main()
