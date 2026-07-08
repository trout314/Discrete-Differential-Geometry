#!/usr/bin/env python3
"""Single-chain re-burn worker: load one seed, equilibrate under the current
build for a fixed number of sweeps, and write a new seed with the same basename.

Intended to be spawned as a subprocess by ``reburn_batch.py`` (memory-aware,
one chain per process for OOM isolation).  For interactive/whole-family use see
``reburn_family.py``.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from seed_utils import build_metadata_comments, get_free_memory_gb

EXIT_LOW_MEMORY = 42


def main():
    p = argparse.ArgumentParser(description="Re-burn a single seed for a fixed sweep count.")
    p.add_argument("--seed-file", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--n-burn", type=int, required=True)
    p.add_argument("--topology", default="S3")
    p.add_argument("--dimension", type=int, default=3)
    p.add_argument("--num-facets-target", type=int, required=True)
    p.add_argument("--num-facets-coef", type=float, required=True)
    p.add_argument("--hinge-degree-target", type=float, required=True)
    p.add_argument("--num-hinges-coef", type=float, required=True)
    p.add_argument("--hinge-degree-variance-coef", type=float, required=True)
    p.add_argument("--codim3-degree-variance-coef", type=float, required=True)
    p.add_argument("--min-free-memory-gb", type=float, default=4.0)
    p.add_argument("--skip-existing", action="store_true")
    args = p.parse_args()

    out_path = os.path.join(args.output_dir, os.path.basename(args.seed_file))
    if args.skip_existing and os.path.exists(out_path):
        return

    if get_free_memory_gb() < args.min_free_memory_gb:
        sys.exit(EXIT_LOW_MEMORY)

    params = SamplerParams(
        num_facets_target=args.num_facets_target,
        num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_degree_target,
        num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hinge_degree_variance_coef,
        codim3_degree_variance_coef=args.codim3_degree_variance_coef,
    )

    mfd = Manifold.load(args.seed_file, args.dimension)
    sampler = ManifoldSampler(mfd, params)
    sampler.reset_stats()

    if args.n_burn > 0:
        def mem_guard(done, total):
            return get_free_memory_gb() < args.min_free_memory_gb
        # callback returning True stops the run; check memory periodically.
        # Fixed move interval (NOT tied to num_facets_target): a large chain
        # otherwise checked only ~once per sweep, far too coarse to catch a
        # runaway before it exhausts RAM.  50k moves is a fraction of a sweep
        # for large N and still cheap (one /proc/meminfo read) for small N.
        sampler.run(sweeps=args.n_burn, callback=mem_guard,
                    callback_interval=50_000)

    if get_free_memory_gb() < args.min_free_memory_gb:
        sys.exit(EXIT_LOW_MEMORY)

    stats = sampler.get_stats()
    comments = build_metadata_comments(
        topology=args.topology,
        dimension=args.dimension,
        initial_triangulation=os.path.basename(args.seed_file),
        num_facets_target=params.num_facets_target,
        num_facets_coef=params.num_facets_coef,
        hinge_degree_target=params.hinge_degree_target,
        num_hinges_coef=params.num_hinges_coef,
        hinge_degree_variance_coef=params.hinge_degree_variance_coef,
        codim3_degree_variance_coef=params.codim3_degree_variance_coef,
        growth_step_size=0,
        eq_sweeps_per_step=0,
        equilibration_sweeps=args.n_burn,
        manifold_view=sampler.manifold,
        objective=sampler.current_objective,
        sampler_stats=stats,
    )
    os.makedirs(args.output_dir, exist_ok=True)
    # Save directly from the read-only view.  A prior `.dup()` here copied the
    # entire manifold (a ~30 GB spike for the largest chains) for no reason --
    # ManifoldView.save() writes the borrowed handle straight to disk.
    sampler.manifold.save(out_path, comments=comments)


if __name__ == "__main__":
    main()
