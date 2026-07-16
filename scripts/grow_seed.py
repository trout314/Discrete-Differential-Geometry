#!/usr/bin/env python3
"""Grow a triangulation to a large facet count, to seed large-N equilibrium runs.

No triangulations exist above N~3e5, so 1e6/1e7 equilibrium chains need an
initial config. This loads a smaller (ideally equilibrated) seed and grows it via
ramped_grow to --target-facets under the target edge-degree pin, then briefly
equilibrates and saves. The result is a high-VDV large-N config to hand to
`equilibrium_vdv.py --produce` as --seed-file (the "above" start); the production
run's warmup provides the "below" start, and its burn-in equilibrates VDV.

Memory-guarded (aborts below --min-free-gb). For very large targets this is a
long job; run it in the background.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from discrete_differential_geometry import (Manifold, ManifoldSampler,
                                             SamplerParams, vertex_degree_target)
from seed_utils import (build_metadata_comments, get_free_memory_gb,
                        make_leg, obj_of, read_history)

EXIT_LOW_MEMORY = 42
_UNMIGRATED = ("source predates history tracking; prior lineage not inlined "
               "(run the back-fill migration to complete it)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--from", dest="src", required=True, help="Smaller seed to grow from.")
    p.add_argument("--target-facets", type=int, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--dim", type=int, default=3)
    p.add_argument("--topology", default="S3")
    p.add_argument("--hinge-target", type=float, default=5.1043)
    p.add_argument("--num-facets-coef", type=float, default=0.1)
    p.add_argument("--num-hinges-coef", type=float, default=0.1)
    p.add_argument("--hdv-coef", type=float, default=0.0)
    p.add_argument("--vdv-coef", type=float, default=0.0,
                   help="VDV coefficient during growth (0 = fastest; production "
                        "drives VDV down afterward).")
    p.add_argument("--vdq-coef", type=float, default=0.0,
                   help="VDQ: fixed-target vertex penalty coupling c_v (raw, per "
                        "vertex; target derived from --hinge-target).")
    p.add_argument("--edq-coef", type=float, default=0.0,
                   help="EDQ: fixed-target edge penalty coupling c_e (raw, per edge).")
    p.add_argument("--step-size", type=int, default=5000)
    p.add_argument("--eq-sweeps", type=int, default=5, help="Equilibration sweeps per step.")
    p.add_argument("--final-eq-sweeps", type=int, default=200)
    p.add_argument("--min-free-gb", type=float, default=4.0)
    args = p.parse_args()

    if os.path.exists(args.out):
        print(f"{args.out} exists; nothing to do."); return

    mfd = Manifold.load(args.src, args.dim)
    params = SamplerParams(
        num_facets_target=args.target_facets, num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=args.vdv_coef,
        hinge_degree_target_coef=args.edq_coef,
        codim3_degree_target_coef=args.vdq_coef,
        codim3_degree_target=(vertex_degree_target(args.hinge_target)
                              if args.vdq_coef else 0.0))
    s = ManifoldSampler(mfd, params)
    start = s.manifold.num_facets
    print(f"growing {os.path.basename(args.src)} ({start} facets) -> "
          f"{args.target_facets}, hinge_target={args.hinge_target}", flush=True)

    report = max(args.step_size, args.target_facets // 20)

    def cb(cur, tgt):
        if get_free_memory_gb() < args.min_free_gb:
            print(f"LOW MEMORY at {cur} facets; aborting.", file=sys.stderr)
            sys.exit(EXIT_LOW_MEMORY)
        if cur % report < args.step_size:
            print(f"  {cur}/{args.target_facets}  edgeDeg={s.manifold.mean_degree(1):.3f}",
                  flush=True)

    s.ramped_grow(args.target_facets, step_size=args.step_size,
                  eq_sweeps_per_step=args.eq_sweeps, callback=cb)
    if args.final_eq_sweeps > 0:
        s.run(sweeps=args.final_eq_sweeps)

    v = s.manifold
    # Provenance legs: a grow leg (sphere/prev -> target N) then a settle leg.
    prior = read_history(args.src)
    st = s.get_stats()
    growth_steps = max(0, -(-(args.target_facets - start) // args.step_size))  # ceil
    legs = [make_leg("grow", obj_of(params), growth_steps * args.eq_sweeps)]
    if args.final_eq_sweeps > 0:
        legs.append(make_leg("equilibrate", obj_of(params), args.final_eq_sweeps))
    # Stats span the whole invocation (grow+settle, no reset between); attribute
    # the totals to the final leg rather than fabricate a per-leg split.
    legs[-1]["tried"], legs[-1]["accepted"] = st.total_tried, st.total_accepted
    comments = build_metadata_comments(
        topology=args.topology, dimension=args.dim,
        initial_triangulation=os.path.basename(args.src),
        num_facets_target=args.target_facets, num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=args.vdv_coef,
        hinge_degree_target_coef=args.edq_coef,
        codim3_degree_target_coef=args.vdq_coef,
        codim3_degree_target=params.codim3_degree_target,
        growth_step_size=args.step_size, eq_sweeps_per_step=args.eq_sweeps,
        equilibration_sweeps=args.final_eq_sweeps,
        manifold_view=v, objective=s.current_objective, sampler_stats=st,
        legs=legs, prior_history=prior,
        history_note=None if prior else _UNMIGRATED)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    v.save(args.out, comments=comments)
    print(f"done: {v.num_facets} facets, edgeDeg={v.mean_degree(1):.4f}, "
          f"VDV={v.degree_variance(0):.1f} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
