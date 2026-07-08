#!/usr/bin/env python3
"""Re-burn one family of seed triangulations under the CURRENT objective build.

Every existing seed in ``seeds/`` was equilibrated under an *older* objective
(the pre-fix codim-3 VDV floor, which was ~raw VDV and pulled d-bar toward even
integers).  The current build uses the corrected, d-bar-agnostic floor, whose
equilibrium is slightly different.  This script re-equilibrates a whole family
of seeds under the current build so the library becomes a valid set of
independent equilibrium samples for the *new* objective.

Because the equilibrium distribution is start-independent (ergodicity), starting
from the old seeds rather than a fresh sphere only changes the *burn-in length*
required -- not the final distribution.  So the whole job reduces to: burn long
enough, gated by a convergence check rather than a fixed sweep count.

Two phases per family:

  A. CALIBRATE.  Advance a subset of the family's chains in lockstep, chunk by
     chunk, tracking the ensemble mean of the edge (hinge) degree -- the most
     sensitive observable, since the codim-3 term is what changed.  When the
     ensemble mean stops drifting (OLS slope over a trailing window is
     statistically indistinguishable from zero), declare convergence.  The
     sweep count at that point is the per-family burn length N_burn.

  B. PRODUCE.  Re-burn every member of the family for N_burn sweeps and write a
     new seed carrying updated metadata (initial_triangulation = the old seed).

The convergence gate is deliberately conservative: N_burn is the *full* number
of sweeps the calibration ran before it declared a plateau, so every production
chain gets at least as much burn as the calibration showed sufficient.
"""

import argparse
import glob
import math
import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from seed_utils import build_metadata_comments, load_seed_metadata


def family_stem(path: str) -> str:
    """Strip the ``_s{NNN}.mfd`` replica suffix to get the family stem."""
    base = os.path.basename(path)
    return re.sub(r"_s\d+\.mfd$", "", base)


def params_from_metadata(meta: dict) -> tuple[SamplerParams, str, int]:
    """Build SamplerParams (and topology, dimension) from an .mfd header."""
    params = SamplerParams(
        num_facets_target=int(meta["num_facets_target"]),
        num_facets_coef=float(meta["num_facets_coef"]),
        hinge_degree_target=float(meta["hinge_degree_target"]),
        num_hinges_coef=float(meta["num_hinges_coef"]),
        hinge_degree_variance_coef=float(meta["hinge_degree_variance_coef"]),
        codim3_degree_variance_coef=float(meta["codim3_degree_variance_coef"]),
    )
    topology = meta.get("topology", "S3")
    dimension = int(meta.get("dimension", 3))
    return params, topology, dimension


def observables(mfd_view):
    """Return (edge_degree, mean_vertex_degree, vertex_degree_variance)."""
    dim = mfd_view.dimension
    return (
        mfd_view.mean_degree(dim - 2),   # hinge / edge degree
        mfd_view.mean_degree(0),         # mean vertex degree (d-bar)
        mfd_view.degree_variance(0),     # VDV
    )


def ols_slope_tstat(xs, ys):
    """Least-squares slope of ys vs xs and |slope|/se(slope).

    Returns (slope, tstat).  tstat is inf when there is no residual scatter and
    a nonzero slope; 0 when the slope is exactly flat.
    """
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return 0.0, 0.0
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    # residual variance -> standard error of the slope
    resid = [y - (my + slope * (x - mx)) for x, y in zip(xs, ys)]
    if n <= 2:
        return slope, float("inf")
    s2 = sum(r * r for r in resid) / (n - 2)
    se = math.sqrt(s2 / sxx) if s2 > 0 else 0.0
    if se == 0:
        return slope, (0.0 if slope == 0 else float("inf"))
    return slope, abs(slope) / se


def calibrate(members, params, dimension, *, subset, chunk_sweeps,
              window, patience, t_thresh, min_sweeps, max_sweeps, log):
    """Phase A: find the burn length by watching ensemble-mean edge degree.

    Returns (n_burn, trajectory) where trajectory is a list of
    (sweeps, edge_mean, edge_sem, dbar_mean, vdv_mean).
    """
    subset = min(subset, len(members))
    log(f"  [calibrate] loading {subset} chains for lockstep drift test")
    samplers = []
    for path in members[:subset]:
        mfd = Manifold.load(path, dimension)
        samplers.append(ManifoldSampler(mfd, params))

    def ensemble():
        rows = [observables(s.manifold) for s in samplers]
        k = len(rows)
        edge = [r[0] for r in rows]
        edge_mean = sum(edge) / k
        edge_var = sum((e - edge_mean) ** 2 for e in edge) / k if k > 1 else 0.0
        edge_sem = math.sqrt(edge_var / k) if k > 1 else 0.0
        dbar_mean = sum(r[1] for r in rows) / k
        vdv_mean = sum(r[2] for r in rows) / k
        return edge_mean, edge_sem, dbar_mean, vdv_mean

    traj = []
    sweeps = 0
    em, es, db, vd = ensemble()
    traj.append((sweeps, em, es, db, vd))
    log(f"    sweeps={sweeps:>6}  edgeDeg={em:.5f}+/-{es:.5f}  dbar={db:.3f}  VDV={vd:.2f}  <- as saved")

    good = 0
    n_burn = None
    while sweeps < max_sweeps:
        for s in samplers:
            s.run(sweeps=chunk_sweeps)
        sweeps += chunk_sweeps
        em, es, db, vd = ensemble()
        traj.append((sweeps, em, es, db, vd))

        # OLS slope of edge degree over the trailing window
        win = traj[-window:]
        tstat = float("inf")
        slope = 0.0
        if len(win) >= 3:
            xs = [w[0] for w in win]
            ys = [w[1] for w in win]
            slope, tstat = ols_slope_tstat(xs, ys)

        flat = (sweeps >= min_sweeps) and (tstat < t_thresh)
        good = good + 1 if flat else 0
        flag = f" flat({good}/{patience})" if flat else ""
        log(f"    sweeps={sweeps:>6}  edgeDeg={em:.5f}+/-{es:.5f}  dbar={db:.3f}  "
            f"VDV={vd:.2f}  slope*W={slope * (window * chunk_sweeps):+.4f} t={tstat:.2f}{flag}")

        if good >= patience:
            n_burn = sweeps
            break

    if n_burn is None:
        n_burn = sweeps
        log(f"  [calibrate] WARNING: no plateau within max_sweeps={max_sweeps}; "
            f"using N_burn={n_burn}")
    else:
        log(f"  [calibrate] converged: N_burn={n_burn} sweeps")
    return n_burn, traj


def reburn_one(path, params, topology, dimension, n_burn, output_dir, log):
    """Phase B worker: re-burn a single seed for n_burn sweeps and save."""
    mfd = Manifold.load(path, dimension)
    sampler = ManifoldSampler(mfd, params)
    sampler.reset_stats()
    if n_burn > 0:
        sampler.run(sweeps=n_burn)

    stats = sampler.get_stats()
    comments = build_metadata_comments(
        topology=topology,
        dimension=dimension,
        initial_triangulation=os.path.basename(path),
        num_facets_target=params.num_facets_target,
        num_facets_coef=params.num_facets_coef,
        hinge_degree_target=params.hinge_degree_target,
        num_hinges_coef=params.num_hinges_coef,
        hinge_degree_variance_coef=params.hinge_degree_variance_coef,
        codim3_degree_variance_coef=params.codim3_degree_variance_coef,
        growth_step_size=0,           # no growth: started from an existing seed
        eq_sweeps_per_step=0,
        equilibration_sweeps=n_burn,
        manifold_view=sampler.manifold,
        objective=sampler.current_objective,
        sampler_stats=stats,
    )
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, os.path.basename(path))
    # Save straight from the read-only view; no need to dup the whole manifold.
    sampler.manifold.save(out_path, comments=comments)
    return out_path


def main():
    p = argparse.ArgumentParser(description="Re-burn one seed family under the current objective.")
    p.add_argument("--family", required=True,
                   help="Family stem (filename without _sNNN.mfd), e.g. "
                        "S3_N1e2_1e-1_ED5p1043_1e-1_VDV_5e-1")
    p.add_argument("--seed-dir", default="seeds")
    p.add_argument("--output-dir", default="seeds_reburned")
    p.add_argument("--subset", type=int, default=16,
                   help="Chains held in lockstep during calibration (default 16)")
    p.add_argument("--chunk-sweeps", type=int, default=100)
    p.add_argument("--window", type=int, default=5,
                   help="Trailing chunks used for the slope test (default 5)")
    p.add_argument("--patience", type=int, default=2,
                   help="Consecutive flat evaluations required (default 2)")
    p.add_argument("--t-thresh", type=float, default=1.5,
                   help="|slope|/se below this counts as flat (default 1.5)")
    p.add_argument("--min-sweeps", type=int, default=200)
    p.add_argument("--max-sweeps", type=int, default=5000)
    p.add_argument("--calibrate-only", action="store_true",
                   help="Run only phase A (report N_burn, write nothing).")
    p.add_argument("--n-burn", type=int, default=None,
                   help="Skip calibration; use this burn length for phase B.")
    p.add_argument("--limit", type=int, default=None,
                   help="Re-burn at most this many members (for trials).")
    args = p.parse_args()

    members = sorted(glob.glob(os.path.join(args.seed_dir, args.family + "_s*.mfd")))
    if not members:
        print(f"No members found for family {args.family!r} in {args.seed_dir}", file=sys.stderr)
        sys.exit(1)

    meta = load_seed_metadata(members[0])
    params, topology, dimension = params_from_metadata(meta)

    def log(msg):
        print(msg, flush=True)

    log(f"Family {args.family}: {len(members)} members")
    log(f"  params: N_target={params.num_facets_target} ED_target={params.hinge_degree_target} "
        f"facets_coef={params.num_facets_coef} hinges_coef={params.num_hinges_coef} "
        f"VDV_coef={params.codim3_degree_variance_coef} HDV_coef={params.hinge_degree_variance_coef}")

    t0 = time.monotonic()
    if args.n_burn is not None:
        n_burn = args.n_burn
        log(f"  [calibrate] skipped; using N_burn={n_burn}")
    else:
        n_burn, _ = calibrate(
            members, params, dimension,
            subset=args.subset, chunk_sweeps=args.chunk_sweeps,
            window=args.window, patience=args.patience, t_thresh=args.t_thresh,
            min_sweeps=args.min_sweeps, max_sweeps=args.max_sweeps, log=log,
        )
    log(f"  [calibrate] elapsed {time.monotonic() - t0:.0f}s")

    if args.calibrate_only:
        log(f"DONE (calibrate-only): N_burn={n_burn}")
        return

    to_do = members if args.limit is None else members[:args.limit]
    log(f"  [produce] re-burning {len(to_do)} members for {n_burn} sweeps each -> {args.output_dir}")
    tprod = time.monotonic()
    for i, path in enumerate(to_do, 1):
        reburn_one(path, params, topology, dimension, n_burn, args.output_dir, log)
        if i % 8 == 0 or i == len(to_do):
            rate = (time.monotonic() - tprod) / i
            log(f"    {i}/{len(to_do)} done ({rate:.1f}s/chain)")
    log(f"  [produce] elapsed {time.monotonic() - tprod:.0f}s")
    log(f"DONE: family {args.family}, N_burn={n_burn}, {len(to_do)} seeds written to {args.output_dir}")


if __name__ == "__main__":
    main()
