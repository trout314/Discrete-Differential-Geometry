#!/usr/bin/env python3
"""Fixed-coefficient VDV equilibrium sampling with two-sided convergence checks.

Goal: find, for a given N, the largest VDV coefficient beta at which we can still
*equilibrate* (the "edge of equilibratability"), and measure the equilibrium
vertex-degree variance <VDV>_beta there -- so we can compare it against what
(non-equilibrium) annealing achieves.

Method. At each beta we run chains from BOTH sides of the equilibrium VDV:
  * "above": start from a plain high-VDV seed and hold at beta -> VDV falls.
  * "below": start from the same seed, over-squeeze at a higher warmup
             coefficient, then hold at beta -> VDV rises.
If the two sides converge to a common value (split-R-hat -> 1, small gap) and the
integrated autocorrelation time tau is small enough for a usable ESS, beta is
equilibratable and <VDV>_beta is bracketed. Every other objective parameter
(facet-count pin, edge-degree pin below the flat value, HDV) is held fixed, so
only VDV is being equilibrated.

Two modes:
  --chain   run ONE chain, write its per-sample time series to --out (a worker).
  (default) driver: spawn chains over a beta-grid x {above,below} x replicas,
            then report R-hat / tau / ESS / <VDV> and the two-sided gap per beta.

Example
-------
    python scripts/equilibrium_vdv.py --n-target 1000 \
        --seed-file seeds/S3_N1e3_1e-1_ED5p0043_1e-1_s000.mfd \
        --coefs 1 5 20 50 100 200 400 --replicas 2 \
        --thin 5 --n-samples 500 --burnin 100 --output-dir data/equil_vdv/N1e3
"""

import argparse
import csv
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import (
    Manifold, ManifoldSampler, SamplerParams,
    integrated_autocorrelation_time, split_rhat,
)
from seed_utils import build_metadata_comments, build_seed_filename

FIELDS = ["side", "replica", "coef", "sample", "sweeps",
          "num_facets", "vdv", "edge_deg", "hdv", "accept"]


# ----------------------------------------------------------------------------
# Worker: one fixed-beta chain
# ----------------------------------------------------------------------------

def run_chain(args):
    mfd = Manifold.load(args.init, args.dim)
    params = SamplerParams(
        num_facets_target=args.n_target,
        num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target,
        num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=(args.warmup_coef if args.warmup_sweeps > 0
                                     else args.coef),
    )
    s = ManifoldSampler(mfd, params)

    # "below" start: over-squeeze VDV at a higher coefficient, then relax at beta.
    if args.warmup_sweeps > 0:
        s.run(sweeps=args.warmup_sweeps)
        s.set_codim3_degree_variance_coef(args.coef)

    # Burn-in at beta before recording (production: relax from the init to
    # equilibrium and decorrelate from the start).
    if args.burnin_sweeps > 0:
        s.run(sweeps=args.burnin_sweeps)

    v = s.manifold
    s.reset_stats()
    rows, total = [], 0
    for i in range(args.n_samples):
        s.reset_stats()
        s.run(sweeps=args.thin)
        total += args.thin
        st = s.get_stats()
        acc = st.total_accepted / st.total_tried if st.total_tried else 0.0
        rows.append([args.side, args.replica, args.coef, i, total,
                     v.num_facets, v.degree_variance(0), v.mean_degree(1),
                     v.degree_variance(1), acc])

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f); w.writerow(FIELDS); w.writerows(rows)

    if args.save_config:
        comments = build_metadata_comments(
            topology=args.topology, dimension=args.dim,
            initial_triangulation=os.path.basename(args.init),
            num_facets_target=args.n_target, num_facets_coef=args.num_facets_coef,
            hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
            hinge_degree_variance_coef=args.hdv_coef,
            codim3_degree_variance_coef=args.coef,
            growth_step_size=0, eq_sweeps_per_step=0,
            equilibration_sweeps=args.warmup_sweeps + args.burnin_sweeps
                                 + args.n_samples * args.thin,
            manifold_view=v, objective=s.current_objective,
            sampler_stats=s.get_stats())
        os.makedirs(os.path.dirname(args.save_config) or ".", exist_ok=True)
        v.save(args.save_config, comments=comments)


# ----------------------------------------------------------------------------
# Driver: grid + diagnostics
# ----------------------------------------------------------------------------

def spawn(args, coef, side, rep, out):
    # "below" start: either a supplied low-VDV config (e.g. an annealed one) held
    # at beta so VDV rises to equilibrium, or the plain seed over-squeezed by a
    # warmup at a higher coefficient. "above" always starts from the plain seed.
    if side == "below" and args.below_init:
        init, warm_sweeps = args.below_init, 0
    elif side == "below":
        init, warm_sweeps = args.seed_file, args.warmup_sweeps
    else:
        init, warm_sweeps = args.seed_file, 0
    cmd = [sys.executable, os.path.abspath(__file__), "--chain",
           "--init", init, "--dim", str(args.dim),
           "--n-target", str(args.n_target),
           "--num-facets-coef", str(args.num_facets_coef),
           "--hinge-target", str(args.hinge_target),
           "--num-hinges-coef", str(args.num_hinges_coef),
           "--hdv-coef", str(args.hdv_coef),
           "--coef", str(coef),
           "--warmup-coef", str(coef * args.warmup_factor),
           "--warmup-sweeps", str(warm_sweeps),
           "--thin", str(args.thin), "--n-samples", str(args.n_samples),
           "--side", side, "--replica", str(rep), "--out", out]
    return subprocess.run(cmd).returncode


def load_series(path, burnin):
    """Return the post-burnin (vdv, edge_deg, hdv, accept) arrays from a chain CSV."""
    v, e, h, a = [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            if int(row["sample"]) >= burnin:
                v.append(float(row["vdv"])); e.append(float(row["edge_deg"]))
                h.append(float(row["hdv"])); a.append(float(row["accept"]))
    return np.array(v), np.array(e), np.array(h), np.array(a)


def analyze(coef, out_dir, burnin, thin):
    chains = {"above": [], "below": []}
    edge_all, acc_all = [], []
    for path in sorted(glob.glob(os.path.join(out_dir, f"coef{coef:g}_*.csv"))):
        side = "above" if "_above_" in path else "below"
        v, e, h, a = load_series(path, burnin)
        if len(v) >= 4:
            chains[side].append(v); edge_all.append(e.mean()); acc_all.append(a.mean())
    allv = chains["above"] + chains["below"]
    if not allv:
        return None
    rhat = split_rhat(allv) if len(allv) >= 2 else float("nan")
    taus = [integrated_autocorrelation_time(v) for v in allv]
    tau = float(np.nanmedian(taus))
    n_post = sum(len(v) for v in allv)
    ess = n_post / tau if tau and not np.isnan(tau) else float("nan")
    mean_above = np.mean([v.mean() for v in chains["above"]]) if chains["above"] else np.nan
    mean_below = np.mean([v.mean() for v in chains["below"]]) if chains["below"] else np.nan
    vdv_all = np.concatenate(allv)
    return dict(coef=coef, rhat=rhat, tau_samples=tau, tau_sweeps=tau * thin,
                ess=ess, vdv=vdv_all.mean(), vdv_std=vdv_all.std(),
                above=mean_above, below=mean_below,
                gap=abs(mean_above - mean_below), edge=np.mean(edge_all),
                accept=np.mean(acc_all), n_chains=len(allv))


def produce(args):
    """Run K equilibrium chains at args.beta; if split-R-hat passes, copy the
    configs into seeds/ with library names. Half the chains start from the plain
    seed (above) and half from --below-init (below) so R-hat tests convergence
    from both directions."""
    import shutil
    if not (args.seed_file and args.beta):
        sys.exit("--produce needs --seed-file and --beta")
    params = SamplerParams(
        num_facets_target=args.n_target, num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=args.beta)
    stem = build_seed_filename(args.topology, params, seed_index=0).rsplit("_s", 1)[0]
    if glob.glob(os.path.join(args.seeds_dir, stem + "_s*.mfd")):
        print(f"seeds already exist for these params (stem {stem}); not duplicating.",
              file=sys.stderr)
        return

    K = args.replicas
    stage = os.path.join(args.output_dir, "staging")
    os.makedirs(stage, exist_ok=True)
    print(f"Producing {K} equilibrium chains: N={args.n_target}, beta={args.beta:g} "
          f"(beta/N={args.beta/args.n_target:.4g}), burnin={args.production_burnin} + "
          f"{args.n_samples}x{args.thin} measured", flush=True)

    def pspawn(i):
        below = (i % 2 == 1) and bool(args.below_init)
        init = args.below_init if below else args.seed_file
        cfg = os.path.join(stage, f"chain_{i}.mfd")
        out = os.path.join(stage, f"chain_{i}.csv")
        cmd = [sys.executable, os.path.abspath(__file__), "--chain",
               "--init", init, "--dim", str(args.dim), "--n-target", str(args.n_target),
               "--num-facets-coef", str(args.num_facets_coef),
               "--hinge-target", str(args.hinge_target),
               "--num-hinges-coef", str(args.num_hinges_coef),
               "--hdv-coef", str(args.hdv_coef), "--coef", str(args.beta),
               "--burnin-sweeps", str(args.production_burnin),
               "--thin", str(args.thin), "--n-samples", str(args.n_samples),
               "--side", "below" if below else "above", "--replica", str(i),
               "--topology", args.topology, "--out", out, "--save-config", cfg]
        return i, subprocess.run(cmd).returncode, cfg, out

    res = {}
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        for fut in as_completed([pool.submit(pspawn, i) for i in range(K)]):
            i, rc, cfg, out = fut.result()
            res[i] = (rc, cfg, out)
            if rc != 0:
                print(f"  chain {i} FAILED rc={rc}", file=sys.stderr)

    series, edges = [], []
    for i in range(K):
        rc, cfg, out = res.get(i, (1, None, None))
        if rc == 0 and out and os.path.exists(out):
            v, e, _, _ = load_series(out, burnin=0)  # in-chain burn-in already done
            if len(v) >= 4:
                series.append(v); edges.append(e.mean())
    rhat = split_rhat(series) if len(series) >= 2 else float("nan")
    allv = np.concatenate(series) if series else np.array([0.0])
    print(f"\nEnsemble: {len(series)} chains  <VDV>={allv.mean():.3f}+/-{allv.std():.3f}  "
          f"edgeDeg={np.mean(edges):.4f}  split-R-hat={rhat:.3f}", flush=True)

    if not (rhat < args.rhat_max):
        print(f"R-hat {rhat:.3f} >= {args.rhat_max}: NOT confidently equilibrated; "
              f"configs left in {stage}, not copied.", file=sys.stderr)
        return
    os.makedirs(args.seeds_dir, exist_ok=True)
    copied = 0
    for i in range(K):
        rc, cfg, _ = res.get(i, (1, None, None))
        if rc == 0 and cfg and os.path.exists(cfg):
            shutil.copy2(cfg, os.path.join(args.seeds_dir,
                         build_seed_filename(args.topology, params, seed_index=copied)))
            copied += 1
    print(f"R-hat OK -> copied {copied} equilibrium seeds to {args.seeds_dir}/ "
          f"({stem}_s000..s{copied - 1:03d})", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chain", action="store_true", help="Worker mode: run one chain.")
    p.add_argument("--init", default=None, help="Initial .mfd config (worker).")
    p.add_argument("--dim", type=int, default=3)
    # objective template
    p.add_argument("--n-target", type=int, required=True)
    p.add_argument("--num-facets-coef", type=float, default=0.1)
    p.add_argument("--hinge-target", type=float, default=5.0043)
    p.add_argument("--num-hinges-coef", type=float, default=0.1)
    p.add_argument("--hdv-coef", type=float, default=0.0)
    # chain controls
    p.add_argument("--coef", type=float, help="Study VDV coefficient beta (worker).")
    p.add_argument("--thin", type=int, default=5, help="Sweeps per recorded sample.")
    p.add_argument("--n-samples", type=int, default=500)
    p.add_argument("--warmup-coef", type=float, default=0.0)
    p.add_argument("--warmup-sweeps", type=int, default=0)
    p.add_argument("--side", default="above")
    p.add_argument("--replica", type=int, default=0)
    p.add_argument("--out", default=None)
    p.add_argument("--burnin-sweeps", type=int, default=0,
                   help="Sweeps at beta before recording (worker).")
    p.add_argument("--save-config", default=None,
                   help="Save the final manifold here (worker).")
    p.add_argument("--topology", default="S3")
    # driver controls
    p.add_argument("--seed-file", help="Plain high-VDV init (driver).")
    p.add_argument("--below-init", default=None,
                   help="Low-VDV config to start 'below' chains from (e.g. an "
                        "annealed one), held at beta so VDV rises to equilibrium. "
                        "If unset, 'below' uses the seed over-squeezed by a warmup.")
    p.add_argument("--coefs", type=float, nargs="+", help="Beta grid (driver).")
    p.add_argument("--replicas", type=int, default=2)
    p.add_argument("--burnin", type=int, default=100, help="Samples discarded as burn-in.")
    p.add_argument("--warmup-factor", type=float, default=8.0,
                   help="'below' warmup coefficient = factor * beta.")
    p.add_argument("--warmup-sweeps-driver", type=int, default=400,
                   help="Sweeps of over-squeeze warmup for 'below' chains.")
    p.add_argument("--max-workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    p.add_argument("--output-dir", default="data/equil_vdv")
    # production mode
    p.add_argument("--produce", action="store_true",
                   help="Produce a validated equilibrium ensemble and, if R-hat<"
                        "--rhat-max, copy the configs into --seeds-dir.")
    p.add_argument("--beta", type=float, help="Study coefficient for --produce.")
    p.add_argument("--production-burnin", type=int, default=1500,
                   help="Burn-in sweeps per production chain.")
    p.add_argument("--rhat-max", type=float, default=1.05)
    p.add_argument("--seeds-dir", default="seeds")
    args = p.parse_args()

    if args.chain:
        run_chain(args)
        return
    if args.produce:
        produce(args)
        return

    # ---- driver ----
    if not (args.seed_file and args.coefs):
        p.error("driver mode needs --seed-file and --coefs")
    args.warmup_sweeps = args.warmup_sweeps_driver
    os.makedirs(args.output_dir, exist_ok=True)

    jobs = []
    for coef in args.coefs:
        for side in ("above", "below"):
            for rep in range(args.replicas):
                out = os.path.join(args.output_dir, f"coef{coef:g}_{side}_{rep}.csv")
                jobs.append((coef, side, rep, out))
    print(f"N={args.n_target}: {len(jobs)} chains "
          f"({len(args.coefs)} coefs x 2 sides x {args.replicas} reps), "
          f"{args.n_samples} samples x {args.thin} sweeps each", flush=True)

    done = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futs = {pool.submit(spawn, args, c, s, r, o): (c, s, r) for c, s, r, o in jobs}
        for fut in as_completed(futs):
            rc = fut.result(); done += 1
            if rc != 0:
                print(f"  FAILED {futs[fut]} rc={rc}", file=sys.stderr)
            if done % max(1, len(jobs) // 10) == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} chains done", flush=True)

    print(f"\n{'coef':>8} {'accept':>7} {'R-hat':>7} {'tau_swp':>8} {'ESS':>6} "
          f"{'<VDV>':>8} {'above':>8} {'below':>8} {'gap':>7} {'edgeDeg':>7}")
    rows = []
    for coef in args.coefs:
        r = analyze(coef, args.output_dir, args.burnin, args.thin)
        if not r:
            continue
        rows.append(r)
        eq = "OK" if (r["rhat"] < 1.05 and r["gap"] < 0.1 * r["vdv"] + 0.5) else "  "
        print(f"{r['coef']:>8g} {r['accept']:>7.3f} {r['rhat']:>7.3f} "
              f"{r['tau_sweeps']:>8.1f} {r['ess']:>6.0f} {r['vdv']:>8.3f} "
              f"{r['above']:>8.3f} {r['below']:>8.3f} {r['gap']:>7.3f} "
              f"{r['edge']:>7.4f} {eq}")

    summ = os.path.join(args.output_dir, "summary.csv")
    with open(summ, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())) if rows else None
        if w:
            w.writeheader(); w.writerows(rows)
    print(f"\nwrote {summ}")
    print("edge-of-equilibratability = largest coef with R-hat<1.05, small gap, usable ESS")


if __name__ == "__main__":
    main()
