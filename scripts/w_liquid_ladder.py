#!/usr/bin/env python3
"""Version-(b) hyperuniformity test: the UNDOPED liquid under the
frustration-free TCP chemistry (EDQ + Z-legality + impurity m^2, NO tilts).

All previous equilibrium hyperuniformity nulls were measured under the
VDV/HDV variance objectives — on-site self-energy terms with no charge-charge
coupling. This driver equilibrates liquid ensembles under the W-family
objective at a ladder of couplings lambda (below the freeze), then runs the
curvature-charge estimators (window variance + low-pass spectral) on each
final state. Question: does the string-net chemistry bend S(k->0) down with
coupling where the variance objective never did?

Per (lambda, seed): equilibrate in segments (census logged each segment as a
drift check), save final state, measure. Output: CSV + JSON + states.

    python scripts/w_liquid_ladder.py --lams 0.05 0.1 0.2 0.3 0.5 0.7 --reps 2
"""
import argparse
import csv
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_one(init_path, lam, rep, cfg, out_dir):
    import discrete_differential_geometry as ddg
    from fk_skeleton import edges_from_facets, vertex_class_census

    m = ddg.Manifold.load(init_path, 3)
    params = ddg.SamplerParams(
        num_facets_target=cfg["n"], num_facets_coef=cfg["nfc"],
        hinge_degree_target=cfg["edge"], num_hinges_coef=cfg["k"],
        hinge_degree_target_coef=lam * cfg["edge"] / 6.0)
    s = ddg.ManifoldSampler(m, params)
    s.set_n6_potential(cfg["w_scale"] * lam, cfg["w_scale"] * lam)
    v = s.manifold

    rows = []
    for seg in range(cfg["segments"]):
        s.reset_stats()
        s.run(sweeps=cfg["seg_sweeps"])
        st = s.get_stats()
        eu, edeg, V = edges_from_facets(v.facets())
        fz, _ = vertex_class_census(eu, edeg, V)
        rows.append(dict(
            seg=seg, sweeps=(seg + 1) * cfg["seg_sweeps"],
            acc=st.total_accepted / max(1, st.total_tried),
            mean_edeg=float(edeg.mean()), edv=float(np.var(edeg)),
            pure56=1.0 - fz["impure"],
            fZ12=fz["Z12"], fZ14=fz["Z14"], fZ15=fz["Z15"], fZ16=fz["Z16"]))
    path = os.path.join(out_dir, f"lam{lam:g}_rep{rep}.mfd")
    v.save(path)

    # estimators on the final state
    from curvature_hyperuniformity_g import (load_fields, bfs_order,
                                             window_variances, lowpass_ratio)
    rng = np.random.default_rng(1000 + rep)
    qR, qV, adj, eu2, V2 = load_fields(path)
    mmax = V2 // 6
    mgrid = np.unique(np.geomspace(8, mmax, 12).astype(int))
    orders = [bfs_order(adj, V2, int(rng.integers(V2)), mmax)
              for _ in range(cfg["ncent"])]
    vr = window_variances(qR, orders, mgrid)
    vs = np.mean([window_variances(rng.permutation(qR), orders, mgrid)
                  for _ in range(4)], axis=0)
    ratio = (vr / vs)
    lp = lowpass_ratio(eu2, V2, qR - qR.mean(), rng=rng)
    return dict(lam=lam, rep=rep, trace=rows, mgrid=mgrid.tolist(),
                win_ratio=ratio.tolist(), win_max=float(ratio[-1]),
                lowpass=float(lp), state=path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--init-glob",
                    default="seeds/S3_N31623_1e-1_ED5p1043_2_VDVs_2e-3_HDVs_2e-1_s0*.mfd")
    ap.add_argument("--lams", type=float, nargs="+",
                    default=[0.05, 0.1, 0.2, 0.3, 0.5, 0.7])
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--n", type=int, default=31623)
    ap.add_argument("--edge", type=float, default=5.1043)
    ap.add_argument("--k", type=float, default=2.0)
    ap.add_argument("--nfc", type=float, default=0.1)
    ap.add_argument("--w-scale", type=float, default=0.3,
                    help="zleg and cimp coefficients = this * lambda")
    ap.add_argument("--segments", type=int, default=3)
    ap.add_argument("--seg-sweeps", type=int, default=2000)
    ap.add_argument("--ncent", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out-dir", default="data/w_liquid_ladder")
    args = ap.parse_args()

    inits = sorted(glob.glob(os.path.join(_ROOT, args.init_glob)))
    if len(inits) < args.reps:
        sys.exit(f"need {args.reps} seeds for {args.init_glob}, found {len(inits)}")
    out_dir = os.path.join(_ROOT, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    cfg = dict(n=args.n, edge=args.edge, k=args.k, nfc=args.nfc,
               w_scale=args.w_scale, segments=args.segments,
               seg_sweeps=args.seg_sweeps, ncent=args.ncent)

    jobs = [(inits[r], lam, r, cfg, out_dir)
            for lam in args.lams for r in range(args.reps)]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_one, *j): j for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            tr = r["trace"][-1]
            print(f"[lam={r['lam']:g} rep{r['rep']}] acc={tr['acc']:.4f} "
                  f"mean_edeg={tr['mean_edeg']:.4f} pure56={tr['pure56']:.4f} "
                  f"winR(max)={r['win_max']:.4f} lowpassR={r['lowpass']:.6f}",
                  flush=True)

    with open(os.path.join(out_dir, "results.json"), "w") as f:
        json.dump(dict(cfg=cfg, lams=args.lams, results=results), f, indent=1)
    print(f"wrote {out_dir}/results.json")


if __name__ == "__main__":
    main()
