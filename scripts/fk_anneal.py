#!/usr/bin/env python3
"""Rate-controlled annealing into the TCP/quasicrystal regime, with FK census.

Companion to scripts/fk_skeleton.py (see memory: QC/phason/TCP program). The
static census showed the certified library is a dense disclination fluid with
zero Frank-Kasper order: EDV freezes ~200x above the TCP scale at the glass
wall. This script asks the dynamic question: when we ramp the degree-variance
couplings BEYOND the wall, does FK/TCP order develop, and does annealing more
SLOWLY get further before arrest?

Protocol: start from a certified equilibrium seed (default: the strongest
combo family VDVs_2e-3_HDVs_0p405, our lowest-EDV equilibrium). Ramp BOTH the
vertex (VDV) and hinge (HDV) degree-variance coefficients geometrically by the
same total factor along the same coupling path; different --rates values give
different dwell times (sweeps per step) on that identical path, so any
difference in the census at matched coupling is pure annealing-rate physics
(coarsening-limited order = slower gets further; ideal-glass wall = rate
independent arrest). After the ramp, hold at max coupling. FK census (from
fk_skeleton) is recorded inline at every step; .mfd snapshots are written at
each decade of the ramp factor and at the end of the hold.

Output per (rate, replica): <out>/rate<R>_rep<K>.csv + snapshots
<out>/rate<R>_rep<K>.snap_f<factor>.mfd + manifest.json.

    python scripts/fk_anneal.py --rates 3 30 300 3000 --reps 2
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
from fk_skeleton import edges_from_facets, vertex_class_census, skeleton_stats

CSV_COLS = ["step", "cum_sweeps", "vdv_coef", "hdv_coef", "acceptance",
            "edv", "vdv_obs", "mean_edeg", "p_le4", "p5", "p6", "p_ge7",
            "pure56", "fFK", "frac_val1", "largest_frac", "cyclomatic",
            "n_sk_edges", "num_facets"]


def census_row(mfd_view):
    eu, edeg, V = edges_from_facets(mfd_view.facets())
    fz, _ = vertex_class_census(eu, edeg, V)
    sk = skeleton_stats(eu, edeg, V)
    return dict(
        edv=float(np.var(edeg)), mean_edeg=float(edeg.mean()),
        p_le4=float(np.mean(edeg <= 4)), p5=float(np.mean(edeg == 5)),
        p6=float(np.mean(edeg == 6)), p_ge7=float(np.mean(edeg >= 7)),
        pure56=1.0 - fz["impure"],
        fFK=fz["Z12"] + fz["Z14"] + fz["Z15"] + fz["Z16"],
        frac_val1=sk["frac_val1"], largest_frac=sk["largest_frac"],
        cyclomatic=sk["cyclomatic"], n_sk_edges=sk["n_edges"])


def run_one(init_path, cfg, rate, rep, out_dir):
    from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
    m = Manifold.load(init_path, 3)
    params = SamplerParams(
        num_facets_target=cfg["n"], hinge_degree_target=cfg["edge"],
        num_facets_coef=cfg["nfc"], num_hinges_coef=cfg["k"],
        hinge_degree_variance_coef=cfg["hdv0"],
        codim3_degree_variance_coef=cfg["vdv0"])
    s = ManifoldSampler(m, params)
    v = s.manifold

    tag = f"rate{rate}_rep{rep}"
    rows = []
    factor = cfg["total_factor"] ** (1.0 / cfg["steps"])
    next_snap_decade = 1  # snapshot when cumulative factor crosses 10**k

    def record(step, cum, acc, vdv_c, hdv_c):
        r = census_row(v)
        r.update(step=step, cum_sweeps=cum, vdv_coef=vdv_c, hdv_coef=hdv_c,
                 acceptance=acc, vdv_obs=v.degree_variance(0),
                 num_facets=v.num_facets)
        rows.append([r[c] for c in CSV_COLS])

    record(0, 0, np.nan, cfg["vdv0"], cfg["hdv0"])
    cum = 0
    for step in range(1, cfg["steps"] + 1):
        f_cum = factor ** step
        vdv_c = cfg["vdv0"] * f_cum
        hdv_c = cfg["hdv0"] * f_cum
        s.set_codim3_degree_variance_coef(vdv_c)
        s.set_hinge_degree_variance_coef(hdv_c)
        if cfg.get("zleg0") or cfg.get("cimp0"):
            s.set_n6_potential(cfg.get("zleg0", 0.0) * f_cum,
                               cfg.get("cimp0", 0.0) * f_cum)
        s.reset_stats()
        s.run(sweeps=rate)
        cum += rate
        st = s.get_stats()
        acc = st.total_accepted / max(1, st.total_tried)
        record(step, cum, acc, vdv_c, hdv_c)
        if f_cum >= 10 ** next_snap_decade - 1e-9:
            v.save(os.path.join(out_dir, f"{tag}.snap_f1e{next_snap_decade}.mfd"))
            next_snap_decade += 1

    # Hold at max coupling, census every hold/16.
    chunk = max(1, cfg["hold"] // 16)
    held = 0
    while held < cfg["hold"]:
        s.reset_stats()
        s.run(sweeps=chunk)
        held += chunk
        cum += chunk
        st = s.get_stats()
        record(cfg["steps"], cum, st.total_accepted / max(1, st.total_tried),
               vdv_c, hdv_c)
    v.save(os.path.join(out_dir, f"{tag}.final.mfd"))

    with open(os.path.join(out_dir, f"{tag}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(CSV_COLS)
        w.writerows(rows)
    return tag, rows[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--init-glob",
                    default="seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p405_s0*.mfd")
    ap.add_argument("--rates", type=int, nargs="+", default=[3, 30, 300, 3000],
                    help="sweeps per ramp step (dwell time; one run per rate)")
    ap.add_argument("--reps", type=int, default=2, help="replicas per rate")
    ap.add_argument("--steps", type=int, default=32, help="ramp steps")
    ap.add_argument("--total-factor", type=float, default=100.0,
                    help="total coupling multiplication over the ramp")
    ap.add_argument("--hold", type=int, default=10000,
                    help="sweeps held at max coupling after the ramp")
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--edge", type=float, default=5.0043)
    ap.add_argument("--k", type=float, default=2.0)
    ap.add_argument("--nfc", type=float, default=0.1)
    ap.add_argument("--vdv0", type=float, default=20.0)
    ap.add_argument("--hdv0", type=float, default=4050.0)
    ap.add_argument("--zleg0", type=float, default=0.0,
                    help="Z-legality potential start coef (ramped by the same "
                         "factor; frustration-free TCP pressure)")
    ap.add_argument("--cimp0", type=float, default=0.0,
                    help="impurity-valence m^2 start coef (ramped)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out-dir", default="data/fk_anneal/N1e4_x100")
    args = ap.parse_args()

    inits = sorted(glob.glob(os.path.join(_ROOT, args.init_glob)))
    if len(inits) < args.reps:
        sys.exit(f"need {args.reps} seeds, found {len(inits)} for {args.init_glob}")
    out_dir = os.path.join(_ROOT, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    cfg = dict(n=args.n, edge=args.edge, k=args.k, nfc=args.nfc,
               vdv0=args.vdv0, hdv0=args.hdv0, steps=args.steps,
               total_factor=args.total_factor, hold=args.hold,
               zleg0=args.zleg0, cimp0=args.cimp0)
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(dict(cfg, rates=args.rates, reps=args.reps,
                       inits=inits[:args.reps]), f, indent=1)

    jobs = [(inits[rep], cfg, rate, rep, out_dir)
            for rate in args.rates for rep in range(args.reps)]
    # Longest jobs first so the pool tail is short.
    jobs.sort(key=lambda j: -j[2])
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_one, *j): j for j in jobs}
        for fut in as_completed(futs):
            tag, last = fut.result()
            d = dict(zip(CSV_COLS, last))
            print(f"[{tag}] done: EDV={d['edv']:.3f} pure56={d['pure56']:.4f} "
                  f"fFK={d['fFK']:.4f} p_le4={d['p_le4']:.3f} "
                  f"acc={d['acceptance']:.4f}", flush=True)
    print(f"wrote {out_dir}")


if __name__ == "__main__":
    main()
