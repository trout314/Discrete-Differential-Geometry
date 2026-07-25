#!/usr/bin/env python3
"""Equilibrium statics of the defect gas: interactions and hyperuniformity.

Two questions, both answerable from existing snapshots without any dynamics:

1. DO DEFECT COMPLEXES INTERACT?  Interaction is an equilibrium observable:
   if complexes are an ideal (Poisson) gas their pair statistics factorise, so
   any structure in g(r), in the charge-charge correlator, or in the
   chord-director correlator is interaction, regardless of how (or whether)
   the complexes move. Reported:
     g(r)      pair correlation of complex centroids (torus minimal image)
     C_QQ(r)   <dQ_i dQ_j>(r) / <dQ^2>, dQ = Q_c - n_verts * qbar
     P2(r)     nematic correlator of chord directors (the degree-3 edge axis),
               for complex pairs that each carry exactly one chord
2. IS THE PERSISTENT COMPONENT HYPERUNIFORM UNDER THE FLICKER FLOOR?
   Ball-variance decomposition: sigma^2(R) of the summed defect charge in
   Euclidean balls of the lift chart, split by flicker attribution
   (crystal_flicker keys, an upper bound as always), with the same-charges
   uniform-resprinkling Poisson reference. Also the F/P cross-covariance:
   negative means the vacuum screens the strings.

Positions are the cocycle lift (tree_positions on the saved cocycle), i.e.
crystal-chart coordinates -- REAL Euclidean balls, which the graph-distance
proxy is blind to (curvature length-scale lesson). Ball radii stop at half
the box so shell/ball volumes are exact on the torus.
"""
import argparse
import glob as globmod
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
import defect_state as ds
import crystal_flicker as cf

CELL = 1.0e6


def minimg(d, box):
    return d - np.round(d / box) * box


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--glob", required=True, action="append")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--flicker", required=True)
    ap.add_argument("--lam", type=float, default=None)
    ap.add_argument("--centers", type=int, default=3000,
                    help="ball centers per snapshot")
    ap.add_argument("--poisson-reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    files = sorted({f for g in args.glob for f in globmod.glob(g)})
    if not files:
        sys.exit(f"no snapshots matched {args.glob}")
    box = float(args.mcell)                       # cells
    flick, fmeta = cf.load_flicker(args.flicker, "full")
    rng = np.random.default_rng(args.seed)

    snaps = []
    nill_per = []
    for f in files:
        cocf = f.replace(".mfd", ".cocycle.npz")
        if not os.path.exists(cocf):
            print(f"  SKIP {f} (no cocycle)")
            continue
        m = ddg.Manifold.load(f, 3)
        st = ds.DefectState(m)
        e1, w1, _ = coc.load_cocycle(cocf)
        e1 = np.asarray(e1)
        pos = coc.tree_positions(e1, np.asarray(w1),
                                 int(e1.max()) + 1)[0].astype(float) / CELL
        q = st.vertex_charges()
        qbar = st.qbar(q)
        nill_per.append(ds.census(st)["n_illegal"])
        comps = []
        for cx in st.components(broad=True):
            verts = cx.verts
            p = pos[verts]
            rel = minimg(p - p[0], box)
            cen = (p[0] + rel.mean(0)) % box
            Q = st.complex_charge(verts, q)
            dQ = Q - len(verts) * qbar
            fac = st.induced_facets(verts)
            vc, ec = st.decorations(verts)
            isf = cf.flicker_key("full", facets=fac, vcolor=vc,
                                 ecolor=ec) in flick
            vset = set(verts)
            chords = [e for e in st.ill_edges
                      if st.edeg[e] == 3 and e[0] in vset and e[1] in vset]
            axis = None
            if len(chords) == 1:
                d = minimg(pos[chords[0][1]] - pos[chords[0][0]], box)
                n = np.linalg.norm(d)
                if n > 0:
                    axis = (d / n).tolist()
            comps.append({"cen": cen.tolist(), "n": len(verts),
                          "Q": Q, "dQ": dQ, "flick": isf, "axis": axis})
        snaps.append({"file": os.path.basename(f), "comps": comps,
                      "qbar": qbar})
        nf = sum(c["flick"] for c in comps)
        print(f"  {os.path.basename(f)}: {len(comps)} complexes "
              f"({nf} flicker-reachable)", flush=True)
        del st, m

    ns = len(snaps)
    print(f"\n{ns} snapshots, {sum(len(s['comps']) for s in snaps)} complexes "
          f"({sum(c['flick'] for s in snaps for c in s['comps'])} "
          f"flicker-reachable)")

    # ---------------- pair statistics, pooled over snapshots ---------------
    edges_r = np.linspace(0.0, box / 2, 13)       # shells to half box
    mid = 0.5 * (edges_r[1:] + edges_r[:-1])
    vshell = 4 * np.pi / 3 * (edges_r[1:] ** 3 - edges_r[:-1] ** 3)
    vbox = box ** 3

    def pair_arrays(sel):
        """(r, dQi*dQj, p2 or nan) over all within-snapshot pairs."""
        rr, qq, pp = [], [], []
        for s in sel:
            cs = s["comps"]
            for i in range(len(cs)):
                for j in range(i + 1, len(cs)):
                    d = minimg(np.array(cs[i]["cen"]) - np.array(cs[j]["cen"]),
                               box)
                    r = float(np.linalg.norm(d))
                    rr.append(r)
                    qq.append(cs[i]["dQ"] * cs[j]["dQ"])
                    if cs[i]["axis"] and cs[j]["axis"]:
                        c2 = float(np.dot(cs[i]["axis"], cs[j]["axis"])) ** 2
                        pp.append((3 * c2 - 1) / 2)
                    else:
                        pp.append(np.nan)
        return (np.array(rr), np.array(qq), np.array(pp))

    rr, qq, pp = pair_arrays(snaps)
    npairs_exp = sum(len(s["comps"]) * (len(s["comps"]) - 1) / 2
                     for s in snaps)
    allQ = np.array([c["dQ"] for s in snaps for c in s["comps"]])
    # CONNECTED correlator: <dQ> != 0 (defect complexes carry net negative
    # curvature vs the background by construction), so the raw product is
    # positive at every r no matter what. Subtract the disconnected part.
    mQ, vQ = allQ.mean(), allQ.var()
    dq2 = float(np.mean(allQ ** 2))
    # P2 NULL is not zero: chords lie along crystallographic axes, so even
    # uncorrelated pairs share a discrete axis set. The null is the
    # cross-snapshot pair average (same axis distribution, no possible
    # interaction).
    axes = [np.array(c["axis"]) for s in snaps for c in s["comps"]
            if c["axis"]]
    snap_of = [i for i, s in enumerate(snaps) for c in s["comps"]
               if c["axis"]]
    cross = [(3 * float(np.dot(a, b)) ** 2 - 1) / 2
             for i in range(len(axes)) for j in range(i + 1, len(axes))
             if snap_of[i] != snap_of[j]
             for a, b in ((axes[i], axes[j]),)]
    p2_null = float(np.mean(cross)) if cross else 0.0
    p2_null_se = float(np.std(cross) / np.sqrt(len(cross))) if cross else 0.0

    def bincurve(rr, w, edges):
        idx = np.digitize(rr, edges) - 1
        out, cnt = [], []
        for b in range(len(edges) - 1):
            m_ = idx == b
            cnt.append(int(m_.sum()))
            out.append(float(np.nanmean(w[m_])) if m_.sum() else np.nan)
        return np.array(out), np.array(cnt)

    # g(r): observed pair count per shell / uniform expectation
    obs, _ = np.histogram(rr, edges_r)
    exp = npairs_exp * vshell / vbox
    g = obs / np.where(exp > 0, exp, 1)
    # bootstrap over snapshots
    def boot(stat, reps=400):
        vals = []
        for _ in range(reps):
            pick = [snaps[i] for i in rng.integers(0, ns, ns)]
            vals.append(stat(pick))
        return np.nanstd(np.array(vals, dtype=float), axis=0)

    def g_of(sel):
        r2, _, _ = pair_arrays(sel)
        o, _ = np.histogram(r2, edges_r)
        e2 = sum(len(s["comps"]) * (len(s["comps"]) - 1) / 2
                 for s in sel) * vshell / vbox
        return o / np.where(e2 > 0, e2, 1)

    gerr = boot(g_of)
    cqq, cqq_n = bincurve(rr, (qq - mQ ** 2) / vQ, edges_r)
    p2, p2_n = bincurve(rr, pp, edges_r)

    print(f"\n<dQ> = {mQ:+.3f}, Var(dQ) = {vQ:.3f} -- C_QQ below is the "
          f"CONNECTED correlator (0 = no interaction)")
    print(f"P2 null from cross-snapshot pairs: {p2_null:+.3f} +- "
          f"{p2_null_se:.3f} (crystal-axis anisotropy; subtract before "
          f"reading alignment)")
    print("NOTE the g(r) hole below r ~ 0.4 is partly DEFINITIONAL: "
          "vertex-adjacent complexes merge into one component.")
    print("\npair statistics (r in cells; box/2 = %.1f):" % (box / 2))
    print(f"   {'r':>5} {'g(r)':>7} {'+-':>5} {'C_QQ conn':>9} {'nq':>4} "
          f"{'P2-null':>8} {'np':>4}")
    for i in range(len(mid)):
        print(f"   {mid[i]:5.2f} {g[i]:7.3f} {gerr[i]:5.3f} "
              f"{cqq[i]:9.3f} {cqq_n[i]:4d} "
              f"{(p2[i]-p2_null):8.3f} {p2_n[i]:4d}")

    # ---------------- ball variance decomposition --------------------------
    Rgrid = np.array([0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0])

    def ball_stats(get_pts):
        """Var of ball charge for [all, flicker, persistent] + F/P cov."""
        tot = {R: [] for R in Rgrid}
        fpart = {R: [] for R in Rgrid}
        ppart = {R: [] for R in Rgrid}
        for s in snaps:
            pts, chg, isf = get_pts(s)
            if len(pts) == 0:
                continue
            cen = rng.uniform(0, box, (args.centers, 3))
            d = minimg(cen[:, None, :] - pts[None, :, :], box)
            dist = np.sqrt((d ** 2).sum(-1))
            for R in Rgrid:
                w = dist < R
                tot[R].extend((w * chg).sum(1))
                fpart[R].extend((w * (chg * isf)).sum(1))
                ppart[R].extend((w * (chg * ~isf)).sum(1))
        out = []
        for R in Rgrid:
            t = np.array(tot[R]); f_ = np.array(fpart[R]); p_ = np.array(ppart[R])
            out.append({"R": float(R), "var_tot": float(np.var(t)),
                        "var_F": float(np.var(f_)),
                        "var_P": float(np.var(p_)),
                        "cov_FP": float(np.cov(f_, p_)[0, 1])})
        return out

    def real_pts(s):
        cs = s["comps"]
        pts = np.array([c["cen"] for c in cs]) if cs else np.empty((0, 3))
        chg = np.array([c["dQ"] for c in cs])
        isf = np.array([c["flick"] for c in cs], bool)
        return pts, chg, isf

    real = ball_stats(real_pts)

    # Poisson reference: same charges + flags per snapshot, uniform positions
    pois_acc = None
    for _ in range(args.poisson_reps):
        def pois_pts(s):
            cs = s["comps"]
            pts = rng.uniform(0, box, (len(cs), 3))
            chg = np.array([c["dQ"] for c in cs])
            isf = np.array([c["flick"] for c in cs], bool)
            return pts, chg, isf
        p = ball_stats(pois_pts)
        if pois_acc is None:
            pois_acc = [{k: v for k, v in row.items()} for row in p]
        else:
            for a, b in zip(pois_acc, p):
                for k in ("var_tot", "var_F", "var_P", "cov_FP"):
                    a[k] += b[k]
    for a in pois_acc:
        for k in ("var_tot", "var_F", "var_P", "cov_FP"):
            a[k] /= args.poisson_reps

    print("\nball variance of defect charge (rad^2), real vs Poisson "
          "resprinkle (same charges):")
    print(f"   {'R':>4} {'var_tot':>9} {'pois':>9} {'ratio':>6}  "
          f"{'var_F':>9} {'var_P':>9} {'cov_FP':>9} {'pois_cov':>9}")
    for r_, p_ in zip(real, pois_acc):
        print(f"   {r_['R']:4.1f} {r_['var_tot']:9.4f} {p_['var_tot']:9.4f} "
              f"{r_['var_tot']/max(p_['var_tot'],1e-12):6.3f}  "
              f"{r_['var_F']:9.4f} {r_['var_P']:9.4f} {r_['cov_FP']:9.4f} "
              f"{p_['cov_FP']:9.4f}")

    with open(args.out, "w") as fh:
        json.dump({"files": [s["file"] for s in snaps], "box": box,
                   "lam": args.lam, "n_ill_mean": float(np.mean(nill_per)),
                   "flicker_ref": fmeta,
                   "snaps": snaps,
                   "pairs": {"edges_r": edges_r.tolist(), "g": g.tolist(),
                             "g_err": gerr.tolist(), "cqq": cqq.tolist(),
                             "cqq_n": cqq_n.tolist(), "p2": p2.tolist(),
                             "p2_n": p2_n.tolist(), "dq2": float(dq2),
                             "meanQ": float(mQ), "varQ": float(vQ),
                             "p2_null": p2_null, "p2_null_se": p2_null_se},
                   "balls": {"real": real, "poisson": pois_acc}}, fh, indent=1)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
