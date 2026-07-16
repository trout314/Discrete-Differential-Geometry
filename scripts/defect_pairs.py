#!/usr/bin/env python3
"""Defect pair-correlations: the effective interaction of "matter" defects.

Matter-as-defects reading (see memory: QC/phason/TCP program): relative to the
degree-5 reference at the 5.0043 pin, edges with deg<=4 are positive
disclination defects (P, charge q=6-deg >= 2) and edges with deg>=6 are
negative defects (N, q <= 0); deg-5 edges are the vacuum. The annealed TCP
states hold a balanced ~9%/9% plasma of P and N. This script measures the
class-conditional pair correlation g_AB(r) in 1-skeleton graph distance:

    g_AB(r) = [ fraction of class-B edges among all edges at distance r
                from a class-A source edge ] / [ global fraction of class B ]

so g == 1 is "no interaction", g > 1 attraction, g < 1 repulsion, and
V_AB(r) = -ln g_AB(r) is the effective potential (in units of the sampler's
temperature). Edge-edge distance = min over endpoint vertex distances (so
r = 0 means sharing a vertex).

CAVEAT: the exact link sum rule (sum of q over a vertex's edges = 12) forces
mechanical anticorrelation between same-sign charges at r <= 1; interaction
physics should be read from r >= 2.

Output: printed table + out/defect_pairs.json (+ optional --plot PNG).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fk_skeleton import load_edges

# label -> path pattern (first match used)
STATES = [
    ("equilibrium", "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p405_s000.mfd"),
    ("dwell 3",     "data/fk_anneal/N1e4_x100/rate3_rep0.final.mfd"),
    ("dwell 30",    "data/fk_anneal/N1e4_x100/rate30_rep0.final.mfd"),
    ("dwell 300",   "data/fk_anneal/N1e4_x100/rate300_rep0.final.mfd"),
    ("dwell 3000",  "data/fk_anneal/N1e4_x100/rate3000_rep0.final.mfd"),
    ("dwell 300 N31623", "data/fk_anneal/N31623_x100/rate300_rep0.final.mfd"),
]
PAIRS = [("P", "N"), ("P", "P"), ("N", "N")]


def edge_classes(edeg):
    """P (deg<=4), V (deg==5), N (deg>=6) as boolean masks."""
    return {"P": edeg <= 4, "V": edeg == 5, "N": edeg >= 6}


def pair_correlations(eu, edeg, V, rmax=10):
    """g_AB(r) for AB in PAIRS. BFS vertex distances from defect endpoints,
    edge-edge distance = min over the 2x2 endpoint combinations."""
    cls = edge_classes(edeg)
    n_edges = len(eu)
    rows = np.concatenate([eu[:, 0], eu[:, 1]])
    cols = np.concatenate([eu[:, 1], eu[:, 0]])
    A = coo_matrix((np.ones(len(rows), np.int8), (rows, cols)),
                   shape=(V, V)).tocsr()

    out = {}
    for a, b in PAIRS:
        src_edges = np.flatnonzero(cls[a])
        src_verts = np.unique(eu[src_edges])
        # vertex distances from every vertex touched by a class-a edge
        vidx = {v: i for i, v in enumerate(src_verts)}
        D = shortest_path(A, method="D", unweighted=True, directed=False,
                          indices=src_verts)          # (nsrc_v, V)
        # per source edge: distance to every edge
        nA_r = np.zeros(rmax + 1)   # class-b targets at distance r
        nT_r = np.zeros(rmax + 1)   # all targets at distance r
        tgt_b = cls[b]
        e0, e1 = eu[:, 0], eu[:, 1]
        for se in src_edges:
            d0 = D[vidx[eu[se, 0]]]
            d1 = D[vidx[eu[se, 1]]]
            dv = np.minimum(d0, d1)                   # (V,)
            de = np.minimum(dv[e0], dv[e1])           # (n_edges,)
            de[se] = rmax + 99                        # exclude self
            de = np.clip(de, 0, rmax + 1).astype(int)
            keep = de <= rmax
            nT_r += np.bincount(de[keep], minlength=rmax + 1)
            nb = np.bincount(de[keep & tgt_b], minlength=rmax + 1)
            nA_r += nb
        pB = tgt_b.sum() / n_edges
        with np.errstate(divide="ignore", invalid="ignore"):
            g = (nA_r / nT_r) / pB
        out[f"{a}{b}"] = dict(g=g.tolist(), n_pairs=nT_r.tolist(),
                              pB=float(pB), n_src=int(len(src_edges)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rmax", type=int, default=10)
    ap.add_argument("--out", default=os.path.join(_ROOT, "out", "defect_pairs.json"))
    ap.add_argument("--plot", default=None, help="optional PNG path")
    args = ap.parse_args()

    results = {}
    for label, pat in STATES:
        paths = sorted(glob.glob(os.path.join(_ROOT, pat)))
        if not paths:
            print(f"[{label}] missing: {pat}", file=sys.stderr)
            continue
        eu, edeg, V = load_edges(paths[0])
        cls = edge_classes(edeg)
        fr = {k: float(m.mean()) for k, m in cls.items()}
        res = pair_correlations(eu, edeg, V, args.rmax)
        res["fractions"] = fr
        res["edv"] = float(np.var(edeg))
        results[label] = res
        gPN = res["PN"]["g"]
        gPP = res["PP"]["g"]
        gNN = res["NN"]["g"]
        print(f"[{label:>18}] P={fr['P']:.3f} N={fr['N']:.3f} EDV={res['edv']:.3f}")
        print(f"   r:      " + " ".join(f"{r:>6d}" for r in range(args.rmax + 1)))
        for name, g in (("g_PN", gPN), ("g_PP", gPP), ("g_NN", gNN)):
            print(f"   {name}: " + " ".join(f"{x:6.3f}" for x in g))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        # sequential blues for the dwell ladder, red for equilibrium
        colors = {"equilibrium": "#e34948", "dwell 3": "#9ec5f4",
                  "dwell 30": "#5598e7", "dwell 300": "#256abf",
                  "dwell 3000": "#104281", "dwell 300 N31623": "#1baf7a"}
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), dpi=150, sharey=True)
        for ax, (a, b) in zip(axes, PAIRS):
            for label in results:
                g = np.array(results[label][f"{a}{b}"]["g"])
                r = np.arange(len(g))
                ax.plot(r, g, "-o", ms=4, lw=1.6, color=colors.get(label, "#888"),
                        label=label)
            ax.axhline(1.0, color="#999", lw=0.8, ls="--")
            ax.set_title(f"g$_{{{a}{b}}}$(r)", fontsize=10)
            ax.set_xlabel("graph distance r", fontsize=9)
            ax.grid(alpha=0.25, lw=0.5)
            ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        axes[0].set_ylabel("pair correlation g(r)", fontsize=9)
        axes[0].legend(fontsize=7.5, frameon=False)
        fig.suptitle("Defect pair correlations: P = deg$\\leq$4, N = deg$\\geq$6"
                     " (sum rule dominates r$\\leq$1)", fontsize=10, y=1.02)
        fig.tight_layout()
        fig.savefig(args.plot, bbox_inches="tight", facecolor="white")
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
