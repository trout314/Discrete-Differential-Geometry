#!/usr/bin/env python3
"""Roundness test for equilibrium seed families.

For each family (glob of replica .mfd files) build the ensemble distance
distribution and compare it to a round S^3:

  * Hausdorff dimension d_H from the small-r volume-growth slope (round S^3 -> 3).
  * KS distance between the mean-rescaled empirical CDF and the round-S^3 CDF
    F(u) = (u - sin u cos u)/pi on u in [0, pi]  (Delta_S3; 0 = perfect round).

Three distance metrics (``--graph``):

  * ``dual``    -- facet-adjacency BFS (not a length; historical default).
  * ``edge``    -- 1-skeleton BFS vertex-to-vertex (a quasi-isometric PL proxy).
  * ``steiner`` -- the TRUE PL geodesic vertex-to-vertex distance (Steiner-Dijkstra,
                   steiner_geodesic.py): removes the direction-dependent lattice
                   distortion of ``edge``. The honest metric for the roundness SHAPE.

Distances are rescaled by matching the mean to the round-S^3 mean pi/2, so only
SHAPE is compared -- which makes the comparison independent of the (constant)
graph-vs-geodesic lattice factor and of the Steiner bin width. Prints a table and
writes a CDF-overlay plot. Use the same ``--vdv``/``--ns``/``--ed`` for edge and
steiner to see whether correcting the metric changes the roundness verdict.
"""

import glob
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.dirname(__file__))
from discrete_differential_geometry import Manifold
from distance_distribution import (read_edge_list, build_csr, distance_histogram,
                                    steiner_distance_histogram)

# Round S^3 (unit): shell density (2/pi) sin^2(u), CDF below, mean pi/2.
def round_s3_cdf(u):
    return (u - np.sin(u) * np.cos(u)) / np.pi


def family_hist(pattern, sources, rng, graph="dual", dim=3,
                order=3, binw=0.25, max_seeds=None):
    """Ensemble distance histogram over the replicas matching `pattern`.

    Returns (hist, n_seeds_used, width) where width is the physical size of one
    histogram bin (1.0 for edge/dual integer shells, `binw` for steiner)."""
    hist = np.zeros(1, dtype=np.int64)
    seeds = sorted(glob.glob(pattern))
    if max_seeds is not None:
        seeds = seeds[:max_seeds]
    width = binw if graph == "steiner" else 1.0
    with tempfile.TemporaryDirectory() as tmp:
        for s in seeds:
            m = Manifold.load(s, dim)
            if graph == "steiner":
                h, _, _ = steiner_distance_histogram(m, order, binw, sources, rng)
            else:
                fv = list(m.f_vector)
                gf = os.path.join(tmp, "d.txt")
                if graph == "edge":
                    m.save_edge_graph(gf); nhint = fv[0]
                else:
                    m.save_dual_graph(gf); nhint = fv[dim]
                csr, n, _ = build_csr(read_edge_list(gf), nhint)
                h, _ = distance_histogram(csr, n, sources, rng)
            if len(h) > len(hist):
                hist = np.pad(hist, (0, len(h) - len(hist)))
            hist[:len(h)] += h
    return hist, len(seeds), width


def analyze(hist, width=1.0):
    r = np.arange(len(hist)); tot = hist.sum()
    cdf = np.cumsum(hist) / tot
    mean = (r * hist).sum() / tot
    diam = int(r[hist > 0][-1])

    def dh(lo, hi):
        m = (cdf > lo) & (cdf < hi) & (r > 0)
        return float(np.polyfit(np.log(r[m]), np.log(cdf[m]), 1)[0]) if m.sum() >= 3 else np.nan

    u = r * (np.pi / 2) / mean               # rescale: match mean -> pi/2 (width cancels)
    m = (u <= np.pi) & (r > 0)
    ks = float(np.max(np.abs(cdf[m] - round_s3_cdf(u[m])))) if m.any() else np.nan
    return dict(mean=mean * width, diam=diam * width, dh1=dh(.02, .3),
                dh2=dh(.05, .5), ks=ks, u=u, cdf=cdf)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", choices=["edge", "dual", "steiner"], default="dual")
    ap.add_argument("--preset", choices=["eq", "hi"], default="eq",
                    help="eq = low-VDV equilibrium families (round candidate); "
                         "hi = uncontrolled beta=0 high-VDV families (crumpled).")
    ap.add_argument("--ed", default="ED5p0043",
                    help="Edge-degree family tag (ED5p0043 / ED5p1043 / ED5p2043).")
    ap.add_argument("--vdv", default="8e-3",
                    help="VDV token beta/N for the eq preset (e.g. 8e-3, 2e-3).")
    ap.add_argument("--ns", default="N1e3,N1e4,N1e5",
                    help="Comma list of N tokens for the eq preset, with the N prefix "
                         "as in the filenames (e.g. N1e3,N1e4,N31623,N56234).")
    ap.add_argument("--sources", default="600",
                    help="Sources per replica, or 'all' (default 600).")
    ap.add_argument("--max-seeds", type=int, default=8,
                    help="Cap replicas per family (default 8; steiner is the pricier metric).")
    ap.add_argument("--steiner-order", type=int, default=3,
                    help="Steiner subdivision order n (default 3).")
    ap.add_argument("--steiner-bin", type=float, default=0.25,
                    help="Steiner histogram bin width, edge-length units (default 0.25).")
    args = ap.parse_args()
    graph = args.graph
    sources = None if args.sources == "all" else int(args.sources)
    if args.preset == "eq":
        ed = args.ed
        families = [(tok, f"seeds/S3_{tok}_1e-1_{ed}_2_VDVs_{args.vdv}_s*.mfd", sources)
                    for tok in args.ns.split(",")]
        tag = f"{ed}_VDVs_{args.vdv}"
    else:
        families = [
            ("N=1e3", "seeds/S3_N1e3_1e-1_ED5p0043_1e-1_s0[01]*.mfd", sources),
            ("N=1e4", "seeds/S3_N1e4_1e-1_ED5p0043_1e-1_s0[01]*.mfd", sources),
            ("N=1e5", "seeds/S3_N1e5_1e-1_ED5p0043_1e-1_s*.mfd", sources),
        ]
        tag = "hi"
    rng = np.random.default_rng(0)
    extra = f" (n={args.steiner_order}, bin={args.steiner_bin})" if graph == "steiner" else ""
    print(f"{graph}-distance roundness{extra}")
    print(f"{'family':>8} {'seeds':>5} {'diam':>7} {'mean':>7} "
          f"{'d_H[.02-.3]':>11} {'d_H[.05-.5]':>11} {'Delta_S3(KS)':>12}")
    results = []
    for name, pat, src in families:
        hist, nseed, width = family_hist(pat, src, rng, graph=graph,
                                         order=args.steiner_order, binw=args.steiner_bin,
                                         max_seeds=args.max_seeds)
        if nseed == 0:
            print(f"{name:>8} : no seeds ({pat})"); continue
        a = analyze(hist, width); a["name"] = name; results.append(a)
        print(f"{name:>8} {nseed:>5} {a['diam']:>7.1f} {a['mean']:>7.1f} "
              f"{a['dh1']:>11.2f} {a['dh2']:>11.2f} {a['ks']:>12.3f}", flush=True)

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except ImportError:
        print("(no matplotlib; skipping plot)"); return
    fig, ax = plt.subplots(figsize=(7, 5))
    uu = np.linspace(0, np.pi, 200)
    ax.plot(uu, round_s3_cdf(uu), "k--", lw=2, label="round S^3: (u - sin u cos u)/pi")
    for a in results:
        m = a["u"] <= np.pi
        ax.plot(a["u"][m], a["cdf"][m], marker=".", ms=4, lw=1,
                label=f"{a['name']}  (Delta_S3={a['ks']:.2f}, d_H={a['dh1']:.1f})")
    ax.set_xlabel("rescaled distance u  (mean -> pi/2)")
    ax.set_ylabel("CDF  P(d <= u)")
    ax.set_title(f"{graph}-distance CDF vs round S^3  [{tag}]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    os.makedirs("out", exist_ok=True)
    out = f"out/roundness_{tag}_{graph}.png"
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
