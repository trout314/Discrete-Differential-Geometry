#!/usr/bin/env python3
"""Roundness test for equilibrium seed families.

For each family (glob of replica .mfd files) build the ensemble dual-graph
distance distribution (multi-source BFS over facet adjacency, averaged over
replicas) and compare it to a round S^3:

  * Hausdorff dimension d_H from the small-r volume-growth slope (round S^3 -> 3).
  * KS distance between the mean-rescaled empirical CDF and the round-S^3 CDF
    F(u) = (u - sin u cos u)/pi on u in [0, pi]  (Delta_S3; 0 = perfect round).

Distances rescaled by matching mean graph-distance to the round-S^3 mean pi/2,
so only SHAPE is compared (graph vs geodesic distance differ by a lattice
constant). Prints a table and writes a CDF-overlay plot.
"""

import glob
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.dirname(__file__))
from discrete_differential_geometry import Manifold
from distance_distribution import read_edge_list, build_csr, distance_histogram

# Round S^3 (unit): shell density (2/pi) sin^2(u), CDF below, mean pi/2.
def round_s3_cdf(u):
    return (u - np.sin(u) * np.cos(u)) / np.pi


def family_hist(pattern, sources, rng, graph="dual", dim=3):
    hist = np.zeros(1, dtype=np.int64)
    seeds = sorted(glob.glob(pattern))
    with tempfile.TemporaryDirectory() as tmp:
        for s in seeds:
            m = Manifold.load(s, dim); fv = list(m.f_vector)
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
    return hist, len(seeds)


def analyze(hist):
    r = np.arange(len(hist)); tot = hist.sum()
    cdf = np.cumsum(hist) / tot
    mean = (r * hist).sum() / tot
    diam = int(r[hist > 0][-1])

    def dh(lo, hi):
        m = (cdf > lo) & (cdf < hi) & (r > 0)
        return float(np.polyfit(np.log(r[m]), np.log(cdf[m]), 1)[0]) if m.sum() >= 3 else np.nan

    u = r * (np.pi / 2) / mean               # rescale: match mean -> pi/2
    m = (u <= np.pi) & (r > 0)
    ks = float(np.max(np.abs(cdf[m] - round_s3_cdf(u[m])))) if m.any() else np.nan
    return dict(mean=mean, diam=diam, dh1=dh(.02, .3), dh2=dh(.05, .5),
                ks=ks, u=u, cdf=cdf)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", choices=["edge", "dual"], default="dual")
    ap.add_argument("--preset", choices=["eq", "hi"], default="eq",
                    help="eq = VDV~13 equilibrium families; hi = uncontrolled "
                         "beta=0 high-VDV families (VDV 570/1540/3079).")
    ap.add_argument("--ed", default="ED5p0043",
                    help="Edge-degree family tag for the eq preset (ED5p0043 / "
                         "ED5p1043 / ED5p2043).")
    args = ap.parse_args()
    graph = args.graph
    if args.preset == "eq":
        ed = args.ed
        families = [
            # beta/N=0.008 (VDVs_8e-3) at the k=2 edge pin; new library naming
            ("N=1e3", f"seeds/S3_N1e3_1e-1_{ed}_2_VDVs_8e-3_s*.mfd", None),
            ("N=1e4", f"seeds/S3_N1e4_1e-1_{ed}_2_VDVs_8e-3_s*.mfd", 1000),
            ("N=1e5", f"seeds/S3_N1e5_1e-1_{ed}_2_VDVs_8e-3_s*.mfd", 400),
        ]
    else:
        families = [
            ("N=1e3", "seeds/S3_N1e3_1e-1_ED5p0043_1e-1_s0[01]*.mfd", None),
            ("N=1e4", "seeds/S3_N1e4_1e-1_ED5p0043_1e-1_s0[01]*.mfd", 1000),
            ("N=1e5", "seeds/S3_N1e5_1e-1_ED5p0043_1e-1_s*.mfd", 400),
        ]
    rng = np.random.default_rng(0)
    print(f"{graph}-distance roundness")
    print(f"{'family':>7} {'seeds':>5} {'diam':>5} {'mean':>6} "
          f"{'d_H[.02-.3]':>11} {'d_H[.05-.5]':>11} {'Delta_S3(KS)':>12}")
    results = []
    for name, pat, src in families:
        hist, nseed = family_hist(pat, src, rng, graph=graph)
        if nseed == 0:
            print(f"{name:>7} : no seeds ({pat})"); continue
        a = analyze(hist); a["name"] = name; results.append(a)
        print(f"{name:>7} {nseed:>5} {a['diam']:>5} {a['mean']:>6.1f} "
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
    ax.set_title(f"{graph}-distance CDF of equilibrium seeds (VDV~13) vs round S^3")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    tag = args.ed if args.preset == "eq" else "hi"
    out = f"/tmp/claude-1000/roundness_{tag}_{graph}.png"
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
