#!/usr/bin/env python3
"""Scale-aware roundness: recover the curvature signal the mean-rescaling throws away.

The roundness test rescales mean graph-distance -> pi/2, which normalizes away the
overall scale -- exactly where a small mean-curvature difference lives. Here we
instead FIT the scale: for each equilibrium family, find kappa minimizing the KS
distance between the empirical dual-distance CDF and the round-S^3 CDF evaluated
at u = kappa*r (so the family is compared to sin^2 of its own best-fit radius).

Report, per (edge-degree target, N): mean dual distance, fitted effective diameter
D_eff = pi/kappa, and the fitted-scale KS (does sin^2 still fit at the best scale?).
Then plot the scale (mean distance) against the Weyl curvature parameter
kappa_c = dbar_flat/dbar - 1 (dbar_flat = 2*pi/arccos(1/3) ~ 5.1043), one line per N.
"""

import glob
import os
import sys
import math

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))               # discrete_differential_geometry
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling: roundness_analysis
from roundness_analysis import family_hist, round_s3_cdf

DBAR_FLAT = 2 * math.pi / math.acos(1 / 3)   # ~5.10430

EDS = [("ED5p0043", 5.0043), ("ED5p1043", 5.1043), ("ED5p2043", 5.2043)]
NS = [("1e3", "VDVs_8e-3", None), ("1e4", "VDVs_8e-3", 1000), ("1e5", "VDVs_8e-3", 400)]
COLOR = {"1e3": "#4477aa", "1e4": "#ee6677", "1e5": "#228833"}


def fit_scale(hist):
    r = np.arange(len(hist)); tot = hist.sum()
    cdf = np.cumsum(hist) / tot
    mean = (r * hist).sum() / tot
    k0 = np.pi / (2 * mean)                    # exact for a perfect round fit
    best_k, best_ks = k0, 1e9
    for k in np.linspace(0.7 * k0, 1.5 * k0, 200):
        u = k * r; m = (u <= np.pi) & (r > 0)
        if m.sum() < 3:
            continue
        ks = float(np.max(np.abs(cdf[m] - round_s3_cdf(u[m]))))
        if ks < best_ks:
            best_k, best_ks = k, ks
    return mean, np.pi / best_k, best_ks


def main():
    rng = np.random.default_rng(0)
    rows = []
    print(f"{'ED':>9} {'curv kappa_c':>12} {'N':>5} {'mean':>7} "
          f"{'D_eff':>7} {'fit_KS':>7}")
    for ed, dbar in EDS:
        kc = DBAR_FLAT / dbar - 1.0
        for nt, vtag, src in NS:
            pat = f"seeds/S3_N{nt}_1e-1_{ed}_2_{vtag}_s*.mfd"
            if not glob.glob(pat):
                continue
            hist, nseed = family_hist(pat, src, rng, graph="dual")
            mean, deff, ks = fit_scale(hist)
            rows.append(dict(ed=ed, kc=kc, N=nt, mean=mean, deff=deff, ks=ks))
            print(f"{ed:>9} {kc:>+12.4f} {nt:>5} {mean:>7.1f} {deff:>7.1f} "
                  f"{ks:>7.3f}", flush=True)

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    for nt, _, _ in NS:
        pts = sorted([r for r in rows if r["N"] == nt], key=lambda r: r["kc"])
        if len(pts) >= 2:
            ax.plot([p["kc"] for p in pts], [p["mean"] for p in pts],
                    "o-", color=COLOR[nt], label=f"N={nt}")
    ax.axvline(0, color="grey", ls=":", lw=1)
    ax.text(0.001, ax.get_ylim()[0], "flat", color="grey", fontsize=8)
    ax.set_xlabel("Weyl curvature parameter  kappa_c = dbar_flat/dbar - 1  "
                  "(>0 positive, <0 negative)")
    ax.set_ylabel("mean dual graph-distance (scale)")
    ax.set_title("Effective scale vs curvature at fixed VDV~10 (equilibrium)")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("/tmp/claude-1000/scale_curvature.png", dpi=130)
    print("wrote /tmp/claude-1000/scale_curvature.png")


if __name__ == "__main__":
    main()
