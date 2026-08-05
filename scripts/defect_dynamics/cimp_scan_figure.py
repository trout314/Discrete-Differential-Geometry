#!/usr/bin/env python3
"""Figure for cimp_scan.py: where the defect gas stops percolating, and what
the contact barrier costs there.

Reads the per-cimp .json files that cimp_scan.py --out writes.

Usage:
    python scripts/defect_dynamics/cimp_scan_figure.py data/ecmc_ab/cimp_*.json \
        --out data/figs/cimp_scan.png
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# dS_m2 = cimp * Delta(sum m^2) is LINEAR in cimp, and the measured median
# barrier at cimp = 0.5 was 16.0 (8,084 rejected contacts, with pins+geometry
# contributing exactly 0) -- so the barrier at any cimp is just this slope.
BARRIER_AT_HALF = 16.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--out", default="data/figs/cimp_scan.png")
    args = ap.parse_args()

    rows = {}
    for f in args.jsons:
        d = json.load(open(f))["scan"]
        for cimp, g in d.items():
            for start, r in g.items():
                if not isinstance(r, dict) or "late" not in r:
                    continue
                rows.setdefault(float(cimp), {})[start] = r["late"]

    cs = sorted(rows)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.2, 4.8))

    for start, mk, ls in (("clumped", "o", "-"), ("dispersed", "s", "--")):
        x, gf, me = [], [], []
        for c in cs:
            r = rows[c].get(start)
            if not r:
                continue
            x.append(c)
            gf.append(100 * r["top1"][0] / max(r["n_ill"][0], 1e-9))
            me.append(r["max_edeg"][0])
        a1.plot(x, gf, mk + ls, ms=6, lw=1.4, label=f"{start} start")
        a2.plot(x, me, mk + ls, ms=6, lw=1.4, label=f"{start} start")

    a1.axhline(50, color="0.6", lw=0.8, ls=":")
    a1.axvspan(0.28, 0.33, color="tab:red", alpha=0.12)
    a1.text(0.305, 60, "percolation\nthreshold", ha="center", fontsize=8,
            color="tab:red")
    a1.set_xlabel("m² coefficient  c_imp")
    a1.set_ylabel("giant component:  largest complex / all illegal edges  (%)")
    a1.set_title("A single clumped defect survives up to c_imp ≈ 0.30",
                 fontsize=10)
    a1.set_ylim(0, 108)
    a1.legend(fontsize=8)
    a1.grid(alpha=0.3)

    a2.axhline(6, color="0.6", lw=0.8, ls=":")
    a2.text(0.42, 6.25, "FK-legal ceiling (deg ≤ 6)", fontsize=7.5, color="0.4")
    a2.set_xlabel("m² coefficient  c_imp")
    a2.set_ylabel("max edge degree anywhere")
    a2.set_title("High-degree edges persist well past the percolation point",
                 fontsize=10)
    a2.legend(fontsize=8)
    a2.grid(alpha=0.3)

    # the cost of dispersal: the contact barrier at each cimp
    ax = a2.twinx()
    bx = np.array(cs)
    ax.plot(bx, BARRIER_AT_HALF * bx / 0.5, color="tab:green", lw=1.2,
            ls="-.", marker="^", ms=4)
    ax.axhline(1.0, color="tab:green", lw=0.8, ls=":")
    ax.text(0.02, 1.6, "dS = 1 (collisions become thermally accessible)",
            fontsize=7, color="tab:green")
    ax.set_ylabel("median contact barrier dS  (green)", color="tab:green")
    ax.tick_params(axis="y", labelcolor="tab:green")

    fig.suptitle(
        "Minimum m² coefficient that keeps the defect gas dispersed   |   "
        "C15 m3 N3672, 12 fliers, nfc=30, nonlocal slide on, 4000 sweeps, "
        "2 seeds, TWO-SIDED (clumped + dispersed starts)\n"
        "PROVISIONAL — uncertified", fontsize=8.2)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
