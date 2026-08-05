#!/usr/bin/env python3
"""Four-arm summary figure for the blob-dispersal transport A/B.

    python scripts/defect_dynamics/ecmc_ab_figure.py \
        data/ecmc_ab/main_nfc30 data/ecmc_ab/matchedN
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecmc_ab_analyze import load_dir, mean_curve  # noqa: E402

ARMS = [
    ("A", "main", "A  plain (local Pachner + 4-4 hinge)", "#4053d3"),
    ("B", "main", "B  + lifted ECMC flight  (~8 transport steps/sweep)",
     "#dd3f74"),
    ("N", "matched",
     "M  + diffusive nonlocal slide (~8/sweep, MATCHED budget)", "#f4a11c"),
    ("N", "main", "N  + diffusive nonlocal slide (~184/sweep)", "#00b25d"),
]


def main():
    main_dir, matched_dir = sys.argv[1], sys.argv[2]
    src = {"main": load_dir(main_dir)[0], "matched": load_dir(matched_dir)[0]}
    cfg = load_dir(main_dir)[1].get("cfg", {})

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    wmax = min(c["work"][-1] for a, s, _, _ in ARMS for c in src[s][a])
    grid = np.geomspace(5e4, wmax, 220)
    for key, a, ylab in (("spread", ax[0],
                          "mean graph distance of defect vertices\n"
                          "from the blob center  (dispersal)"),
                         ("m2", ax[1], r"$\sum_v m(v)^2$  (crowding)")):
        for arm, s, label, col in ARMS:
            cs = src[s][arm]
            mu, se = mean_curve(cs, key, grid)
            a.plot(grid, mu, color=col, lw=1.8,
                   label=f"{label}   [{len(cs)} chains]")
            a.fill_between(grid, mu - se, mu + se, color=col, alpha=0.18,
                           lw=0)
        a.set_xscale("log")
        a.set_xlabel("work  (attempted moves, all channels)")
        a.set_ylabel(ylab)
        a.grid(alpha=0.15, lw=0.5)
    ax[0].axhline(3.98, color="gray", ls=":", lw=1)
    ax[0].annotate("fully dispersed (3.98)", (grid[0], 3.985), fontsize=7.5,
                   color="gray", va="bottom")
    ax[0].annotate("packed blob", (grid[0], 3.07), fontsize=7.5,
                   color="gray", va="top")
    ax[0].legend(fontsize=7.6, loc="lower right", framealpha=0.95)

    fig.suptitle(
        "Blob dispersal on C15 m3 N3672 (T3, 648 v) -- "
        f"S = {cfg.get('nfc')}(f3-f3*)^2 + {cfg.get('k1')}(f1-6f3/e*)^2 "
        f"+ {cfg.get('cimp')}Σm^2,  k={cfg.get('k_fliers')} fliers, "
        f"{cfg.get('sweeps')} sweeps/chain\n"
        "the volume pin closes the kill/rebirth teleport channel, so "
        "relaxation is transport-limited  --  PROVISIONAL (uncertified)",
        fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = os.path.join(os.path.dirname(main_dir.rstrip("/")),
                       "ecmc_ab_four_arm.png")
    fig.savefig(out, dpi=150)
    print(f"Saved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
