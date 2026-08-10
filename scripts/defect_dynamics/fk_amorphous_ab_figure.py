#!/usr/bin/env python3
"""A/B figure for the contract/split channel as a defect annealer.

Three single-chain arms from the SAME near-legal start (a15 debt-canceling
volume, nh=30 final state: 128 illegal edges, f_FK 81.2%), identical action
except for the contract/split probability and the volume-pin stiffness.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARMS = [
    ("ab_cs000", "cs OFF, fixed zleg 2 / $\\mu$ 2", "#888888"),
    ("ab_softvol", "cs 0.05, fixed, soft volume pin", "#2ca02c"),
    ("ab_cs005", "cs 0.05, fixed zleg 2 / $\\mu$ 2", "#1f77b4"),
    ("stiffen", "cs 0.05, RAMP zleg 2$\\to$8, $\\mu$ 2$\\to$12", "#d62728"),
]


def load(path):
    rows = [json.loads(l) for l in open(path)]
    return [r for r in rows if "sw" in r]


def main():
    base = sys.argv[1]
    out = sys.argv[2]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

    for tag, label, col in ARMS:
        p = os.path.join(base, tag + ".obs.jsonl")
        if not os.path.exists(p):
            continue
        r = load(p)
        sw = [x["sw"] for x in r]
        ax[0].plot(sw, [x["n_ill_e"] for x in r], "o-", color=col, label=label)
        ax[1].plot(sw, [100.0 * x["n_fk"] / x["f0"] for x in r], "o-",
                   color=col, label=label)

    ax[0].axhline(128, ls=":", color="k", lw=1)
    ax[0].text(0.02, 0.06, "start: 128 illegal edges", transform=ax[0].transAxes,
               fontsize=8)
    ax[0].set_xlabel("sweeps")
    ax[0].set_ylabel("illegal edges  (deg $\\notin\\{5,6\\}$)")
    ax[0].set_title("defect population")
    ax[1].set_xlabel("sweeps")
    ax[1].set_ylabel("$f_{FK}$  (%)")
    ax[1].set_title("Frank-Kasper vertex fraction")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=8)

    fig.suptitle(
        "Driving a near-legal FK state toward n_ill = 0 -- PROVISIONAL "
        "(1 chain/arm, uncertified)\n"
        "start a15 m5 at debt-canceling volume, nh30 final state\n"
        "(f0 = 1000, f3* = 5700, 128 illegal edges, $f_{FK}$ = 81.2%)\n"
        "S = c_N(f3-f3*)$^2$ + 30(f1-6f3/e*)$^2$ + "
        "zleg$\\sum_v$dist$^2$(n6, \\{0,2,3,4\\}) + 1.0$\\sum_v m^2$ + "
        "$\\mu\\sum_v m$;   e* = 5.1042993;   seeds 4001 / 4002",
        fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.82])
    fig.savefig(out, dpi=150)
    print(f"Saved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
