#!/usr/bin/env python3
"""Figure for defect_density_hu.py: the k-resolved whole-defect-shuffle
structure factor + the estimator comparison.

Reads the .json files that `defect_density_hu.py --out` writes (one per
ensemble, plus the --control run) so the figure is free -- no re-measurement.

Usage:
    python scripts/defect_dynamics/defect_density_hu_figure.py \
        data/defect_hu/lam40_pooled.json data/defect_hu/lam35_pooled.json \
        data/defect_hu/lam40_control.json \
        --out data/figs/defect_density_hu.png
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BARS = [
    ("charge_permnull_lowk", "charge vs\nPERMUTATION\n(OLD)"),
    ("charge_relocnull_lowk", "charge, per-VERTEX\nrelocation\n(= NEUTRALITY)"),
    ("vertex_lowk", "defect vertices\n(form-factor\nloaded)"),
    ("centroid_lowk", "centroids\nunweighted"),
    ("rigid_charge_lowk", "WHOLE DEFECTS\nshuffled, charge-wtd\n(= ARRANGEMENT)"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--out", default="data/figs/defect_density_hu.png")
    args = ap.parse_args()

    ens = []
    for f in args.jsons:
        d = json.load(open(f))
        for label, s in d["summary"].items():
            ens.append((label, s, d["per_snapshot"][label],
                        "CONTROL" in label.upper()))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.0, 5.0))

    # --- left: k-resolved rigid-complex S(k) --------------------------------
    for label, s, rows, is_ctl in ens:
        sh = [r["rigid_charge_shells"] for r in rows if r.get("rigid_charge_shells")]
        if not sh:
            continue
        k = [x["k"] for x in sh[0]]
        M = np.array([[x["ratio"] for x in row] for row in sh])
        a1.errorbar(k, M.mean(0), yerr=M.std(0, ddof=1) / np.sqrt(len(M)),
                    marker="s" if is_ctl else "o", ms=5, lw=1.3, capsize=3,
                    ls=":" if is_ctl else "-", label=label)
        if not is_ctl:
            nc = np.mean([r["n_complex"] for r in rows])
            a1.axvline(nc ** (1 / 3), color="grey", lw=0.7, ls="-.", alpha=0.6)
    a1.axhline(1.0, color="k", ls="--", lw=1.0)
    a1.text(0.02, 0.05,
            "1.0 = complexes placed at random\nhyperuniform would fall toward 0",
            transform=a1.transAxes, fontsize=7.5, color="0.3")
    a1.text(0.98, 0.96, "grey line: inter-complex spacing\n"
                        "(only left of it is hydrodynamic)",
            transform=a1.transAxes, fontsize=7, color="0.45",
            ha="right", va="top")
    a1.set_xlabel("|k|   (units of 2π/L,  L = 4 unit cells)")
    a1.set_ylabel(r"$\left|\sum_i F_i(k)\right|^2\ /\ \sum_i |F_i(k)|^2$")
    a1.set_ylim(0, 2.0)
    a1.set_title("Whole defects shuffled rigidly: the ARRANGEMENT\n"
                 "is Poisson at every accessible wavevector", fontsize=10)
    a1.legend(fontsize=8, loc="lower right")
    a1.grid(alpha=0.3)

    # --- right: every estimator on the same snapshots -----------------------
    x = np.arange(len(BARS))
    w = 0.8 / len(ens)
    for i, (label, s, rows, is_ctl) in enumerate(ens):
        v = [s[k] if s.get(k) is not None else np.nan for k, _ in BARS]
        e = [s.get(k + "_sem") or 0 for k, _ in BARS]
        a2.bar(x + i * w, v, w, yerr=e, capsize=3, label=label,
               hatch="///" if is_ctl else None,
               edgecolor="k" if is_ctl else "none", lw=0.6,
               color="0.75" if is_ctl else None)
    a2.axhline(1.0, color="k", ls="--", lw=1.0)
    a2.set_yscale("log")
    a2.set_ylim(8e-4, 90)
    a2.set_xticks(x + 0.4 - w / 2)
    a2.set_xticklabels([n for _, n in BARS], fontsize=7.2)
    a2.set_ylabel("low-k ratio   (|k| ≤ 2,  1 = no effect)")
    a2.set_title("Same snapshots, five estimators.  The old one calls\n"
                 "RANDOMISED defects hyperuniform too.", fontsize=10)
    a2.legend(fontsize=8, loc="upper left")
    a2.grid(alpha=0.3, axis="y")

    real = [e for e in ens if not e[3]]
    if real:
        txt = "\n".join(
            f"{l}:  Σdq² share {s['s2_frac_defect']*100:.2f}%,   "
            f"Q_i = {s['Q_mean']:+.2f} ({100*s['Q_frac_negative']:.0f}% neg),  "
            f"|Q|/Σ|dq| = {s['neutrality_frac']:.2f}"
            for l, s, _, _ in real)
        fig.text(0.5, 0.005, txt, fontsize=7.2, color="0.25",
                 ha="center", va="bottom")

    ns = ",  ".join(f"{l}: {s['n_snapshots']} snaps, "
                    f"⟨n_complex⟩={s['n_complex']:.1f}, ⟨m⟩={s['mean_size']:.1f}"
                    for l, s, _, _ in real)
    fig.suptitle("Defect arrangement vs curvature-charge hyperuniformity   |   "
                 "R m4, T³, V=10176, box 4 cells   |   " + ns +
                 "\nharmonic (Tutte) frame — ratios robust, k gauge   |   "
                 "PROVISIONAL: uncertified snapshots", fontsize=8.0)
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
