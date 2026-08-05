#!/usr/bin/env python3
"""Figure for defect_density_hu.py: the estimator comparison + the k-resolved
complex-centroid structure factor.

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

BARS = [("charge_permnull_lowk", "curvature charge\nvs PERMUTATION null\n(the old estimator)"),
        ("charge_relocnull_lowk", "defect source charge\nvs relocation null"),
        ("vertex_lowk", "defect vertices\n(form-factor loaded)"),
        ("centroid_lowk", "complex centroids\n(form-factor FREE)")]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--out", default="data/figs/defect_density_hu.png")
    args = ap.parse_args()

    ens = []
    for f in args.jsons:
        d = json.load(open(f))
        for label, s in d["summary"].items():
            rows = d["per_snapshot"][label]
            ens.append((label, s, rows, "CONTROL" in label.upper()))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.8))

    # --- left: k-resolved centroid S(k) ------------------------------------
    for label, s, rows, is_ctl in ens:
        sh = [r["centroid_shells"] for r in rows if r.get("centroid_shells")]
        if not sh:
            continue
        k = [x["k"] for x in sh[0]]
        M = np.array([[x["ratio"] for x in row] for row in sh])
        a1.errorbar(k, M.mean(0), yerr=M.std(0, ddof=1) / np.sqrt(len(M)),
                    marker="s" if is_ctl else "o", ms=5, lw=1.3, capsize=3,
                    ls=":" if is_ctl else "-", label=label)
        if not is_ctl:
            nc = np.mean([r["n_complex"] for r in rows])
            # k (box units) of the mean inter-complex spacing L/n_c^(1/3)
            a1.axvline(nc ** (1 / 3), color="grey", lw=0.7, ls="-.", alpha=0.6)
    a1.axhline(1.0, color="k", ls="--", lw=1.0)
    a1.text(0.02, 0.06, "1.0 = Poisson (random placement)\nhyperuniform would fall to 0",
            transform=a1.transAxes, fontsize=7.5, color="0.3")
    a1.text(0.60, 0.93, "grey: inter-complex spacing\n(left of it = hydrodynamic)",
            transform=a1.transAxes, fontsize=7, color="0.45", va="top")
    a1.set_xlabel("|k|   (units of 2π/L,  L = 4 unit cells)")
    a1.set_ylabel(r"$S_{\rm centroid}(k)\ /\ S_{\rm random}$")
    a1.set_ylim(0, 1.6)
    a1.set_title("Defect COMPLEX arrangement is Poisson at every\n"
                 "accessible wavevector", fontsize=10)
    a1.legend(fontsize=8, loc="lower right")
    a1.grid(alpha=0.3)

    # --- right: the four estimators ----------------------------------------
    x = np.arange(len(BARS))
    w = 0.8 / len(ens)
    for i, (label, s, rows, is_ctl) in enumerate(ens):
        v = [s[k] if s[k] is not None else np.nan for k, _ in BARS]
        e = [s[k + "_sem"] or 0 for k, _ in BARS]
        a2.bar(x + i * w, v, w, yerr=e, capsize=3, label=label,
               hatch="///" if is_ctl else None,
               edgecolor="k" if is_ctl else "none", lw=0.6,
               color="0.75" if is_ctl else None)
    a2.axhline(1.0, color="k", ls="--", lw=1.0)
    a2.set_yscale("log")
    a2.set_xticks(x + 0.4 - w / 2)
    a2.set_xticklabels([n for _, n in BARS], fontsize=7.5)
    a2.set_ylabel("low-k ratio   (|k| ≤ 2,  1 = no effect)")
    a2.set_title("Same snapshots, four estimators.  The old one calls\n"
                 "RANDOMISED defects hyperuniform too.", fontsize=10)
    a2.set_ylim(8e-4, 60)
    a2.legend(fontsize=8, loc="upper left")
    a2.grid(alpha=0.3, axis="y")

    real = [e for e in ens if not e[3]]
    if real:
        s2 = ", ".join(f"{l}: {s['s2_frac_defect']*100:.2f}%"
                       for l, s, _, _ in real)
        a2.text(0.98, 0.97,
                f"defect share of Σdq²\n{s2.replace(', ', chr(10))}\n"
                "(crystal carries the rest,\nand only at Bragg k)",
                transform=a2.transAxes, fontsize=7, color="0.35",
                ha="right", va="top",
                bbox=dict(fc="white", ec="0.8", lw=0.5, pad=2.5))

    ns = ", ".join(f"{l}: {s['n_snapshots']} snapshots, "
                   f"⟨n_complex⟩={s['n_complex']:.1f}" for l, s, _, _ in real)
    fig.suptitle("Defect-density vs curvature-charge hyperuniformity   |   "
                 "R m4 crystal, T³, V=10176, box 4 cells   |   " + ns +
                 "\nharmonic (Tutte) frame — ratios robust, k gauge   |   "
                 "PROVISIONAL: uncertified snapshots", fontsize=8.2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
