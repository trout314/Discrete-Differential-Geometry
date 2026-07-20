"""Fit and plot R dopant-solubility isotherms from r_solubility.py output.

For each (species, mu) it pools replicas, takes the equilibrium (last-third)
average of the dopant count, subtracts the native crystal baseline to get the
EXCESS n(mu), and also averages the curvature q-bar and the web-graft/companion
observables. Emits out/r_isotherm.csv + a figure: n(mu), q-bar(mu), and the
curvature lever q-bar vs excess dopant.

Usage: python scripts/r_isotherm.py
"""
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import discrete_differential_geometry as ddg
from tcp_melt import CRYSTALS
from fk_skeleton import edges_from_facets, vertex_class_census
from dopant_pairs import CLASS_N6

OUT = "data/r_solubility"


def native_counts(struct):
    m = ddg.Manifold.load(CRYSTALS[struct], 3)
    eu, ed, V = edges_from_facets(m.facets())
    fz, _ = vertex_class_census(eu, ed, V)
    V0 = len(np.unique(m.facets()))
    return {c: fz[c] * V0 for c in ("Z12", "Z14", "Z15", "Z16")}


def read_csv(path):
    with open(path) as f:
        rows = [r for r in csv.DictReader(x for x in f
                                          if not x.startswith("#"))]
    return rows


def equil(rows, col):
    v = [float(r[col]) for r in rows if float(r["sweeps"]) >= 0.6 *
         float(rows[-1]["sweeps"])]
    return np.mean(v) if v else float("nan")


def main():
    native = native_counts("r")
    # group files by (species, mu)
    groups = defaultdict(list)
    for p in sorted(glob.glob(f"{OUT}/r_sol_*_r*.csv")):
        m = re.search(r"r_sol_(Z\d+)_mu([0-9p]+)_r\d+\.csv$", p)
        if not m:
            continue
        sp = m.group(1)
        mu = float(m.group(2).replace("p", "."))
        groups[(sp, mu)].append(p)

    rows = []
    for (sp, mu), files in sorted(groups.items()):
        n_native = native[sp]
        exc, qb, comp, graft = [], [], [], []
        for p in files:
            r = read_csv(p)
            if not r:
                continue
            exc.append(equil(r, "n_dop") - n_native)
            qb.append(equil(r, "mean_edeg"))
            comp.append(equil(r, "n_companion"))
            # graft fraction: dopant-host six-edges / all dopant six-edges
            edh = equil(r, "e_dop_host"); edd = equil(r, "e_dop_dop")
            graft.append(edh / (edh + edd) if (edh + edd) > 0 else np.nan)
        rows.append(dict(
            species=sp, mu=mu, nrep=len(exc),
            n_excess=np.mean(exc), n_excess_sem=np.std(exc, ddof=1) /
            np.sqrt(len(exc)) if len(exc) > 1 else 0.0,
            qbar=np.mean(qb), companions=np.mean(comp),
            graft_frac=np.nanmean(graft)))
        print(f"  {sp} mu={mu:<4g} nrep={len(exc)} "
              f"n_excess={np.mean(exc):8.1f} qbar={np.mean(qb):.5f} "
              f"companions={np.mean(comp):6.1f} graft={np.nanmean(graft):.2f}")

    os.makedirs("out", exist_ok=True)
    with open("out/r_isotherm.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote out/r_isotherm.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    SP = sorted({r["species"] for r in rows})
    COL = {"Z14": "#b2182b", "Z15": "#2166ac", "Z16": "#1a9850"}
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    for sp in SP:
        s = sorted([r for r in rows if r["species"] == sp],
                   key=lambda r: r["mu"])
        mu = [r["mu"] for r in s]
        ax[0].errorbar(mu, [r["n_excess"] for r in s],
                       yerr=[r["n_excess_sem"] for r in s], marker="o",
                       ms=5, capsize=2, color=COL.get(sp), label=sp)
        ax[1].plot(mu, [r["qbar"] for r in s], "o-", color=COL.get(sp),
                   label=sp)
        ax[2].plot([max(r["n_excess"], 0) for r in s],
                   [r["qbar"] for r in s], "o-", color=COL.get(sp), label=sp)
    ax[0].set_xlabel("chemical potential μ"); ax[0].set_ylabel("excess dopants n(μ)")
    ax[0].set_title("Solubility isotherm"); ax[0].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("chemical potential μ"); ax[1].set_ylabel("mean edge degree q̄")
    ax[1].axhline(5.10430, color="green", ls="--", lw=1); ax[1].set_title("Curvature response")
    ax[2].set_xlabel("excess dopants n"); ax[2].set_ylabel("q̄")
    ax[2].set_title("Curvature lever dq̄/dn"); ax[2].axhline(5.10430, color="green", ls="--", lw=1)
    for a in ax:
        a.legend(fontsize=8); a.grid(alpha=.3)
    fig.suptitle("R dopant-solubility isotherms (own-q̄ pin)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("out/r_isotherm.png", dpi=125)
    print("wrote out/r_isotherm.png")


if __name__ == "__main__":
    main()
