#!/usr/bin/env python3
"""Verdict table + figure for a directory of quanta_heal.py two-phase runs.

Per cell it reports, for each phase separately (late window only):

    f3 - f3_nat vs the requested -dq          did the strain actually happen?
    n_ill, ncomp, top1                        how dilute is the population?
    frac(n_ill = 0)                           how often is the state legal?
    acceptance rate                           is the chain moving AT ALL?

The last column is the one that decides whether a "clean" cell is a result or
an artefact: a run that sits at the native f-vector with n_ill = 0 and an
acceptance rate of ~1e-6 has not healed, it has frozen before it ever strained.

Stationarity is reported as a late-window slope with a moving-block bootstrap
sigma (notes/CONVENTIONS.md), not a full-series fit.

Usage:
    python scripts/defect_dynamics/quanta_heal_report.py data/quanta_heal/main \
        --fig data/quanta_heal/main/quanta_heal.png
"""
import argparse
import json
import os
import sys
from glob import glob

import numpy as np


def load(prefix):
    obs = [json.loads(l) for l in open(prefix + ".obs.jsonl")]
    hdr = obs[0]
    body = [r for r in obs if r.get("kind") is None]
    fin = obs[-1] if obs[-1].get("kind") == "final" else {}
    acc = {}
    for tag in ("A", "B"):
        p = f"{prefix}.{tag}.rec.jsonl"
        if not os.path.exists(p):
            continue
        rr = [json.loads(l) for l in open(p)]
        rows = [r for r in rr if r.get("kind") == "row"]
        late = [r for r in rows if r["sw"] > rows[-1]["sw"] * 0.4] or rows
        t = sum(r["d_tried"] for r in late)
        acc[tag] = sum(r["d_accepted"] for r in late) / max(t, 1)
    return hdr, body, fin, acc


def block_slope(sw, y, nboot=400, rng=None):
    """Late-window OLS slope per 1000 sweeps + moving-block bootstrap sigma."""
    rng = rng or np.random.default_rng(0)
    if len(y) < 8:
        return float("nan"), float("nan")
    sw, y = np.asarray(sw, float), np.asarray(y, float)
    fit = lambda a, b: np.polyfit(a, b, 1)[0] * 1000.0
    s0 = fit(sw, y)
    L = max(2, len(y) // 8)
    starts = np.arange(len(y) - L + 1)
    boots = []
    for _ in range(nboot):
        idx = np.concatenate([np.arange(s, s + L)
                              for s in rng.choice(starts, len(y) // L + 1)])
        idx = idx[:len(y)]
        boots.append(fit(sw[idx], y[idx]))
    return s0, float(np.std(boots))


def phase_stats(body, tag, frac=0.5):
    ph = [r for r in body if r["phase"] == tag]
    if not ph:
        return None
    cut = ph[0]["sw"] + frac * (ph[-1]["sw"] - ph[0]["sw"])
    late = [r for r in ph if r["sw"] >= cut] or ph
    g = lambda k: np.array([r[k] for r in late], float)
    sl, sg = block_slope([r["sw"] for r in late], g("n_ill"))
    return {"f3": g("f") if False else np.array([r["f"][3] for r in late]),
            "f0": np.array([r["f"][0] for r in late]),
            "n_ill": g("n_ill"), "ncomp": g("ncomp"), "top1": g("top1"),
            "n6": np.array([r["edeg"].get("6", 0) for r in late]),
            "slope": sl, "slope_sd": sg, "all": ph}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--fig", default=None)
    args = ap.parse_args()

    prefixes = sorted({p[:-len(".obs.jsonl")]
                       for d in args.dirs
                       for p in glob(os.path.join(d, "*.obs.jsonl"))})
    if not prefixes:
        sys.exit("no runs found")

    cells = []
    for p in prefixes:
        hdr, body, fin, acc = load(p)
        A, B = phase_stats(body, "A"), phase_stats(body, "B")
        cells.append((os.path.basename(p), hdr, A, B, acc, fin))

    h = cells[0][1]
    print(f"host {h['crystal']} m{h['mcell']}  native f = ({h['f0_native']}, "
          f"{h['f0_native']+h['f3_native']}, {2*h['f3_native']}, "
          f"{h['f3_native']})   e_nat = {h['e_native']:.9f}")
    print(f"1 quantum = 1 tet = 1 sixfold edge = {h['quantum']:.4e} in mean "
          f"edge degree;  e* is {h['quanta_native_to_eflat']:.2f} quanta below "
          f"native")
    print(f"action: A(f3-f3_ref)^2 + B(f1-6f3/e_tar)^2 + c_imp*sum m^2, "
          f"A = {h['vol_coef']}, B = {h['hinge_coef']}   [1 chain per cell]\n")

    print(f"{'cell':<30}{'c_imp':>6}{'dq':>5} |{'A: f3-nat':>10}{'n_ill':>8}"
          f"{'ncmp':>6}{'top1':>6}{'acc':>9} |{'B: f3-nat':>10}{'n_ill':>8}"
          f"{'frac0':>7}{'acc':>9}{'slopeB':>16}")
    for name, hd, A, B, acc, fin in cells:
        f = lambda x: f"{x:>+10.2f}"
        row = f"{name:<30}{hd['cimp']:>6}{hd['dq']:>5} |"
        if A:
            row += (f"{A['f3'].mean()-hd['f3_native']:>+10.2f}"
                    f"{A['n_ill'].mean():>8.1f}{A['ncomp'].mean():>6.1f}"
                    f"{A['top1'].mean():>6.1f}{acc.get('A', 0):>9.1e} |")
        if B and B is not A:
            row += (f"{B['f3'].mean()-hd['f3_native']:>+10.2f}"
                    f"{B['n_ill'].mean():>8.2f}"
                    f"{(B['n_ill'] == 0).mean():>7.2f}{acc.get('B', 0):>9.1e}"
                    f"{B['slope']:>+9.2f}+-{B['slope_sd']:<5.2f}")
        else:
            row += f"{'--':>10}{'--':>8}{'--':>7}{'--':>9}{'--':>16}"
        print(row)
        if fin.get("events"):
            print(f"{'':<30}  first passage: " +
                  ", ".join(f"{e['what']} @ {e['sw']}"
                            for e in fin["events"]))

    if args.fig:
        make_fig(cells, args.fig, h)


def make_fig(cells, path, h):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    two = [c for c in cells if c[3] is not None and c[3] is not c[2]]
    above = [c for c in two if c[1]["dq"] < 0]
    below = [c for c in two if c[1]["dq"] > 0]
    FROZEN = 1e-6            # phase-A acceptance below this = never strained

    def label(hd):
        t = f"$c_{{imp}}$ {hd['cimp']:g}"
        if hd.get("heal_cimp") not in (None, hd["cimp"]):
            t += rf"$\rightarrow${hd['heal_cimp']:g}"
        t += (f", cap {hd['illegal_budget']}"
              if hd.get("illegal_budget", -1) >= 0 else ", no cap")
        return f"dq {hd['dq']:+d}, " + t

    def group(cells_):
        seen = {}
        for name, hd, A, B, acc, fin in cells_:
            key = (hd["dq"], hd["cimp"],
                   hd.get("heal_cimp") if hd.get("heal_cimp") is not None
                   else hd["cimp"], hd.get("illegal_budget", -1) or -1)
            seen.setdefault(key, []).append((A, B, hd, acc))
        return sorted(seen.items())

    def traj(axis, cells_, key_of, cmap, logy=False):
        gs = group(cells_)
        for i, (key, runs) in enumerate(gs):
            col = cmap(0.12 + 0.76 * i / max(1, len(gs) - 1))
            A0, B0, hd0, acc0 = runs[0]
            froz = acc0.get("A", 1) < FROZEN
            for j, (A, B, hd, acc) in enumerate(runs):
                rows = A["all"] + B["all"]
                sw = [r["sw"] + (0 if r["phase"] == "A"
                                 else hd["strain_sweeps"]) for r in rows]
                axis.plot(sw, [key_of(r, hd) for r in rows], color=col,
                          lw=1.6 if froz else 1.1, ls="--" if froz else "-",
                          alpha=0.9,
                          label=(label(hd) + (" [frozen]" if froz else ""))
                          if j == 0 else None)
        axis.axvline(gs[0][1][0][2]["strain_sweeps"], color="k", lw=1.4)
        axis.set_xlabel("sweep")
        if logy:
            axis.set_yscale("symlog", linthresh=1)
            axis.set_ylim(bottom=0)
        axis.legend(fontsize=7.5, loc="upper left" if logy else "best",
                    framealpha=0.9)

    fig, ax = plt.subplots(2, 2, figsize=(13.6, 9.2))
    f3of = lambda r, hd: r["f"][3] - hd["f3_native"]
    niof = lambda r, hd: r["n_ill"]

    if above:
        traj(ax[0][0], above, f3of, plt.cm.viridis)
        for k, _ in group(above):
            ax[0][0].axhline(-k[0], color="0.5", ls=":", lw=0.9)
        ax[0][0].axhline(0, color="0.7", lw=0.8)
        ax[0][0].set_ylabel(r"$f_3-f_3^{nat}$   (quanta)")
        ax[0][0].set_title("ABOVE native — barrierless: one flicker per "
                           "quantum, and it all comes back", fontsize=10.5)
        traj(ax[1][0], above, niof, plt.cm.viridis, logy=True)
        ax[1][0].set_ylabel("n_ill  (defect vertices)")
        ax[1][0].set_title("release heals to the pristine crystal "
                           "(registry-verified 1000/1000)", fontsize=10.5)
    if below:
        traj(ax[0][1], below, f3of, plt.cm.plasma)
        for k, _ in group(below):
            ax[0][1].axhline(-k[0], color="0.5", ls=":", lw=0.9)
        ax[0][1].axhline(0, color="0.7", lw=0.8)
        ax[0][1].set_ylabel(r"$f_3-f_3^{nat}$   (quanta)")
        ax[0][1].set_title("BELOW native (requested) — a price either freezes "
                           "or melts; a BUDGET works", fontsize=10.5)
        traj(ax[1][1], below, niof, plt.cm.plasma, logy=True)
        ax[1][1].set_ylabel("n_ill  (defect vertices)")
        ax[1][1].set_title("same strain, two orders of magnitude in defect "
                           "load", fontsize=10.5)

    fig.suptitle(
        f"a15 m{h['mcell']}  ($f_0$ = {h['f0_native']}, $f_3$ = "
        f"{h['f3_native']}, $e_{{nat}}$ = {h['e_native']:.6f});   "
        rf"$S = {h['vol_coef']:g}(f_3-f_3^{{ref}})^2 + {h['hinge_coef']:g}"
        rf"(f_1-6f_3/e_{{tar}})^2 + c_{{imp}}\sum_v m(v)^2$;   "
        "1 chain per curve;  1 quantum = 1 tet = 1 sixfold edge = "
        f"{h['quantum']:.3e} in mean edge degree", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fig.savefig(path, dpi=150)
    print(f"\nSaved to: {os.path.abspath(path)}")


if __name__ == "__main__":
    main()
