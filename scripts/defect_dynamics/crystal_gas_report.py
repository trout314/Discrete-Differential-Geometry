#!/usr/bin/env python3
"""Verdict table for the crystal-library x c_imp defect-gas sweep.

Reads every `*.obs.jsonl` written by crystal_gas_scan.py and answers, per cell:
is the post-burn state a STABLE DILUTE DISPERSED gas?

The three words are tested separately, because they fail separately:

  DISPERSED  frac_top1 = (largest complex)/(all defect vertices). 1.0 = the
             single percolating clump that is the melt; ~0 = a gas. The
             companion cut is top1 itself in vertices -- frac_top1 can look
             small merely because n_ill is huge.
  DILUTE     n_ill / f0, the defect vertex fraction. A dispersed state with
             30% of vertices defective is a dense fluid, not a gas.
  STABLE     late-window OLS slope of n_ill vs sweep, in units of its
             BLOCK-BOOTSTRAP sigma (CONVENTIONS sec 3: never full-series OLS
             alone). |z| < 3 = no detectable drift over the window. This is a
             one-sided test -- it cannot see kinetic arrest, so `turnover` is
             reported alongside: a glass is stationary AND frozen.

`turnover` = 1 - Jaccard(defect vertex set, one chunk earlier). ~1 means the
defect population is completely re-drawn every chunk (a gas in the strong
sense); ~0 with a stationary n_ill is the glassy-arrest signature that no
single-start stationarity test can distinguish from equilibrium
(notes/memory/ecmc-blob-ab.md).

Usage:
    python scripts/defect_dynamics/crystal_gas_report.py [--dir data/crystal_gas]
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from crystal_gas_scan import E_FLAT, LIBRARY, pin_gap

# Verdict cuts. Deliberately loose on DILUTE (the pin FORCES a defect debt
# that scales with the host's e_native mismatch, so a15 must carry ~5x c15's
# population and still counts as a gas if it stays dispersed and mobile).
MAX_FRAC_TOP1 = 0.25     # largest complex may hold at most a quarter of defects
MAX_TOP1 = 12            # ...and at most this many vertices outright
MAX_DEFECT_FRAC = 0.30   # defect vertices / f0
MAX_M = 4                # deep burial
MIN_TURNOVER = 0.5       # per chunk
MAX_DRIFT_Z = 3.0
MIN_NCOMP = 3            # a "gas" needs several objects, not one lone defect
MIN_PIN_PAID = 0.25      # below this the state simply ignores the flat target


def pin_paid(f1gap_now, gap0):
    """Fraction of the crystal's native pin debt actually discharged.

    0 = the state sits at the pristine crystal's mismatch and is eating the
    full (f1 - 6f3/e*)^2 penalty rather than paying it in defects; 1 = the
    flat target is exactly met. This is the axis that separates "no gas
    because it melted" from "no gas because defects cost more than the pin".
    """
    return float("nan") if not gap0 else 1.0 - abs(f1gap_now) / abs(gap0)


def block_bootstrap_slope(x, y, nblock=8, ndraw=2000, rng=None):
    """(slope, sigma) with sigma from a moving-block bootstrap over the series.

    Blocks preserve the short-range autocorrelation that makes a naive OLS
    standard error far too small on an MCMC time series."""
    rng = rng or np.random.default_rng(0)
    n = len(x)
    if n < 4:
        return float("nan"), float("nan")
    slope = np.polyfit(x, y, 1)[0]
    blen = max(2, n // nblock)
    starts = np.arange(0, n - blen + 1)
    if len(starts) == 0:
        return slope, float("nan")
    nb = int(math.ceil(n / blen))
    out = np.empty(ndraw)
    for i in range(ndraw):
        idx = np.concatenate([np.arange(s, s + blen)
                              for s in rng.choice(starts, nb)])[:n]
        out[i] = np.polyfit(x, y[idx], 1)[0]
    return float(slope), float(out.std(ddof=1))


def load(stem):
    rows, head = [], None
    with open(stem) as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "header":
                head = d
            elif "kind" not in d:
                rows.append(d)
    return head, rows


def analyse(stem):
    head, rows = load(stem)
    late = [r for r in rows if r.get("burned")]
    if not late:
        return None
    f0 = late[-1]["f"][0]
    sw = np.array([r["sw"] for r in late], float)
    n_ill = np.array([r["n_ill"] for r in late], float)
    slope, sig = block_bootstrap_slope(sw, n_ill)
    # per-1000-sweeps drift, and its significance
    z = abs(slope) / sig if sig and sig > 0 else float("inf")
    top1 = np.array([r["top1"] for r in late], float)
    frac = np.array([r["frac_top1"] for r in late], float)
    turn = np.array([r["turnover"] for r in late if r["turnover"] is not None],
                    float)
    e_mean = np.array([r["e_mean"] for r in late], float)
    f1gap = np.array([r["f"][1] - 6.0 * r["f"][3] / E_FLAT for r in late])
    o = dict(
        crystal=head["crystal"], cimp=head["cimp"], f0=f0,
        e_native=head["e_native"], pin_gap=head["pin_gap_f1"],
        n_forced=head["n_forced_moves"],
        n_ill=n_ill.mean(), n_ill_sd=n_ill.std(),
        dfrac=n_ill.mean() / f0,
        ncomp=np.mean([r["ncomp"] for r in late]),
        top1=top1.mean(), top1_max=top1.max(),
        frac_top1=frac.mean(),
        max_m=max(r["max_m"] for r in late),
        sum_m2=np.mean([r["sum_m2"] for r in late]),
        turnover=turn.mean() if len(turn) else float("nan"),
        e_mean=e_mean.mean(), f1gap=f1gap.mean(),
        paid=pin_paid(f1gap.mean(), head["pin_gap_f1"]),
        slope_k=slope * 1000.0, sigma_k=sig * 1000.0, drift_z=z,
        nsamp=len(late))
    fails = []
    # NOT a gas criterion: `paid` answers a different question (does the state
    # discharge the pin debt?) and is reported as an annotation. A crystal can
    # perfectly well hold a stable dilute gas while still eating most of the
    # pin penalty -- that is the physics, not a failed run.
    #
    # Dispersal is judged on top1 in ABSOLUTE vertices, not on frac_top1:
    # frac_top1 -> 1 whenever the gas is dilute enough to hold only a couple
    # of complexes, which flags the most dilute gases as the most clumped.
    # frac_top1 stays in the table as a percolation readout.
    if o["top1"] > MAX_TOP1:
        fails.append("clumped")
    if o["dfrac"] > MAX_DEFECT_FRAC:
        fails.append("dense")
    if o["max_m"] > MAX_M:
        fails.append("buried")
    if o["ncomp"] < MIN_NCOMP:
        fails.append("sparse")
    if not (o["turnover"] >= MIN_TURNOVER):
        fails.append("frozen")
    if o["drift_z"] > MAX_DRIFT_Z:
        fails.append("drifting")
    o["fails"] = fails
    o["verdict"] = "GAS" if not fails else "/".join(fails)
    # Annotation, not a verdict. Meaningless where the native gap is already
    # ~0 (r: gap = +0.12 f1, so any thermal wobble divides by nothing).
    o["pin_note"] = ("n/a" if abs(o["pin_gap"]) < 1.0
                     else ("unpaid" if o["paid"] < MIN_PIN_PAID else "paid"))
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=os.path.join(_ROOT, "data", "crystal_gas"))
    ap.add_argument("--json", default=None, help="also dump rows as JSON")
    args = ap.parse_args()

    res = [r for r in (analyse(p) for p in
                       sorted(glob.glob(os.path.join(args.dir, "*.obs.jsonl"))))
           if r]
    order = ["a15", "sigma", "z", "p", "delta", "r", "mu", "c15", "c14", "c36"]
    res.sort(key=lambda r: (order.index(r["crystal"])
                            if r["crystal"] in order else 99, r["cimp"]))

    print(f"e_flat = {E_FLAT:.7f}   action = 0.1(f3-f3ref)^2 + "
          f"1.0(f1-6f3/e*)^2 + c_imp*sum m^2   (no EDQ, no U(n6), no VDV/HDV)\n")
    hdr = (f"{'crystal':7s} {'c_imp':>5s} | {'n_ill':>7s} {'dfrac':>6s} "
           f"{'ncomp':>6s} {'top1':>5s} {'ftop1':>6s} {'maxm':>4s} "
           f"{'turn':>5s} | {'drift/kSw':>10s} {'z':>5s} | "
           f"{'<edeg>':>8s} {'paid':>6s} {'pin':6s} | verdict")
    print(hdr)
    print("-" * len(hdr))
    cur = None
    for r in res:
        if cur != r["crystal"]:
            if cur is not None:
                print()
            cur = r["crystal"]
            print(f"# {r['crystal']}  f0={r['f0']}  e_nat={r['e_native']:.7f}"
                  f"  (e*-e_nat)={E_FLAT - r['e_native']:+.6f}"
                  f"  pin gap={r['pin_gap']:+.2f} f1"
                  f"  -> {r['n_forced']:.0f} forced moves")
        print(f"{r['crystal']:7s} {r['cimp']:5.2f} | {r['n_ill']:7.1f} "
              f"{r['dfrac']:6.3f} {r['ncomp']:6.1f} {r['top1']:5.1f} "
              f"{r['frac_top1']:6.3f} {r['max_m']:4d} {r['turnover']:5.2f} | "
              f"{r['slope_k']:+7.2f}+-{r['sigma_k']:<5.2f}{r['drift_z']:5.1f} | "
              f"{r['e_mean']:8.6f} {r['paid']:6.2f} {r['pin_note']:6s} | "
              f"{r['verdict']}")

    print("\n=== recommended c_imp* per crystal (smallest c_imp with GAS) ===")
    best = {}
    for r in res:
        if r["verdict"] == "GAS" and r["crystal"] not in best:
            best[r["crystal"]] = r
    for c in order:
        r = best.get(c)
        if r:
            print(f"  {c:7s} c_imp* = {r['cimp']:.2f}   n_ill={r['n_ill']:.0f} "
                  f"({100*r['dfrac']:.1f}% of vertices) in {r['ncomp']:.0f} "
                  f"complexes, top1={r['top1']:.1f}, turnover={r['turnover']:.2f}")
        elif any(x["crystal"] == c for x in res):
            f = sorted({w for x in res if x["crystal"] == c for w in x["fails"]})
            print(f"  {c:7s} NO GAS at any scanned c_imp  (failure modes: "
                  f"{', '.join(f)})")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=1)
        print(f"\nrows -> {args.json}")


if __name__ == "__main__":
    main()
