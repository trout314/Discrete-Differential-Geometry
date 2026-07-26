#!/usr/bin/env python3
"""Aggregate the FP production run (fp_encounter / fp_transport JSONs)
into tables + a four-panel figure.

Panels: (1) encounter time vs initial separation (dock times + censored
marks, medians via Kaplan-Meier when censoring < 50%); (2) dock census
(on/off-chain, additive vs contact docks by |V_dock|); (3) pooled MSD
vs lag with power-law fit and D_slide from the linear window, lam 0.40
vs 1.0; (4) recombination: outcome fractions + t_unbind / t_return
distributions. All kinetics are SLIDE CHANNEL, FROZEN BACKGROUND
(v1 FROZEN); time unit = attempted moves (1 sweep = N3 attempts).
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")


def load_eps(pattern):
    eps = []
    for f in sorted(glob.glob(os.path.join(DATA, pattern))):
        d = json.load(open(f))
        for e in d["episodes"]:
            e["_sep"] = d["sep_sites"]
            e["_maxfl"] = d.get("max_flights", 60)
        eps += d["episodes"]
    return eps


def km_median(times, censored):
    """Kaplan-Meier median with right censoring (times, censor times)."""
    pts = sorted([(t, 1) for t in times] + [(t, 0) for t in censored])
    n = len(pts)
    surv, at_risk = 1.0, n
    for t, ev in pts:
        if ev:
            surv *= (at_risk - 1) / at_risk
        at_risk -= 1
        if surv <= 0.5:
            return t
    return None       # median not reached


def msd_curve(trajs, nbins=24):
    lags, sq = [], []
    for tr in trajs:
        t = np.array([p["t"] for p in tr])
        x = np.array([p["x"] for p in tr])
        n = len(t)
        for i in range(n - 1):
            dt = t[i + 1:] - t[i]
            d2 = ((x[i + 1:] - x[i]) ** 2).sum(1)
            lags.append(dt)
            sq.append(d2)
    lags = np.concatenate(lags).astype(float)
    sq = np.concatenate(sq)
    lo, hi = np.log10(max(lags.min(), 1)), np.log10(lags.max())
    edges = np.logspace(lo, hi, nbins + 1)
    mid, msd, cnt = [], [], []
    for k in range(nbins):
        m = (lags >= edges[k]) & (lags < edges[k + 1])
        if m.sum() < 20:
            continue
        mid.append(np.sqrt(edges[k] * edges[k + 1]))
        msd.append(sq[m].mean())
        cnt.append(int(m.sum()))
    return np.array(mid), np.array(msd), np.array(cnt)


def main():
    # ---------------- encounter ----------------
    eps = load_eps("prodA_s*_p*.json")
    seps = sorted({e["_sep"] for e in eps})
    print("== Phase A: encounter kinetics (attempted-move units) ==")
    A = {}
    for s0 in seps:
        es = [e for e in eps if e["_sep"] == s0]
        t = [e["t"] for e in es if e["reason"] == "dock"]
        cen = [e["t"] for e in es if e["reason"] == "censored"]
        onc = sum(1 for e in es if e.get("on_chain_site") is not None)
        vd = [e["V_dock"] for e in es if e.get("V_dock") is not None]
        add = sum(1 for v in vd if abs(v) < 1e-6)
        contact = sum(1 for v in vd if abs(v) >= 1e-6)
        med = km_median(t, cen)
        A[s0] = {"t": t, "cen": cen, "med": med, "onc": onc,
                 "add": add, "contact": contact}
        ms = f"{med:.3e}" if med else f">{np.median(cen):.1e} (censored)"
        print(f"  s0={s0:2d}: {len(t):3d} dock / {len(cen):3d} censored; "
              f"KM median {ms}; on-chain {onc}; docks additive/contact "
              f"{add}/{contact}")

    # ---------------- transport ----------------
    print("\n== Phase C: transport (slide channel, frozen background) ==")
    curves = {}
    for lam, pat in ((0.40, "prodC_p*.json"), (1.0, "prodClam1_p*.json")):
        trajs = []
        for f in sorted(glob.glob(os.path.join(DATA, pat))):
            d = json.load(open(f))
            if abs(d["lam"] - lam) < 1e-9 and len(d["traj"]) > 50:
                trajs.append(d["traj"])
        if not trajs:
            continue
        mid, msd, cnt = msd_curve(trajs)
        m = (mid > 1e5) & (msd > 0)
        alpha, lnA = np.polyfit(np.log(mid[m]), np.log(msd[m]), 1)
        D = float(np.exp(lnA)) / 6.0 if abs(alpha - 1) < 0.25 else None
        # D from the largest well-populated decade regardless, for scale
        Dt = msd[m][-1] / (6.0 * mid[m][-1])
        curves[lam] = (mid, msd, alpha, Dt, len(trajs))
        print(f"  lam={lam:.2f}: {len(trajs)} trajectories; MSD exponent "
              f"{alpha:.3f}; MSD/6t at t={mid[m][-1]:.1e} -> "
              f"D = {Dt:.3e} cells^2/attempt "
              f"({Dt * 7249:.3e} cells^2/sweep)")

    # ---------------- recombination ----------------
    print("\n== Phase B: contact-pair recombination ==")
    eps = load_eps("prodB_p*.json") + load_eps("prodB2_p*.json")
    for s0 in sorted({e["_sep"] for e in eps}):
        es = [e for e in eps if e["_sep"] == s0]
        fr = [e for e in es if e.get("freed")]
        o = {}
        for e in fr:
            o[e.get("outcome", "?")] = o.get(e.get("outcome", "?"), 0) + 1
        nr = o.get("recombine", 0)
        print(f"  start sep {s0}: {len(es)} episodes, {len(fr)} freed; "
              f"{o}; P(rec|freed) = {nr / max(1, len(fr)):.2f}")
    freed = [e for e in eps if e.get("freed")]
    out = {}
    for e in freed:
        out[e.get("outcome", "?")] = out.get(e.get("outcome", "?"), 0) + 1
    tub = np.array([e["t_unbind"] for e in freed])
    tre = np.array([e["t_return"] for e in freed
                    if e.get("outcome") == "recombine"])
    nrec = out.get("recombine", 0)
    if len(tub):
        print(f"  pooled: t_unbind median {np.median(tub):.3e}; "
              f"t_return median {np.median(tre):.3e}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax = axes[0, 0]
    for s0 in seps:
        t = A[s0]["t"]
        ax.scatter([s0] * len(t), t, s=14, alpha=0.5, color="tab:blue")
        ax.scatter([s0] * len(A[s0]["cen"]), A[s0]["cen"], s=18, alpha=0.6,
                   marker="^", color="tab:red")
    med_s = [s for s in seps if A[s]["med"]]
    ax.plot(med_s, [A[s]["med"] for s in med_s], "k-o", lw=2,
            label="KM median")
    ax.set_yscale("log")
    ax.set_xlabel("initial separation s0 (chain sites)")
    ax.set_ylabel("t_dock (attempted moves)")
    ax.set_title("encounter time vs separation (^ = censored)")
    ax.legend()

    ax = axes[0, 1]
    labels = [f"s0={s}" for s in seps]
    onc = [A[s]["onc"] for s in seps]
    total = [len(A[s]["t"]) for s in seps]
    contact = [A[s]["contact"] for s in seps]
    x = np.arange(len(seps))
    ax.bar(x - 0.2, total, 0.4, label="docks")
    ax.bar(x + 0.2, onc, 0.2, label="on-chain")
    ax.bar(x + 0.4, contact, 0.2, label="contact (V>0)")
    ax.set_xticks(x, labels)
    ax.set_ylabel("count")
    ax.set_title("dock census")
    ax.legend()

    ax = axes[1, 0]
    for lam, (mid, msd, alpha, Dt, ntr) in curves.items():
        ax.loglog(mid, msd, "o-", label=f"lam={lam:.2f} (a={alpha:.2f}, "
                  f"{ntr} traj)")
    if 0.40 in curves:
        mid = curves[0.40][0]
        ax.loglog(mid, mid * curves[0.40][3] * 6, "k--", alpha=0.5,
                  label="slope 1")
    ax.set_xlabel("lag (attempted moves)")
    ax.set_ylabel("MSD (cells^2)")
    ax.set_title("lone-knot MSD, slide channel, frozen background")
    ax.legend()

    ax = axes[1, 1]
    if len(tub):
        bins = np.logspace(np.log10(max(tub.min(), 1)),
                           np.log10(max(tub.max(), tre.max() if len(tre)
                                        else 1) + 1), 18)
        ax.hist(tub, bins=bins, alpha=0.6, label="t_unbind")
        if len(tre):
            ax.hist(tre, bins=bins, alpha=0.6, label="t_return (recombine)")
        ax.set_xscale("log")
    ax.set_xlabel("time (attempted moves)")
    ax.set_ylabel("count")
    ax.set_title(f"recombination: P(rec|freed) = "
                 f"{nrec}/{len(freed)}")
    ax.legend()

    fig.suptitle(
        "FP production (M3, v1 FROZEN): slide-channel kinetics -- "
        "host R m2 N7248 pristine + knots; "
        "S = 0.1(N3-N*)^2 + lam[(e*/6)SUM(deg-e*)^2 + 0.6 U(n6) + m^2], "
        "e*=5.105025; time = attempted moves", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(DATA, "fp_production_report.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
