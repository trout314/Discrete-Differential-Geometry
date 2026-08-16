#!/usr/bin/env python3
"""Analyze ecmc_ab.py output: blob-dispersal observables vs WORK, arm by arm.

PRIMARY observable is ``spread_mean`` -- the mean graph distance of defect
vertices from the (fixed) blob center. It rises from ~3.1 (packed) toward
the fully-dispersed value (mean graph distance from a fixed vertex, 3.98 on
C15 m3 N3672), and it is the SLOW mode: sum_v m^2 relieves the worst
crowding locally within ~10 sweeps and then saturates, while actual spatial
dispersal is transport-limited. Both are reported; first-passage statistics
use spread (rising thresholds) and m2 (falling thresholds).

WORK is the common currency: plain attempted moves + 2 x flight-kernel
nonlocal_slide_at calls (the driver already doubles them). The D-side
nonlocal channel of arm N is counted by totalTried natively.

Arms: A = plain, B = plain + lifted ECMC flights, N = plain + D-side
diffusive nonlocal slide (the control that separates "ballistic lifting
helps" from "any f-preserving nonlocal move helps").

    python scripts/defect_dynamics/ecmc_ab_analyze.py data/ecmc_ab/main_nfc30
"""
import glob
import json
import os
import sys

import numpy as np

ARM_LABEL = {"A": "A: plain", "B": "B: + lifted ECMC flight",
             "N": "N: + diffusive nonlocal slide"}
ARM_COLOR = {"A": "#4053d3", "B": "#dd3f74", "N": "#00b25d"}


def load_dir(d):
    chains = {}
    for p in sorted(glob.glob(os.path.join(d, "arm?_s??.ab.jsonl"))):
        tag = os.path.basename(p).split(".")[0]
        rows = [json.loads(l) for l in open(p)]
        if len(rows) < 3:
            continue
        work = np.array([r["work_plain"] + r["work_ep"] for r in rows], float)
        chains.setdefault(tag[3], []).append({
            "tag": tag,
            "work": work,
            "sw": np.array([r["sw"] for r in rows], float),
            "m2": np.array([r["m2"] for r in rows], float),
            "spread": np.array([r.get("spread_mean") or np.nan
                                for r in rows], float),
            "nill": np.array([r["n_ill"] for r in rows], float),
            "n3": np.array([r["n3"] for r in rows], float),
            "ep_frac": rows[-1]["work_ep"] / max(1.0, work[-1]),
        })
    mp = os.path.join(d, "manifest.json")
    manifest = json.load(open(mp)) if os.path.exists(mp) else {}
    return chains, manifest


def smooth(y, w=9):
    """CAUSAL (trailing) running mean, NaN-safe.

    Must be trailing, not centered: a centered window leaks later values
    backwards, so a fast-relaxing chain appears to have crossed the
    threshold before it did -- which biases exactly the arms we are trying
    to rank, and hardest for the fastest one.
    """
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 2:
        return y
    yy = np.interp(np.arange(len(y)), np.nonzero(ok)[0], y[ok])
    out = np.empty_like(yy)
    for i in range(len(yy)):
        out[i] = yy[max(0, i - w + 1):i + 1].mean()
    return out


def first_passage(work, series, thr, rising):
    """First work at which the trailing-smoothed series crosses thr.

    Index 0 is excluded: it is the t = 0 state at zero work, shared by every
    arm by construction, so a "crossing" there is the initial condition, not
    a relaxation event (and log(0) would poison the ratio test).
    """
    s = smooth(series)[1:]
    w = np.asarray(work, float)[1:]
    idx = np.nonzero(s >= thr if rising else s <= thr)[0]
    return w[idx[0]] if len(idx) else np.nan


def mean_curve(chains, key, grid):
    vals = []
    for c in chains:
        y, ok = c[key], np.isfinite(c[key])
        if ok.sum() < 2:
            continue
        vals.append(np.interp(grid, c["work"][ok], y[ok]))
    v = np.array(vals)
    return v.mean(0), v.std(0, ddof=1) / np.sqrt(max(1, len(v)))


def endpoints(chains, frac=0.25):
    """Late-window levels + late-window drift per arm.

    THE control for any speedup claim: a biased channel relaxes 'faster'
    only by drifting to a DIFFERENT stationary state, which is
    indistinguishable from acceleration on a first-passage plot. Same
    late-time level across arms + faster approach = real acceleration.
    Different level = a different stationary law.

    Drift is the late-window OLS slope per 1e7 work with a chain-level
    bootstrap sigma (chains are independent, so resampling chains is the
    correct block); |slope| < 2 sigma is consistent with stationarity.
    """
    print(f"\n  -- late-window levels (last {frac:.0%} of work) --")
    print(f"    {'arm':3} {'observable':10} {'mean +- se':>18} "
          f"{'drift/1e7 work +- sigma':>26}  stationary?")
    rng = np.random.default_rng(0)
    for arm in sorted(chains):
        for key in ("spread", "m2", "nill", "n3"):
            per_chain, slopes = [], []
            for c in chains[arm]:
                n0 = int(len(c["work"]) * (1 - frac))
                w, y = c["work"][n0:], c[key][n0:]
                ok = np.isfinite(y)
                if ok.sum() < 3:
                    continue
                per_chain.append(y[ok].mean())
                slopes.append(np.polyfit(w[ok] / 1e7, y[ok], 1)[0])
            if not per_chain:
                continue
            v = np.array(per_chain)
            s = np.array(slopes)
            se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
            bs = np.array([s[rng.integers(0, len(s), len(s))].mean()
                           for _ in range(4000)])
            sig = bs.std(ddof=1)
            flag = "yes" if abs(s.mean()) < 2 * sig else "NO -- still moving"
            print(f"    {arm:3} {key:10} {v.mean():10.2f} +- {se:5.2f} "
                  f"{s.mean():16.3f} +- {sig:6.3f}  {flag}")


def compare(chains, key, thrs, rising, label):
    """First-passage table + log-ratio two-sample test vs arm A."""
    print(f"\n  -- first passage in WORK, {label} --")
    fps = {}
    for arm in sorted(chains):
        fps[arm] = {t: np.array([first_passage(c["work"], c[key], t, rising)
                                 for c in chains[arm]]) for t in thrs}
    for t in thrs:
        line = f"    thr {t:6.2f}: "
        for arm in sorted(chains):
            fp = fps[arm][t]
            ok = np.isfinite(fp)
            line += (f"{arm} {ok.sum()}/{len(fp)}"
                     + (f" {fp[ok].mean():11,.0f}" if ok.sum() else
                        " " * 12) + "   ")
        print(line)
        if "A" not in fps:
            continue
        fa = fps["A"][t]
        for arm in sorted(chains):
            if arm == "A":
                continue
            fb = fps[arm][t]
            ok = np.isfinite(fa) & np.isfinite(fb)
            # paired: same seed index => same initial blob
            if ok.sum() < 2:
                nb, na = np.isfinite(fb), np.isfinite(fa)
                if nb.sum() < na.sum():
                    print(f"        {arm} vs A: {arm} reached it in "
                          f"{nb.sum()}/{len(fb)} chains vs A "
                          f"{na.sum()}/{len(fa)} -- no paired test")
                continue
            d = np.log(fa[ok]) - np.log(fb[ok])
            se = d.std(ddof=1) / np.sqrt(ok.sum())
            z = d.mean() / se if se > 0 else np.nan
            print(f"        {arm}/A work ratio {np.exp(-d.mean()):.3f}"
                  f"   speedup {np.exp(d.mean()):.3f}x"
                  f"   (paired z = {z:+.2f}, n = {ok.sum()};"
                  f" +z => {arm} faster)")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: ecmc_ab_analyze.py <run-dir>")
    d = sys.argv[1].rstrip("/")
    chains, manifest = load_dir(d)
    cfg = manifest.get("cfg", {})
    if not chains:
        sys.exit(f"no chains in {d}")

    print(f"=== {d} ===")
    print(f"  C15 m3 N3672, k={cfg.get('k_fliers')} fliers in a radius-"
          f"{cfg.get('radius')} ball; action = {cfg.get('nfc')}*(f3-f3*)^2 "
          f"+ {cfg.get('k1')}*(f1-6f3/e*)^2 + {cfg.get('cimp')}*sum m^2")
    print(f"  {cfg.get('sweeps')} sweeps, chunk {cfg.get('chunk')}; "
          f"arm B: {cfg.get('eps_per_chunk')} x {cfg.get('ep_steps')} kernel "
          f"steps/chunk, kscan {cfg.get('kscan')}, beta {cfg.get('beta')}")
    for arm in sorted(chains):
        cs = chains[arm]
        wf = np.mean([c["ep_frac"] for c in cs])
        w = np.mean([c["work"][-1] for c in cs])
        s0 = np.nanmean([c["spread"][0] for c in cs])
        s1 = np.nanmean([c["spread"][-1] for c in cs])
        m0 = np.mean([c["m2"][0] for c in cs])
        m1 = np.mean([c["m2"][-1] for c in cs])
        print(f"  arm {arm} ({ARM_LABEL[arm]}): {len(cs)} chains, "
              f"total work {w:,.0f}, flight-kernel fraction {wf:.1%}; "
              f"spread {s0:.2f} -> {s1:.2f}, m2 {m0:.0f} -> {m1:.0f}")
    if "B" in chains:
        neutral = 1.0 / (1.0 - np.mean([c["ep_frac"] for c in chains["B"]]))
        print(f"  NEUTRAL prediction if flights do nothing: "
              f"B/A total-work ratio = {neutral:.3f}")

    endpoints(chains)
    compare(chains, "spread", [3.3, 3.5, 3.7, 3.8], True,
            "spread_mean (rising; dispersed ~ 3.98)")
    compare(chains, "m2", [200, 160, 140], False,
            "sum m^2 (falling)")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    wmax = min(c["work"][-1] for cs in chains.values() for c in cs)
    grid = np.geomspace(max(1e4, wmax * 2e-4), wmax, 200)
    for key, a, ylab in ((("spread"), ax[0],
                          "mean graph distance of defect\n"
                          "vertices from blob center"),
                         (("m2"), ax[1], r"$\sum_v m(v)^2$")):
        for arm in sorted(chains):
            mu, se = mean_curve(chains[arm], key, grid)
            a.plot(grid, mu, color=ARM_COLOR[arm],
                   label=f"{ARM_LABEL[arm]} ({len(chains[arm])})")
            a.fill_between(grid, mu - se, mu + se, color=ARM_COLOR[arm],
                           alpha=0.2, lw=0)
        a.set_xscale("log")
        a.set_xlabel("work (attempted moves, all channels)")
        a.set_ylabel(ylab)
    ax[0].axhline(3.98, color="gray", ls=":", lw=1)
    ax[0].annotate("fully dispersed (3.98)", (grid[0], 3.98), fontsize=7,
                   color="gray", va="bottom")
    ax[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(
        f"C15 m3 N3672 blob dispersal -- {cfg.get('nfc')}(f3-f3*)^2 + "
        f"{cfg.get('k1')}(f1-6f3/e*)^2 + {cfg.get('cimp')}sum m^2, "
        f"k={cfg.get('k_fliers')}, {len(chains.get('A', []))} chains/arm "
        "-- PROVISIONAL (uncertified)", fontsize=9)
    fig.tight_layout()
    out = os.path.join(d, "ecmc_ab_dispersal.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
