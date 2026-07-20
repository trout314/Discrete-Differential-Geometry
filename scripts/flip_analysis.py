#!/usr/bin/env python3
"""Dynamics of the disclination network from the six-edge flip stream.

Consumes ``<base>_flips.bin`` files (dtype disclination.SIX_FLIP_DTYPE: one
(clock, u, v, dir) record per edge crossing degree 5<->6 in an accepted
move) and reports the dynamical structure of the web:

  * gross vs net flux: total birth+death events per sweep vs the net drift
    of E6 -- how much rewiring underlies the slow observable motion;
  * turnover: six-edge LIFETIME distribution (birth -> death on the same
    edge) and the web turnover time  <n_six> / death rate;
  * blinking: dead-time (death -> rebirth) distribution -- fast recurrence
    means local double-well flickering, not diffusing lines;
  * dynamic heterogeneity: flips-per-edge distribution and the share of all
    activity carried by the top 1% most active edges (glass physics: is
    motion localized or spread?).

Clock (attempted moves) is converted to sweeps via the sibling CSV's sweep
span. Vertex-label REUSE after 1->4/4->1 churn can alias distinct physical
edges under one key; the reported alternation-violation fraction (same-edge
consecutive events with equal dir, which is impossible for one physical
edge) bounds that contamination.

Usage:
    python scripts/flip_analysis.py data/replicas/*_flips.bin \
        --out out/flip_stats.csv --plot out/flip_dynamics.png
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry.disclination import SIX_FLIP_DTYPE

DT_BINS = np.logspace(-5, 4.5, 96)     # sweeps
FPE_BINS = np.logspace(0, 5, 51)       # flips per edge


def csv_span(path):
    """(segment_start, last_sweep) of the sibling dope_hold CSV."""
    meta, header, last = {}, None, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        meta[k] = v
            elif header is None:
                header = line
            elif line:
                last = line
    return float(meta.get("segment_start", 0)), float(last.split(",")[0])


def analyze(path, sweeps_span, time_bin=500.0):
    ev = np.fromfile(path, dtype=SIX_FLIP_DTYPE)
    clock = ev["clock"].astype(np.int64)
    key = (ev["u"].astype(np.int64) << 32) | ev["v"].astype(np.int64)
    d = ev["dir"].astype(np.int64)
    aps = clock[-1] / max(sweeps_span, 1.0)          # attempts per sweep

    order = np.lexsort((clock, key))
    k2, c2, d2 = key[order], clock[order], d[order]
    same = k2[1:] == k2[:-1]
    dt = (c2[1:] - c2[:-1]) / aps                    # sweeps
    viol = float(np.mean(d2[1:][same] == d2[:-1][same])) if same.any() else 0.0
    lifetimes = dt[same & (d2[:-1] == 1) & (d2[1:] == -1)]
    deadtimes = dt[same & (d2[:-1] == -1) & (d2[1:] == 1)]

    uk, ck = np.unique(key, return_counts=True)
    ck_desc = np.sort(ck)[::-1]
    top1 = float(ck_desc[:max(1, len(ck) // 100)].sum() / ck.sum())

    births, deaths = int((d == 1).sum()), int((d == -1).sum())
    res = dict(
        records=len(ev), births=births, deaths=deaths,
        net=births - deaths, sweeps=sweeps_span, attempts_per_sweep=aps,
        gross_per_sweep=len(ev) / sweeps_span,
        edges_touched=len(uk), max_flips_one_edge=int(ck.max()),
        mean_flips_per_edge=float(ck.mean()), top1pct_share=top1,
        alt_violation=viol, n_lifetimes=int(len(lifetimes)),
        life_median=float(np.median(lifetimes)) if len(lifetimes) else 0.0,
        life_mean=float(lifetimes.mean()) if len(lifetimes) else 0.0,
        dead_median=float(np.median(deadtimes)) if len(deadtimes) else 0.0)
    hists = dict(
        life=np.histogram(lifetimes, DT_BINS)[0],
        dead=np.histogram(deadtimes, DT_BINS)[0],
        fpe=np.histogram(ck, FPE_BINS)[0])
    # gross activity per time bin (sweeps)
    tb = np.arange(0, sweeps_span + time_bin, time_bin)
    hists["rate_t"] = np.histogram(clock / aps, tb)[0] / time_bin
    hists["rate_edges"] = tb
    return res, hists


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("flips", nargs="+", help="*_flips.bin files")
    ap.add_argument("--group", default=r"_r\d+.*$",
                    help="regex stripped from basename to pool replicas")
    ap.add_argument("--out", default=None, help="per-file scalar CSV")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    rows, pooled = [], defaultdict(lambda: defaultdict(lambda: 0))
    for path in args.flips:
        base = re.sub(r"_flips\.bin$", "", path)
        if not os.path.exists(base + ".csv"):
            print(f"  SKIP {os.path.basename(path)}: no sibling CSV")
            continue
        s0, s1 = csv_span(base + ".csv")
        label = re.sub(args.group, "",
                       os.path.basename(re.sub(r"_flips\.bin$", "", path)))
        res, hists = analyze(path, s1 - s0)
        res = {"file": os.path.basename(path), "label": label, **res}
        rows.append(res)
        for hk in ("life", "dead", "fpe"):
            cur = pooled[label][hk]
            pooled[label][hk] = (cur + hists[hk] if not np.isscalar(cur)
                                 else hists[hk].astype(np.float64).copy())
        print(f"  {res['file']}: {res['records']} ev over "
              f"{res['sweeps']:.0f} sw | gross {res['gross_per_sweep']:.1f}/sw"
              f" net {res['net']:+d} | median life "
              f"{res['life_median']:.2f} sw | top1% share "
              f"{res['top1pct_share']:.2f} | viol {res['alt_violation']:.4f}")

    print()
    for label in sorted({r["label"] for r in rows}):
        sel = [r for r in rows if r["label"] == label]
        g = np.mean([r["gross_per_sweep"] for r in sel])
        net = np.mean([abs(r["net"]) / r["sweeps"] for r in sel])
        lm = np.mean([r["life_median"] for r in sel])
        dm = np.mean([r["dead_median"] for r in sel])
        t1 = np.mean([r["top1pct_share"] for r in sel])
        print(f"{label} ({len(sel)} replicas): gross {g:.1f}/sweep, "
              f"|net| {net:.3f}/sweep (ratio {g/max(net,1e-12):.0f}:1), "
              f"median life {lm:.2f} sw, median dead {dm:.3f} sw, "
              f"top-1% edges carry {t1:.0%}")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        mid_dt = np.sqrt(DT_BINS[1:] * DT_BINS[:-1])
        wid_dt = np.diff(DT_BINS)
        mid_f = np.sqrt(FPE_BINS[1:] * FPE_BINS[:-1])
        wid_f = np.diff(FPE_BINS)
        for label in sorted(pooled):
            h = pooled[label]
            m = h["life"] > 0
            axes[0, 0].loglog(mid_dt[m], (h["life"] / wid_dt)[m], ".-",
                              ms=3, lw=0.8, label=label)
            m = h["dead"] > 0
            axes[0, 1].loglog(mid_dt[m], (h["dead"] / wid_dt)[m], ".-",
                              ms=3, lw=0.8, label=label)
            m = h["fpe"] > 0
            axes[1, 0].loglog(mid_f[m], (h["fpe"] / wid_f)[m], ".-",
                              ms=3, lw=0.8, label=label)
        axes[0, 0].set_title("six-edge lifetime (birth→death)", fontsize=10)
        axes[0, 0].set_xlabel("lifetime (sweeps)")
        axes[0, 1].set_title("dead time (death→rebirth)", fontsize=10)
        axes[0, 1].set_xlabel("dead time (sweeps)")
        axes[1, 0].set_title("flips per edge", fontsize=10)
        axes[1, 0].set_xlabel("total flips of one edge")
        for r in rows:
            if r["file"].endswith("r0_flips.bin"):
                pass
        for label in sorted(pooled):
            first = next(r for r in rows if r["label"] == label)
        for label in sorted(pooled):
            sel = [r for r in rows if r["label"] == label]
            axes[1, 1].bar([f"{label}\ngross" ], [np.mean(
                [r["gross_per_sweep"] for r in sel])], alpha=0.7)
            axes[1, 1].bar([f"{label}\n|net|"], [np.mean(
                [abs(r['net']) / r['sweeps'] for r in sel])], alpha=0.7)
        axes[1, 1].set_yscale("log")
        axes[1, 1].set_title("gross vs |net| flux (per sweep)", fontsize=10)
        axes[1, 1].tick_params(axis="x", labelsize=7)
        for ax in axes.flat[:3]:
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        fig.suptitle("disclination-web rewiring dynamics (pooled replicas)")
        fig.tight_layout()
        os.makedirs(os.path.dirname(args.plot) or ".", exist_ok=True)
        fig.savefig(args.plot, dpi=120)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
