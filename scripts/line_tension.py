#!/usr/bin/env python3
"""Line (string) tension from disclination segment-length statistics.

For a living-polymer network in equilibrium the open-segment length
distribution is exponential,  n(l) ~ exp(-f1 * l),  with f1 the segment free
energy per unit length (edge) at kT = 1 -- the "string tension" of the
disclination lines. This script pools the D-side census ``seg_len_hist``
over a set of .mfd snapshots (replicas and/or checkpoint times of the same
state point) and fits  log n(l) = a - f1 * l  by Poisson-weighted least
squares, with a tail-only refit (--lmin) to drop the lattice-stiff l=1 bin.

Snapshots are grouped into state points by regex: --group strips the
replica/time part of the filename so e.g. flatpin_mu3_r3_final.mfd and
flatpin_mu3_r7_final.mfd pool together.

Usage:
    python scripts/line_tension.py data/replicas/*_final.mfd \
        --group '_r\\d+.*$' --out out/line_tension.csv --plot out/line_tension.png
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg


def fit_exp(hist, lmin=1, lmax=None, min_count=5):
    """Weighted LSQ fit of log n(l) = a - f1*l. Returns (f1, err, n_bins,
    (l, n) arrays used) or None if fewer than 3 usable bins."""
    l = np.arange(len(hist))
    m = (l >= lmin) & (hist >= min_count)
    if lmax is not None:
        m &= l <= lmax
    l, n = l[m], hist[m].astype(float)
    if len(l) < 3:
        return None
    w = n                       # var(log n) ~ 1/n for Poisson counts
    A = np.stack([np.ones_like(l, dtype=float), -l.astype(float)], axis=1)
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    coef = cov @ A.T @ W @ np.log(n)
    # chi^2-scaled parameter error
    resid = np.log(n) - A @ coef
    chi2 = float(resid @ W @ resid) / max(1, len(l) - 2)
    return coef[1], float(np.sqrt(cov[1, 1] * max(chi2, 1.0))), len(l), (l, n)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+", help=".mfd files to pool")
    ap.add_argument("--group", default=r"_r\d+.*$",
                    help="regex stripped from the basename to form the "
                         "state-point label (default pools over replicas)")
    ap.add_argument("--lmin", type=int, default=2,
                    help="first segment length of the tail fit (default 2; "
                         "the full fit from l=1 is always also reported)")
    ap.add_argument("--min-count", type=int, default=5)
    ap.add_argument("--out", default=None, help="results CSV")
    ap.add_argument("--plot", default=None, help="histogram + fit PNG")
    args = ap.parse_args()

    groups = defaultdict(lambda: np.zeros(64, dtype=np.int64))
    nsnap = defaultdict(int)
    for path in args.snapshots:
        label = re.sub(args.group, "", os.path.splitext(os.path.basename(path))[0])
        v = ddg.Manifold.load(path, 3)
        c = v.disclination_census()
        groups[label] += c["seg_len_hist"]
        nsnap[label] += 1
        print(f"  {os.path.basename(path)} -> {label}: "
              f"{c['n_segments']} segs, mean {c['mean_seg_len']:.2f}, "
              f"loops {c['n_pure_loops']}")

    results = []
    for label, hist in sorted(groups.items()):
        row = {"label": label, "n_snapshots": nsnap[label],
               "n_segments": int(hist.sum()),
               "mean_seg": (np.arange(64) * hist).sum() / max(1, hist.sum())}
        full = fit_exp(hist, lmin=1, min_count=args.min_count)
        tail = fit_exp(hist, lmin=args.lmin, min_count=args.min_count)
        for name, ft in (("f1_full", full), ("f1_tail", tail)):
            row[name] = ft[0] if ft else float("nan")
            row[name + "_err"] = ft[1] if ft else float("nan")
            row[name + "_bins"] = ft[2] if ft else 0
        results.append((row, hist, tail or full))
        print(f"{label}: {row['n_segments']} segments over "
              f"{nsnap[label]} snapshots | f1(full) = "
              f"{row['f1_full']:.3f} ± {row['f1_full_err']:.3f} | "
              f"f1(l>={args.lmin}) = {row['f1_tail']:.3f} ± "
              f"{row['f1_tail_err']:.3f}  [{row['f1_tail_bins']} bins]")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        keys = list(results[0][0])
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row, _, _ in results:
                w.writerow(row)
        print(f"wrote {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4.2))
        for (row, hist, ft) in results:
            l = np.arange(64)
            m = hist > 0
            (line,) = ax.semilogy(l[m], hist[m], "o", ms=4,
                                  label=f"{row['label']}  "
                                        f"f1={row['f1_tail']:.2f}")
            if ft:
                lf = np.array(ft[3][0], dtype=float)
                a = np.exp(np.log(ft[3][1]).mean() + ft[0] * lf.mean())
                ax.semilogy(lf, a * np.exp(-ft[0] * lf), "-", lw=1,
                            color=line.get_color())
        ax.set_xlabel("segment length l (edges)")
        ax.set_ylabel("pooled count n(l)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs(os.path.dirname(args.plot) or ".", exist_ok=True)
        fig.savefig(args.plot, dpi=120)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
