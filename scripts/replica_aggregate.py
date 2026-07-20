#!/usr/bin/env python3
"""Aggregate a replica fleet of dope_hold CSVs into mean +/- SEM time series.

Handles kill/resume splicing: a replica with both  <fam>_rK.csv  and
<fam>_rK_resume.csv  is merged by truncating the primary at the snapshot the
resume started from (the largest multiple of --snap-quantum at or below the
primary's last row) and appending the resume rows with that offset added.
The resume's t=0 row (a re-measurement of the snapshot state) is dropped.

Outputs <out>_agg.csv with columns  sweeps, n_reps, <col>_mean, <col>_sem ...
for every numeric column, on the time grid common to all replicas; optionally
a 4-panel PNG of the headline disclination observables.

Usage:
    python scripts/replica_aggregate.py data/replicas/ownpin_mu3 \
        --out out/ownpin_mu3 --plot
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np


def read_rows(path):
    with open(path, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        rows = [[float(x) for x in r] for r in rd if r]
    return header, np.array(rows)


def merge_replica(primary, resume, quantum):
    header, rows = read_rows(primary)
    if resume is None:
        return header, rows
    _, rrows = read_rows(resume)
    offset = quantum * int(rows[-1, 0] // quantum)
    keep = rows[rows[:, 0] <= offset]
    tail = rrows[rrows[:, 0] > 0].copy()
    tail[:, 0] += offset
    return header, np.vstack([keep, tail])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("family", help="fleet basename, e.g. data/replicas/ownpin_mu3"
                                   " (globs <family>_r*.csv)")
    ap.add_argument("--snap-quantum", type=int, default=5000,
                    help="snapshot cadence in sweeps (resume-offset rule)")
    ap.add_argument("--out", default=None,
                    help="output basename (default: out/<family basename>)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    primaries = sorted(p for p in glob.glob(f"{args.family}_r*.csv")
                       if "_resume" not in p and "_agg" not in p)
    if not primaries:
        raise SystemExit(f"no replica CSVs match {args.family}_r*.csv")
    out = args.out or os.path.join("out", os.path.basename(args.family))
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    header = None
    series = []
    for p in primaries:
        res = re.sub(r"\.csv$", "_resume.csv", p)
        h, rows = merge_replica(p, res if os.path.exists(res) else None,
                                args.snap_quantum)
        header = header or h
        if h != header:
            raise SystemExit(f"column mismatch in {p}")
        # de-duplicate any repeated sweep values (keep first)
        _, idx = np.unique(rows[:, 0], return_index=True)
        series.append(rows[np.sort(idx)])
        print(f"  {os.path.basename(p)}: {len(rows)} rows, "
              f"t = {rows[0, 0]:.0f}..{rows[-1, 0]:.0f}")

    common = sorted(set.intersection(*(set(s[:, 0]) for s in series)))
    grid = np.array(common)
    print(f"{len(series)} replicas, common grid {len(grid)} points "
          f"(t <= {grid[-1]:.0f})")

    aligned = np.stack([s[np.isin(s[:, 0], grid)][:, 1:] for s in series])
    mean = aligned.mean(axis=0)
    sem = (aligned.std(axis=0, ddof=1) / np.sqrt(len(series))
           if len(series) > 1 else np.zeros_like(mean))

    cols = header[1:]
    with open(f"{out}_agg.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sweeps", "n_reps"]
                   + [f"{c}_{s}" for c in cols for s in ("mean", "sem")])
        for i, t in enumerate(grid):
            w.writerow([int(t), len(series)]
                       + [x for j in range(len(cols))
                          for x in (mean[i, j], sem[i, j])])
    print(f"wrote {out}_agg.csv")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        panels = [c for c in ("n_dop", "n_six", "net_giant_frac", "mean_seg")
                  if c in cols]
        fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
        for ax, c in zip(axes.flat, panels):
            j = cols.index(c)
            ax.fill_between(grid, mean[:, j] - sem[:, j],
                            mean[:, j] + sem[:, j], alpha=0.3)
            ax.plot(grid, mean[:, j], lw=1.2)
            ax.set_title(c, fontsize=10)
            ax.grid(alpha=0.3)
        for ax in axes[-1]:
            ax.set_xlabel("sweeps")
        fig.suptitle(f"{os.path.basename(args.family)}  "
                     f"({len(series)} replicas, mean ± SEM)")
        fig.tight_layout()
        fig.savefig(f"{out}_agg.png", dpi=120)
        print(f"wrote {out}_agg.png")


if __name__ == "__main__":
    main()
