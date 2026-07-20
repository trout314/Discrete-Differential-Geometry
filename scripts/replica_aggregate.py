#!/usr/bin/env python3
"""Aggregate a replica fleet of dope_hold CSVs into mean +/- SEM time series.

Segment splicing (new-style, absolute sweep clocks): a replica's trajectory
is  <fam>_rK.csv  plus any  <fam>_rK_seg*.csv  continuations. Each segment
CSV carries an ``# ensemble=<hash> segment_start=<n>`` comment written by
dope_hold; splicing truncates the accumulated rows at the continuation's
segment_start and concatenates. Ensemble hashes are verified when both sides
carry one -- a mismatch is a hard error (--allow-mixed downgrades it), and
a splice with a hashless legacy side is warned about. When _seg* files
exist, a legacy _resume.csv for the same replica is IGNORED (superseded).

Legacy kill/resume splicing (relative clocks): with only  <fam>_rK_resume.csv
present, the primary is truncated at the snapshot the resume started from
(largest multiple of --snap-quantum at or below the primary's last row) and
the resume rows are appended with that offset added; the resume's t=0 row
(a re-measurement of the snapshot state) is dropped. No hash protection.

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
    """Returns (meta, header, rows). meta parsed from leading '#' comment
    lines of the form 'key=value ...' (empty for legacy CSVs)."""
    meta = {}
    with open(path, newline="") as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines) and lines[i].startswith("#"):
        for tok in lines[i][1:].split():
            if "=" in tok:
                k, val = tok.split("=", 1)
                meta[k] = val
        i += 1
    rd = csv.reader(lines[i:])
    header = next(rd)
    rows = [[float(x) for x in r] for r in rd if r]
    return meta, header, np.array(rows)


def splice(base_rows, seg_rows, seg_start):
    return np.vstack([base_rows[base_rows[:, 0] <= seg_start],
                      seg_rows[seg_rows[:, 0] > seg_start]])


def merge_replica(primary, allow_mixed, quantum):
    meta, header, rows = read_rows(primary)
    stem = re.sub(r"\.csv$", "", primary)
    segs = sorted(glob.glob(f"{stem}_seg*.csv"))
    resume = f"{stem}_resume.csv"
    notes = []
    if segs:
        if os.path.exists(resume):
            notes.append(f"ignoring {os.path.basename(resume)} "
                         "(superseded by _seg*)")
        for sp in segs:
            smeta, sheader, srows = read_rows(sp)
            if sheader != header:
                raise SystemExit(f"column mismatch in {sp}")
            if "segment_start" not in smeta:
                raise SystemExit(f"{sp}: no segment_start comment; "
                                 "not a dope_hold segment CSV")
            hp, hs = meta.get("ensemble"), smeta.get("ensemble")
            if hp and hs and hp != hs:
                msg = (f"ensemble hash mismatch splicing {sp}: "
                       f"{hp} vs {hs}")
                if not allow_mixed:
                    raise SystemExit(msg + " (use --allow-mixed to force)")
                notes.append("MIXED-ENSEMBLE SPLICE: " + msg)
            elif not (hp and hs):
                notes.append(f"unverified splice with {os.path.basename(sp)} "
                             "(missing ensemble hash on one side)")
            rows = splice(rows, srows, float(smeta["segment_start"]))
    elif os.path.exists(resume):
        _, rheader, rrows = read_rows(resume)
        if rheader != header:
            raise SystemExit(f"column mismatch in {resume}")
        offset = quantum * int(rows[-1, 0] // quantum)
        tail = rrows[rrows[:, 0] > 0].copy()
        tail[:, 0] += offset
        rows = np.vstack([rows[rows[:, 0] <= offset], tail])
        notes.append(f"legacy resume splice at {offset} (no hash protection)")
    return header, rows, notes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("family", help="fleet basename, e.g. data/replicas/ownpin_mu3"
                                   " (globs <family>_r*.csv)")
    ap.add_argument("--snap-quantum", type=int, default=5000,
                    help="snapshot cadence in sweeps (legacy resume-offset "
                         "rule)")
    ap.add_argument("--allow-mixed", action="store_true",
                    help="downgrade ensemble-hash mismatches to a warning")
    ap.add_argument("--out", default=None,
                    help="output basename (default: out/<family basename>)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    primaries = sorted(p for p in glob.glob(f"{args.family}_r*.csv")
                       if not re.search(r"_(resume|seg\d*|agg)", p))
    if not primaries:
        raise SystemExit(f"no replica CSVs match {args.family}_r*.csv")
    out = args.out or os.path.join("out", os.path.basename(args.family))
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    header = None
    series = []
    for p in primaries:
        h, rows, notes = merge_replica(p, args.allow_mixed, args.snap_quantum)
        header = header or h
        if h != header:
            raise SystemExit(f"column mismatch in {p}")
        # de-duplicate any repeated sweep values (keep first)
        _, idx = np.unique(rows[:, 0], return_index=True)
        series.append(rows[np.sort(idx)])
        print(f"  {os.path.basename(p)}: {len(rows)} rows, "
              f"t = {rows[0, 0]:.0f}..{rows[-1, 0]:.0f}")
        for n in notes:
            print(f"    NOTE: {n}")

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
