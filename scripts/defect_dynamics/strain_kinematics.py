#!/usr/bin/env python3
"""Defect kinematics from a --track ts.jsonl series (strain/channel gases).

Protocols follow pass1_kinematics.py (run5h) and mgas_analyze.py part 1,
adapted to the crystal_gas_scan --track producer:
  * frames may include a header row (skipped);
  * centroids come from the MAINTAINED vertex lift, so the tree-rebuild
    gauge glitches pass1 filters (single-frame steps > 1 cell) cannot occur
    -- any large step is real (or a genuine complex merge/split);
  * worldlines: complex identity linked frame-to-frame by member overlap
    (best match with Jaccard >= --jmin, pass1's linkage rule);
  * outputs: population Jaccard vs lag (mgas part 1), worldline lifetime
    census, per-worldline net displacement vs path length (glide vs jitter),
    pooled centroid MSD vs lag by size class, blink recurrence (death ->
    rebirth overlap, pass1's recurrence check).

Usage:
    python scripts/defect_dynamics/strain_kinematics.py PRE [--jmin 0.3]
    (PRE such that PRE.ts.jsonl exists; box read from the header row or
     --mcell)
"""
import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np

CELL = 1e6


def minimg(d, box):
    return (d + box / 2) % box - box / 2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pre")
    ap.add_argument("--jmin", type=float, default=0.3,
                    help="min member-Jaccard to link a worldline")
    ap.add_argument("--mcell", type=int, default=0,
                    help="override torus cells/axis (else from header)")
    ap.add_argument("--burn-frames", type=int, default=0,
                    help="drop the first k frames (quench transient)")
    args = ap.parse_args()

    raw = [json.loads(l) for l in open(args.pre + ".ts.jsonl")]
    hdr = next((r for r in raw if r.get("kind") == "header"), {})
    frames = [r for r in raw if "sweep" in r][args.burn_frames:]
    mcell = args.mcell or int(hdr.get("mcell", 0))
    assert mcell > 0, "pass --mcell (no header with mcell found)"
    box = CELL * mcell
    sw = [r["sweep"] for r in frames]
    dt = sw[1] - sw[0]
    print(f"{os.path.basename(args.pre)}: {len(frames)} frames, "
          f"dt={dt} sweeps, box={mcell} cells")

    # ---- population mobility: Jaccard of the defect-vertex set vs lag ----
    sets = [frozenset(x for c in r["members"] for x in c) for r in frames]
    print("  population Jaccard vs lag:", end=" ")
    for lag_sw in (dt, 5 * dt, 20 * dt, 100 * dt):
        st = lag_sw // dt
        vals = [len(sets[i] & sets[i + st]) / max(len(sets[i] | sets[i + st]), 1)
                for i in range(0, len(sets) - st, max(1, st // 2))]
        if vals:
            print(f"J({lag_sw})={np.mean(vals):.2f}", end="  ")
    print()

    # ---- worldline linking by member overlap (pass1 rule) ----
    # worldline: list of (frame idx, complex idx)
    worldlines = []
    active = {}          # complex idx in previous frame -> worldline idx
    for fi, r in enumerate(frames):
        cur = [frozenset(c) for c in r["members"]]
        nxt_active = {}
        if fi == 0:
            for ci in range(len(cur)):
                worldlines.append([(fi, ci)])
                nxt_active[ci] = len(worldlines) - 1
        else:
            prev = [frozenset(c) for c in frames[fi - 1]["members"]]
            used_prev = set()
            for ci, cs in enumerate(cur):
                best, bj = -1, args.jmin
                for pi, ps in enumerate(prev):
                    if pi in used_prev or pi not in active:
                        continue
                    j = len(cs & ps) / len(cs | ps)
                    if j > bj:
                        best, bj = pi, j
                if best >= 0:
                    wi = active[best]
                    worldlines[wi].append((fi, ci))
                    nxt_active[ci] = wi
                    used_prev.add(best)
                else:
                    worldlines.append([(fi, ci)])
                    nxt_active[ci] = len(worldlines) - 1
        active = nxt_active

    # ---- per-worldline kinematics ----
    lifetimes = []
    rows = []
    for wl in worldlines:
        lifetimes.append(len(wl))
        if len(wl) < 3:
            continue
        cents = np.array([frames[fi]["cents"][ci] for fi, ci in wl], float)
        steps = minimg(np.diff(cents, axis=0), box)
        traj = np.vstack([[0, 0, 0], np.cumsum(steps, axis=0)])
        net = np.linalg.norm(traj[-1]) / CELL
        path = np.sum(np.linalg.norm(steps, axis=1)) / CELL
        size = np.mean([len(frames[fi]["members"][ci]) for fi, ci in wl])
        rows.append((len(wl), size, net, path))

    lifetimes = np.array(lifetimes)
    print(f"  worldlines: {len(worldlines)} total; lifetime frames "
          f"median {np.median(lifetimes):.0f}, max {lifetimes.max()}, "
          f"single-frame (blink-fast) {np.mean(lifetimes == 1):.0%}")

    if rows:
        rows.sort(key=lambda r: -r[0])
        print("  longest worldlines (frames, <verts>, net cells, path cells,"
              " net/path):")
        for L, size, net, path in rows[:8]:
            print(f"    {L:4d}  {size:5.1f}  {net:6.2f}  {path:7.2f}  "
                  f"{net / path if path > 0 else 0:.2f}")

    # ---- pooled MSD vs lag for worldlines >= 10 frames ----
    print("  pooled centroid MSD (cells^2) vs lag (sweeps):")
    trajs = []
    for wl in worldlines:
        if len(wl) < 10:
            continue
        cents = np.array([frames[fi]["cents"][ci] for fi, ci in wl], float)
        steps = minimg(np.diff(cents, axis=0), box)
        trajs.append(np.vstack([[0, 0, 0], np.cumsum(steps, axis=0)]) / CELL)
    for st in (1, 2, 5, 10, 25):
        v = [np.sum((tr[st:] - tr[:-st]) ** 2, axis=1)
             for tr in trajs if len(tr) > st]
        if v:
            v = np.concatenate(v)
            print(f"    lag {st * dt:5d}: MSD = {np.mean(v):.4f} "
                  f"(n={len(v)})")

    # ---- blink recurrence (pass1): dead worldline reborn at same site ----
    ends = {}
    reborn = 0
    births = 0
    for wl in worldlines:
        fi0, ci0 = wl[0]
        if fi0 > 0:
            births += 1
            m0 = frozenset(frames[fi0]["members"][ci0])
            for (fe, ce), wend in list(ends.items()):
                if fe < fi0 and len(m0 & wend) / len(m0 | wend) >= args.jmin:
                    reborn += 1
                    break
        fiE, ciE = wl[-1]
        if fiE < len(frames) - 1:
            ends[(fiE, ciE)] = frozenset(frames[fiE]["members"][ciE])
    if births:
        print(f"  blink recurrence: {reborn}/{births} births overlap a "
              f"prior death site ({reborn / births:.0%})")


if __name__ == "__main__":
    main()
