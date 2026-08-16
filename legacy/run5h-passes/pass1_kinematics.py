#!/usr/bin/env python3
"""Pass 1 -- defect kinematics from the run5h ts.jsonl trajectories.

Tracks complex identity through the 150-sweep frames by member-overlap linkage,
then classifies each worldline and measures lifetimes, churn, excursions, MSD,
and blinking (death->rebirth recurrence).

Position protocol (validated): centroids are raw cocycle tree-lift coords,
internally consistent in time (99.4% of unchanged-membership steps are exactly
0); displacements use min-image steps in cell units (CELL=1e6 raw, period 4
cells); single-frame steps > 1 cell are spanning-tree gauge glitches (0.09%)
and are dropped from the unwrapped trajectory.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SP = sys.argv[1]
CELL, P4 = 1e6, 4e6
DT = 150                       # sweeps per frame
GLITCH = 1.0                   # cells; single-frame step above this = gauge jump
SIZE_CLASSES = [("5", lambda s: s == 5), ("6-9", lambda s: 6 <= s <= 9),
                ("10-14", lambda s: 10 <= s <= 14), (">=15", lambda s: s >= 15)]


def size_class(s):
    for name, f in SIZE_CLASSES:
        if f(s):
            return name
    return "<5"


def load_frames(f):
    by_sweep = {}
    for line in open(f):
        r = json.loads(line)
        by_sweep.setdefault(r["sweep"], r)
    frames = []
    for sw in sorted(by_sweep):
        comps = [(set(c["members"]), np.asarray(c["centroid"], float))
                 for c in by_sweep[sw]["comps"] if c["centroid"] is not None]
        frames.append((sw, comps))
    return frames


def mi(d):
    return d - np.round(d / P4) * P4


class Track:
    def __init__(self, tid, sw, mem, cen):
        self.tid = tid
        self.sweeps = [sw]
        self.mems = [mem]
        self.sizes = [len(mem)]
        self.raw_cen = [cen]
        self.unwrap = [np.zeros(3)]        # displacement from birth, de-glitched
        self.jac = []
        self.glitches = 0

    def extend(self, sw, mem, cen):
        a, b = self.mems[-1], mem
        self.jac.append(len(a & b) / len(a | b))
        step = mi(cen - self.raw_cen[-1]) / CELL
        if np.linalg.norm(step) > GLITCH:
            self.glitches += 1
            step = np.zeros(3)
        self.unwrap.append(self.unwrap[-1] + step)
        self.sweeps.append(sw)
        self.mems.append(mem)
        self.sizes.append(len(mem))
        self.raw_cen.append(cen)

    def finish(self, last_sweep_of_chain):
        self.life = self.sweeps[-1] - self.sweeps[0] + DT
        self.censored = (self.sweeps[-1] >= last_sweep_of_chain)
        self.founding = False  # set by caller
        U = np.array(self.unwrap)
        self.excursion = float(np.linalg.norm(U, axis=1).max())
        self.net = float(np.linalg.norm(U[-1]))
        self.msize = float(np.median(self.sizes))
        self.mjac = float(np.mean(self.jac)) if self.jac else 1.0


def build_tracks(frames):
    tracks, active, tid = [], {}, 0     # active: idx-in-frame -> Track
    prev = []
    for fi, (sw, comps) in enumerate(frames):
        used_new = set()
        new_active = {}
        # link: greedy by shared-member count, require >=2 shared or J>=0.3
        cand = []
        for ai, tr in active.items():
            a = tr.mems[-1]
            for ni, (mem, cen) in enumerate(comps):
                sh = len(a & mem)
                if sh >= 2 or (sh and sh / len(a | mem) >= 0.3):
                    cand.append((sh, ai, ni))
        used_old = set()
        for sh, ai, ni in sorted(cand, reverse=True):
            if ai in used_old or ni in used_new:
                continue
            used_old.add(ai); used_new.add(ni)
            tr = active[ai]
            tr.extend(sw, comps[ni][0], comps[ni][1])
            new_active[ni] = tr
        for ai, tr in active.items():
            if ai not in used_old:
                tracks.append(tr)                     # death
        for ni, (mem, cen) in enumerate(comps):
            if ni not in used_new:
                t = Track(tid, sw, mem, cen); tid += 1
                if fi == 0:
                    t.founding0 = True
                new_active[ni] = t
        active = new_active
        prev = comps
    tracks += list(active.values())
    return tracks


# ------------------------------------------------------------------ run
all_tracks = []
recur = []          # (deadtime_sweeps, dist_cells, size) for death->rebirth matches
chain_of = {}
for tsf in sorted(glob.glob(f"{SP}/run5h/*.ts.jsonl")):
    chain = os.path.basename(tsf)[:-9]
    frames = load_frames(tsf)
    last_sw = frames[-1][0]
    first_sw = frames[0][0]
    trs = build_tracks(frames)
    for t in trs:
        t.finish(last_sw)
        t.founding = (t.sweeps[0] == first_sw)
        chain_of[id(t)] = chain
    # blinking: match each death to the nearest subsequent birth
    deaths = [(t.sweeps[-1], t.raw_cen[-1], t.msize) for t in trs if not t.censored]
    births = [(t.sweeps[0], t.raw_cen[0]) for t in trs if not t.founding]
    for dsw, dc, dsz in deaths:
        best = None
        for bsw, bc in births:
            if bsw <= dsw:
                continue
            dist = np.linalg.norm(mi(bc - dc)) / CELL
            if dist < 0.7 and (best is None or bsw < best[0]):
                best = (bsw, dist)
        if best:
            recur.append((best[0] - dsw, best[1], dsz))
    all_tracks += trs

T = all_tracks
print(f"tracks: {len(T)} total over 8 chains "
      f"(founding {sum(t.founding for t in T)}, censored-at-end {sum(t.censored for t in T)}, "
      f"gauge glitches dropped {sum(t.glitches for t in T)})\n")

# ---- classification ----
def classify(t):
    if t.excursion < 0.35:
        return "frozen" if t.mjac >= 0.8 else "pinned-churning"
    return "mobile" if t.excursion >= 1.0 else "jiggling"

cls = defaultdict(list)
for t in T:
    cls[classify(t)].append(t)

print(f"{'class':16s} {'n':>4s} {'med size':>8s} {'med life(sw)':>12s} "
      f"{'med excur':>9s} {'med Jac':>8s} {'%censored':>9s}")
for name in ("frozen", "pinned-churning", "jiggling", "mobile"):
    g = cls.get(name, [])
    if not g:
        continue
    print(f"{name:16s} {len(g):4d} {np.median([t.msize for t in g]):8.1f} "
          f"{np.median([t.life for t in g]):12.0f} "
          f"{np.median([t.excursion for t in g]):9.2f} "
          f"{np.median([t.mjac for t in g]):8.2f} "
          f"{100*np.mean([t.censored for t in g]):8.0f}%")

# ---- by size class ----
print(f"\n{'size class':10s} {'n':>4s} {'med life':>9s} {'%immortal':>9s} "
      f"{'med excur':>9s} {'med Jac':>8s} {'births/10ksw':>12s}")
span = 17700 * 8            # total observed sweeps
for name, sel in SIZE_CLASSES:
    g = [t for t in T if sel(t.msize)]
    if not g:
        continue
    births = sum(1 for t in g if not t.founding)
    print(f"{name:10s} {len(g):4d} {np.median([t.life for t in g]):9.0f} "
          f"{100*np.mean([t.censored and t.founding for t in g]):8.0f}% "
          f"{np.median([t.excursion for t in g]):9.2f} "
          f"{np.median([t.mjac for t in g]):8.2f} "
          f"{1e4*births/span:12.2f}")

# ---- founding-core survival ----
print("\nfounding cores (alive at first frame):")
fnd = [t for t in T if t.founding]
for name, sel in SIZE_CLASSES:
    g = [t for t in fnd if sel(t.msize)]
    if g:
        alive = np.mean([t.censored for t in g])
        print(f"  size {name:6s}: n={len(g):3d}  survive-to-end {100*alive:3.0f}%  "
              f"med life {np.median([t.life for t in g]):6.0f} sw")

# ---- blinking ----
if recur:
    dt_, dist_, sz_ = map(np.array, zip(*recur))
    print(f"\nblinking (death -> rebirth within 0.7 cell): {len(recur)} events")
    print(f"  dead time: med {np.median(dt_):.0f} sw  p90 {np.percentile(dt_,90):.0f} sw")
    print(f"  rebirth distance: med {np.median(dist_):.2f} cells")
    print(f"  size of blinker: med {np.median(sz_):.0f}")
ndeath = sum(1 for t in T if not t.censored)
print(f"  recurrence fraction of all deaths: {len(recur)}/{ndeath} = {len(recur)/ndeath:.2f}")

# ---- MSD by size class ----
lags = np.array([1, 2, 4, 8, 16, 32, 64])
msd = {}
for name, sel in SIZE_CLASSES:
    acc = {L: [] for L in lags}
    for t in T:
        if not sel(t.msize) or len(t.unwrap) < 3:
            continue
        U = np.array(t.unwrap)
        for L in lags:
            if L < len(U):
                d = U[L:] - U[:-L]
                acc[L] += list((d ** 2).sum(1))
    msd[name] = np.array([np.mean(acc[L]) if acc[L] else np.nan for L in lags])

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))

ax = axes[0, 0]                       # excursion vs time, all tracks
for t in T:
    U = np.array(t.unwrap)
    sw = np.array(t.sweeps) - t.sweeps[0]
    c = {"5": "C1", "6-9": "C2", "10-14": "C0", ">=15": "C3"}.get(size_class(t.msize), "gray")
    ax.plot(sw, np.linalg.norm(U, axis=1), lw=0.7, alpha=0.45, color=c)
for name, col in [("5", "C1"), ("6-9", "C2"), ("10-14", "C0"), (">=15", "C3")]:
    ax.plot([], [], color=col, lw=2, label=f"size {name}")
ax.axhline(0.35, color="k", lw=0.7, ls=":")
ax.axhline(1.0, color="k", lw=0.7, ls="--")
ax.set_xlabel("sweeps since birth"); ax.set_ylabel("|displacement| (cells)")
ax.set_title("worldline excursions (all 8 chains)")
ax.legend(fontsize=7); ax.grid(alpha=0.25)

ax = axes[0, 1]                       # lifetime vs size
for t in T:
    ax.plot(t.msize, t.life, "^" if t.censored else "o",
            ms=5 if t.censored else 3.5,
            color="C3" if t.censored else "C0", alpha=0.55)
ax.plot([], [], "^", color="C3", label="censored (alive at end)")
ax.plot([], [], "o", color="C0", label="completed")
ax.set_xlabel("median complex size"); ax.set_ylabel("lifetime (sweeps)")
ax.set_yscale("log"); ax.set_title("lifetime vs size")
ax.legend(fontsize=7); ax.grid(alpha=0.25)

ax = axes[1, 0]                       # MSD
for name, col in [("5", "C1"), ("6-9", "C2"), ("10-14", "C0"), (">=15", "C3")]:
    ax.loglog(lags * DT, msd[name], "o-", color=col, ms=4, label=f"size {name}")
ll = np.array([lags[0], lags[-1]]) * DT
ax.loglog(ll, 0.02 * ll / ll[0], "k--", lw=0.9, label="~t (diffusive)")
ax.set_xlabel("lag (sweeps)"); ax.set_ylabel(r"MSD (cells$^2$)")
ax.set_title("MSD by size class"); ax.legend(fontsize=7)
ax.grid(alpha=0.25, which="both")

ax = axes[1, 1]                       # size spectrum (pooled frame occupancy)
allsizes = []
for tsf in sorted(glob.glob(f"{SP}/run5h/*.ts.jsonl")):
    for sw, comps in load_frames(tsf):
        allsizes += [len(m) for m, _ in comps]
ax.hist(allsizes, bins=np.arange(4.5, 30.5), color="C0", alpha=0.8)
ax.set_xlabel("complex size (illegal vertices)"); ax.set_ylabel("frame-occupancy count")
ax.set_title("size spectrum (pooled)"); ax.grid(alpha=0.25)

fig.suptitle("Pass 1: defect-complex kinematics (run5h, 8 chains, 150-sweep frames)",
             y=0.995)
fig.tight_layout()
out = f"{SP}/pass1_kinematics.png"
fig.savefig(out, dpi=130)
print(f"\nSaved to: {out}")
