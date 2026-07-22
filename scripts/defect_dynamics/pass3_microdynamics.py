#!/usr/bin/env python3
"""Pass 3 -- microdynamics from the accepted-move event stream (.events.bin).

Per chain: map each event's clock -> sweep (linear calibration over the sampled
span), classify each event by the vertices it touches at the nearest ts frame
(priority: multimer core [illegal complex >=10] > blinker [<=9] > registry halo
[union of snapshot halo memberships] > bulk), then measure:
  * per-class activity rate: events / vertex / 1000 sweeps (+ share of total)
  * move-type composition per class (1->4, 2->3, 3->2, 4->1, 4-4)
  * gross vs net at multimer cores: vertex creations - destructions, and the
    time series of core-touching activity (steady flicker vs bursts)
  * dynamic heterogeneity: per-vertex event counts, top-1% share
"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
for p in ("python", "scripts"):
    sys.path.insert(0, os.path.join(_ROOT, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.move_geometry import (EVENT_DTYPE,
                                                          MOVE_NAMES,
                                                          SUPPORT_SPLIT)
from tcp_melt import CRYSTALS
from defect_census import defect_cores

SP = sys.argv[1]
V_TOT = 10176
DT = 150

REF_FACETS = np.asarray(ddg.Manifold.load(CRYSTALS["r"], 3).facets())


def load_frames(f):
    by = {}
    for line in open(f):
        r = json.loads(line)
        by.setdefault(r["sweep"], r)
    out = []
    for sw in sorted(by):
        big, small = set(), set()
        for c in by[sw]["comps"]:
            (big if c["size"] >= 10 else small).update(c["members"])
        out.append((sw, big, small))
    return out


def halo_labels(chain):
    """Union of registry-defect vertices (manifold labels) over the chain's
    snapshots, minus illegal members (those are classed by the frames)."""
    halo = set()
    for snap in sorted(glob.glob(f"{SP}/run5h/{chain}_snap*.mfd")):
        fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
        lab = np.unique(fac)
        dv, vclass, adj, V, _ = defect_cores(fac, REF_FACETS, "r", 30)
        halo.update(int(lab[v]) for v in dv)
    return halo


CLS = ["multimer", "blinker", "halo", "bulk"]
rate_rows = []
type_comp = {c: Counter() for c in CLS}
per_vertex = Counter()
core_ts = {}                     # chain -> (sweep_bins, core event counts)
netflux = Counter()              # creations - destructions touching cores
gross = Counter()
class_events = Counter()
class_vsw = Counter()            # vertex-sweeps per class (for rates)

chains = sorted(os.path.basename(f)[:-9]
                for f in glob.glob(f"{SP}/run5h/*.ts.jsonl"))
for chain in chains:
    frames = load_frames(f"{SP}/run5h/{chain}.ts.jsonl")
    fsw = np.array([f[0] for f in frames])
    halo = halo_labels(chain) - set().union(*[b | s for _, b, s in frames])
    ev = np.fromfile(f"{SP}/run5h/{chain}.events.bin", dtype=EVENT_DTYPE)
    # clock -> sweep calibration over the sampled span
    c0, c1 = ev["clock"].min(), ev["clock"].max()
    s0, s1 = 2500, fsw[-1]
    sweeps = s0 + (ev["clock"] - c0) * (s1 - s0) / max(c1 - c0, 1)
    fidx = np.clip(np.searchsorted(fsw, sweeps), 0, len(frames) - 1)

    nbins = int((s1 - s0) // DT) + 1
    cts = np.zeros(nbins)
    for e, sw, fi in zip(ev, sweeps, fidx):
        _, big, small = frames[fi]
        nA, nB = SUPPORT_SPLIT[int(e["type"])]
        labs = [int(x) for x in e["labels"][:nA + nB] if x >= 0]
        if any(l in big for l in labs):
            c = "multimer"
        elif any(l in small for l in labs):
            c = "blinker"
        elif any(l in halo for l in labs):
            c = "halo"
        else:
            c = "bulk"
        class_events[c] += 1
        type_comp[c][MOVE_NAMES[int(e["type"])]] += 1
        for l in labs:
            per_vertex[l] += 1
        if c == "multimer":
            cts[min(int((sw - s0) // DT), nbins - 1)] += 1
            gross[chain] += 1
            if int(e["type"]) == 0:
                netflux[chain] += 1
            elif int(e["type"]) == 3:
                netflux[chain] -= 1
    core_ts[chain] = cts
    # vertex-sweep denominators (avg set sizes over frames)
    nbig = np.mean([len(b) for _, b, _ in frames])
    nsml = np.mean([len(s) for _, _, s in frames])
    span = s1 - s0
    class_vsw["multimer"] += nbig * span
    class_vsw["blinker"] += nsml * span
    class_vsw["halo"] += len(halo) * span
    class_vsw["bulk"] += (V_TOT - nbig - nsml - len(halo)) * span
    print(f"  {chain}: {len(ev)} events  core-share "
          f"{class_events['multimer'] and 0 or 0}", flush=True)

tot = sum(class_events.values())
print(f"\n=== {tot} accepted moves over {len(chains)} chains ===")
print(f"{'class':10s} {'events':>8s} {'share':>7s} {'ev/vtx/1000sw':>14s} {'enrich':>7s}")
bulk_rate = 1000 * class_events["bulk"] / class_vsw["bulk"]
for c in CLS:
    r = 1000 * class_events[c] / max(class_vsw[c], 1)
    print(f"{c:10s} {class_events[c]:8d} {class_events[c]/tot:7.1%} {r:14.2f} "
          f"{r/bulk_rate:7.1f}x")

print("\nmove-type composition by class (%):")
print(f"{'class':10s} " + " ".join(f"{m:>6s}" for m in MOVE_NAMES))
for c in CLS:
    n = sum(type_comp[c].values())
    if n:
        print(f"{c:10s} " + " ".join(f"{100*type_comp[c][m]/n:6.1f}"
                                     for m in MOVE_NAMES))

print("\nmultimer cores, gross vs net (vertex creations - destructions):")
for chain in chains:
    print(f"  {chain}: gross {gross[chain]:6d}   net {netflux[chain]:+5d}")

# heterogeneity
counts = np.array(sorted(per_vertex.values(), reverse=True))
active_verts = len(counts)
top1 = max(1, int(0.01 * V_TOT * len(chains)))
share = counts[:top1].sum() / counts.sum()
print(f"\ndynamic heterogeneity: {active_verts} vertices ever touched "
      f"({active_verts/(V_TOT*len(chains)):.1%} of all); "
      f"top 1% of all vertices carry {share:.1%} of activity")

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.4))
ax = axes[0, 0]
rates = [1000 * class_events[c] / max(class_vsw[c], 1) for c in CLS]
ax.bar(CLS, rates, color=["C3", "C1", "C2", "C0"])
ax.set_yscale("log")
ax.set_ylabel("events / vertex / 1000 sweeps")
ax.set_title("activity rate by class")
for i, r in enumerate(rates):
    ax.text(i, r * 1.1, f"{r:.2f}", ha="center", fontsize=8)
ax.grid(alpha=0.25, axis="y")

ax = axes[0, 1]
bot = np.zeros(len(CLS))
for m, col in zip(MOVE_NAMES, ["C0", "C1", "C2", "C3", "C4"]):
    v = np.array([100 * type_comp[c][m] / max(sum(type_comp[c].values()), 1)
                  for c in CLS])
    ax.bar(CLS, v, bottom=bot, label=m, color=col)
    bot += v
ax.set_ylabel("% of class events"); ax.set_title("move-type composition")
ax.legend(fontsize=7, ncol=5, loc="upper center")
ax.grid(alpha=0.25, axis="y")

ax = axes[1, 0]
for chain, cts in core_ts.items():
    col = "C2" if chain.startswith("above") else "C3"
    ax.plot(2500 + np.arange(len(cts)) * DT, cts / DT * 1000, lw=0.8,
            alpha=0.6, color=col)
ax.plot([], [], color="C2", label="above chains")
ax.plot([], [], color="C3", label="below chains")
ax.set_xlabel("sweep"); ax.set_ylabel("core events / 1000 sw")
ax.set_title("multimer-core activity time series (steady flicker?)")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

ax = axes[1, 1]
ax.loglog(np.arange(1, len(counts) + 1), counts, lw=1.2)
ax.set_xlabel("vertex rank"); ax.set_ylabel("events touching vertex")
ax.set_title(f"activity Zipf plot (top 1% carry {share:.0%})")
ax.grid(alpha=0.25, which="both")

fig.suptitle("Pass 3: microdynamics of the accepted-move stream (8 chains)", y=0.995)
fig.tight_layout()
out = f"{SP}/pass3_microdynamics.png"
fig.savefig(out, dpi=130)
print(f"\nSaved to: {out}")
