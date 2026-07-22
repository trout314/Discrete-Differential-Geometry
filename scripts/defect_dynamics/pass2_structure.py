#!/usr/bin/env python3
"""Pass 2 -- structural anatomy of the defect complexes (all 32 snapshots).

Per illegal complex (vertex_classes imp>0, connected):
  * registry halo: size of the crystal_grains registry-defect component
    containing it (covering map vs R reference, min grain 30);
  * edge-degree motif: intra-complex edge degrees (deg<=4 = positive
    disclination, deg>=7 = negative) -> knot content vs the (3,4,4) 5-knot;
  * net curvature charge Q = sum qR - n * qbar_crystal;
  * gyration radius Rg from harmonic torus coords (min-image) -> Rg vs N.
Plus: blinker-birth-site test -- distance from small-complex (<=9) sites to the
nearest big core (>=10) vs a uniform random null.
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
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS, edges_and_degrees
import crystal_grains as cg
from tcp_melt import CRYSTALS
from tcp_reference import STRUCTURES
from dopant_pairs import vertex_classes
from defect_census import defect_cores

SP = sys.argv[1]
CELL = 1e6

# reference (build once)
REF_FACETS = np.asarray(ddg.Manifold.load(CRYSTALS["r"], 3).facets())

rows = []          # one per illegal complex
blink_d, null_d = [], []
rng = np.random.default_rng(0)

snaps = sorted(glob.glob(f"{SP}/run5h/*_snap*.mfd"))
for snap in snaps:
    tag = os.path.basename(snap)[:-4]
    fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
    eu, ecnt, deg, V = edges_and_degrees(fac)
    qR = FIELDS["curvature_charge"](fac)
    n6, imp, adj = vertex_classes(fac)

    # illegal complexes (dense labels)
    illv = [v for v in range(V) if imp[v] > 0]
    seen, comps = set(), []
    for s0 in illv:
        if s0 in seen:
            continue
        st, comp = [s0], []
        seen.add(s0)
        while st:
            u = st.pop(); comp.append(u)
            for w in adj[u]:
                if imp[w] > 0 and w not in seen:
                    seen.add(w); st.append(w)
        comps.append(sorted(comp))

    # registry cores + their components
    dv, vclass, adj2, V2, gsizes = defect_cores(fac, REF_FACETS, "r", 30)
    dvset = set(dv)
    rcomp_of = {}
    seenr = set()
    rcomps = []
    for s0 in dv:
        if s0 in seenr:
            continue
        st, comp = [s0], []
        seenr.add(s0)
        while st:
            u = st.pop(); comp.append(u)
            for w in adj2[u]:
                if w in dvset and w not in seenr:
                    seenr.add(w); st.append(w)
        rid = len(rcomps)
        rcomps.append(comp)
        for u in comp:
            rcomp_of[u] = rid

    # harmonic coords for geometry
    edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    X = frac * np.diag(basis)         # Cartesian raw units
    P = np.abs(np.diag(basis))

    qbar = np.median(qR)              # crystal background (defects are rare)
    # edge lookup for intra-complex motifs
    edge_deg = {(int(a), int(b)): int(c) for (a, b), c in zip(eu, ecnt)}

    cens_big, cens_small = [], []
    for comp in comps:
        n = len(comp)
        cset = set(comp)
        # registry halo
        rids = {rcomp_of.get(u) for u in comp} - {None}
        halo = sum(len(rcomps[r]) for r in rids) if rids else n
        # intra-complex edge motif
        degs = [edge_deg[(min(u, w), max(u, w))]
                for u in comp for w in adj[u] if w > u and w in cset]
        mot = Counter(degs)
        n3 = mot.get(3, 0); n4 = mot.get(4, 0)
        n7p = sum(v for k, v in mot.items() if k >= 7)
        # charge + gyration (min-image about first member)
        Q = float(qR[comp].sum() - n * qbar)
        d = X[comp] - X[comp[0]]
        d -= np.round(d / P) * P
        cen = X[comp[0]] + d.mean(0)
        Rg = float(np.sqrt(((d - d.mean(0)) ** 2).sum(1).mean())) / CELL
        rows.append(dict(tag=tag, n=n, halo=halo, n3=n3, n4=n4, n7p=n7p,
                         Q=Q, Rg=Rg))
        (cens_big if n >= 10 else cens_small).append(cen)

    # blinker-site test within this snapshot
    if cens_big and cens_small:
        B = np.array(cens_big)
        def dmin(pt):
            dd = B - pt
            dd -= np.round(dd / P) * P
            return np.linalg.norm(dd, axis=1).min() / CELL
        for c in cens_small:
            blink_d.append(dmin(c))
        for c in rng.uniform(0, 1, (40, 3)) * P:
            null_d.append(dmin(c))
    print(f"  {tag}: {len(comps)} complexes "
          f"({sum(1 for c in comps if len(c) >= 10)} big)", flush=True)

# ------------------------------------------------------------------ report
R = rows
print(f"\n=== {len(R)} complexes over {len(snaps)} snapshots ===")
bins = [("5", lambda n: n == 5), ("6-9", lambda n: 6 <= n <= 9),
        ("10-14", lambda n: 10 <= n <= 14), (">=15", lambda n: n >= 15)]
print(f"{'size':6s} {'n':>4s} {'halo/n':>7s} {'<n3>':>5s} {'<n4>':>5s} {'<n7+>':>6s} "
      f"{'(n3:n4)':>8s} {'<Q>':>7s} {'sd(Q)':>6s} {'<Rg>':>6s}")
for name, sel in bins:
    g = [r for r in R if sel(r["n"])]
    if not g:
        continue
    h = np.mean([r["halo"] / r["n"] for r in g])
    a3 = np.mean([r["n3"] for r in g]); a4 = np.mean([r["n4"] for r in g])
    a7 = np.mean([r["n7p"] for r in g])
    Qm = np.mean([r["Q"] for r in g]); Qs = np.std([r["Q"] for r in g])
    Rgm = np.mean([r["Rg"] for r in g])
    print(f"{name:6s} {len(g):4d} {h:7.2f} {a3:5.2f} {a4:5.2f} {a7:6.2f} "
          f"{a3/max(a4,1e-9):8.2f} {Qm:7.2f} {Qs:6.2f} {Rgm:6.2f}")

# Rg vs N exponent (big complexes)
big = [r for r in R if r["n"] >= 5 and r["Rg"] > 0]
lN = np.log([r["n"] for r in big]); lR = np.log([r["Rg"] for r in big])
nu = np.polyfit(lN, lR, 1)[0]
print(f"\nRg ~ N^nu:  nu = {nu:.2f}   (1/3 blob, ~0.59 SAW, 1 straight string)")

# blinker-site null
blink_d, null_d = np.array(blink_d), np.array(null_d)
print(f"\nblinker-site distance to nearest big core (cells): "
      f"med {np.median(blink_d):.2f} (n={len(blink_d)})  "
      f"vs null med {np.median(null_d):.2f} (n={len(null_d)})")
from scipy.stats import ks_2samp
ks = ks_2samp(blink_d, null_d)
print(f"KS test: D={ks.statistic:.3f}  p={ks.pvalue:.2e}  "
      f"-> {'blinkers cluster near cores' if np.median(blink_d) < np.median(null_d) and ks.pvalue < 0.01 else 'no significant templating' if ks.pvalue >= 0.01 else 'blinkers AVOID cores'}")

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
ax = axes[0, 0]
ax.scatter([r["n"] for r in R], [r["halo"] / r["n"] for r in R],
           s=14, alpha=0.5, c=["C1" if r["n"] < 10 else "C0" for r in R])
ax.set_xlabel("illegal-core size N"); ax.set_ylabel("registry halo / N")
ax.set_title("off-registry dressing vs core size"); ax.grid(alpha=0.25)

ax = axes[0, 1]
nn = np.array([r["n"] for r in big]); rg = np.array([r["Rg"] for r in big])
ax.loglog(nn, rg, "o", ms=4, alpha=0.5)
xs = np.array([nn.min(), nn.max()])
b = np.exp(np.mean(lR) - nu * np.mean(lN))
ax.loglog(xs, b * xs ** nu, "k-", lw=1.2, label=fr"$R_g \sim N^{{{nu:.2f}}}$")
for e, ls in [(1 / 3, ":"), (1.0, "--")]:
    ax.loglog(xs, rg.mean() * (xs / np.exp(np.mean(lN))) ** e, ls, color="gray",
              lw=0.9, label=f"$N^{{{e:.2f}}}$")
ax.set_xlabel("N"); ax.set_ylabel(r"$R_g$ (cells)")
ax.set_title("gyration radius vs size"); ax.legend(fontsize=8); ax.grid(alpha=0.25)

ax = axes[1, 0]
for name, sel, col in [("5", lambda n: n == 5, "C1"),
                       ("10-14", lambda n: 10 <= n <= 14, "C0"),
                       (">=15", lambda n: n >= 15, "C3")]:
    qs = [r["Q"] for r in R if sel(r["n"])]
    if qs:
        ax.hist(qs, bins=24, alpha=0.55, color=col, label=f"size {name}", density=True)
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("net curvature charge Q (vs crystal background)")
ax.set_ylabel("density"); ax.set_title("complex net charge"); ax.legend(fontsize=8)
ax.grid(alpha=0.25)

ax = axes[1, 1]
ax.hist(null_d, bins=24, alpha=0.5, color="gray", density=True, label="uniform null")
ax.hist(blink_d, bins=24, alpha=0.6, color="C1", density=True, label="blinker sites")
ax.set_xlabel("distance to nearest big core (cells)"); ax.set_ylabel("density")
ax.set_title(f"blinker templating (KS p={ks.pvalue:.1e})")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

fig.suptitle("Pass 2: defect-complex structure (32 snapshots, 8 chains)", y=0.995)
fig.tight_layout()
out = f"{SP}/pass2_structure.png"
fig.savefig(out, dpi=130)
json.dump(rows, open(f"{SP}/pass2_rows.json", "w"))
print(f"\nSaved to: {out}")
