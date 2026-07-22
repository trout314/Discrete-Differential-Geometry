#!/usr/bin/env python3
"""Validate the ts.jsonl raw-tree-lift centroids against exact min-image
centroids recomputed from the saved snapshot cocycles.

For each snapshot (chain x 4): match ts.jsonl complexes at the same sweep by
member set; recompute the complex centroid properly on the torus (min-image
about one member, harmonic-gauged fractional coords -> Cartesian); report
  * intra-cluster spread in the raw lift vs min-image (lift-splitting test:
    a split cluster has raw spread ~ a full period, min-image spread ~ 1 cell)
  * the offset between stored and exact centroid (mod lattice).
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, "python")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import cocycle as coc

SP = sys.argv[1]
CELL = 1e6                      # raw units per unit cell (m4: period 4e6)


def load_ts(f):
    out = {}
    for line in open(f):
        r = json.loads(line)
        out.setdefault(r["sweep"], r)     # driver may duplicate; keep first
    return out


nbad = ntot = 0
worst = []
for tsf in sorted(glob.glob(f"{SP}/run5h/*.ts.jsonl")):
    chain = os.path.basename(tsf)[:-9]
    ts = load_ts(tsf)
    for snap in sorted(glob.glob(f"{SP}/run5h/{chain}_snap*.mfd")):
        sweep = int(snap.split("snap")[1].split(".")[0])
        if sweep not in ts:
            continue
        fac = np.asarray(Manifold.load(snap, 3).facets())
        edges, omega, _ = coc.load_cocycle(os.path.splitext(snap)[0] + ".cocycle.npz")
        frac, basis = coc.torus_positions(fac, edges, omega)
        lab = np.unique(fac)
        d_of = {int(v): i for i, v in enumerate(lab)}    # manifold -> dense
        P = np.abs(np.diag(basis))                       # periods (raw units)
        for c in ts[sweep]["comps"]:
            mem = [d_of[v] for v in c["members"] if v in d_of]
            if len(mem) < len(c["members"]) or c["centroid"] is None:
                continue
            X = frac[mem] * np.diag(basis)               # Cartesian raw units
            # min-image about first member
            ref = X[0]
            d = X - ref
            d -= np.round(d / P) * P
            cen_exact = ref + d.mean(0)
            spread_mi = np.sqrt((d - d.mean(0)).var(0).sum())
            # stored centroid comparison, mod lattice
            off = np.asarray(c["centroid"]) - cen_exact
            off -= np.round(off / P) * P
            ntot += 1
            bad = spread_mi > 2.0 * CELL or np.linalg.norm(off) > 1.0 * CELL
            if bad:
                nbad += 1
            worst.append((np.linalg.norm(off) / CELL, spread_mi / CELL,
                          c["size"], chain, sweep))

worst.sort(reverse=True)
offs = np.array([w[0] for w in worst])
spreads = np.array([w[1] for w in worst])
print(f"complexes checked: {ntot}  (chains x snapshots x complexes)")
print(f"stored-vs-exact centroid offset (cells): median {np.median(offs):.3f}  "
      f"p90 {np.percentile(offs, 90):.3f}  max {offs.max():.3f}")
print(f"min-image intra-cluster RMS spread (cells): median {np.median(spreads):.3f}  "
      f"max {spreads.max():.3f}")
print(f"flagged (offset>1 cell or spread>2 cells): {nbad}/{ntot}")
print("\nworst 8 by offset:")
for off, sp, size, chain, sweep in worst[:8]:
    print(f"  off={off:6.3f} cells  spread={sp:5.2f}  size={size:3d}  {chain}@{sweep}")
