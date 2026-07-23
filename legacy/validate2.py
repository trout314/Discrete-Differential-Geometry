#!/usr/bin/env python3
"""Two sharper centroid-reliability checks.

A. Per-snapshot residual: remove each snapshot's GLOBAL offset (median stored-
   minus-exact over its complexes, mod lattice) -> residual scatter is the real
   per-complex inconsistency of the stored centroids.
B. Empirical: within ts.jsonl, |d centroid| between consecutive frames for
   complexes with IDENTICAL member sets (physically unmoved by construction).
   Distribution should be ~0; any large values = tree-lift gauge jumps.
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
CELL = 1e6
P4 = 4e6


def load_ts(f):
    out = {}
    for line in open(f):
        r = json.loads(line)
        out.setdefault(r["sweep"], r)
    return out


# ---- A: residual after per-snapshot global offset removal ----
resid = []
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
        d_of = {int(v): i for i, v in enumerate(lab)}
        P = np.abs(np.diag(basis))
        offs = []
        for c in ts[sweep]["comps"]:
            mem = [d_of[v] for v in c["members"] if v in d_of]
            if len(mem) != len(c["members"]) or c["centroid"] is None:
                continue
            X = frac[mem] * np.diag(basis)
            ref = X[0]; d = X - ref; d -= np.round(d / P) * P
            cen = ref + d.mean(0)
            off = np.asarray(c["centroid"]) - cen
            off -= np.round(off / P) * P
            offs.append(off)
        if len(offs) >= 3:
            offs = np.array(offs)
            g = np.median(offs, axis=0)              # global gauge offset
            r = offs - g
            r -= np.round(r / P) * P
            resid += list(np.linalg.norm(r, axis=1) / CELL)
resid = np.array(resid)
print("A. per-complex residual after removing snapshot-global offset (cells):")
print(f"   median {np.median(resid):.3f}  p90 {np.percentile(resid,90):.3f}  "
      f"max {resid.max():.3f}   n={len(resid)}")

# ---- B: frame-to-frame |dc| for identical member sets ----
dstill = []
for tsf in sorted(glob.glob(f"{SP}/run5h/*.ts.jsonl")):
    recs = sorted(load_ts(tsf).items())
    for (s1, r1), (s2, r2) in zip(recs[:-1], recs[1:]):
        m2 = {frozenset(c["members"]): c["centroid"] for c in r2["comps"]}
        for c in r1["comps"]:
            key = frozenset(c["members"])
            if key in m2 and c["centroid"] and m2[key]:
                d = np.asarray(m2[key]) - np.asarray(c["centroid"])
                d -= np.round(d / P4) * P4
                dstill.append(np.linalg.norm(d) / CELL)
dstill = np.array(dstill)
print("\nB. frame-to-frame |dcentroid| for IDENTICAL member sets (cells):")
print(f"   n={len(dstill)}  zero-frac={np.mean(dstill < 1e-9):.3f}  "
      f"median {np.median(dstill):.4f}  p99 {np.percentile(dstill,99):.4f}  "
      f"max {dstill.max():.3f}")
print(f"   frac > 0.25 cell = {np.mean(dstill > 0.25):.4f}   "
      f"frac > 1 cell = {np.mean(dstill > 1.0):.4f}")
