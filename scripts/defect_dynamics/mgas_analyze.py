#!/usr/bin/env python3
"""Analyze a mobile-gas campaign chain (mobile_gas.py output):

  1. mobility at m4: illegal-set Jaccard vs lag (decay -> gas rearranges);
  2. carrier census: complex-size spectrum (still quantized 5/11/...?),
     monopole per complex vs CELL-MEAN reference (the correct far-field ref);
  3. S_knot(k): number/charge structure factor of complex centroids, pooled
     over snapshots -- the Stillinger-Lovett test (0.98 = frozen Poisson;
     sub-Poisson at small k = collective screening has set in);
  4. charge S(k) plateau vs the bare monopole floor (mean-referenced);
  5. charge budget: achieved <edeg> vs target; carrier monopole total vs the
     pin demand theta*E*(et - native).

Usage: mgas_analyze.py OUTPREFIX [kcut]
"""
import glob
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

THETA = float(np.arccos(1.0 / 3.0))
PRE = sys.argv[1]
KCUT = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
NMAX = 8

# ---- 1. mobility from ts.jsonl ----
rows = [json.loads(l) for l in open(PRE + ".ts.jsonl")]
sets = [frozenset(x for c in r["members"] for x in c) for r in rows]
sw = [r["sweep"] for r in rows]
DT = sw[1] - sw[0]
print(f"{os.path.basename(PRE)}: {len(rows)} frames, sweeps {sw[0]}..{sw[-1]}")
print(f"  <nill>={np.mean([r['n_illegal'] for r in rows]):.1f}  "
      f"<ncomp>={np.mean([len(r['sizes']) for r in rows]):.1f}  "
      f"legaledge={np.mean([r['legaledge'] for r in rows]):.4f}  "
      f"legalvert={np.mean([r['legalvert'] for r in rows]):.4f}")
Vtot = 10176
base = np.mean([r["n_illegal"] for r in rows]) / (2 * Vtot)
print("  Jaccard vs lag: ", end="")
for lag in (150, 600, 1500, 3000, 6000, 12000):
    st = lag // DT
    vals = [len(sets[i] & sets[i + st]) / max(len(sets[i] | sets[i + st]), 1)
            for i in range(0, len(sets) - st, max(1, st // 2))]
    if vals:
        print(f"J{lag}={np.mean(vals):.2f}", end=" ")
print(f"(rand~{base:.3f})")

# ---- 2-4. snapshot-based ----
snaps = sorted(glob.glob(PRE + "_snap*.mfd"))
size_hist = Counter()
mono = []
plate_meas, floor_mean, sknotN, sknotQ = [], [], [], []
edeg_ach = []
for snap in snaps:
    fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
    qR = FIELDS["curvature_charge"](fac)
    dq = qR - qR.mean()
    S2 = float(dq @ dq)
    eu, edeg, V = edges_from_facets(fac)
    edeg_ach.append(float(edeg.mean()))
    n6, imp, adj = vertex_classes(fac)
    edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    P = np.abs(np.diag(basis))
    illv = [i for i in range(V) if imp[i] > 0]
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
        comps.append(comp)
    if not comps:
        continue
    for c in comps:
        size_hist[len(c)] += 1
    Qs = np.array([float(qR[c].sum() - len(c) * qR.mean()) for c in comps])
    mono += [(len(c), q) for c, q in zip(comps, Qs)]
    cens = []
    for c in comps:
        X = frac[c] * P
        d = X - X[0]; d -= np.round(d / P) * P
        cens.append(((X[0] + d.mean(0)) / P) % 1.0)
    cens = np.array(cens)
    r = np.arange(-NMAX, NMAX + 1)
    nv = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    nv = nv[np.any(nv != 0, 1)]
    keep = ((nv[:, 0] > 0) | ((nv[:, 0] == 0) & (nv[:, 1] > 0))
            | ((nv[:, 0] == 0) & (nv[:, 1] == 0) & (nv[:, 2] > 0)))
    nv = nv[keep]
    km = np.linalg.norm(nv, axis=1)
    nvl = nv[km <= 2.0]                     # LOW-k modes for S_knot headline
    nva = nv[km <= KCUT]
    ph = np.exp(2j * np.pi * (frac @ nva.T))
    plate_meas.append(float(np.mean(np.abs(dq @ ph) ** 2) / S2))
    floor_mean.append(float((Qs ** 2).sum()) / S2)
    for nvx, acc in ((nvl, None),):
        AQ = (Qs[:, None] * np.exp(2j * np.pi * cens @ nvx.T)).sum(0)
        AN = np.exp(2j * np.pi * cens @ nvx.T).sum(0)
        n = len(comps)
        sknotQ.append(float(np.mean(np.abs(AQ) ** 2 / (Qs ** 2).sum())))
        sknotN.append(float(np.mean((np.abs(AN) ** 2 - 0) / n)))

print(f"\n  snapshots analyzed: {len(snaps)}")
print(f"  carrier size spectrum: "
      + " ".join(f"{k}x{v}" for k, v in sorted(size_hist.items())))
if mono:
    m5 = [q for n, q in mono if n == 5]
    print(f"  monopole (cell-mean ref): all sizes <Q/N>="
          f"{np.mean([q/n for n, q in mono]):+.3f}; "
          f"5-knots <Q>={np.mean(m5):+.2f} (n={len(m5)})" if m5 else "")
print(f"\n  S_knot NUMBER (|k|<=2): {np.mean(sknotN):.3f} +/- "
      f"{np.std(sknotN)/max(np.sqrt(len(sknotN)),1):.3f}   (1=Poisson; <1=HU onset)")
print(f"  S_knot CHARGE (|k|<=2): {np.mean(sknotQ):.3f} +/- "
      f"{np.std(sknotQ)/max(np.sqrt(len(sknotQ)),1):.3f}")
print(f"  charge plateau (k<={KCUT:.0f}): {np.mean(plate_meas):.4f}   "
      f"bare monopole floor: {np.mean(floor_mean):.4f}   "
      f"ratio {np.mean(plate_meas)/max(np.mean(floor_mean),1e-12):.2f}")

# ---- 5. charge budget ----
ref = ddg.Manifold.load(os.path.join(os.path.dirname(PRE), "..", "r_m4.mfd"), 3)
native = float(edges_from_facets(ref.facets())[1].mean())
E = len(eu)
et_row = json.loads(open(PRE + ".ts.jsonl").readline())
ach = np.mean(edeg_ach)
demand = -THETA * E * (ach - native)          # deficit change actually realized
mono_tot = np.mean([sum(q for _, q in mono[i::len(snaps)]) for i in range(1)]) \
    if mono else 0.0
per_snap_mono = (sum(q for _, q in mono) / max(len(snaps), 1))
print(f"\n  charge budget: native={native:.6f} achieved <edeg>={ach:.6f} "
      f"(target {native + 8e-4:.6f})")
print(f"    realized deficit change  : {demand:+.1f} rad")
print(f"    carrier monopole total   : {per_snap_mono:+.1f} rad/snapshot")
print(f"    -> background carries    : {demand - per_snap_mono:+.1f} rad")
