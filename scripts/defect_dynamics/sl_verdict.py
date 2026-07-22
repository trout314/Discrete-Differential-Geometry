#!/usr/bin/env python3
"""The Stillinger-Lovett verdict: shell-resolved carrier structure factors over
all lam=0.40 snapshots (5 chains x 6). For each snapshot: carrier complexes,
harmonic centroids, monopoles (cell-mean ref); pooled per k-shell:
  S_N(k) = <|sum_j e^{ik x_j}|^2>/n  - number structure factor (1 = Poisson)
  S_Q(k) = <|sum_j Q_j e^{ik x_j}|^2>/sum Q^2  - charge-weighted
Jackknife over chains for errors. Also the charge plateau / bare floor ratio
per shell (medium screening), pooled.
"""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
for p in ("python", "scripts"):
    sys.path.insert(0, os.path.join(_R, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from dopant_pairs import vertex_classes

SP = sys.argv[1]
NMAX = 8
SHELLS = [(1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0)]

r = np.arange(-NMAX, NMAX + 1)
NV = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
NV = NV[np.any(NV != 0, 1)]
keep = ((NV[:, 0] > 0) | ((NV[:, 0] == 0) & (NV[:, 1] > 0))
        | ((NV[:, 0] == 0) & (NV[:, 1] == 0) & (NV[:, 2] > 0)))
NV = NV[keep]
KM = np.linalg.norm(NV, axis=1)
MASKS = [(KM >= lo) & (KM < hi) for lo, hi in SHELLS]

CHAINS = ["lam40", "l40s201", "l40s202", "l40s203", "l40s204"]
per_chain = {c: defaultdict(list) for c in CHAINS}   # metric -> per-snapshot vals

for chain in CHAINS:
    for snap in sorted(glob.glob(f"{SP}/mgas/{chain}_snap*.mfd")):
        fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
        qR = FIELDS["curvature_charge"](fac)
        dq = qR - qR.mean()
        S2 = float(dq @ dq)
        n6, imp, adj = vertex_classes(fac)
        edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
        frac, basis = coc.torus_positions(fac, edges, omega)
        P = np.abs(np.diag(basis))
        illv = np.nonzero(imp > 0)[0]
        seen, comps = set(), []
        for s0 in illv:
            s0 = int(s0)
            if s0 in seen:
                continue
            st, comp = [s0], []
            seen.add(s0)
            while st:
                u = st.pop(); comp.append(u)
                for w in adj[u]:
                    if imp[w] > 0 and w not in seen:
                        seen.add(int(w)); st.append(int(w))
            comps.append(comp)
        if len(comps) < 3:
            continue
        Qs = np.array([float(qR[c].sum() - len(c) * qR.mean()) for c in comps])
        cens = []
        for c in comps:
            X = frac[c] * P
            d = X - X[0]; d -= np.round(d / P) * P
            cens.append(((X[0] + d.mean(0)) / P) % 1.0)
        cens = np.array(cens)
        n = len(comps)
        ph = np.exp(2j * np.pi * (cens @ NV.T))
        AN = ph.sum(0)
        AQ = (Qs[:, None] * ph).sum(0)
        phv = np.exp(2j * np.pi * (frac @ NV.T))
        Adq = dq @ phv
        for si, m in enumerate(MASKS):
            per_chain[chain][("SN", si)].append(float(np.mean(np.abs(AN[m])**2)) / n)
            per_chain[chain][("SQ", si)].append(
                float(np.mean(np.abs(AQ[m])**2)) / float((Qs**2).sum()))
            per_chain[chain][("plate", si)].append(
                float(np.mean(np.abs(Adq[m])**2)) / S2)
            per_chain[chain][("floor", si)].append(float((Qs**2).sum()) / S2)

def jack(metric, si):
    ch_means = [np.mean(per_chain[c][(metric, si)]) for c in CHAINS
                if per_chain[c][(metric, si)]]
    nc = len(ch_means)
    full = np.mean(ch_means)
    loo = [np.mean([m for j, m in enumerate(ch_means) if j != i])
           for i in range(nc)]
    se = np.sqrt((nc - 1) / nc * np.sum((np.array(loo) - np.mean(loo))**2))
    return full, se

nsnap = sum(len(per_chain[c][("SN", 0)]) for c in CHAINS)
print(f"pooled snapshots: {nsnap} over {len(CHAINS)} chains "
      f"(jackknife over chains)\n")
print(f"{'shell':>8s} {'S_N (number)':>16s} {'S_Q (charge)':>16s} "
      f"{'plateau/floor':>14s}")
for si, (lo, hi) in enumerate(SHELLS):
    sn, esn = jack("SN", si)
    sq, esq = jack("SQ", si)
    pl, _ = jack("plate", si)
    fl, _ = jack("floor", si)
    print(f"{lo:.1f}-{hi:.1f} {sn:9.3f}+/-{esn:.3f} {sq:9.3f}+/-{esq:.3f} "
          f"{pl/fl:14.2f}")
print("\nS=1 Poisson; S<1 = correlation holes (Stillinger-Lovett onset); "
      "plateau/floor <1 = medium screening")
