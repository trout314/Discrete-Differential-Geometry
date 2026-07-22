#!/usr/bin/env python3
"""Is the real-space charge surface law (R^-2) genuine, or a boundary artifact?

Decisive control: for the mean-subtracted charge dq, compare OBSERVED vs its
own SHUFFLE (dq permuted over the fixed positions = Poisson floor). If observed
~ R^-2 while shuffle ~ R^-1.5, the surface law is a real anticorrelation
(hyperuniform charge), not an artifact. Error bars = SEM over the 4 chains.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, "python")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS

SP = sys.argv[1]
RNG = np.random.default_rng(0)


def var_Q(frac, basis, q, Rgrid, ncent):
    V = len(q)
    cen = RNG.choice(V, size=min(ncent, V), replace=False)
    Q = np.empty((len(cen), len(Rgrid)))
    for j, c in enumerate(cen):
        df = frac - frac[c]; df -= np.round(df)
        dist = np.linalg.norm(df @ basis, axis=1)
        o = np.argsort(dist); qs = np.cumsum(q[o])
        idx = np.searchsorted(dist[o], Rgrid, side="right")
        Q[j] = np.where(idx > 0, qs[np.clip(idx - 1, 0, V - 1)], 0.0)
    return Q.var(axis=0)


def rms_slope(R, vq):
    s = vq > 0
    p = np.polyfit(np.log(R[s]), np.log(vq[s]), 1)[0]
    return p / 2 - 3               # RMS-density ~ R^(p/2 - 3)


def chain_slopes(paths):
    raw, ms, sh = [], [], []
    for p in paths:
        fac = np.asarray(Manifold.load(p, 3).facets())
        edges, omega, _ = coc.load_cocycle(os.path.splitext(p)[0] + ".cocycle.npz")
        frac, basis = coc.torus_positions(fac, edges, omega)
        q = FIELDS["curvature_charge"](fac)
        dq = q - q.mean()
        periods = np.linalg.norm(basis, axis=1)
        rho = len(q) / abs(np.linalg.det(basis))
        Rg = np.geomspace((25 / (rho * 4 * np.pi / 3)) ** (1 / 3.0),
                          0.45 * periods.min(), 18)
        raw.append(rms_slope(Rg, var_Q(frac, basis, q, Rg, 120)))
        ms.append(rms_slope(Rg, var_Q(frac, basis, dq, Rg, 120)))
        sh.append(rms_slope(Rg, var_Q(frac, basis, RNG.permutation(dq), Rg, 120)))
    return np.mean(raw), np.mean(ms), np.mean(sh)


chains = {i: sorted(glob.glob(f"{SP}/run5h/above1{i}_snap*.mfd")) for i in range(4)}
raw, ms, sh = [], [], []
for i, paths in chains.items():
    r, m, s = chain_slopes(paths)
    raw.append(r); ms.append(m); sh.append(s)

for name, arr, pred in [("raw Q (incl. point pattern)", raw, "-2 surface"),
                        ("mean-subtracted OBSERVED   ", ms, "-2 if HU charge"),
                        ("mean-subtracted SHUFFLE null", sh, "-1.5 Poisson")]:
    a = np.array(arr)
    print(f"  {name}: RMS-density slope = {a.mean():+.3f} +/- {a.std(ddof=1)/np.sqrt(len(a)):.3f}"
          f"   (pred {pred})")
