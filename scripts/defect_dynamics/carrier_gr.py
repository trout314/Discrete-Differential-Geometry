#!/usr/bin/env python3
"""Carrier pair correlation g_cc(r) and effective potential u_eff(r) = -ln g
from the lam=0.40 snapshot set. Euclidean min-image distances between complex
centroids (harmonic coords); null = uniform random points in the same torus
(Monte Carlo, exact for the min-image metric). Decides: short-ranged u_eff
(no constraint-induced long-range interaction) vs slow tail (weak plasma,
kappa below the box window).
"""
import glob
import os
import sys

import numpy as np

_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
for p in ("python", "scripts"):
    sys.path.insert(0, os.path.join(_R, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from dopant_pairs import vertex_classes

SP = sys.argv[1]
CELL = 1e6
RMAX, DR = 2.0, 0.2          # cells; box is 4 cells/side
rng = np.random.default_rng(0)

bins = np.arange(0.0, RMAX + DR, DR)
obs = np.zeros(len(bins) - 1)
obs_pp = np.zeros(len(bins) - 1)     # same-sign (all negative) weighting: Q_i Q_j > 0 pairs
nsnap = 0
npairs = 0
per_snap_n = []

for chain in ("lam40", "l40s201", "l40s202", "l40s203", "l40s204"):
    for snap in sorted(glob.glob(f"{SP}/mgas/{chain}_snap*.mfd")):
        fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
        qR = FIELDS["curvature_charge"](fac)
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
        if len(comps) < 2:
            continue
        cens = []
        for c in comps:
            X = frac[c] * P
            d = X - X[0]; d -= np.round(d / P) * P
            cens.append((X[0] + d.mean(0)) / CELL)          # cells
        cens = np.array(cens)
        L = P / CELL
        n = len(cens)
        per_snap_n.append(n)
        for i in range(n):
            for j in range(i + 1, n):
                d = cens[i] - cens[j]
                d -= np.round(d / L) * L
                r = np.linalg.norm(d)
                if r < RMAX:
                    k = int(r / DR)
                    obs[k] += 1
                    npairs += 1
        nsnap += 1

# MC null: same per-snapshot carrier counts, uniform positions
null = np.zeros(len(bins) - 1)
NMC = 400
L = np.array([4.0, 4.0, 4.0])
for n in per_snap_n:
    for _ in range(NMC):
        pts = rng.uniform(0, 1, (n, 3)) * L
        d = pts[:, None, :] - pts[None, :, :]
        d -= np.round(d / L) * L
        r = np.linalg.norm(d, axis=2)
        iu = np.triu_indices(n, 1)
        rr = r[iu]
        h, _ = np.histogram(rr[rr < RMAX], bins=bins)
        null += h
null /= NMC

print(f"snapshots {nsnap}, carrier pairs (r<{RMAX}) {npairs}, "
      f"<n_carriers> {np.mean(per_snap_n):.1f}")
print(f"\n{'r (cells)':>10s} {'g(r)':>7s} {'+/-':>6s} {'u_eff=-ln g':>12s}")
for k in range(len(bins) - 1):
    rc = (bins[k] + bins[k + 1]) / 2
    if null[k] > 0:
        g = obs[k] / null[k]
        eg = np.sqrt(max(obs[k], 1)) / null[k]
        u = -np.log(g) if g > 0 else np.inf
        print(f"{rc:10.2f} {g:7.3f} {eg:6.3f} {u:12.3f}")
