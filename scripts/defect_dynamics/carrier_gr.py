#!/usr/bin/env python3
"""Carrier pair correlation g_cc(r) and effective potential u_eff(r) = -ln g
from the lam=0.40 snapshot set. Euclidean min-image distances between complex
centroids (harmonic coords); null = uniform random points in the same torus
(Monte Carlo, exact for the min-image metric). Decides: short-ranged u_eff
(no constraint-induced long-range interaction) vs slow tail (weak plasma,
kappa below the box window).

FRAME (fixed 2026-08-13): distances go through the FULL lattice basis via
`coc.min_image`. The previous version used `frac * np.diag(basis)` and
min-imaged per axis, which on the R m4 host (basis [[4,4,0],[0,4,0],[0,0,4]],
not diagonal) mis-stated 68% of pair distances by >10% and 18% by >30%, range
0.62-1.62x -- i.e. it smeared g(r) rather than shifting it. The MC null also
hardcoded a 4x4x4 cube; it now samples the true cell.
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
        B = basis / CELL                       # rows = lattice vectors, cells
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
        # centroids stay FRACTIONAL; the metric enters only in min_image
        cens = []
        for c in comps:
            df = frac[c] - frac[c][0]
            df -= np.round(df)
            cens.append((frac[c][0] + df.mean(0)) % 1.0)
        cens = np.array(cens)
        n = len(cens)
        per_snap_n.append((n, B))
        iu = np.triu_indices(n, 1)
        _, rr = coc.min_image(cens[iu[0]] - cens[iu[1]], B)
        rr = rr[rr < RMAX]
        obs += np.histogram(rr, bins=bins)[0]
        npairs += len(rr)
        nsnap += 1

# Beyond half the shortest lattice vector the min-image distance stops being a
# faithful separation (the periodic images interfere), so bins past it are not
# a pair correlation at all -- flagged rather than silently plotted.
B0 = per_snap_n[0][1]
offs = np.array([[i, j, k] for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2)
                 for k in (-2, -1, 0, 1, 2) if (i, j, k) != (0, 0, 0)], float)
r_valid = 0.5 * np.linalg.norm(offs @ B0, axis=1).min()

# MC null: same per-snapshot carrier counts, uniform in the TRUE cell
null = np.zeros(len(bins) - 1)
NMC = 400
for n, Bs in per_snap_n:
    iu = np.triu_indices(n, 1)
    for _ in range(NMC):
        pts = rng.uniform(0, 1, (n, 3))
        _, rr = coc.min_image(pts[iu[0]] - pts[iu[1]], Bs)
        null += np.histogram(rr[rr < RMAX], bins=bins)[0]
null /= NMC

print(f"snapshots {nsnap}, carrier pairs (r<{RMAX}) {npairs}, "
      f"<n_carriers> {np.mean([n for n, _ in per_snap_n]):.1f}")
print(f"cell (cells) {np.round(B0, 3).tolist()}, volume "
      f"{abs(np.linalg.det(B0)):.2f}; min-image valid to r < {r_valid:.2f}")
print(f"\n{'r (cells)':>10s} {'g(r)':>7s} {'+/-':>6s} {'u_eff=-ln g':>12s}")
for k in range(len(bins) - 1):
    rc = (bins[k] + bins[k + 1]) / 2
    if null[k] > 0:
        g = obs[k] / null[k]
        eg = np.sqrt(max(obs[k], 1)) / null[k]
        u = -np.log(g) if g > 0 else np.inf
        flag = "" if bins[k + 1] <= r_valid else "   (beyond min-image radius)"
        print(f"{rc:10.2f} {g:7.3f} {eg:6.3f} {u:12.3f}{flag}")
