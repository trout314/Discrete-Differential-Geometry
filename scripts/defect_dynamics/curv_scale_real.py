#!/usr/bin/env python3
"""Curvature-charge fluctuation vs length scale, with REAL Euclidean balls on the
metric torus (cocycle harmonic coordinates) -- the graph-BFS proxy fails here
because graph balls have rough (volume-scaling) boundaries.

For random centers and radius R (min-image Euclidean), Q(R) = sum_{|x-c|<R} q_v.
  Var_Q(R)  : extensive charge variance.   class-I HU ~ R^2 (surface); Poisson ~ R^3.
  RMS_dens  = sqrt(Var_Q)/((4/3)pi R^3).   class-I HU ~ R^-2;        Poisson ~ R^-1.5.
Charge-shuffle (permute q over fixed positions) is the built-in Poisson null:
it keeps the point pattern but randomizes the charge arrangement -> slope 3 / -1.5.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "python")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS

SP = sys.argv[1]
RNG = np.random.default_rng(0)


def var_Q(frac, basis, q, Rgrid, ncent):
    """Var over centers of the charge in a min-image Euclidean ball of radius R."""
    V = len(q)
    cen = RNG.choice(V, size=min(ncent, V), replace=False)
    Q = np.empty((len(cen), len(Rgrid)))
    for j, c in enumerate(cen):
        df = frac - frac[c]
        df -= np.round(df)                    # min image (fractional)
        dist = np.linalg.norm(df @ basis, axis=1)
        o = np.argsort(dist)
        qs = np.cumsum(q[o])
        idx = np.searchsorted(dist[o], Rgrid, side="right")
        Q[j] = np.where(idx > 0, qs[np.clip(idx - 1, 0, V - 1)], 0.0)
    return Q.var(axis=0)


def state_curves(paths, ncent):
    vobs, vshuf, Rg = [], [], None
    for p in paths:
        fac = np.asarray(Manifold.load(p, 3).facets())
        edges, omega, _ = coc.load_cocycle(os.path.splitext(p)[0] + ".cocycle.npz")
        frac, basis = coc.torus_positions(fac, edges, omega)
        q = FIELDS["curvature_charge"](fac)
        periods = np.linalg.norm(basis, axis=1)
        rho = len(q) / abs(np.linalg.det(basis))
        Rmin = (25 / (rho * 4 * np.pi / 3)) ** (1 / 3.0)      # ~25 verts in smallest ball
        Rmax = 0.45 * periods.min()
        Rg = np.geomspace(Rmin, Rmax, 18)
        vobs.append(var_Q(frac, basis, q, Rg, ncent))
        vshuf.append(var_Q(frac, basis, RNG.permutation(q), Rg, ncent))
    return Rg, np.mean(vobs, axis=0), np.mean(vshuf, axis=0)


def slopes(R, vq):
    """Var_Q~R^p (p: 2=surface/HU, 3=volume/Poisson); RMS-density slope = p/2 - 3."""
    s = vq > 0
    p = np.polyfit(np.log(R[s]), np.log(vq[s]), 1)[0]
    return p, p / 2 - 3


STATES = [
    ("perfect crystal",            [f"{SP}/r_m4.mfd"],                                "C0", 300),
    ("constrained (above-native)", sorted(glob.glob(f"{SP}/run5h/above1*_snap*.mfd")), "C2", 120),
    ("constrained (below-native)", sorted(glob.glob(f"{SP}/run5h/below*_snap*.mfd")),   "C3", 120),
    # melt dropped: its cocycle harmonic embedding collapses (no crystalline
    # connectivity), so real-space balls are meaningless. Non-HU shown via S(k)=1.95.
]

fig, ax = plt.subplots(figsize=(7.4, 5.6))
print(f"{'state':28s} {'Var_Q~R^p':>10s} {'RMS~R^s':>8s}   null(shuffle): p, s")
anchor = None
for name, paths, col, nc in STATES:
    if not paths:
        continue
    R, vobs, vshuf = state_curves(paths, nc)
    vol = (4 / 3) * np.pi * R ** 3
    rms, rms_sh = np.sqrt(vobs) / vol, np.sqrt(vshuf) / vol
    p, s = slopes(R, vobs)
    ps, ss = slopes(R, vshuf)
    print(f"{name:28s} {p:10.2f} {s:8.2f}   p={ps:.2f}, s={ss:.2f}")
    ax.loglog(R, rms, "o-", color=col, ms=4, lw=1.5, label=f"{name}  (slope {s:+.2f})")
    ax.loglog(R, rms_sh, ":", color=col, lw=1.0, alpha=0.7)
    if name.startswith("perfect"):
        anchor = (R, rms_sh)

# reference slopes
Rr, yr = anchor
Rref = np.array([Rr[3], Rr[-1]])
for slope, lab, ls in [(-2.0, r"$R^{-2}$ (class-I / Hamiltonian constraint)", "--"),
                       (-1.5, r"$R^{-3/2}$ (Poisson / CLT)", "-.")]:
    y0 = 1.3 * yr[3]
    ax.loglog(Rref, y0 * (Rref / Rref[0]) ** slope, ls, color="k", lw=1.2, label=lab)

ax.set_xlabel(r"ball radius $R$  (metric-torus units)")
ax.set_ylabel(r"RMS of mean curvature in ball  $\sqrt{\mathrm{Var}\,Q(R)}\,/\,V(R)$")
ax.set_title("Curvature-charge fluctuation vs length scale (real Euclidean balls)\n"
             "dotted = charge-shuffle null (positions kept, charges randomized)")
ax.legend(fontsize=8, loc="lower left")
ax.grid(alpha=0.25, which="both")
fig.tight_layout()
out = f"{SP}/curv_fluctuation_scale_real.png"
fig.savefig(out, dpi=130)
print(f"\nSaved to: {out}")
