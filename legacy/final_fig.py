#!/usr/bin/env python3
"""Headline figure: mean-subtracted curvature-charge fluctuation vs scale, with
the shuffle (Poisson) null and slope fits +/- SEM. Real Euclidean balls on the
metric torus; constrained above-native ensemble (16 snapshots, 4 chains)."""
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


def var_Q(frac, basis, q, Rgrid, ncent=150):
    V = len(q)
    cen = RNG.choice(V, size=ncent, replace=False)
    Q = np.empty((ncent, len(Rgrid)))
    for j, c in enumerate(cen):
        df = frac - frac[c]; df -= np.round(df)
        dist = np.linalg.norm(df @ basis, axis=1)
        o = np.argsort(dist); qs = np.cumsum(q[o])
        idx = np.searchsorted(dist[o], Rgrid, side="right")
        Q[j] = np.where(idx > 0, qs[np.clip(idx - 1, 0, V - 1)], 0.0)
    return Q.var(axis=0)


paths = sorted(glob.glob(f"{SP}/run5h/above1*_snap*.mfd"))
# fixed absolute Rgrid from the first snapshot
fac0 = np.asarray(Manifold.load(paths[0], 3).facets())
e0, o0, _ = coc.load_cocycle(os.path.splitext(paths[0])[0] + ".cocycle.npz")
_, b0 = coc.torus_positions(fac0, e0, o0)
rho = len(np.unique(fac0)) / abs(np.linalg.det(b0))
Rg = np.geomspace((25 / (rho * 4 * np.pi / 3)) ** (1 / 3.0),
                  0.45 * np.linalg.norm(b0, axis=1).min(), 18)
vol = (4 / 3) * np.pi * Rg ** 3

obs, shf = [], []
for p in paths:
    fac = np.asarray(Manifold.load(p, 3).facets())
    edges, omega, _ = coc.load_cocycle(os.path.splitext(p)[0] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    dq = FIELDS["curvature_charge"](fac); dq = dq - dq.mean()
    obs.append(var_Q(frac, basis, dq, Rg))
    shf.append(var_Q(frac, basis, RNG.permutation(dq), Rg))
obs, shf = np.array(obs), np.array(shf)

rms_o = np.sqrt(obs.mean(0)) / vol
rms_s = np.sqrt(shf.mean(0)) / vol
# per-chain slopes for error bar
def chain_slopes(V4):
    sl = []
    for i in range(4):
        v = V4[4 * i:4 * i + 4].mean(0)
        p = np.polyfit(np.log(Rg), np.log(v), 1)[0]
        sl.append(p / 2 - 3)
    return np.mean(sl), np.std(sl, ddof=1) / 2.0
so, seo = chain_slopes(obs)
ss, ses = chain_slopes(shf)

fig, ax = plt.subplots(figsize=(7.2, 5.6))
ax.loglog(Rg, rms_o, "o-", color="C0", ms=5, lw=1.7,
          label=fr"observed charge  (slope ${so:.2f}\pm{seo:.2f}$)")
ax.loglog(Rg, rms_s, "s--", color="C3", ms=4, lw=1.3,
          label=fr"shuffle null  (slope ${ss:.2f}\pm{ses:.2f}$)")
for al, ls, lab in [(-2.0, ":", r"$R^{-2}$ (class-I / surface law)"),
                    (-1.5, "-.", r"$R^{-3/2}$ (Poisson)")]:
    ax.loglog(Rg[[2, -1]], rms_o[2] * (Rg[[2, -1]] / Rg[2]) ** al * (1.4 if al == -2 else 3.2),
              ls, color="k", lw=1.1, label=lab)
ax.set_xlabel(r"ball radius $R$ (metric-torus units)")
ax.set_ylabel(r"RMS mean curvature $\sqrt{\mathrm{Var}\,\delta Q(R)}/V(R)$")
ax.set_title("Curvature-charge (mean-subtracted) fluctuation vs scale\n"
             "constrained above-native ensemble — genuine surface law vs Poisson null")
ax.legend(fontsize=8.5, loc="lower left")
ax.grid(alpha=0.25, which="both")
fig.tight_layout()
out = f"{SP}/curv_charge_surfacelaw.png"
fig.savefig(out, dpi=130)
print(f"observed slope {so:+.3f} +/- {seo:.3f} ;  shuffle {ss:+.3f} +/- {ses:.3f}")
print(f"number-variance exponent p = {2*(so+3):.2f} +/- {2*seo:.2f}  (2 = surface law)")
print(f"Saved to: {out}")
