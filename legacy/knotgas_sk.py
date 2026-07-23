#!/usr/bin/env python3
"""Is the sub-Bragg S(k) plateau the Poisson floor of the frozen knot gas?

Point-charge model: each complex c -> charge Q_c at its harmonic centroid.
  * predicted plateau if knots were Poisson-placed: sum_c Q_c^2 / S2,
    S2 = sum_v dq_v^2 (the same normalization as the measured S(k));
  * knot-gas structure factor S_knot(k) = |sum_c Q_c e^{ik.x_c}|^2 / sum Q_c^2
    at the commensurate sub-Bragg modes (k <= 5): ~1 = Poisson arrangement,
    < 1 = the frozen arrangement is itself anticorrelated.
Compare predicted plateau * S_knot against the measured 0.013-0.025.
"""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

_ROOT = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
for p in ("python", "scripts"):
    sys.path.insert(0, os.path.join(_ROOT, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from dopant_pairs import vertex_classes

SP = sys.argv[1]
NMAX, KCUT = 8, 5.0

plate_pred, sknot_low, plate_meas = [], [], []
for snap in sorted(glob.glob(f"{SP}/run5h/*_snap*.mfd")):
    fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
    qR = FIELDS["curvature_charge"](fac)
    dq = qR - qR.mean()
    S2 = float(dq @ dq)
    n6, imp, adj = vertex_classes(fac)
    qbar = float(np.median(qR))
    edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    P = np.abs(np.diag(basis))

    illv = [v for v in range(len(n6)) if imp[v] > 0]
    seen, Qs, cens = set(), [], []
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
        X = frac[comp] * P
        d = X - X[0]; d -= np.round(d / P) * P
        cens.append(((X[0] + d.mean(0)) / P) % 1.0)      # fractional centroid
        Qs.append(float(qR[comp].sum() - len(comp) * qbar))
    Qs = np.array(Qs); cens = np.array(cens)

    plate_pred.append(float((Qs ** 2).sum()) / S2)

    # knot-gas S(k) on the commensurate sub-Bragg grid
    rng = np.arange(-NMAX, NMAX + 1)
    nvec = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"), -1).reshape(-1, 3)
    nvec = nvec[np.any(nvec != 0, 1)]
    keep = ((nvec[:, 0] > 0) | ((nvec[:, 0] == 0) & (nvec[:, 1] > 0))
            | ((nvec[:, 0] == 0) & (nvec[:, 1] == 0) & (nvec[:, 2] > 0)))
    nvec = nvec[keep]
    kmag = np.linalg.norm(nvec, axis=1)
    nvec = nvec[kmag <= KCUT]
    A = (Qs[:, None] * np.exp(2j * np.pi * cens @ nvec.T)).sum(0)
    sknot_low.append(float(np.mean(np.abs(A) ** 2 / (Qs ** 2).sum())))

    # measured charge S(k) plateau on the same modes (for like-for-like)
    ph = np.exp(2j * np.pi * (frac @ nvec.T))
    plate_meas.append(float(np.mean(np.abs(dq @ ph) ** 2) / S2))

plate_pred, sknot_low, plate_meas = map(np.array, (plate_pred, sknot_low, plate_meas))
print(f"snapshots: {len(plate_pred)}")
print(f"Poisson-floor prediction  sum Q^2 / S2 : {plate_pred.mean():.4f} "
      f"+/- {plate_pred.std()/np.sqrt(len(plate_pred)):.4f}")
print(f"knot-gas S(k) sub-Bragg  (1=Poisson)   : {sknot_low.mean():.3f} "
      f"+/- {sknot_low.std()/np.sqrt(len(sknot_low)):.3f}")
print(f"point-model plateau (pred x S_knot)     : "
      f"{(plate_pred*sknot_low).mean():.4f}")
print(f"measured charge plateau (same modes)    : {plate_meas.mean():.4f} "
      f"+/- {plate_meas.std()/np.sqrt(len(plate_meas)):.4f}")
print(f"ratio measured / point-model            : "
      f"{plate_meas.mean() / (plate_pred*sknot_low).mean():.2f}")
