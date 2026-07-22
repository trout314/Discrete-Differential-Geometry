#!/usr/bin/env python3
"""S_66(k): structure factor of the six-edge (disclination) web -- the
emergent-gauge / Coulomb-phase diagnostic.

Fields on edge midpoints (harmonic torus coords, min-image):
  scalar : s_e = 1 for deg-6 edges (six-density channel)
  tensor : Q_e^{ab} = nhat_e^a nhat_e^b - delta/3 on six-edges (director field;
           sign-free, correct for UNORIENTED lines)
Fourier at sub-Bragg commensurate k; helicity-decompose the tensor channel into
TT / vector / longitudinal exactly as in tt_channel. Null = random edge subset
of the same size (occupation shuffle over the fixed edge skeleton).

Coulomb-phase (closed-web) prediction: web closure = lattice Gauss law =>
coarse-grained six-flux is divergence-free => LONGITUDINAL channel suppressed
at small k while transverse stays thermal (the pinch-point projector), with the
suppression smeared over a width ~ sqrt(monopole density). Monopole-dense
states -> L/T ratio ~ 1.

Usage: web_s66.py LABEL:snapshot.mfd [LABEL:...]   (needs sibling .cocycle.npz)
"""
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

NMAX, KCUT = 8, 5.0
NSHUF = 8
rng = np.random.default_rng(0)


def kgrid():
    r = np.arange(-NMAX, NMAX + 1)
    nv = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    nv = nv[np.any(nv != 0, 1)]
    keep = ((nv[:, 0] > 0) | ((nv[:, 0] == 0) & (nv[:, 1] > 0))
            | ((nv[:, 0] == 0) & (nv[:, 1] == 0) & (nv[:, 2] > 0)))
    nv = nv[keep]
    km = np.linalg.norm(nv, axis=1)
    m = km <= KCUT
    return nv[m], km[m]


NV, KM = kgrid()
KHAT = None
SHELLS = [(1.0, 2.0), (2.0, 3.0), (3.0, 5.0)]


def channels(mid_frac, nhat, sel):
    """(scalar, TT, V, L) powers of the six-web fields for edge subset sel."""
    ph = np.exp(2j * np.pi * (mid_frac[sel] @ NV.T))          # (E6, M)
    S = np.abs(ph.sum(0)) ** 2
    t = nhat[sel][:, :, None] * nhat[sel][:, None, :] - np.eye(3) / 3
    Qk = np.einsum("eab,em->mab", t, ph)
    P = np.eye(3)[None] - KHAT[:, :, None] * KHAT[:, None, :]
    PQP = np.einsum("mac,mcd,mdb->mab", P, Qk, P)
    PtrQ = np.einsum("mcd,mcd->m", P, Qk)
    LQ = PQP - 0.5 * P * PtrQ[:, None, None]
    tot = np.einsum("mab,mab->m", np.conj(Qk), Qk).real
    TT = np.einsum("mab,mab->m", np.conj(Qk), LQ).real
    Lg = 1.5 * np.abs(np.einsum("ma,mab,mb->m", KHAT, Qk, KHAT)) ** 2
    Vc = tot - TT - Lg
    return S, TT, Vc, Lg


for arg in sys.argv[1:]:
    label, path = arg.split(":", 1)
    fac = np.asarray(ddg.Manifold.load(path, 3).facets())
    eu, edeg, Vn = edges_from_facets(fac)          # DENSE labels
    lab = np.unique(fac)
    n6v, imp, adj = vertex_classes(fac)
    nmono = int((imp > 0).sum())
    edges, omega, _ = coc.load_cocycle(path[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    # edge midpoints + directions in FRACTIONAL coords (min image);
    # NOTE eu is dense-labelled and frac rows are np.unique(fac)-ordered = dense ✓
    d = frac[eu[:, 1]] - frac[eu[:, 0]]
    d -= np.round(d)
    mid = (frac[eu[:, 0]] + 0.5 * d) % 1.0
    dc = d @ basis
    nhat = dc / np.linalg.norm(dc, axis=1, keepdims=True)
    KHAT_local = (NV @ np.linalg.inv(basis).T)
    KHAT = KHAT_local / np.linalg.norm(KHAT_local, axis=1, keepdims=True)

    six = np.nonzero(edeg == 6)[0]
    f6 = len(six) / len(eu)
    So, TTo, Vo, Lo = channels(mid, nhat, six)
    acc = None
    for _ in range(NSHUF):
        r_ = rng.choice(len(eu), size=len(six), replace=False)
        vals = channels(mid, nhat, r_)
        acc = vals if acc is None else tuple(a + b for a, b in zip(acc, vals))
    Ss, TTs, Vs, Ls = (a / NSHUF for a in acc)

    print(f"\n═══ {label}: E={len(eu)} six-frac={f6:.3f} monopoles(nill)={nmono} ═══")
    print(f"{'shell':>7s} {'scalar':>8s} {'TT':>8s} {'vector':>8s} {'longit':>8s}"
          f" {'L/T':>7s}")
    for lo, hi in SHELLS:
        m = (KM >= lo) & (KM < hi)
        sc = So[m].mean() / Ss[m].mean()
        tt = TTo[m].mean() / TTs[m].mean()
        vv = Vo[m].mean() / Vs[m].mean()
        ll = Lo[m].mean() / Ls[m].mean()
        print(f"{lo:.0f}-{hi:.0f} {sc:8.3f} {tt:8.3f} {vv:8.3f} {ll:8.3f}"
              f" {ll/max(tt,1e-12):7.2f}")
    print("  (ratios to occupation-shuffle null; Coulomb phase => longit "
          "suppressed, L/T << 1 at small k)")
