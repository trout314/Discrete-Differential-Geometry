#!/usr/bin/env python3
"""Vector-channel (oriented-flux) structure factor of the six-edge web via the
EXACT SEGMENT INTEGRAL -- the identity-grade Gauss-law observable.

web_s66.py measures the NEMATIC (director) channel: six-edges are unoriented,
so the sign-free field is Q_e = nn - I/3, and web closure only *suppresses* its
helicity-0 power (bending leaks in; frozen L/T = 0.60, not 0).  The Gauss law
proper lives in the VECTOR channel: an oriented unit flux J with div J = 0 has
k.J(k) = 0 identically.  Two obstacles, two fixes:

1. ORIENTATION.  Unit +/-1 fluxes cannot be divergence-free at odd-valence
   junctions (Z15 Y-junctions) -- the abelian shadow of the non-abelian
   disclination algebra.  Fix: real-valued fluxes phi_e from a CYCLE-SPACE
   PROJECTION: draw random signs sigma0, then solve

       phi = argmin ||phi - sigma0||^2   s.t.   (B phi)_v = 0  at every
       legal vertex (imp == 0);  divergence is permitted only at the
       edge-illegal (monopole) vertices.

   Solved via LSMR on the normal equations (B B^T) mu = B sigma0,
   phi = sigma0 - B^T mu; rank deficiency (fully closed components) is
   handled gracefully.  Averaged over NSIG sign draws (Dirac-string gauge
   average).

2. EXACT SEGMENT KERNEL.  Fourier the flux as a line integral along each
   segment, not a midpoint point-charge:

       I_e(k) = int_0^1 dt e^{ik.(x0 + t d)}
              = e^{ik.x0} e^{i theta/2} sinc(theta/2pi),  theta = k.d,

   so that the longitudinal projection TELESCOPES exactly:

       k.J(k) = (1/i) sum_v D_v e^{ik.x_v},   D_v = net divergence at v.

   Every legal vertex drops out identically =>

       L_vec(k) = |khat.J|^2 = S_mono(k)/|k|^2

   -- the longitudinal power IS the monopole-charge structure factor over k^2.
   A closed web gives L_vec = 0 to machine precision at every k; the script
   verifies the telescoping identity numerically per state.

   Small-k reading of L*k^2 = S_mono(k): -> n_mono (flat) for a Poisson
   monopole gas; -> 0 as k->0 iff the monopole plasma itself screens
   (Stillinger-Lovett for the emergent charges).

Null baseline: the UNPROJECTED random signs sigma0 (same web geometry, no
constraint) -- its L/T is the no-Gauss-law reference, O(1).

Usage: web_vector.py LABEL:snapshot.mfd [LABEL:...]  (needs sibling .cocycle.npz)
"""
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsmr

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

NMAX, KCUT = 8, 5.0
NSIG = 4
SHELLS = [(1.0, 2.0), (2.0, 3.0), (3.0, 5.0)]
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


def segment_kernel(x0, d):
    """Exact per-segment Fourier kernel, shape (E6, M).
    I_e(k) = e^{2pi i n.x0} e^{i theta/2} sinc(theta/2pi), theta = 2pi n.d."""
    theta = 2 * np.pi * (d @ NV.T)                       # = k_cart . d_cart
    return (np.exp(2j * np.pi * (x0 @ NV.T))
            * np.exp(0.5j * theta) * np.sinc(theta / (2 * np.pi)))


def flux_channels(phi, dc, kern, khat):
    """(total, T, L) vector powers of J(k) = sum_e phi_e dc_e I_e(k)."""
    J = np.einsum("e,ea,em->ma", phi, dc, kern, optimize=True)
    L = np.abs(np.einsum("ma,ma->m", khat, J)) ** 2
    tot = np.einsum("ma,ma->m", np.conj(J), J).real
    return J, tot, tot - L, L


for arg in sys.argv[1:]:
    label, path = arg.split(":", 1)
    fac = np.asarray(ddg.Manifold.load(path, 3).facets())
    eu, edeg, Vn = edges_from_facets(fac)                # DENSE labels
    n6v, imp, adj = vertex_classes(fac)
    edges, omega, _ = coc.load_cocycle(path[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)

    six = np.nonzero(edeg == 6)[0]
    eu6 = eu[six]
    E6 = len(six)
    x0 = frac[eu6[:, 0]]
    d = frac[eu6[:, 1]] - frac[eu6[:, 0]]
    d -= np.round(d)                                     # min image
    dc = d @ basis

    kcart = 2 * np.pi * (NV @ np.linalg.inv(basis).T)
    kmag = np.linalg.norm(kcart, axis=1)
    khat = kcart / kmag[:, None]

    # incidence over ALL vertices: +1 at head, -1 at tail (D = B_full phi)
    rows = np.concatenate([eu6[:, 1], eu6[:, 0]])
    cols = np.concatenate([np.arange(E6), np.arange(E6)])
    vals = np.concatenate([np.ones(E6), -np.ones(E6)])
    Bfull = sp.csr_matrix((vals, (rows, cols)), shape=(Vn, E6))
    legal = np.nonzero(imp == 0)[0]                      # constraint rows
    B = Bfull[legal]
    nmono = Vn - len(legal)

    kern = segment_kernel(x0, d)

    acc = np.zeros((3, len(NV)))                         # tot, T, L (projected)
    accN = np.zeros((3, len(NV)))                        # unprojected null
    surv, resid, ident, sumD2 = [], [], [], []
    for _ in range(NSIG):
        sig = rng.choice([-1.0, 1.0], size=E6)
        mu = lsmr(B.T, sig, atol=1e-14, btol=1e-14, maxiter=20000)[0]
        phi = sig - B.T @ mu
        resid.append(np.abs(B @ phi).max())              # legal-vertex div ~ 0
        surv.append(phi @ phi / E6)
        J, tot, T, L = flux_channels(phi, dc, kern, khat)
        acc += np.stack([tot, T, L])
        # telescoping identity: k.J == (1/i) sum_v D_v e^{2pi i n.f_v}
        D = Bfull @ phi
        sumD2.append(D @ D)                              # Poisson level of S_mono
        hot = np.nonzero(np.abs(D) > 1e-9)[0]
        rhs = (np.exp(2j * np.pi * (frac[hot] @ NV.T)).T @ D[hot]) / 1j
        lhs = np.einsum("ma,ma->m", kcart, J)
        ident.append(np.abs(lhs - rhs).max())
        _, totN, TN, LN = flux_channels(sig, dc, kern, khat)
        accN += np.stack([totN, TN, LN])
    tot, T, L = acc / NSIG
    totN, TN, LN = accN / NSIG

    D2 = np.mean(sumD2)                                  # uncorrelated reference
    print(f"\n=== {label}: E6={E6} monopoles(imp>0)={nmono} "
          f"flux-survival={np.mean(surv):.3f} sum(D^2)={D2:.2f} "
          f"max|div_legal|={max(resid):.1e} identity-err={max(ident):.1e} ===")
    print(f"{'shell':>7s} {'T/edge':>9s} {'L/edge':>10s} {'L/T':>9s} "
          f"{'L/T null':>9s} {'S_mono=L*k^2':>13s} {'S_mono/D2':>10s}")
    for lo, hi in SHELLS:
        m = (KM >= lo) & (KM < hi)
        t, l = T[m].mean() / E6, L[m].mean() / E6
        tn, ln = TN[m].mean() / E6, LN[m].mean() / E6
        smono = (L[m] * kmag[m] ** 2).mean()
        print(f"{lo:.0f}-{hi:.0f} {t:9.4f} {l:10.2e} {l/max(t,1e-300):9.2e} "
              f"{ln/max(tn,1e-300):9.3f} {smono:13.4f} "
              f"{smono/max(D2,1e-300):10.3f}")
    print("  (projected flux: L = monopole structure factor / k^2, exactly; "
          "closed web => L = 0 at machine precision.\n   null = unprojected "
          "random signs on the same web. S_mono -> n_mono: Poisson monopole "
          "gas; -> 0 at small k: SL-screening plasma.)")
