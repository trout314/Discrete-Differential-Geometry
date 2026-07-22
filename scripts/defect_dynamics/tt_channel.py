#!/usr/bin/env python3
"""Radiative-channel structure factor: does the deficit quadrupole carry a free
transverse-traceless (graviton-like) channel while the scalar sector is gapped?

Deficit quadrupole per vertex (traceless symmetric rank-2):
    Q_ab(v) = 1/2 sum_{e ni v} delta_e (nhat_e^a nhat_e^b - 1/3 delta_ab),
delta_e = 2pi - theta deg(e), nhat_e = min-image edge direction (harmonic coords).

Q_ab(k) = sum_v Q_ab(v) e^{ik.x_v}; decompose vs khat into the SO(2) helicity
channels of a symmetric traceless 3x3 (5 dof = 2 TT + 2 vector + 1 longitudinal):
    P_ab = delta_ab - khat_a khat_b
    TT   : (P Q P)_ab - 1/2 P_ab (P:Q)          [helicity +-2, radiative]
    L    : 3/2 |khat.Q.khat|^2                    [helicity 0, scalar-like]
    V    : total - TT - L                         [helicity +-1]
Ratio S_c(k) = P_c^obs / P_c^shuffle (tensors permuted over fixed positions):
<1 suppressed/rigid, ~1 free/Poisson. Scalar q_R channel computed alongside for
reference. Shuffle split TT:V:L -> 2:2:1 is the built-in sanity check.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, "python"); sys.path.insert(0, "scripts")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS, edges_and_degrees

THETA = float(np.arccos(1.0 / 3.0))
SP = sys.argv[1]
NMAX, KCUT = 8, 5.0
NSHUF = 6


def deficit_quadrupole(fac, frac, basis):
    eu, ecnt, deg, V = edges_and_degrees(fac)
    delta = 2 * np.pi - THETA * ecnt
    df = frac[eu[:, 1]] - frac[eu[:, 0]]
    df -= np.round(df)                          # min image
    d = df @ basis
    nhat = d / np.linalg.norm(d, axis=1, keepdims=True)
    t = nhat[:, :, None] * nhat[:, None, :] - np.eye(3) / 3
    w = 0.5 * delta
    Q = np.zeros((V, 3, 3))
    np.add.at(Q, eu[:, 0], w[:, None, None] * t)
    np.add.at(Q, eu[:, 1], w[:, None, None] * t)
    Q -= Q.mean(0, keepdims=True)
    qsc = FIELDS["curvature_charge"](fac)
    return Q, qsc - qsc.mean()


def kgrid(basis):
    r = np.arange(-NMAX, NMAX + 1)
    nvec = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    nvec = nvec[np.any(nvec != 0, 1)]
    keep = ((nvec[:, 0] > 0) | ((nvec[:, 0] == 0) & (nvec[:, 1] > 0))
            | ((nvec[:, 0] == 0) & (nvec[:, 1] == 0) & (nvec[:, 2] > 0)))
    nvec = nvec[keep]
    Binv = np.linalg.inv(basis)
    kphys = nvec @ Binv.T
    kmag = np.linalg.norm(kphys, axis=1) / np.abs(Binv).max()
    m = kmag <= KCUT
    kp = kphys[m]
    khat = kp / np.linalg.norm(kp, axis=1, keepdims=True)
    return nvec[m], kmag[m], khat


def powers(Q, qsc, frac, nvec, khat):
    phase = np.exp(2j * np.pi * (frac @ nvec.T))       # (V, M)
    Qk = np.einsum("vab,vm->mab", Q, phase)            # (M,3,3) complex
    P = np.eye(3)[None] - khat[:, :, None] * khat[:, None, :]
    PQP = np.einsum("mac,mcd,mdb->mab", P, Qk, P)
    PtrQ = np.einsum("mcd,mcd->m", P, Qk)
    LQ = PQP - 0.5 * P * PtrQ[:, None, None]
    tot = np.einsum("mab,mab->m", np.conj(Qk), Qk).real
    TT = np.einsum("mab,mab->m", np.conj(Qk), LQ).real
    QL = np.einsum("ma,mab,mb->m", khat, Qk, khat)
    L = 1.5 * np.abs(QL) ** 2
    V = tot - TT - L
    sc = np.abs(qsc @ phase) ** 2
    return dict(TT=TT, V=V, L=L, scal=sc, tot=tot)


def null(Q, qsc, frac, nvec, khat, rng):
    acc = None
    for _ in range(NSHUF):
        p = rng.permutation(len(Q))
        pw = powers(Q[p], qsc[p], frac, nvec, khat)
        acc = pw if acc is None else {k: acc[k] + pw[k] for k in acc}
    return {k: v / NSHUF for k, v in acc.items()}


PATTERN = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].isdigit() \
    else "run5h/*_snap*.mfd"
snaps = sorted(glob.glob(f"{SP}/{PATTERN}"))
if len(sys.argv) > 2 and sys.argv[2].isdigit():
    snaps = snaps[:int(sys.argv[2])]
if len(sys.argv) > 3 and sys.argv[3].isdigit():
    snaps = snaps[:int(sys.argv[3])]
rng = np.random.default_rng(0)

CH = ["scal", "TT", "V", "L"]
edg = np.geomspace(1.0, KCUT, 6)
binn = {c: [[] for _ in range(len(edg) - 1)] for c in CH}
split = {c: [] for c in ("TT", "V", "L")}     # shuffle dof split sanity
for snap in snaps:
    fac = np.asarray(Manifold.load(snap, 3).facets())
    edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    Q, qsc = deficit_quadrupole(fac, frac, basis)
    nvec, kmag, khat = kgrid(basis)
    obs = powers(Q, qsc, frac, nvec, khat)
    nul = null(Q, qsc, frac, nvec, khat, rng)
    for c in ("TT", "V", "L"):
        split[c].append(nul[c].sum() / nul["tot"].sum())
    for bi, (lo, hi) in enumerate(zip(edg[:-1], edg[1:])):
        sel = (kmag >= lo) & (kmag < hi)
        if sel.sum() >= 3:
            for c in CH:
                binn[c][bi].append((obs[c][sel] / nul[c][sel]).mean())
    print(f"  {os.path.basename(snap)[:-4]}: {len(kmag)} sub-Bragg modes", flush=True)

print(f"\nshuffle-null dof split TT:V:L = "
      f"{np.mean(split['TT']):.2f} : {np.mean(split['V']):.2f} : "
      f"{np.mean(split['L']):.2f}   (expect 0.40 : 0.40 : 0.20)")
kc = np.sqrt(edg[:-1] * edg[1:])
print(f"\n{'k':>5s} " + " ".join(f"{c:>8s}" for c in CH))
for bi in range(len(kc)):
    if all(binn[c][bi] for c in CH):
        print(f"{kc[bi]:5.2f} " + " ".join(f"{np.mean(binn[c][bi]):8.4f}" for c in CH))

# headline: lowest-k shell
print("\nlowest-k shell (k~1-1.4) ratio to shuffle null:")
for c in CH:
    v = binn[c][0]
    if v:
        print(f"  {c:5s}: {np.mean(v):.4f} +/- {np.std(v)/np.sqrt(len(v)):.4f}"
              f"   ({'suppressed' if np.mean(v) < 0.5 else 'FREE/Poisson' if np.mean(v) > 0.7 else 'partial'})")

# figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.4, 5.4))
col = {"scal": "C0", "TT": "C3", "V": "C2", "L": "C1"}
lab = {"scal": "scalar q_R", "TT": "TT (radiative)", "V": "vector", "L": "longitudinal"}
for c in CH:
    y = [np.mean(binn[c][bi]) if binn[c][bi] else np.nan for bi in range(len(kc))]
    e = [np.std(binn[c][bi]) / np.sqrt(len(binn[c][bi])) if binn[c][bi] else 0
         for bi in range(len(kc))]
    ax.errorbar(kc, y, yerr=e, fmt="o-", color=col[c], ms=5, capsize=2, label=lab[c])
ax.axhline(1.0, color="k", lw=0.8, ls="--", label="Poisson (free)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$|k|$ (sub-Bragg, box-mode units)")
ax.set_ylabel(r"channel structure factor $S_c(k)$ / shuffle null")
ax.set_title("Deficit-quadrupole channels: is the TT (radiative) sector free\n"
             "while the scalar sector is gapped?  (32 run5h snapshots)")
ax.legend(fontsize=8.5); ax.grid(alpha=0.25, which="both")
fig.tight_layout()
out = f"{SP}/tt_channel.png"
fig.savefig(out, dpi=130)
print(f"\nSaved to: {out}")
