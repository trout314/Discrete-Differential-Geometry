#!/usr/bin/env python3
"""Curvature structure factor S(k) on T^3 via cocycle harmonic coordinates.

The flagship constraint test: in a state satisfying the (discrete) Hamiltonian
constraint at long wavelength, curvature-charge fluctuations are surface-law,
so the charge-weighted structure factor is suppressed as k -> 0
(hyperuniformity). This measures it with REAL positions, not graph proxies:

  1. load snapshot .mfd + its .cocycle.npz (canonical labels),
  2. integrate the integer cocycle over a spanning tree (lift), read the
     winding lattice diag(M1, M2, M3) off the fundamental cycles,
  3. harmonic gauge (periodic Tutte embedding): X_v = lift_v - phi_v gives
     every vertex a position on the metric torus,
  4. charge field q_v (default: n6, the six-edge count at each vertex = the
     local disclination-line density; note the naive vertex deficit
     sum_{e at v} (6 - deg e) is IDENTICALLY 12 on every legal FK vertex --
     the link sum rule, i.e. the vertex-level Gauss constraint is exact --
     so it only sees impurities),
  5. S(k) = |sum_v dq_v exp(i k.X_v)|^2 / sum dq^2 at torus-commensurate
     k = 2 pi (n1/M1, n2/M2, n3/M3), against the EXACT permutation null
     E[S]/S2 = 1 - (|F(k)|^2 - N) / (N(N-1)),  F(k) = sum_v exp(i k.X_v),
     (charges shuffled over fixed positions -- flat marginal-preserving null).

Ratio < 1 at small |n| = hyperuniform. Pools replicas by --group regex.

Usage:
    python scripts/sk_torus.py data/replicas/flatpin_mu3_r*_final.mfd \
        --out out/sk_flatpin.csv --plot out/sk_flatpin.png
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from fk_skeleton import edges_from_facets


def charge_field(facets, n_verts, kind):
    eu, edeg, _ = edges_from_facets(facets)
    q = np.zeros(n_verts)
    if kind == "deficit":
        w = 6.0 - edeg
    elif kind == "n6":
        w = (edeg >= 6).astype(float)
    else:
        raise ValueError(kind)
    np.add.at(q, eu[:, 0], w)
    np.add.at(q, eu[:, 1], w)
    return q


def sk_one(mfd_path, coc_path, field, nmax):
    v = ddg.Manifold.load(mfd_path, 3)
    facets = np.asarray(v.facets())
    edges, omega, _ = coc.load_cocycle(coc_path)
    # legacy npz may carry raw sampler labels (holes); the .mfd writer
    # renumbers to rank order, which is exactly canonicalize (idempotent)
    edges, omega = coc.canonicalize_labels(edges, omega)
    n_verts = int(facets.max()) + 1
    if edges.max() + 1 != n_verts:
        raise RuntimeError(f"vertex count mismatch {mfd_path}: "
                           f"mfd {n_verts} vs cocycle {edges.max() + 1}")

    pos, cyc, _ = coc.tree_positions(edges, omega, n_verts)
    basis = coc.lattice_basis(cyc).astype(float)   # rows = lattice vectors
    if basis.shape != (3, 3) or abs(np.linalg.det(basis)) < 1:
        raise RuntimeError(f"winding lattice degenerate for {mfd_path}: "
                           f"{basis}")
    # General (possibly oblique/hexagonal, e.g. R phase) torus: fractional
    # coordinates s = X @ B^{-1} in [0,1)^3; a commensurate mode has integer
    # reciprocal index n with phase exp(2 pi i s.n), physical |k| ~ |B^{-T} n|.
    Binv = np.linalg.inv(basis)
    _, phi, _ = coc.harmonic_gauge(edges, omega, n_verts)
    s = ((pos - phi) @ Binv) % 1.0

    q = charge_field(facets, n_verts, field)
    dq = q - q.mean()
    S2 = float(dq @ dq)
    N = n_verts

    rng = np.arange(-nmax, nmax + 1)
    nvec = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"),
                    axis=-1).reshape(-1, 3)
    nvec = nvec[np.any(nvec != 0, axis=1)]
    # half-space only (S(-k) = S(k) for a real field)
    keep = ((nvec[:, 0] > 0)
            | ((nvec[:, 0] == 0) & (nvec[:, 1] > 0))
            | ((nvec[:, 0] == 0) & (nvec[:, 1] == 0) & (nvec[:, 2] > 0)))
    nvec = nvec[keep]

    s_obs = np.empty(len(nvec))
    s_null = np.empty(len(nvec))
    for i0 in range(0, len(nvec), 256):
        blk = nvec[i0:i0 + 256]
        phase = np.exp(2j * np.pi * s @ blk.T)         # (N, blk)
        A = dq @ phase
        F = phase.sum(axis=0)
        s_obs[i0:i0 + len(blk)] = np.abs(A) ** 2 / S2
        s_null[i0:i0 + len(blk)] = 1.0 - (np.abs(F) ** 2 - N) / (N * (N - 1))
    kphys = nvec @ Binv.T                               # reciprocal vectors
    kmag = np.linalg.norm(kphys, axis=1) / np.abs(Binv).max()
    return kmag, s_obs, s_null, np.abs(np.diag(basis)), N


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+",
                    help=".mfd files; each needs a sibling .cocycle.npz")
    ap.add_argument("--field", choices=["n6", "deficit"], default="n6")
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--group", default=r"_r\d+.*$")
    ap.add_argument("--out", default=None)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    pooled = defaultdict(list)
    for path in args.snapshots:
        cpath = os.path.splitext(path)[0] + ".cocycle.npz"
        if not os.path.exists(cpath):
            print(f"  SKIP {os.path.basename(path)}: no cocycle")
            continue
        label = re.sub(args.group, "", os.path.splitext(os.path.basename(path))[0])
        kmag, s_obs, s_null, M, N = sk_one(path, cpath, args.field, args.nmax)
        pooled[label].append((kmag, s_obs, s_null))
        low = kmag <= 2.0 + 1e-9
        print(f"  {os.path.basename(path)}: N={N} M={M.astype(int).tolist()} "
              f"low-k ratio (|n|<=2): "
              f"{np.mean(s_obs[low] / s_null[low]):.4f}")

    rows = []
    for label, runs in sorted(pooled.items()):
        kmag = runs[0][0]
        s_obs = np.mean([r[1] for r in runs], axis=0)
        s_null = np.mean([r[2] for r in runs], axis=0)
        ratio = s_obs / s_null
        edges = np.arange(0.5, args.nmax + 0.5)
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (kmag > lo) & (kmag <= hi)
            if m.any():
                rows.append({"label": label, "k_shell": (lo + hi) / 2,
                             "n_modes": int(m.sum()),
                             "s_obs": float(s_obs[m].mean()),
                             "s_null": float(s_null[m].mean()),
                             "ratio": float(ratio[m].mean()),
                             "ratio_sem": float(ratio[m].std(ddof=1)
                                                / np.sqrt(m.sum()))
                             if m.sum() > 1 else 0.0})
        low = kmag <= 2.0 + 1e-9
        print(f"{label} ({len(runs)} snapshots): low-k (|n|<=2) ratio = "
              f"{ratio[low].mean():.4f} ± "
              f"{ratio[low].std(ddof=1) / np.sqrt(low.sum()):.4f}")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4.2))
        for label in sorted(pooled):
            sel = [r for r in rows if r["label"] == label]
            k = [r["k_shell"] for r in sel]
            ax.errorbar(k, [r["ratio"] for r in sel],
                        yerr=[r["ratio_sem"] for r in sel],
                        marker="o", ms=4, lw=1, capsize=2, label=label)
        ax.axhline(1.0, color="k", lw=0.8, ls="--", label="shuffled null")
        ax.set_xlabel("|k| (units of 2π/L)")
        ax.set_ylabel(f"S(k) / S_null  ({args.field})")
        ax.set_yscale("log")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs(os.path.dirname(args.plot) or ".", exist_ok=True)
        fig.savefig(args.plot, dpi=120)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
