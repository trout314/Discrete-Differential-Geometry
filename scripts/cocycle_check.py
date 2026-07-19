#!/usr/bin/env python3
"""Validation ladder for T^3 integer-cocycle tracking (Tier 3 instrumentation).

Builds an A15 reference crystal with known coordinates, initializes the three
winding cocycles from them, then climbs:

  0. negative control — a tampered cocycle must be REJECTED at enable
  1. enable on the pristine crystal (D-side closedness + edge-cover audit)
  2. zero-churn round trip — harmonic gauge must recover the crystal's own
     coordinates (RMS << edge length)
  3. drift audit after melting-strength churn (closedness is exact or dead)
  4. class conservation — the fundamental-cycle winding lattice must still be
     exactly (scale * m) Z^3 after churn
  5. harmonic gauge after churn: Gram (cell-shape) matrix + position RMS +
     six-edge nematic order, as a smoke test of the consumers

Usage:
    python scripts/cocycle_check.py [--structure a15] [--m 3] [--sweeps 10]
"""
import argparse
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from tcp_reference import STRUCTURES, build_t3_triangulation

SCALE = 10**6


def reference_frac_positions(name, m, perturb=1e-6):
    """Fractional torus coordinates (period m) per canonical vertex id,
    reproducing build_t3_triangulation's site perturbation and id scheme."""
    L, sites, cn, _ = STRUCTURES[name]
    ns = len(sites)
    rng = np.random.default_rng(12345)
    sites = sites + perturb * rng.standard_normal(sites.shape)
    n = ns * m**3
    v = np.arange(n)
    s = v % ns
    c = v // ns
    cz = c % m
    cy = (c // m) % m
    cx = c // (m * m)
    return (sites[s] + np.stack([cx, cy, cz], axis=1)) % m


def wrapped_rms(a, b, M):
    """RMS deviation between torus position sets after removing the global
    translation (wrap-aware)."""
    d = a - b
    d = d - d[0]
    d -= M * np.round(d / M)
    d -= d.mean(axis=0)
    d -= M * np.round(d / M)
    return float(np.sqrt((d**2).sum(axis=1).mean()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--structure", default="a15")
    ap.add_argument("--m", type=int, default=3)
    ap.add_argument("--sweeps", type=int, default=10)
    ap.add_argument("--lam", type=float, default=0.55)
    args = ap.parse_args()

    fac, n_verts = build_t3_triangulation(args.structure, args.m)
    frac = reference_frac_positions(args.structure, args.m)
    M = SCALE * args.m
    mfd = ddg.Manifold(3, fac.tolist())
    qbar = mfd.mean_degree(1)
    print(f"[{args.structure} m={args.m}] f0={n_verts} f3={len(fac)} "
          f"qbar={qbar:.4f} M={M}")

    params = ddg.SamplerParams(
        num_facets_target=mfd.num_facets, num_facets_coef=0.1,
        hinge_degree_target=qbar, num_hinges_coef=2.0,
        hinge_degree_target_coef=args.lam * qbar / 6.0)
    s = ddg.ManifoldSampler(mfd, params)
    v = s.manifold

    edges = np.asarray(v.simplices(1))
    omega = coc.build_from_positions(edges, frac, args.m, scale=SCALE)
    edge_len = float(np.linalg.norm(omega, axis=1).mean())

    # --- rung 0: tampered cocycle must be rejected ---
    bad = omega.copy()
    bad[0, 0] += 1
    try:
        s.enable_cocycle(edges, bad)
        print("RUNG 0 FAIL: tampered cocycle accepted")
        sys.exit(1)
    except RuntimeError as e:
        print(f"rung 0 PASS: tampered cocycle rejected ({e})")

    # --- rung 1: pristine cocycle accepted ---
    s.enable_cocycle(edges, omega)
    print("rung 1 PASS: crystal cocycle validated at enable "
          f"(mean |omega| = {edge_len:.0f} = {edge_len / SCALE:.3f} cells)")

    # --- rung 2: zero-churn harmonic round trip ---
    _, phi, gram = coc.harmonic_gauge(edges, omega, n_verts)
    pos, cyc, _ = coc.tree_positions(edges, omega, n_verts)
    xh = pos - phi
    rms0 = wrapped_rms(xh, SCALE * frac, M)
    print(f"rung 2 {'PASS' if rms0 < 0.3 * edge_len else 'FAIL'}: "
          f"harmonic round trip RMS = {rms0:.0f} = {rms0 / edge_len:.3f} edges")

    # --- rung 3: churn + drift audit ---
    s.run(sweeps=args.sweeps)
    st = s.get_stats()
    acc = st.total_accepted
    s.check_cocycle()
    print(f"rung 3 PASS: closedness exact after {acc} accepted moves "
          f"({args.sweeps} sweeps at lam={args.lam})")

    # --- rung 4: winding lattice still (SCALE*m) Z^3 ---
    e2, w2 = s.read_cocycle()
    n2 = int(e2.max()) + 1
    pos2, cyc2, _ = coc.tree_positions(e2, w2, n2)
    basis = coc.lattice_basis(cyc2)
    det = 0
    if len(basis) == 3:
        det = abs(int(basis[0, 0]) * int(basis[1, 1]) * int(basis[2, 2]))
    ok4 = det == M**3
    print(f"rung 4 {'PASS' if ok4 else 'FAIL'}: winding lattice det = "
          f"{det} vs M^3 = {M**3} (basis diag = "
          f"{[int(basis[i, i]) for i in range(len(basis))] if len(basis) == 3 else basis})")

    # --- rung 5: harmonic gauge + consumers on the churned state ---
    oh2, phi2, gram2 = coc.harmonic_gauge(e2, w2, n2)
    xh2 = coc.tree_positions(e2, w2, n2)[0] - phi2
    surv = np.isin(np.arange(n2), np.unique(e2))   # labels still alive
    both = surv[:n_verts]
    rms1 = wrapped_rms(xh2[:n_verts][both], SCALE * frac[both], M)
    g = gram2 / np.trace(gram2) * 3
    ev = np.linalg.eigvalsh(g)
    # six-edge nematic order via edge degrees
    degs = np.array([v.degree(e2[i].tolist()) for i in range(len(e2))])
    sixdirs = oh2[degs >= 6]
    _, _, S_nem = coc.nematic_q(sixdirs)
    print(f"rung 5: churned-state position RMS vs crystal = "
          f"{rms1 / edge_len:.2f} edges (memory of the lattice); "
          f"Gram anisotropy eigs = [{ev[0]:.3f}, {ev[1]:.3f}, {ev[2]:.3f}] "
          f"(1 = cubic); nematic S(six-edges) = {S_nem:.3f} over {len(sixdirs)} edges")

    print("\nALL RUNGS PASS" if ok4 and rms0 < 0.3 * edge_len else "\nFAILURES — see above")


if __name__ == "__main__":
    main()
