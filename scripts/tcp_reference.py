#!/usr/bin/env python3
"""Build T^3 triangulations of real tetrahedrally-close-packed (TCP) crystals.

A Frank-Kasper / TCP phase is BY DEFINITION a structure whose Delaunay
decomposition contains only tetrahedra, with every edge shared by 5 or 6 of
them. This script takes published Wyckoff positions of classical TCP phases,
builds an m x m x m supercell, computes the periodic Delaunay triangulation,
and quotients it to a simplicial triangulation of the 3-torus, saved as .mfd.

Structures (atoms/cell, mean coordination CN, mean edge degree q = 6-12/CN):
    a15   Cr3Si type,  8 atoms, CN 13.500, q 5.1111  (Z12_2 Z14_6)
    c15   MgCu2 Laves, 24 atoms, CN 13.333, q 5.1000  (Z12_16 Z16_8)
    c14   MgZn2 Laves, 12 atoms, CN 13.333, q 5.1000  (Z12_8 Z16_4)
    c36   MgNi2 Laves, 24 atoms, CN 13.333, q 5.1000  (Z12_16 Z16_8)
    sigma CrFe sigma,  30 atoms, CN 13.467, q 5.1090  (Z12_10 Z14_16 Z15_4)
          (dodecagonal-QC approximant; internal parameters from refinements)
    z     Zr4Al3,      7 atoms, CN 13.429, q 5.1064  (Z12_3 Z14_2 Z15_2)
    mu    W6Fe7 type, 39 atoms (hex cell), CN 13.385, q 5.1034
          (Z12_21 Z14_6 Z15_6 Z16_6; refined Co7Mo6 isotype parameters)
    r     Co-Cr-Mo R phase, 159 atoms (hex cell), CN 13.396, q 5.1042
          (Z12_81 Z14_36 Z15_18 Z16_24) — q matches the flat pin 5.1043
    p     Cr-Ni-Mo P phase, 56 atoms, CN 13.429, q 5.1064
          (Z12_24 Z14_20 Z15_8 Z16_4)
    delta MoNi delta, 56 atoms, CN 13.429, q 5.1064  (Z12_24 Z14_20 Z15_8 Z16_4)

Coordinates for c36/z/mu/r/p/delta are full expanded cells from published
refinements (COD/AFLOW; see per-block comments), independently validated by
periodic Voronoi CN census before inclusion. Two COD typos corrected: R-phase
A7 z = 0.3969 (COD 2310299 has 0.3696, transposed) and the P-phase cell
setting (COD 2310313 permutes the cell without permuting coordinates; we use
the original Pbnm cell a=9.070, b=16.983, c=4.752 of Shoemaker 1955/57).

Method notes:
  * Site positions get one deterministic ~1e-6 perturbation (fixed RNG seed)
    BEFORE tiling, so every periodic image is perturbed identically: this
    resolves co-spherical Delaunay degeneracies CONSISTENTLY across images
    (a per-point joggle like Qhull QJ could triangulate two images of the
    same degeneracy differently and break the quotient).
  * The supercell is tiled with a 2-unit-cell margin, Delaunay is computed on
    the extended cloud, and exactly one representative per translation class
    of tets is kept (centroid in the fundamental domain).
  * Validation: expected f-vector, no duplicate facets, Euler characteristic 0,
    orientable, the exact link sum rule (sum(6-deg)=12 at every vertex, i.e.
    a closed 3-manifold), edge degrees all in {5,6}, and the literature
    Z-class census. m >= 2 required (m >= 3 recommended) for a simplicial
    quotient.

Usage:
    python scripts/tcp_reference.py a15 -m 3                  # validate+save
    python scripts/tcp_reference.py c15 -m 6 --out data/tcp_reference/
    python scripts/tcp_reference.py all -m 3                  # all structures
"""
import argparse
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold
from discrete_differential_geometry.fk_skeleton import (
    edges_from_facets, vertex_class_census)
# Construction core promoted to the package (2026-08 cleanup, Phase 1a);
# re-exported here so legacy ``from tcp_reference import ...`` keeps working.
from discrete_differential_geometry.tcp_reference import (  # noqa: F401
    STRUCTURES, build_t3_triangulation, reference_frac_positions)

def validate_and_save(name, m, out_dir):
    fac, f0_expected = build_t3_triangulation(name, m)
    _, _, cn, census_per_cell = STRUCTURES[name]

    eu, edeg, V = edges_from_facets(fac)
    assert V == f0_expected, f"vertex count {V} != {f0_expected}"
    hist = np.bincount(edeg)
    ok_56 = hist.sum() == (hist[5] if len(hist) > 5 else 0) + (hist[6] if len(hist) > 6 else 0)
    fz, n_broken = vertex_class_census(eu, edeg, V)   # asserts link sum rule
    m3 = m ** 3

    mfd = Manifold(3, fac.tolist())
    euler = mfd.euler_characteristic
    orient = mfd.is_orientable

    print(f"[{name} m={m}] f0={V} f1={len(eu)} f3={len(fac)} "
          f"q̄={edeg.mean():.4f} (expect {6 - 12 / cn:.4f})")
    print(f"  edge degrees: " + ", ".join(f"deg{d}:{c}" for d, c in enumerate(hist) if c)
          + ("  [all {5,6} OK]" if ok_56 else "  [NON-{5,6} EDGES — positions off]"))
    print(f"  Z census (per cell, expect {census_per_cell}): "
          + ", ".join(f"{k}:{fz[k] * V / m3:.2f}" for k in ("Z12", "Z14", "Z15", "Z16"))
          + f"  fFK={fz['Z12'] + fz['Z14'] + fz['Z15'] + fz['Z16']:.4f}")
    print(f"  euler={euler} orientable={orient} broken56={n_broken}")
    good = (ok_56 and euler == 0 and orient
            and all(abs(fz[k] * V / m3 - census_per_cell[k]) < 1e-6
                    for k in census_per_cell))
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"T3_{name.upper()}_m{m}_N{len(fac)}.mfd")
    from discrete_differential_geometry.seed_utils import (
        get_git_info, make_leg, history_fields)
    commit, dirty = get_git_info()
    qbar = float(edeg.mean())
    # Root leg: this crystal is a from-scratch lineage ORIGIN (built exactly from
    # Wyckoff positions, not sampled). Stamping it means downstream melt/dope/
    # quench runs append to it via read_history() and can never lose provenance.
    # See the history section in ddg.seed_utils (is_root_leg_from).
    root_leg = make_leg(
        "build",
        {"struct": name, "m": m, "cn": cn, "qbar": round(qbar, 6),
         "source": "wyckoff-delaunay"},
        sweeps=0, from_=f"crystal:{name}@wyckoff",
        commit=commit, dirty=dirty, tried=0, accepted=0)
    mfd.save(path, comments=[
        f"topology = T3", f"structure = {name}", f"supercell = {m}",
        f"source = scripts/tcp_reference.py (periodic Delaunay of Wyckoff positions)",
        f"git_commit = {commit}", f"git_dirty = {'true' if dirty else 'false'}",
        f"mean_edge_degree = {qbar:.6f}",
        f"validation = {'PASS' if good else 'PARTIAL (see census)'}",
    ] + history_fields([root_leg]))
    print(f"  {'VALIDATED' if good else 'PARTIAL'} -> {path}")
    return good


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("structure", choices=list(STRUCTURES) + ["all"])
    ap.add_argument("-m", type=int, default=3, help="supercell multiplier (>=2)")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "tcp_reference"))
    args = ap.parse_args()
    names = list(STRUCTURES) if args.structure == "all" else [args.structure]
    results = {n: validate_and_save(n, args.m, args.out) for n in names}
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
