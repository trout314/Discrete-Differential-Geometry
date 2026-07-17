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
    sigma CrFe sigma,  30 atoms, CN 13.467, q 5.1090  (Z12_10 Z14_16 Z15_4)
          (dodecagonal-QC approximant; internal parameters from refinements)

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
from scipy.spatial import Delaunay

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from discrete_differential_geometry import Manifold
from fk_skeleton import edges_from_facets, vertex_class_census

SQ3 = np.sqrt(3.0)


def _fcc(base):
    """Apply +(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0) translations."""
    tr = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    return np.vstack([(np.asarray(base) + t) % 1.0 for t in tr])


def _sigma_sites():
    """CrFe sigma phase, P4_2/mnm (136), 30 atoms; refined internal params."""
    def orbit_4f(x):
        return [(x, x, 0), (-x, -x, 0), (.5 + x, .5 - x, .5), (.5 - x, .5 + x, .5)]

    def orbit_8i(x, y):
        return [(x, y, 0), (-x, -y, 0), (y, x, 0), (-y, -x, 0),
                (.5 + x, .5 - y, .5), (.5 - x, .5 + y, .5),
                (.5 + y, .5 - x, .5), (.5 - y, .5 + x, .5)]

    def orbit_8j(x, z):
        return [(x, x, z), (x, x, -z), (-x, -x, z), (-x, -x, -z),
                (.5 + x, .5 - x, .5 + z), (.5 + x, .5 - x, .5 - z),
                (.5 - x, .5 + x, .5 + z), (.5 - x, .5 + x, .5 - z)]

    sites = [(0, 0, 0), (.5, .5, .5)]                     # 2a
    sites += orbit_4f(0.39864)
    sites += orbit_8i(0.46349, 0.13122)
    sites += orbit_8i(0.73933, 0.06609)
    sites += orbit_8j(0.18267, 0.25202)
    return np.array(sites) % 1.0


def _c14_sites(z_mg=0.06286, x_h=0.83048):
    """MgZn2, P6_3/mmc (194), 12 atoms (hexagonal axes)."""
    mg = [(1/3, 2/3, z_mg), (2/3, 1/3, z_mg + .5),
          (2/3, 1/3, -z_mg), (1/3, 2/3, .5 - z_mg)]
    zn2a = [(0, 0, 0), (0, 0, .5)]
    x = x_h
    zn6h = [(x, 2 * x, .25), (-2 * x, -x, .25), (x, -x, .25),
            (-x, -2 * x, .75), (2 * x, x, .75), (-x, x, .75)]
    return np.array(mg + zn2a + zn6h) % 1.0


STRUCTURES = {
    # name: (lattice rows (unit a=1), fractional sites, CN_mean, census {Z: count})
    "a15": (np.eye(3),
            np.array([(0, 0, 0), (.5, .5, .5),
                      (.25, 0, .5), (.75, 0, .5),
                      (.5, .25, 0), (.5, .75, 0),
                      (0, .5, .25), (0, .5, .75)]),
            13.5, {"Z12": 2, "Z14": 6, "Z15": 0, "Z16": 0}),
    "c15": (np.eye(3),
            np.vstack([_fcc([(0, 0, 0), (.25, .25, .25)]),          # 8a (Mg)
                       _fcc([(5/8, 5/8, 5/8), (5/8, 7/8, 7/8),      # 16d (Cu)
                             (7/8, 5/8, 7/8), (7/8, 7/8, 5/8)])]),
            40 / 3, {"Z12": 16, "Z14": 0, "Z15": 0, "Z16": 8}),
    "c14": (np.array([[1, 0, 0], [-.5, SQ3 / 2, 0], [0, 0, np.sqrt(8 / 3)]]),
            _c14_sites(),
            40 / 3, {"Z12": 8, "Z14": 0, "Z15": 0, "Z16": 4}),
    "sigma": (np.diag([1, 1, 0.51759]),
              _sigma_sites(),
              404 / 30, {"Z12": 10, "Z14": 16, "Z15": 4, "Z16": 0}),
}


def build_t3_triangulation(name, m, perturb=1e-6, pad=2):
    """Return (facets ndarray (f3,4), n_vertices) for an m^3 supercell of the
    named TCP structure, as a triangulation of T^3."""
    if m < 2:
        raise SystemExit("need m >= 2 for a simplicial quotient (m >= 3 safer)")
    L, sites, cn, _ = STRUCTURES[name]
    ns = len(sites)
    rng = np.random.default_rng(12345)
    sites = sites + perturb * rng.standard_normal(sites.shape)  # once, pre-tiling

    # Tile cells in [-pad, m+pad)^3; canonical id = site + m^3-cell (mod m).
    rng_cells = np.arange(-pad, m + pad)
    cells = np.array(np.meshgrid(rng_cells, rng_cells, rng_cells,
                                 indexing="ij")).reshape(3, -1).T   # (nc, 3)
    frac = (sites[None, :, :] + cells[:, None, :]).reshape(-1, 3)   # (nc*ns, 3)
    pts = frac @ L
    cell_mod = (cells[:, None, :] % m).repeat(ns, axis=1).reshape(-1, 3)
    site_idx = np.tile(np.arange(ns), len(cells))
    canon = ((cell_mod[:, 0] * m + cell_mod[:, 1]) * m + cell_mod[:, 2]) * ns + site_idx

    tri = Delaunay(pts)
    tets = tri.simplices                                            # (nt, 4)
    # keep one representative per translation class: centroid in [0, m)^3
    cent = pts[tets].mean(axis=1) @ np.linalg.inv(L)                # unit-cell frac
    keep = np.all((cent >= 0) & (cent < m), axis=1)
    fac = canon[tets[keep]]
    fac.sort(axis=1)
    if np.any(fac[:, :-1] == fac[:, 1:]):
        raise SystemExit(f"{name} m={m}: quotient tet with repeated vertex — increase m")
    fac = fac[np.lexsort(fac.T[::-1])]
    uniq = np.unique(fac, axis=0)
    if len(uniq) != len(fac):
        raise SystemExit(f"{name} m={m}: duplicate facets after quotient "
                         f"({len(fac)} vs {len(uniq)} unique) — degeneracy inconsistency")
    f3_expected = round(ns * m ** 3 * (cn / 2 - 1))
    if len(fac) != f3_expected:
        raise SystemExit(f"{name} m={m}: f3 = {len(fac)} != expected {f3_expected} "
                         f"(bad positions or selection)")
    return fac, ns * m ** 3


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
    from seed_utils import get_git_info
    git = get_git_info()
    mfd.save(path, comments=[
        f"topology = T3", f"structure = {name}", f"supercell = {m}",
        f"source = scripts/tcp_reference.py (periodic Delaunay of Wyckoff positions)",
        f"git_commit = {git}",
        f"mean_edge_degree = {edeg.mean():.6f}",
        f"validation = {'PASS' if good else 'PARTIAL (see census)'}",
    ])
    print(f"  {'VALIDATED' if good else 'PARTIAL'} -> {path}")
    return good


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("structure", choices=list(STRUCTURES) + ["all"])
    ap.add_argument("-m", type=int, default=3, help="supercell multiplier (>=2)")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "tcp_reference"))
    args = ap.parse_args()
    sys.path.insert(0, os.path.join(_ROOT, "tools"))
    names = list(STRUCTURES) if args.structure == "all" else [args.structure]
    results = {n: validate_and_save(n, args.m, args.out) for n in names}
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
