"""Regression + sanity tests for the hyperuniformity refactor: the shared field
library (vertex_fields), the coordinate builder (cocycle.torus_positions), and
the two estimators (structure_factor real-k, graph_hyperuniformity graph proxy).

Built in-memory from Wyckoff positions (tcp_reference), so no data files needed.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry import vertex_fields as vf
from discrete_differential_geometry import graph_hyperuniformity as gh
from discrete_differential_geometry.structure_factor import structure_factor
import tcp_reference as tr
from cocycle_check import reference_frac_positions

THETA = float(np.arccos(1.0 / 3.0))


@pytest.fixture(scope="module")
def crystal():
    """Perfect R crystal (m=2) + its cocycle -> (facets, frac coords, lattice)."""
    fac = np.asarray(tr.build_t3_triangulation("r", 2)[0])
    eu, ecnt, deg, V = vf.edges_and_degrees(fac)
    omega = coc.build_from_positions(eu, reference_frac_positions("r", 2), 2)
    frac, basis = coc.torus_positions(fac, eu, omega)
    return fac, frac, basis


# ---- field library -------------------------------------------------------

def test_field_shapes_and_labeling(crystal):
    fac, _, _ = crystal
    V = len(np.unique(fac))
    for name, fn in vf.FIELDS.items():
        q = fn(fac)
        assert q.shape == (V,), name
        assert np.isfinite(q).all(), name


def test_naive_deficit_is_identically_12(crystal):
    """(6 - deg) deficit is a topological constant 12 (link sum rule) -- the
    reason it is NOT a FIELDS entry (a constant field has no structure factor)."""
    fac, _, _ = crystal
    eu, ecnt, deg, V = vf.edges_and_degrees(fac)
    naive = np.zeros(V)
    np.add.at(naive, eu[:, 0], 6.0 - ecnt)
    np.add.at(naive, eu[:, 1], 6.0 - ecnt)
    assert np.allclose(naive, 12.0)


def test_curvature_charge_is_half_regge_deficit(crystal):
    fac, _, _ = crystal
    eu, ecnt, deg, V = vf.edges_and_degrees(fac)
    dsum = np.zeros(V)
    delta = 2 * np.pi - THETA * ecnt
    np.add.at(dsum, eu[:, 0], delta)
    np.add.at(dsum, eu[:, 1], delta)
    assert np.allclose(vf.curvature_charge(fac), dsum / 2)


def test_curvature_density_is_charge_over_dual_volume(crystal):
    fac, _, _ = crystal
    eu, ecnt, deg, V = vf.edges_and_degrees(fac)
    assert np.allclose(vf.curvature_density(fac),
                       vf.curvature_charge(fac) / (deg / 4.0))


def test_defect_indicator_zero_on_perfect_crystal(crystal):
    fac, _, _ = crystal
    assert vf.defect_indicator(fac).sum() == 0        # no illegal edges


# ---- coordinates ---------------------------------------------------------

def test_torus_positions_shapes(crystal):
    fac, frac, basis = crystal
    V = int(fac.max()) + 1
    assert frac.shape == (V, 3)
    assert basis.shape == (3, 3)
    assert (frac >= 0).all() and (frac < 1).all()
    assert abs(np.linalg.det(basis)) >= 1


def _minima_ratios(B, R=4):
    """Successive-minima ratios of the lattice spanned by the rows of B -- a
    lattice invariant, unlike the shape of B itself (lattice_basis returns a
    Euclidean-reduced, not Minkowski-reduced, basis)."""
    import itertools
    combos = np.array([c for c in itertools.product(range(-R, R + 1), repeat=3)
                       if any(c)])
    V = combos @ B
    chosen = []
    for i in np.argsort(np.linalg.norm(V, axis=1)):
        cand = chosen + [V[i]]
        if np.linalg.matrix_rank(np.array(cand), tol=1e-9) == len(cand):
            chosen.append(V[i])
            if len(chosen) == 3:
                break
    m = np.linalg.norm(np.array(chosen), axis=1)
    return m / m[0]


def test_whitening_is_a_metric_fix_not_a_coordinate_change():
    """Whitening acts on the value space, so fractional coordinates are
    untouched and the fundamental-domain volume is preserved (unimodular);
    only the cell SHAPE carried by ``basis`` changes."""
    fac = np.asarray(tr.build_t3_triangulation("r", 2)[0])
    eu, _, _, _ = vf.edges_and_degrees(fac)
    omega = coc.build_from_positions(eu, reference_frac_positions("r", 2), 2)
    f_raw, b_raw = coc.torus_positions(fac, eu, omega, whiten=False)
    f_w, b_w = coc.torus_positions(fac, eu, omega, whiten=True)
    assert np.array_equal(f_raw, f_w)
    assert np.linalg.det(b_w) == pytest.approx(np.linalg.det(b_raw), rel=1e-9)


def test_min_image_matches_naive_rule_only_on_an_orthogonal_cell():
    """On a diagonal basis the per-axis rule is right; on the Euclidean-reduced
    basis torus_positions actually returns for R m4 it is not, in two separate
    ways -- the primary image is not always nearest, and diag(basis) is not the
    cell. This is the bug that hit carrier_gr / pass2_structure."""
    rng = np.random.default_rng(0)
    dfrac = rng.uniform(-1, 1, (2000, 3))

    B = np.diag([3.0, 5.0, 7.0])
    _, d = coc.min_image(dfrac, B)
    naive = np.linalg.norm((dfrac - np.round(dfrac)) @ B, axis=1)
    assert np.allclose(d, naive)

    B = np.array([[4.0, 4.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]])
    _, d = coc.min_image(dfrac, B)
    # exact against a much wider image search
    import itertools
    offs = np.array(list(itertools.product(range(-3, 4), repeat=3)), float)
    brute = np.min(np.linalg.norm(
        ((dfrac - np.round(dfrac))[:, None, :] + offs) @ B, axis=-1), axis=1)
    assert np.allclose(d, brute)
    # and both naive rules are materially wrong here
    primary = np.linalg.norm((dfrac - np.round(dfrac)) @ B, axis=1)
    assert (primary > d + 1e-9).mean() > 0.2
    diag_only = np.linalg.norm((dfrac - np.round(dfrac)) * np.diag(B), axis=1)
    assert (np.abs(diag_only / d - 1) > 0.1).mean() > 0.5


@pytest.mark.parametrize("name,m", [("c15", 3), ("r", 2), ("c14", 3)])
def test_whitening_recovers_the_true_cell_shape(name, m):
    """The winding lattice in cochain units is cubic for EVERY host (the
    coordinates are fractional), so the unwhitened basis reports a cubic cell
    even for hexagonal R / C14. Whitening by G^{-1/2} recovers the real
    aspect ratio to a few percent without being told the lattice matrix."""
    fac = np.asarray(tr.build_t3_triangulation(name, m)[0])
    eu, _, _, _ = vf.edges_and_degrees(fac)
    omega = coc.build_from_positions(eu, reference_frac_positions(name, m), m)
    _, b_raw = coc.torus_positions(fac, eu, omega, whiten=False)
    _, b_w = coc.torus_positions(fac, eu, omega, whiten=True)
    ideal = _minima_ratios(np.asarray(tr.STRUCTURES[name][0], float) * m)
    assert _minima_ratios(b_raw) == pytest.approx([1.0, 1.0, 1.0], abs=1e-6)
    assert _minima_ratios(b_w) == pytest.approx(ideal, rel=0.03)


# ---- real-k structure factor: crystal is hyperuniform --------------------

@pytest.mark.parametrize("field", ["n6", "curvature_charge"])
def test_structure_factor_crystal_hyperuniform(crystal, field):
    fac, frac, basis = crystal
    kmag, s_obs, s_null = structure_factor(frac, basis, vf.FIELDS[field](fac), 4)
    low = kmag <= 2.0 + 1e-9
    ratio = s_obs[low] / s_null[low]
    # a crystal has S(k) = 0 between Bragg peaks -> low-k ratio ~ 0
    assert np.median(ratio) < 1e-6
    assert ratio.mean() < 0.2


# ---- graph estimator: crystal is hyperuniform ----------------------------

def test_graph_lowpass_crystal_hyperuniform(crystal):
    fac, _, _ = crystal
    eu, ecnt, deg, V = vf.edges_and_degrees(fac)
    q = vf.curvature_charge(fac)
    rng = np.random.default_rng(0)
    r = gh.lowpass_ratio(eu, V, q - q.mean(), rng=rng)
    assert r < 0.05                                   # crystal: sub-BZ power ~ 0
