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
