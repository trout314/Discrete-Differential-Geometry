"""Ground-truth tests for the TCP reference-crystal builder (ddg.tcp_reference).

Every structure in the library is built at m=2 and checked against exact
combinatorics: expected f-vector, pure-{5,6} edge degrees, Euler
characteristic 0, the literature Z-class census -- and, strictly, that every
Z16 site is the T_d Friauf polyhedron (a Z16_D2 in a reference crystal would
mean the builder or the detector is wrong).
"""
import numpy as np
import pytest

from discrete_differential_geometry import Manifold
from discrete_differential_geometry import tcp_reference as tr
from discrete_differential_geometry.fk_skeleton import (
    edges_from_facets, vertex_class_census)

M = 2  # smallest valid supercell; keeps the whole library in the fast suite


@pytest.fixture(scope="module", params=list(tr.STRUCTURES))
def crystal(request):
    name = request.param
    fac, n_verts = tr.build_t3_triangulation(name, M)
    return name, np.asarray(fac), n_verts


def test_f_vector(crystal):
    name, fac, n_verts = crystal
    _, sites, cn, _ = tr.STRUCTURES[name]
    ns = len(sites)
    assert n_verts == ns * M ** 3
    assert len(fac) == round(ns * M ** 3 * (cn / 2 - 1))
    eu, edeg, V = edges_from_facets(fac)
    assert V == n_verts
    # closed 3-manifold with only tets: f2 = 2 f3, so Euler = f0 - f1 + f3
    assert V - len(eu) + len(fac) == 0


def test_euler_characteristic(crystal):
    _, fac, _ = crystal
    assert Manifold(3, fac.tolist()).euler_characteristic == 0


def test_edge_degrees_and_census(crystal):
    name, fac, n_verts = crystal
    _, _, cn, census_per_cell = tr.STRUCTURES[name]
    eu, edeg, V = edges_from_facets(fac)
    # TCP by definition: every edge degree in {5, 6}
    assert set(np.unique(edeg).tolist()) <= {5, 6}
    assert abs(float(edeg.mean()) - (6 - 12 / cn)) < 1e-9
    fz, n_broken = vertex_class_census(eu, edeg, V, facets=fac)
    assert n_broken == 0
    m3 = M ** 3
    for k, expect in census_per_cell.items():
        assert abs(fz[k] * V / m3 - expect) < 1e-9, (name, k)
    # STRICT: no D2 isomer in a reference crystal; FK_strict is exactly 1
    assert fz["Z16_D2"] == 0.0
    assert abs(fz["FK_strict"] - 1.0) < 1e-12


def test_reference_frac_positions(crystal):
    name, fac, n_verts = crystal
    rp = tr.reference_frac_positions(name, M)
    assert rp.shape == (n_verts, 3)
    assert np.all((rp >= 0) & (rp < M))
    # id scheme v = cell*ns + site: the first ns vertices are cell (0,0,0),
    # i.e. the (perturbed) sites themselves
    _, sites, _, _ = tr.STRUCTURES[name]
    ns = len(sites)
    d = (rp[:ns] - sites) % M
    d = np.minimum(d, M - d)
    assert np.abs(d).max() < 1e-4  # only the 1e-6 deterministic perturbation


def test_m1_rejected():
    with pytest.raises(SystemExit):
        tr.build_t3_triangulation("a15", 1)
