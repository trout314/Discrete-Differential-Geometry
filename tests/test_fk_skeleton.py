"""Tests for the FK census primitives (ddg.fk_skeleton).

The heavyweight check is fast-vs-exact cross-validation: on a churned (defect-
carrying) state, the vectorized per-vertex classification of
``vertex_class_census`` + ``link_classes.link_class_census`` must agree with
re-deriving every single vertex link from scratch and classifying it exactly.
"""
import numpy as np
import pytest

import discrete_differential_geometry as ddg
from discrete_differential_geometry import tcp_reference as tr
from discrete_differential_geometry import link_classes as lc
from discrete_differential_geometry.fk_skeleton import (
    edges_from_facets, vertex_classes, vertex_class_census, skeleton_stats)


@pytest.fixture(scope="module")
def a15():
    fac, n = tr.build_t3_triangulation("a15", 2)
    return np.asarray(fac), n


def test_edges_from_facets_counts(a15):
    fac, n = a15
    eu, edeg, V = edges_from_facets(fac)
    assert V == 64
    # f1 = V * CN / 2 with CN = 13.5 for A15
    assert len(eu) == round(64 * 13.5 / 2)
    # each tet contributes 6 edge slots; degrees sum to 6 f3
    assert edeg.sum() == 6 * len(fac)


def test_vertex_classes_matches_census(a15):
    fac, _ = a15
    eu, edeg, V = edges_from_facets(fac)
    n6, imp, adj = vertex_classes(fac)
    fz, _ = vertex_class_census(eu, edeg, V)
    for name, k in (("Z12", 0), ("Z14", 2), ("Z15", 3), ("Z16", 4)):
        assert abs(float(np.mean((imp == 0) & (n6 == k))) - fz[name]) < 1e-12
    # adjacency degrees consistent with the vertex degree law deg v = 12 + n6
    assert all(len(adj[v]) == 12 + n6[v] for v in range(V) if imp[v] == 0)


def test_skeleton_stats_a15(a15):
    fac, _ = a15
    eu, edeg, V = edges_from_facets(fac)
    st = skeleton_stats(eu, edeg, V)
    # A15 six-web: 6 six-edge slots per cell (Z14_6 with n6=2, each edge
    # counted at both ends), all lines closed
    assert st["n_edges"] == 6 * 8
    assert st["frac_val1"] == 0.0


def test_fast_census_vs_exact_links_on_churned_state():
    """Every vertex of a defect-carrying state, classified two ways."""
    fac, n = tr.build_t3_triangulation("c15", 2)
    mfd = ddg.Manifold(3, np.asarray(fac).tolist())
    qbar = mfd.mean_degree(1)
    params = ddg.SamplerParams(
        num_facets_target=mfd.num_facets, num_facets_coef=0.1,
        hinge_degree_target=qbar, num_hinges_coef=2.0,
        hinge_degree_target_coef=0.4 * qbar / 6.0)
    ddg.set_random_seed(7)
    s = ddg.ManifoldSampler(mfd, params)
    s.run(sweeps=2)                      # enough churn to mint defects
    F = np.asarray(s.manifold.facets(), np.int64)

    eu, edeg, V = edges_from_facets(F)
    labels, info = lc.link_class_census(F, eu, edeg, V)
    assert set(np.unique(labels).tolist()) - {0, 1, 2, 3} , \
        "churn produced no defects; increase sweeps/coupling"

    # exact oracle: extract and classify every link from scratch
    Frel, _ = lc.relabel(F)
    lab_to_name = dict(enumerate(lc.CLASS_NAMES))
    for v in range(V):
        name, n6, cof = lc.classify_link_exact(lc.extract_link(Frel, v))
        want = lab_to_name[int(labels[v])]
        got = name
        if got.startswith("Z") and got not in lc.CLASS_NAMES:
            got = "Z17plus"
        assert got == want, (v, name, want)
        if cof is not None:
            assert cof == info["cofacial"][v], (v, cof, info["cofacial"][v])
