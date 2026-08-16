"""Tests for the exact 5/6-link classification (ddg.link_classes).

The two 16-vertex spheres below are the heart of the module: n_6 = 4 admits
exactly two 5/6-triangulated 2-spheres -- the T_d Friauf polyhedron (the
Frank-Kasper Z16, six-edges pairwise non-cofacial) and a D2 isomer with two
cofacial six-pairs that no FK phase admits. They are hard-coded from the
exhaustive enumeration so the fast suite can check the discriminator without
re-running the (slow) search; the slow tests re-derive them.
"""
import numpy as np
import pytest

from discrete_differential_geometry import link_classes as lc

# enumerate_5_6_spheres(4) output, verified against CLASS_COUNTS[4] == 2.
Z16_TD = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5), (1, 2, 6),
          (2, 3, 7), (3, 4, 8), (4, 5, 9), (1, 5, 10), (1, 6, 10), (2, 6, 7),
          (6, 7, 11), (3, 7, 12), (3, 8, 12), (7, 11, 12), (8, 12, 13),
          (11, 12, 13), (11, 13, 14), (8, 9, 13), (5, 10, 15), (5, 9, 15),
          (10, 14, 15), (6, 10, 14), (6, 11, 14), (13, 14, 15), (4, 8, 9),
          (9, 13, 15)]
Z16_D2 = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5), (1, 2, 6),
          (2, 3, 7), (3, 4, 8), (4, 5, 9), (1, 5, 10), (1, 6, 10), (2, 6, 7),
          (6, 7, 11), (3, 7, 8), (7, 8, 11), (8, 11, 12), (6, 11, 13),
          (6, 10, 13), (5, 9, 10), (9, 10, 14), (10, 13, 14), (9, 14, 15),
          (12, 13, 14), (12, 14, 15), (4, 9, 15), (4, 8, 15), (8, 12, 15),
          (11, 12, 13)]


def _check_sphere(tris, n_verts):
    """Closed triangulated surface with V - E + F = 2, all degrees in {5,6}."""
    deg, adj = lc.link_degree_map(tris)
    assert len(deg) == n_verts
    edges = {frozenset((a, b)) for t in tris for a, b in
             ((t[0], t[1]), (t[1], t[2]), (t[0], t[2]))}
    assert len(deg) - len(edges) + len(tris) == 2
    assert set(deg.values()) <= {5, 6}
    return deg


def test_hardcoded_z16_pair():
    for tris, want_name, want_cof in ((Z16_TD, "Z16", 0), (Z16_D2, "Z16_D2", 2)):
        deg = _check_sphere(tris, 16)
        assert sum(1 for d in deg.values() if d == 6) == 4
        name, n6, cof = lc.classify_link_exact(tris)
        assert (name, n6, cof) == (want_name, 4, want_cof)


def test_fk_names_and_caps():
    assert lc.FK_NAMES == {0: "Z12", 2: "Z14", 3: "Z15", 4: "Z16"}
    assert lc.MAX_FK_N6 == 4
    assert lc.CLASS_COUNTS[1] == 0          # n_6 = 1 is IMPOSSIBLE
    assert lc.CLASS_COUNTS[4] == 2          # T_d + D2


def test_summarize_arithmetic():
    labels = np.array([0, 0, 1, 2, 3, 4, 5, 6], dtype=np.int8)
    s = lc.summarize(labels)
    assert s["Z12"] == 2 / 8 and s["Z16_D2"] == 1 / 8
    assert s["n_Z16_D2"] == 1
    assert abs(s["f_FK_strict"] - 5 / 8) < 1e-12
    assert abs(s["f_FK"] - 6 / 8) < 1e-12   # n_6 bucketing counts the D2


@pytest.mark.parametrize("k", [0, 1, 2])
def test_enumeration_small_k(k):
    reps = lc.enumerate_5_6_spheres(k)
    assert len(reps) == lc.CLASS_COUNTS[k]
    for s in reps:
        deg = _check_sphere(s, 12 + k)
        assert sum(1 for d in deg.values() if d == 6) == k


@pytest.mark.slow
@pytest.mark.parametrize("k", [3, 4])
def test_enumeration_larger_k(k):
    reps = lc.enumerate_5_6_spheres(k)
    assert len(reps) == lc.CLASS_COUNTS[k]
    if k == 4:
        got = sorted(lc.classify_link_exact(s)[0] for s in reps)
        assert got == ["Z16", "Z16_D2"]
