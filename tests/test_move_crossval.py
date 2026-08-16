"""Cross-validate ddg_lab.worm_moves (pure-python move arithmetic) against the
D core's doMove/doHingeMove via the targeted-move C API.

Pytest port of scripts/defect_dynamics/crossval_moves.py, on an in-memory
C15 m=2 crystal (the invariant is exact move arithmetic -- size- and
phase-independent):

  * apply a 2-3 both ways; facet sets must agree exactly;
  * on the defected state, apply every valid 3-2 and 4-4 both ways; agree;
  * negative controls: malformed moves must be REJECTED by the D-side
    validation (which cannot rely on asserts -- compiled out in release).
"""
import numpy as np
import pytest

import discrete_differential_geometry as ddg
from discrete_differential_geometry.tcp_reference import build_t3_triangulation
from ddg_lab import worm_moves as wm

NSITES = 5


def facset(F):
    return {tuple(sorted(int(x) for x in t)) for t in F}


def facset_m(m):
    return {tuple(sorted(int(x) for x in t)) for t in np.asarray(m.facets())}


@pytest.fixture(scope="module")
def crystal():
    fac, _ = build_t3_triangulation("c15", 2)
    F0 = np.asarray(fac)
    faces, edeg, vedges = wm.build_tables(F0)
    sites = list(wm.two_three_sites(F0, faces, edeg))
    return F0, faces, sites


def _fresh(F0):
    return ddg.Manifold(3, F0.tolist())


def test_moves_agree_exactly(crystal):
    F0, faces, sites = crystal
    rng = np.random.default_rng(7)
    pick = rng.choice(len(sites), NSITES, replace=False)
    n23 = n32 = n44 = 0
    for si in pick:
        face, d, e, valid = sites[si]
        assert valid
        Fpy = wm.apply_two_three(F0, faces, face, d, e)
        m = _fresh(F0)
        assert m.has_bistellar_move(sorted(face), [d, e])
        m.do_bistellar_move(sorted(face), [d, e])
        assert facset(Fpy) == facset_m(m), "2-3 MISMATCH"
        n23 += 1
        faces2, edeg2, _ = wm.build_tables(Fpy)
        for edge, link, v2 in wm.three_two_sites(Fpy, faces2, edeg2):
            if not v2:
                continue
            Fpy2 = wm.apply_three_two(Fpy, faces2, edge, link)
            m2 = ddg.Manifold(3, Fpy)
            assert m2.has_bistellar_move(sorted(edge), sorted(link))
            m2.do_bistellar_move(sorted(edge), sorted(link))
            assert facset(Fpy2) == facset_m(m2), "3-2 MISMATCH"
            n32 += 1
        for edge, cyc, diag, v4 in wm.four_four_sites(Fpy, faces2, edeg2):
            if not v4:
                continue
            Fpy2 = wm.apply_four_four(Fpy, faces2, edge, cyc, diag)
            m2 = ddg.Manifold(3, Fpy)
            dg = 0 if set([cyc[0], cyc[2]]) == set(diag) else 1
            assert m2.has_hinge_move(sorted(edge), cyc, dg)
            m2.do_hinge_move(sorted(edge), cyc, dg)
            assert facset(Fpy2) == facset_m(m2), "4-4 MISMATCH"
            n44 += 1
    assert n23 == NSITES and n32 > 0 and n44 > 0


def test_invalid_moves_rejected(crystal):
    F0, faces, sites = crystal
    m = _fresh(F0)
    face, d, e, _ = sites[0]
    a, b, c = sorted(face)
    assert not m.has_bistellar_move([a, b, c], [d, a]), "shared vertex accepted"
    assert not m.has_bistellar_move([a, b], [c, d, e]), "star mismatch accepted"
    with pytest.raises(RuntimeError):
        m.do_bistellar_move([a, b], [c, d, e])
