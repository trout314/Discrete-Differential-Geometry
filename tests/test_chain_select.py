"""Tests for BC chain class selection and provenance (tools/chain_select.py).

The critical property is LEGACY EQUIVALENCE: converting a driver from
``bc_orbit(m, F[seed_tet])`` to ``chain_for_run(..., selector=None)`` must
hand back the identical vertex sequence, so nothing measured changes and the
only difference is that the chain's class is recorded. Everything else here
guards the selectors and the refusal paths.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
for _p in ("scripts", "tools", "scripts/defect_dynamics"):
    sys.path.insert(0, os.path.join(_ROOT, _p))

import tcp_reference as tr
from discrete_differential_geometry import Manifold
from chain_select import ChainClasses, chain_for_run


@pytest.fixture(scope="module")
def a15(tmp_path_factory):
    """A perfect A15 m=2 crystal, saved under the tcp_reference naming so the
    position/winding path resolves."""
    facets = np.asarray(tr.build_t3_triangulation("a15", 2)[0])
    d = tmp_path_factory.mktemp("refs")
    path = str(d / f"T3_A15_m2_N{len(facets)}.mfd")
    Manifold(3, facets.tolist()).save(path)
    return path, facets


def test_matches_bc_orbit_for_every_seed_tet(a15):
    """The whole conversion rests on this: reading the sequence off the
    precomputed cycle must equal re-walking it with bc_orbit."""
    from worm_helix import bc_orbit
    path, facets = a15
    m = Manifold.load(path, 3)
    cc = ChainClasses(path)
    for t in (0, 1, 5, len(facets) // 2, len(facets) - 1):
        window = [int(x) for x in facets[t]]
        assert cc.vertices_from_frame(window) == \
            [int(x) for x in bc_orbit(m, window)]


def test_frame_round_trips_through_bc_orbit(a15):
    from worm_helix import bc_orbit
    path, _ = a15
    m = Manifold.load(path, 3)
    cc = ChainClasses(path)
    for k in range(cc.n_classes):
        assert [int(x) for x in bc_orbit(m, list(cc.frame(k)))] == cc.vertices(k)


def test_legacy_selector_changes_nothing_but_records_the_class(a15):
    from worm_helix import bc_orbit
    path, facets = a15
    m = Manifold.load(path, 3)
    for t in (0, 3, 11):
        _, k, seq, prov = chain_for_run(path, facets, None, seed_tet=t,
                                        verbose=False)
        assert seq == [int(x) for x in bc_orbit(m, [int(x) for x in facets[t]])]
        assert prov["chain_class"] == k
        assert prov["chain_length_walked"] == len(seq)
        assert "legacy" in prov["selected_by"]


def test_seed_tet_was_silently_a_class_selector(a15):
    """Different seed tets land in different classes -- which is exactly why
    the class has to be recorded."""
    path, facets = a15
    cc = ChainClasses(path)
    hit = {cc.class_of_frame([int(x) for x in facets[t]])
           for t in range(0, len(facets), 7)}
    assert len(hit) > 1


def test_selectors(a15):
    path, _ = a15
    cc = ChainClasses(path)
    assert cc.select("longest") == max(range(cc.n_classes), key=cc.length)
    assert cc.select("shortest") == min(range(cc.n_classes), key=cc.length)
    assert cc.select(0) == 0 and cc.select("0") == 0
    axis = cc.select("axis")
    assert axis in cc.pure_axis_classes()
    # 'axis' picks the shortest pure-axis class, deterministically
    assert cc.length(axis) == min(cc.length(k) for k in cc.pure_axis_classes())
    w = cc.windings(axis)[0]
    assert sum(1 for x in w if x) == 1
    assert cc.select(f"w={w[0]},{w[1]},{w[2]}") in range(cc.n_classes)
    with pytest.raises(ValueError):
        cc.select(cc.n_classes)
    with pytest.raises(ValueError):
        cc.select("nonsense")


def test_windings_are_integral_and_a_class_property(a15):
    path, _ = a15
    cc = ChainClasses(path)
    assert "certified" in cc.position_note
    for k in range(cc.n_classes):
        ws = cc.windings(k)
        assert ws and all(len(w) == 3 for w in ws)
        # the representative's own winding is one of the class's
        assert tuple(cc.representative_winding(k)) in set(ws)


def test_winding_selectors_refuse_without_certified_positions(tmp_path):
    """A file with no tcp_reference positions must raise on winding-based
    selection, not fall back to an arbitrary chain."""
    facets = np.asarray(tr.build_t3_triangulation("a15", 2)[0])
    path = str(tmp_path / "unnamed.mfd")          # not a tcp_reference name
    Manifold(3, facets.tolist()).save(path)
    cc = ChainClasses(path)
    assert cc.rp is None
    assert cc.select("shortest") in range(cc.n_classes)   # combinatorial: fine
    with pytest.raises(ValueError, match="positions"):
        cc.select("axis")
    with pytest.raises(ValueError, match="positions"):
        cc.windings(0)
