"""Ground-truth tests for crystalline-grain / defect detection (crystal_grains).

The cleanest possible test of the defect finder: take a PERFECT crystal, apply a
single Pachner (bistellar) move, and assert that the defect set is *exactly* the
set of vertices the move touched -- no more (the rest of the crystal stays
crystalline), no fewer (every touched vertex loses its crystallinity).

A perfect crystal has all edge degrees in {5,6}, so only the expansive moves are
applicable: 1->4 (subdivide a tet: touches its 4 vertices + 1 new vertex) and
2->3 (bipyramid flip: touches the shared face's 3 vertices + the 2 apexes). Both
touch exactly 5 vertices, all of which must become defects.

References are generated in-memory from Wyckoff positions (tcp_reference), so the
test needs no data files and no sampler.
"""
import itertools
import os
import sys
from functools import lru_cache

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in ("scripts", "tools"):
    sys.path.insert(0, os.path.join(_ROOT, _p))

import crystal_grains as cg
import tcp_reference as tr

# (phase, supercell) -- smallest quotient that validates, kept tiny for speed.
CRYSTALS = [("a15", 2), ("c15", 2), ("r", 2)]


@lru_cache(maxsize=None)
def _reference(name, m):
    """(facets, refs, idx, ns_of) for a freshly built perfect crystal, matched
    to itself so the perfect crystal develops at identity registry (100%)."""
    facets = np.asarray(tr.build_t3_triangulation(name, m)[0])
    refs = {name: cg.build_struct(facets)}
    return facets, refs, cg.ref_index(refs), {name: len(tr.STRUCTURES[name][1])}


def _defects(facets, refs, idx, ns_of):
    """Set of vertices NOT interior-crystalline (the defect set) for `facets`
    analyzed against `refs`. Vertex ids are the crystal's own dense labels."""
    st = cg.build_struct(np.asarray(facets))
    grain_of_tet, sig_of_tet, phase_of_grain = cg.find_grains(st, refs, idx)
    interior = cg.interior_vertices(st, grain_of_tet, sig_of_tet,
                                    phase_of_grain, ns_of)
    all_v = {int(v) for f in facets for v in f}
    return all_v - set(interior)


def _perturb_1to4(facets, tet_index=0):
    """Subdivide one tet by a fresh central vertex (always a valid 1->4 move).
    Returns (new_facets, touched_vertices)."""
    facets = [tuple(int(x) for x in f) for f in facets]
    a, b, c, d = facets[tet_index]
    nv = max(v for f in facets for v in f) + 1
    new = [f for i, f in enumerate(facets) if i != tet_index]
    new += [(nv, b, c, d), (a, nv, c, d), (a, b, nv, d), (a, b, c, nv)]
    return np.array(new), {a, b, c, d, nv}


def _perturb_2to3(facets):
    """Flip the first face whose two apexes are not already adjacent (a valid
    2->3 move). Returns (new_facets, touched) or (None, None) if none exists."""
    facets = [tuple(int(x) for x in f) for f in facets]
    edges = set()
    face2 = {}
    for i, tet in enumerate(facets):
        for a, b in itertools.combinations(tet, 2):
            edges.add((a, b) if a < b else (b, a))
        for tri in itertools.combinations(tet, 3):
            apex = next(v for v in tet if v not in tri)
            face2.setdefault(tuple(sorted(tri)), []).append((i, apex))
    for tri, tl in face2.items():
        if len(tl) != 2:
            continue
        (i1, p), (i2, q) = tl
        e = (p, q) if p < q else (q, p)
        if p == q or e in edges:
            continue
        f0, f1, f2 = tri
        new = [f for k, f in enumerate(facets) if k not in (i1, i2)]
        new += [(p, q, f0, f1), (p, q, f1, f2), (p, q, f0, f2)]
        return np.array(new), {p, q, f0, f1, f2}
    return None, None


@pytest.mark.parametrize("name,m", CRYSTALS)
def test_perfect_crystal_has_no_defects(name, m):
    facets, refs, idx, ns = _reference(name, m)
    assert _defects(facets, refs, idx, ns) == set()


@pytest.mark.parametrize("name,m", CRYSTALS)
def test_1to4_move_defects_are_exactly_touched(name, m):
    facets, refs, idx, ns = _reference(name, m)
    perturbed, touched = _perturb_1to4(facets)
    assert _defects(perturbed, refs, idx, ns) == touched


@pytest.mark.parametrize("name,m", CRYSTALS)
def test_2to3_move_defects_are_exactly_touched(name, m):
    facets, refs, idx, ns = _reference(name, m)
    perturbed, touched = _perturb_2to3(facets)
    if perturbed is None:
        pytest.skip(f"{name}: no 2->3 move with non-adjacent apexes")
    assert _defects(perturbed, refs, idx, ns) == touched


@pytest.mark.parametrize("name,m", CRYSTALS)
def test_two_disjoint_moves_give_two_defect_clusters(name, m):
    """Two well-separated moves -> the touched vertices of both are defects,
    and nothing else is (the disconnected-crystal-domain scenario)."""
    facets, refs, idx, ns = _reference(name, m)
    p1, t1 = _perturb_1to4(facets, tet_index=0)
    # a second 1->4 on a tet sharing no vertex with the first
    far = next(i for i, f in enumerate(np.asarray(p1))
               if not (set(int(x) for x in f) & t1))
    p2, t2 = _perturb_1to4(p1, tet_index=far)
    assert _defects(p2, refs, idx, ns) == (t1 | t2)


def test_phase_discrimination_correct_phase_dominates():
    """Perfect C15 develops fully against a C15 reference, but only marginally
    against an R reference (finite regions are phase-ambiguous, but the correct
    phase must dominate by a wide margin)."""
    c15, c15_refs, c15_idx, c15_ns = _reference("c15", 2)
    V = len({int(v) for f in c15 for v in f})
    assert _defects(c15, c15_refs, c15_idx, c15_ns) == set()   # 100% c15

    _, r_refs, r_idx, r_ns = _reference("r", 2)
    wrong = V - len(_defects(c15, r_refs, r_idx, r_ns))        # interior vs R
    assert wrong < 0.25 * V
