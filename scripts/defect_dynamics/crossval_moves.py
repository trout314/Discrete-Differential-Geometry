#!/usr/bin/env python3
"""Cross-validate worm_moves.py (pure-python move arithmetic) against the D
core's doMove/doHingeMove via the targeted-move C API (do_bistellar_move /
do_hinge_move, added 2026-07-24).

For NSITES random 2-3 sites of the r reference crystal:
  * apply the 2-3 both ways; facet sets must agree exactly;
  * on the defected state, apply every valid 3-2 and 4-4 both ways; agree.
Plus negative controls: malformed moves must be rejected by the D-side
validation (which cannot rely on asserts -- they are compiled out in release).

Run after touching either implementation.  Exit 0 = all exact.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs

NSITES = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def facset(F):
    return {tuple(sorted(int(x) for x in t)) for t in F}


def facset_m(m):
    return {tuple(sorted(int(x) for x in t)) for t in np.asarray(m.facets())}


path = best_refs(REF_GLOB)["r"]
F0 = np.asarray(ddg.Manifold.load(path, 3).facets())
faces, edeg, vedges = wm.build_tables(F0)
sites = list(wm.two_three_sites(F0, faces, edeg))
rng = np.random.default_rng(7)
pick = rng.choice(len(sites), NSITES, replace=False)

n23 = n32 = n44 = 0
for si in pick:
    face, d, e, valid = sites[si]
    assert valid
    Fpy = wm.apply_two_three(F0, faces, face, d, e)
    m = ddg.Manifold.load(path, 3)
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

# negative controls
m = ddg.Manifold.load(path, 3)
face, d, e, _ = sites[0]
a, b, c = sorted(face)
assert not m.has_bistellar_move([a, b, c], [d, a]), "shared vertex accepted"
assert not m.has_bistellar_move([a, b], [c, d, e]), "star mismatch accepted"
try:
    m.do_bistellar_move([a, b], [c, d, e])
    raise SystemExit("FAIL: invalid move was not rejected")
except RuntimeError:
    pass

print(f"CROSS-VALIDATED: {n23} 2-3, {n32} 3-2, {n44} 4-4 moves exact; "
      "invalid moves rejected")
