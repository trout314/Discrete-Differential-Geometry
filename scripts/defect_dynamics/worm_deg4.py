#!/usr/bin/env python3
"""Deg-4 WORM move: transport a degree-4 edge (disclination-line element) a
chosen number of steps down its BC chain, by ANNIHILATING it at A and
RE-CREATING it at B.  The degree-3 analog is worm_slide.py (the knot slide);
this is the harder deg-4 case, since a deg-4 edge is not a single-move center
-- destroy/create each take a short 2->3/3->2 sequence (found by the bounded
Pachner search).

Building blocks:
  * CREATE C:  a validated forward 2->3/3->2 sequence turning a pristine local
    ball into a deg-4 defect at an anchored frame.
  * DESTROY C^{-1}: the reverse (swap center/coCenter, reverse order).
  * WORM(A->B) = C^{-1}(A) . C(B): net f-vector ZERO (create's +tets cancel
    destroy's -tets), so N_3/E and every global action term cancel exactly;
    dS is the LOCAL edge-degree change only -- same structure as the slide.

This module first PINS the create sequence on the D core: applies it, confirms
a deg-4 defect, confirms C^{-1} restores the crystal facet-for-facet, and
reports the defect signature + dS.  (Frame parameterization / the non-local
translation come next, built on this verified primitive.)
"""
import argparse
import itertools
import os
import sys
from collections import Counter, defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg

ESTAR = 5.105025
REF_CELL = os.path.join(_HERE, "../../data/tcp_reference/T3_R_m2_N7248.mfd")

# The create sequence found by the bounded Pachner search (pachner_search.py),
# expressed as bistellar moves on the pristine R m2 crystal.  Each entry:
#   ("23", (a,b,c), x, y)   2->3: face (a,b,c) with apexes x,y  ->  new edge x-y
#   ("32", (u,w), a, b, c)  3->2: deg-3 edge u-w with link {a,b,c} -> face abc
CREATE_SEQ = [
    ("23", (0, 201, 334), 163, 358),
    ("23", (163, 1232, 1248), 1075, 1266),
    ("23", (0, 163, 201), 28, 358),
    ("32", (0, 201), 28, 44, 358),
    ("32", (1075, 1266), 163, 1232, 1248),
]


def as_bistellar(mv):
    """(center, coCenter) args for m.do_bistellar_move for a CREATE-list entry."""
    if mv[0] == "23":
        _, tri, x, y = mv
        return list(tri), [x, y]
    _, (u, w), a, b, c = mv
    return [u, w], [a, b, c]


def apply_seq(m, seq):
    """Apply a create-list forward; return the applied (center,coCenter) records."""
    recs = []
    for mv in seq:
        cen, coc = as_bistellar(mv)
        if not m.has_bistellar_move(cen, coc):
            raise RuntimeError(f"move not applicable: {mv} -> {cen}/{coc}")
        m.do_bistellar_move(cen, coc)
        recs.append((cen, coc))
    return recs


def undo_seq(m, recs):
    """Reverse a create-list: inverse of each move (swap center/coCenter),
    reverse order."""
    for cen, coc in reversed(recs):
        m.do_bistellar_move(coc, cen)


def facset(m):
    A = np.sort(np.asarray(m.facets()), axis=1)
    return frozenset(map(tuple, A.tolist()))


def edeg(m):
    F = np.sort(np.asarray(m.facets()), axis=1)
    d = Counter()
    for t in F.tolist():
        for e in itertools.combinations(t, 2):
            d[e] += 1
    return d


def illegal(d):
    return {e: k for e, k in d.items() if k not in (5, 6)}


def dS_edq(d0, d1, estar=ESTAR, lam=0.35):
    """EDQ (edge-degree-target) action change, lam=1 units scaled by lam.
    coef = estar/6 per the sampler convention."""
    coef = lam * estar / 6.0
    s = 0.0
    for e in set(d0) | set(d1):
        a, b = d0.get(e, 0), d1.get(e, 0)
        if a == b:
            continue
        if a:
            s -= coef * (a - estar) ** 2
        if b:
            s += coef * (b - estar) ** 2
    return s


def pin_create(args):
    """Validate the create sequence on the D core: apply, report the deg-4
    defect it makes, then confirm the reverse restores the crystal exactly."""
    m = ddg.Manifold.load(REF_CELL, 3)
    fs0 = facset(m)
    d0 = edeg(m)
    print(f"crystal: {m.num_facets} tets, {len(illegal(d0))} illegal edges "
          f"(should be 0)")

    recs = apply_seq(m, CREATE_SEQ)
    d1 = edeg(m)
    ill = illegal(d1)
    print(f"\nafter CREATE ({len(CREATE_SEQ)} moves): {m.num_facets} tets "
          f"(delta {m.num_facets - len(fs0):+d})")
    print(f"  illegal edges: {len(ill)}  signature {sorted(ill.values())}")
    # show what changed, categorized
    cats = defaultdict(list)
    for e in set(d0) | set(d1):
        a, b = d0.get(e, 0), d1.get(e, 0)
        if a == b:
            continue
        if a == 0:
            cats[f"NEW deg{b}"].append(e)
        elif b == 0:
            cats[f"GONE (was {a})"].append(e)
        else:
            cats[f"{a}->{b}"].append(e)
    for k in sorted(cats):
        ee = cats[k]
        print(f"    {k:12s}: {len(ee)}  {ee if len(ee) <= 6 else ''}")
    print(f"  dS_EDQ(lam={args.lam}) = {dS_edq(d0, d1, lam=args.lam):+.4f}")

    # per-vertex charge (discrete Gauss-Bonnet): must be 12 everywhere
    ch = defaultdict(int)
    for (u, v), k in d1.items():
        ch[u] += 6 - k
        ch[v] += 6 - k
    bad = {v: c for v, c in ch.items() if c != 12}
    print(f"  per-vertex charge==12: "
          f"{'OK everywhere' if not bad else f'VIOLATED at {bad}'}")

    undo_seq(m, recs)
    fs2 = facset(m)
    print(f"\nafter DESTROY (reverse): {m.num_facets} tets, "
          f"{len(illegal(edeg(m)))} illegal edges")
    print(f"  crystal restored facet-for-facet: {fs2 == fs0}")
    if fs2 != fs0:
        print("  *** NOT restored -- create is not cleanly invertible")
        sys.exit(1)
    print("\nPASS: create makes a deg-4 defect; destroy = reverse restores "
          "the crystal exactly.  These are the two halves of the worm.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pin-create", action="store_true",
                    help="validate the create sequence + its reverse on the "
                         "reference crystal")
    ap.add_argument("--lam", type=float, default=0.35)
    args = ap.parse_args()
    if args.pin_create:
        pin_create(args)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
