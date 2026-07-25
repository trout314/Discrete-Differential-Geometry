#!/usr/bin/env python3
"""The KNOT-SLIDE move: translate a (3,4,4) knot 4 positions along its local
BC (Boerdijk-Coxeter) chain in 4 bistellar moves -- the elementary transport
operator, packaged as a proposable sampler move.

Frame (from worm_helix motif, verified move-by-move 2026-07-25): the knot's
deg-3 edge is the CHORD (v_b, v_b+4); its link is {v_b+1, v_b+2, v_b+3}. The
slide template, with c_i = chain vertex v_b+i derived locally:

    M1  3-2 on chord (c0, c4), link {l1,l2,l3}
    M2  2-3 on face (c3, c4, c5), apexes (c2, c6)
    M3  2-3 on face (c6, c5, c7), apexes (c4, c8)
    M4  3-2 on chord (c2, c6)

after which the knot chord is (c4, c8): the knot moved +4 chain steps. The
chain continuation derives LOCALLY by the sliding-window rule (face_apexes
with exclusion) -- no global chain needed:

    c5 = apex of (c2,c3,c4) != c0,  c6 = apex of (c3,c4,c5) != c2,
    c7 = apex of (c4,c5,c6) != c3,  c8 = apex of (c5,c6,c7) != c4.

A direction candidate = choice of chord orientation (2) x ordered pair
(c2,c3) from the link (6) = up to 12 per knot; each validates (or not) by
existence + apex checks + has_bistellar_move. Invalid in a thermal
environment => proposal rejected (correct MH behavior). The inverse slide is
the mirrored template (alphabet inverse-closed).

Modes:
  --crystal-test         create a knot on the reference crystal, slide +N,
                         slide -N, close; assert EXACT facet restoration.
  --thermal SNAP.mfd     count, for every deg-3 edge of the state, the valid
                         slide directions + exact dS -- the mobilization
                         measurement.
"""
import argparse
import os
import sys
from collections import defaultdict
from itertools import permutations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs

ESTAR = 5.105025


def edge_link_verts(m, a, b):
    try:
        lk = m.edge_link(a, b).tolist()
    except RuntimeError:
        return []
    return sorted({x for pr in lk for x in pr})


def apex_excluding(m, face, excl):
    """The apex of `face` that is not `excl`; None if face missing or excl
    is not one of its apexes (frame assumption broken)."""
    try:
        a1, a2 = m.face_apexes(*face)
    except RuntimeError:
        return None
    if a1 == excl:
        return a2
    if a2 == excl:
        return a1
    return None


def derive_frame(m, c0, c4, c2, c3):
    """Given chord (c0,c4) and ordered link pick (c2,c3), derive c5..c8 by
    the sliding-window rule. Returns dict or None."""
    c5 = apex_excluding(m, (c2, c3, c4), c0)
    if c5 is None:
        return None
    c6 = apex_excluding(m, (c3, c4, c5), c2)
    if c6 is None:
        return None
    c7 = apex_excluding(m, (c4, c5, c6), c3)
    if c7 is None:
        return None
    c8 = apex_excluding(m, (c5, c6, c7), c4)
    if c8 is None:
        return None
    cs = dict(c0=c0, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8)
    if len(set(cs.values())) != len(cs):
        return None                          # degenerate frame
    return cs


def slide_moves(m, cs):
    """The 4-move template for frame cs, with full validity checking.
    Returns list of ('23'|'32', center, cocenter) or None."""
    mv = []
    l1 = edge_link_verts(m, cs["c0"], cs["c4"])
    if len(l1) != 3 or not m.has_bistellar_move([cs["c0"], cs["c4"]], l1):
        return None
    mv.append(("32", [cs["c0"], cs["c4"]], l1))
    f2 = (cs["c3"], cs["c4"], cs["c5"])
    ap2 = (cs["c2"], cs["c6"])
    f3 = (cs["c6"], cs["c5"], cs["c7"])
    ap3 = (cs["c4"], cs["c8"])
    mv.append(("23", sorted(f2), list(ap2)))
    mv.append(("23", sorted(f3), list(ap3)))
    mv.append(("32", [cs["c2"], cs["c6"]], None))   # link resolved at apply
    return mv


def apply_slide(m, mv, dv=None):
    """Apply the 4 moves, checking validity stepwise; returns applied records
    for undo, or None (with automatic rollback) if any step invalid."""
    recs = []

    def undo_all():
        for kind, cen, coc in reversed(recs):
            m.do_bistellar_move(coc, cen)   # inverse: center/cocenter swap
        return None

    for i, (kind, cen, coc) in enumerate(mv):
        if kind == "32":
            if coc is None:
                coc = edge_link_verts(m, cen[0], cen[1])
                if len(coc) != 3:
                    return undo_all()
            if not m.has_bistellar_move(cen, coc):
                return undo_all()
            m.do_bistellar_move(cen, coc)
            recs.append(("32", cen, coc))
        else:
            try:
                a1, a2 = m.face_apexes(*cen)
            except RuntimeError:
                return undo_all()
            if {a1, a2} != set(coc):
                return undo_all()
            if not m.has_bistellar_move(cen, coc):
                return undo_all()
            m.do_bistellar_move(cen, coc)
            recs.append(("23", cen, coc))
    return recs


def undo_slide(m, recs):
    for kind, cen, coc in reversed(recs):
        m.do_bistellar_move(coc, cen)   # inverse: center/cocenter swap


def knot_chords(m):
    """All deg-3 edges with their complex signature sizes."""
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    e3 = [e for e, d in ill.items() if d == 3]
    return e3, ill


def candidate_slides(m, chord):
    """All valid slide frames from this chord (both orientations x 6 link
    orderings), deduplicated by resulting move list."""
    out = []
    link = edge_link_verts(m, *chord)
    if len(link) != 3:
        return out
    for c0, c4 in (chord, chord[::-1]):
        for c2, c3 in permutations(link, 2):
            cs = derive_frame(m, c0, c4, c2, c3)
            if cs is None:
                continue
            mv = slide_moves(m, cs)
            if mv is not None:
                out.append((cs, mv))
    return out


def facset(m):
    A = np.sort(np.asarray(m.facets()), axis=1)
    return frozenset(map(tuple, A.tolist()))


def dS_between(em_before, em_after):
    x = ESTAR - int(ESTAR)
    dS = 0.0
    keys = set(em_before) | set(em_after)
    for k in keys:
        d0, d1 = em_before.get(k), em_after.get(k)
        if d0 == d1:
            continue
        if d0 is not None:
            dS -= (ESTAR / 6.0) * ((d0 - ESTAR) ** 2 - x * (1 - x))
        if d1 is not None:
            dS += (ESTAR / 6.0) * ((d1 - ESTAR) ** 2 - x * (1 - x))
    return dS


def edeg_dict(m):
    eu, cnt = wm.edge_degrees_np(np.asarray(m.facets()))
    return {(int(a), int(b)): int(c) for (a, b), c in zip(eu, cnt)}


def crystal_test(args):
    path = args.ref or best_refs(REF_GLOB)["r"]
    m = ddg.Manifold.load(path, 3)
    fs_crystal = facset(m)
    F = np.asarray(m.facets())
    faces0, edeg0, vedges0 = wm.build_tables(F)
    # create a knot on a (5,5,6) face (clean (3,4,4) creation)
    site = None
    for face, d, e, valid in wm.two_three_sites(F, faces0, edeg0):
        if not valid:
            continue
        degs = sorted(edeg0[frozenset(p)]
                      for p in [(sorted(face)[0], sorted(face)[1]),
                                (sorted(face)[1], sorted(face)[2]),
                                (sorted(face)[2], sorted(face)[0])])
        if degs == [5, 5, 6]:
            deltas = wm.delta_two_three(face, d, e, edeg0)
            if sorted(nw for _, nw in deltas.values()
                      if nw is not None and nw not in (5, 6)) == [3, 4, 4]:
                site = (face, d, e)
                break
    assert site is not None
    face, d, e = site
    m.do_bistellar_move(sorted(face), [d, e])
    e3, ill = knot_chords(m)
    print(f"created knot: chord {e3[0]}, illegal {sorted(ill.values())}")

    chord = e3[0]
    trail = []
    clean_frac = []
    for step in range(args.steps):
        cands = candidate_slides(m, chord)
        done = None
        nclean = 0
        for cs, mv in cands:
            recs = apply_slide(m, mv)
            if recs is None:
                continue
            want = tuple(sorted((cs["c4"], cs["c8"])))
            e3n, illn = knot_chords(m)
            clean = (sorted(illn.values()) == [3, 4, 4]
                     and tuple(sorted(e3n[0])) == want)
            if clean:
                done = (recs, want)
                nclean += 1
                break                        # keep it applied; stop scanning
            undo_slide(m, recs)
        if done is None:
            print(f"stuck at step {step} ({len(cands)} template-valid, "
                  f"0 clean)")
            return
        # re-apply the chosen one if we undid it while counting
        clean_frac.append((nclean, len(cands)))
        recs, chord = done
        trail.append(recs)
    cf = np.array([c for c, _ in clean_frac])
    tv = np.array([t for _, t in clean_frac])
    print(f"slid +{args.steps} steps (chord now {chord}); species (3,4,4) "
          f"preserved every step")
    print(f"  per step: template-valid {tv.mean():.1f}, "
          f"clean {cf.mean():.1f} directions")
    print("sliding back...")
    for recs in reversed(trail):
        undo_slide(m, recs)
    e3, ill = knot_chords(m)
    print(f"returned: chord {e3[0]}")
    # close and verify exact restoration
    link = edge_link_verts(m, *e3[0])
    m.do_bistellar_move(list(e3[0]), link)
    pairs, _ = m.illegal_edges()
    assert len(pairs) == 0
    assert facset(m) == fs_crystal, "crystal NOT exactly restored"
    print("closed: crystal restored EXACTLY (facet-set equality). "
          f"Slide is invertible transport over {args.steps} steps.")


def thermal_test(args):
    m = ddg.Manifold.load(args.thermal, 3)
    em = edeg_dict(m)
    e3, ill = knot_chords(m)
    print(f"state: {os.path.basename(args.thermal)}  "
          f"illegal edges {len(ill)}  deg-3 chords {len(e3)}")
    n_mobile = 0
    dSs = []
    per = []
    for chord in e3:
        cands = candidate_slides(m, chord)
        ok = 0
        for cs, mv in cands:
            recs = apply_slide(m, mv)
            if recs is None:
                continue
            em2 = edeg_dict(m)
            dSs.append(dS_between(em, em2))
            undo_slide(m, recs)
            ok += 1
        per.append(ok)
        if ok:
            n_mobile += 1
    print(f"chords with >=1 valid slide: {n_mobile}/{len(e3)}")
    if per:
        print(f"valid directions per chord: "
              f"{np.bincount(per, minlength=1).tolist()} (histogram)")
    if dSs:
        dSs = np.array(dSs)
        print(f"slide dS_shape (lam=1 units): mean {dSs.mean():+.3f}  "
              f"med {np.median(dSs):+.3f}  min {dSs.min():+.3f}  "
              f"max {dSs.max():+.3f}  (n={len(dSs)})")
        print(f"  acceptance at lam=0.4 (exp(-0.4 dS)): "
              f"med {np.exp(-0.4 * np.median(dSs)):.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--crystal-test", action="store_true")
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--thermal", default=None)
    args = ap.parse_args()
    if args.crystal_test:
        crystal_test(args)
    if args.thermal:
        thermal_test(args)


if __name__ == "__main__":
    main()
