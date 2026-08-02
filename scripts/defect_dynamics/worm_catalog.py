#!/usr/bin/env python3
"""Worm-move catalog, stage 1: every 2-3 creation site of a reference crystal,
classified by EXACT symmetry orbit, with junction transitions, dS, exact
verification, and the available follow-up (3-2 / 4-4) moves per class.

Driver over the crystal-independent library `worm_moves.py` (which a future
exhaustive, non-crystal search will reuse; keep all move arithmetic THERE).

On a perfect FK crystal only 2-3 sites exist (no deg-3/4 edges), so stage 1 =
the pair-creation alphabet + the first walk steps available after each
creation. A 2-3 site IS a face (its two apexes are determined by the face), so
the inequivalent sites are exactly the FACE ORBITS of Aut(K), computed by
`discrete_differential_geometry.symmetry`.

WHY THIS REPLACED SIGNATURE BUCKETING (2026-08-02). This script used to bucket
sites by a canonical local signature -- per face vertex, its two face-edge
degrees, its (to-d, to-e) degrees and its (n6, m), canonicalised over S3 x S2.
That signature is an automorphism invariant, so each bucket is a UNION of
orbits and the bucket count is a LOWER bound. It is not a tight one: on the R
phase it gives 47 buckets where there are 102 orbits, and 34 of those buckets
each span several orbits. Measured consequences of the fusion:

  * dS_shape is CONSTANT on every fused bucket -- the local signature does
    determine the energetics, so no dS ever reported here was wrong;
  * the FOLLOW-UP SET is NOT. It differs between orbits in 28 of the 34 fused
    buckets. The old code computed follow-ups from `sites[0]`, one arbitrary
    representative, and reported them for the whole bucket -- so the
    "available follow-up moves", which is the entire point of stage 1, came
    from the wrong orbit in 28 of 47 rows.

Follow-ups are now computed PER ORBIT. The signature coarsening is still
reported, because it is real physics: it says the 5-vertex local data does not
determine the site, and names exactly where the exterior embedding takes over.

SUPERCELL INDEPENDENCE. Orbit counts are properties of the crystal, not of the
m x m x m box (R gives 11/61/102/53 at m = 2, 3 and 4 alike), and the per-orbit
data is local, so the alphabet computed on a small supercell is the same
alphabet. The default reference is therefore the SMALLEST available for the
structure rather than the largest -- exact verification rebuilds the whole edge
table per orbit, so this is the difference between seconds and many minutes.
`--cross-check REF` recomputes on a second supercell and asserts the per-orbit
alphabet agrees (verified m=2 vs m=3 for R).

Usage: worm_catalog.py [REF.mfd] [--full] [--json OUT] [--cross-check REF2]
"""
import argparse
import json
import os
import sys
from collections import defaultdict
from itertools import permutations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import CrystalSymmetry
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs

ESTAR = 5.105025


def all_refs(structure="r"):
    """Every reference .mfd for a structure, smallest supercell first."""
    import glob
    import re
    out = []
    for p in glob.glob(REF_GLOB):
        mm = re.match(r"T3_([A-Z0-9]+)_m(\d+)_N\d+\.mfd", os.path.basename(p))
        if mm and mm.group(1).lower() == structure:
            out.append((int(mm.group(2)), p))
    return [p for _, p in sorted(out)]


def canon_sig(face, d, e, edeg, vedges):
    """Canonical LOCAL signature of a 2-3 site under S3 (face) x S2 (apexes).

    An Aut invariant, so it can only fuse orbits, never split one. Retained as
    a reported coarsening of the exact classification -- NOT as the
    classification (see the module docstring).
    """
    abc = sorted(face)
    best = None
    for pa, pb, pc in permutations(abc):
        for ap0, ap1 in ((d, e), (e, d)):
            row = []
            for v in (pa, pb, pc):
                others = [x for x in (pa, pb, pc) if x != v]
                fe = tuple(sorted(edeg[frozenset((v, o))] for o in others))
                se = (edeg[frozenset((v, ap0))], edeg[frozenset((v, ap1))])
                row.append((fe, se, wm.vertex_counters(v, edeg, vedges)))
            apx = tuple(wm.vertex_counters(a, edeg, vedges) for a in (ap0, ap1))
            cand = (tuple(row), apx)
            if best is None or cand < best:
                best = cand
    return best


def ill_signature(deltas):
    """Sorted tuple of post-move illegal edge degrees."""
    return tuple(sorted(new for e, (old, new) in deltas.items()
                        if new is not None and new not in (5, 6)))


def describe_transitions(trans):
    return sorted(f"{wm.class_name(*before)}->{wm.class_name(*after)}"
                  for v, (before, after) in trans.items() if before != after)


def verify_exact(F, faces, edeg, vedges, face, d, e, deltas, trans):
    """Apply the move to the facet array and recheck every prediction."""
    F2 = wm.apply_two_three(F, faces, face, d, e)
    _, edeg2, vedges2 = wm.build_tables(F2)
    for ed, (old, new) in deltas.items():
        got = edeg2.get(ed)
        if new is None:
            assert got is None, f"edge {sorted(ed)} should be gone, deg {got}"
        else:
            assert got == new, f"edge {sorted(ed)}: predicted {new}, got {got}"
    for v, (_, after) in trans.items():
        got = wm.vertex_counters(v, edeg2, vedges2)
        assert got == after, f"vertex {v}: predicted {after}, got {got}"
    return F2, edeg2, vedges2


def followups(F2, edeg2, vedges2):
    """Available 3-2 / 4-4 moves in the post-creation state, with outcomes.

    Enumerates at the deg-3 / deg-4 edges read straight off `edeg2` (the site
    sets are defined by those degrees, so this is exactly the global
    enumeration) instead of rebuilding the whole face/edge table -- the old
    version called build_tables a second time per class, which dominated the
    runtime and now would be paid 102 times rather than 47.
    """
    e3 = [tuple(sorted(e)) for e, d in edeg2.items() if d == 3]
    e4 = [tuple(sorted(e)) for e, d in edeg2.items() if d == 4]
    out = []
    for edge, link, valid in wm.three_two_sites_at(F2, e3):
        dl = wm.delta_three_two(edge, link, edeg2)
        tr = wm.vertex_transitions(dl, edeg2, vedges2)
        out.append(("3-2", valid, ill_signature(dl),
                    round(wm.delta_S_shape(dl, tr, ESTAR), 4)))
    for edge, cyc, diag, valid in wm.four_four_sites_at(F2, e4):
        dl = wm.delta_four_four(edge, cyc, diag, edeg2)
        tr = wm.vertex_transitions(dl, edeg2, vedges2)
        out.append(("4-4", valid, ill_signature(dl),
                    round(wm.delta_S_shape(dl, tr, ESTAR), 4)))
    return out


def followup_key(fu):
    """Fingerprint of an orbit's follow-up alphabet (for comparing orbits)."""
    return (sum(1 for k, v, _, _ in fu if k == "3-2" and v),
            tuple(sorted((s, ds) for k, v, s, ds in fu if k == "4-4" and v)))


def build_catalog(path, verbose=True):
    """Per-orbit catalog rows for one reference crystal."""
    F = np.asarray(ddg.Manifold.load(path, 3).facets())
    faces, edeg, vedges = wm.build_tables(F)
    sym = CrystalSymmetry.for_manifold_path(path)
    fo = sym.orbit_id_map("face")

    if verbose:
        deg_census = {d: sum(1 for x in edeg.values() if x == d)
                      for d in sorted(set(edeg.values()))}
        print(f"reference: {path}")
        print(f"N3={len(F)} E={len(edeg)} V={len(vedges)}  "
              f"deg census: {deg_census}")
        print(f"|Aut| = {sym.order}; face orbits = {sym.n_orbits('face')}")

    orbits = defaultdict(list)
    n_invalid = 0
    for face, d, e, valid in wm.two_three_sites(F, faces, edeg):
        if not valid:
            n_invalid += 1
            continue
        orbits[fo[tuple(sorted(face))]].append((face, d, e))
    if verbose:
        print(f"\n2-3 sites: {sum(len(v) for v in orbits.values())} valid "
              f"({n_invalid} invalid: apex edge already present), "
              f"{len(orbits)} EXACT orbits")

    rows = []
    for oid, sites in sorted(orbits.items(), key=lambda kv: -len(kv[1])):
        face, d, e = sites[0]
        deltas = wm.delta_two_three(face, d, e, edeg)
        trans = wm.vertex_transitions(deltas, edeg, vedges)
        dS = wm.delta_S_shape(deltas, trans, ESTAR)
        F2, edeg2, vedges2 = verify_exact(F, faces, edeg, vedges,
                                          face, d, e, deltas, trans)
        fu = followups(F2, edeg2, vedges2)
        abc = sorted(face)
        rows.append(dict(
            orbit=int(oid), count=len(sites),
            stabilizer=sym.order // len(sites),
            signature=repr(canon_sig(face, d, e, edeg, vedges)),
            face_degs=tuple(sorted(edeg[frozenset(p)] for p in (
                (abc[0], abc[1]), (abc[1], abc[2]), (abc[2], abc[0])))),
            ill=ill_signature(deltas),
            dS=round(dS, 4),
            transitions=describe_transitions(trans),
            followups_32=sum(1 for k, v, _, _ in fu if k == "3-2" and v),
            followups_44=[(x[2], x[3]) for x in fu if x[0] == "4-4" and x[1]],
            fu_key=repr(followup_key(fu)),
        ))
    return sym, rows


def coarsening_report(rows):
    """How the old local-signature bucketing related to the exact orbits."""
    by_sig = defaultdict(list)
    for r in rows:
        by_sig[r["signature"]].append(r)
    fused = {s: rs for s, rs in by_sig.items() if len(rs) > 1}
    dS_split = sum(1 for rs in fused.values()
                   if len({r["dS"] for r in rs}) > 1)
    fu_split = sum(1 for rs in fused.values()
                   if len({r["fu_key"] for r in rs}) > 1)
    return dict(n_signature_classes=len(by_sig), n_fused=len(fused),
                fused_with_differing_dS=dS_split,
                fused_with_differing_followups=fu_split)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None,
                    help="reference .mfd (default: smallest R supercell)")
    ap.add_argument("--full", action="store_true",
                    help="detailed block per orbit instead of a summary table")
    ap.add_argument("--json", default=None)
    ap.add_argument("--cross-check", default=None,
                    help="second reference; assert the alphabet agrees")
    args = ap.parse_args()

    refs = all_refs("r")
    path = args.ref or (refs[0] if refs else best_refs(REF_GLOB)["r"])
    sym, rows = build_catalog(path)

    co = coarsening_report(rows)
    print(f"\nlocal-signature coarsening: the {len(rows)} orbits fall into "
          f"{co['n_signature_classes']} signature classes; "
          f"{co['n_fused']} of those fuse two or more orbits.")
    print(f"  fused classes with differing dS_shape : "
          f"{co['fused_with_differing_dS']}  (0 expected: dS is a function "
          f"of the local signature)")
    print(f"  fused classes with differing FOLLOW-UPS: "
          f"{co['fused_with_differing_followups']}  <-- these are what "
          f"signature bucketing got wrong")

    if args.full:
        for r in rows:
            print(f"\norbit {r['orbit']:>4} x{r['count']:<5} |stab|="
                  f"{r['stabilizer']}  face degs {r['face_degs']}  "
                  f"-> ill {r['ill']}  dS_shape={r['dS']:+.3f}")
            print(f"       transitions: {', '.join(r['transitions'])}")
            print(f"       follow-ups: {r['followups_32']} valid 3-2 (undo); "
                  f"{len(r['followups_44'])} valid 4-4 -> "
                  + (", ".join(f"ill{s} dS{ds:+.3f}"
                               for s, ds in r['followups_44'][:4]) or "none"))
            print(f"       [verified exactly]")
    else:
        print(f"\n  {'orbit':>5} {'n':>6} {'|stab|':>7} {'face degs':<12} "
              f"{'ill':<12} {'dS':>8} {'32':>3} {'44':>3}  transitions")
        for r in rows:
            print(f"  {r['orbit']:>5} {r['count']:>6} {r['stabilizer']:>7} "
                  f"{str(r['face_degs']):<12} {str(r['ill']):<12} "
                  f"{r['dS']:>+8.3f} {r['followups_32']:>3} "
                  f"{len(r['followups_44']):>3}  "
                  f"{', '.join(r['transitions'])[:60]}")
        print(f"  ({len(rows)} orbits; --full for per-orbit detail)")

    if args.cross_check:
        print(f"\ncross-check against {args.cross_check} ...")
        _, rows2 = build_catalog(args.cross_check, verbose=False)
        key = lambda rs: sorted((r["face_degs"], r["ill"], r["dS"],
                                 tuple(r["transitions"]), r["fu_key"])
                                for r in rs)
        a, b = key(rows), key(rows2)
        if a == b:
            print(f"  AGREES: {len(rows)} orbits, identical alphabet "
                  f"(supercell-independent, as expected)")
        else:
            only_a = [x for x in a if x not in b]
            only_b = [x for x in b if x not in a]
            print(f"  MISMATCH: {len(only_a)} rows only in {os.path.basename(path)}, "
                  f"{len(only_b)} only in {os.path.basename(args.cross_check)}")
            for x in (only_a + only_b)[:5]:
                print(f"    {x}")
            return 1

    if args.json:
        with open(args.json, "w") as f:
            json.dump({"crystal": os.path.basename(path),
                       "aut_order": sym.order,
                       "classification": "aut-face-orbit",
                       "n_orbits": len(rows),
                       "coarsening": co,
                       "classes": [{**r, "face_degs": list(r["face_degs"]),
                                    "ill": list(r["ill"])} for r in rows]},
                      f, indent=1)
        print(f"\nwrote {os.path.abspath(args.json)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
