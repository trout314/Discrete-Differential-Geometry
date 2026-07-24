#!/usr/bin/env python3
"""Worm-move catalog, stage 1: every 2-3 creation site of a reference crystal,
classified by exact local signature, with junction transitions, dS, exact
verification, and the available follow-up (3-2 / 4-4) moves per class.

Driver over the crystal-independent library `worm_moves.py` (which a future
exhaustive, non-crystal search will reuse; keep all move arithmetic THERE).

On a perfect FK crystal only 2-3 sites exist (no deg-3/4 edges), so stage 1 =
the pair-creation alphabet + the first walk steps available after each
creation.  Signature (canonical under S3 x S2 relabeling of face x apexes):
  per face vertex: (its two face-edge degrees sorted, its (to-d, to-e) degrees
  sorted, its (n6, m)); plus the apexes' (n6, m) sorted.
Merges space-group-equivalent sites; may also merge distinct orbits with
identical local data (noted -- orbit refinement comes with the registry pass).

Usage: worm_catalog.py [REF.mfd] [--json OUT]  (default: tcp_reference r)
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
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs

ESTAR = 5.105025


def canon_sig(face, d, e, edeg, vedges):
    """Canonical signature of a 2-3 site under S3 (face) x S2 (apexes)."""
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
    out = []
    for e, (old, new) in deltas.items():
        if new is not None and new not in (5, 6):
            out.append(new)
    return tuple(sorted(out))


def describe_transitions(trans):
    out = []
    for v, (before, after) in sorted(trans.items()):
        if before != after:
            out.append(f"{wm.class_name(*before)}->{wm.class_name(*after)}")
    out.sort()
    return out


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


def followups(F2, edeg2):
    """Available 3-2 / 4-4 moves in the post-creation state, with outcomes."""
    faces2, edeg2b, vedges2 = wm.build_tables(F2)
    out = []
    for edge, link, valid in wm.three_two_sites(F2, faces2, edeg2b):
        dl = wm.delta_three_two(edge, link, edeg2b)
        tr = wm.vertex_transitions(dl, edeg2b, vedges2)
        out.append(("3-2", valid, wm.ill_signature(dl) if hasattr(wm, "ill_signature")
                    else ill_signature(dl),
                    round(wm.delta_S_shape(dl, tr, ESTAR), 4)))
    for edge, cyc, diag, valid in wm.four_four_sites(F2, faces2, edeg2b):
        dl = wm.delta_four_four(edge, cyc, diag, edeg2b)
        tr = wm.vertex_transitions(dl, edeg2b, vedges2)
        out.append(("4-4", valid, ill_signature(dl),
                    round(wm.delta_S_shape(dl, tr, ESTAR), 4)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    path = args.ref or best_refs(REF_GLOB)["r"]
    F = np.asarray(ddg.Manifold.load(path, 3).facets())
    faces, edeg, vedges = wm.build_tables(F)
    print(f"reference: {path}")
    print(f"N3={len(F)} E={len(edeg)} V={len(vedges)}  "
          f"deg census: { {d: sum(1 for x in edeg.values() if x == d) for d in sorted(set(edeg.values()))} }")

    buckets = defaultdict(list)
    n_invalid = 0
    for face, d, e, valid in wm.two_three_sites(F, faces, edeg):
        if not valid:
            n_invalid += 1
            continue
        buckets[canon_sig(face, d, e, edeg, vedges)].append((face, d, e))
    print(f"\n2-3 sites: {sum(len(v) for v in buckets.values())} valid "
          f"({n_invalid} invalid: apex edge already present), "
          f"{len(buckets)} signature classes\n")

    rows = []
    for sig, sites in sorted(buckets.items(),
                             key=lambda kv: -len(kv[1])):
        face, d, e = sites[0]
        deltas = wm.delta_two_three(face, d, e, edeg)
        trans = wm.vertex_transitions(deltas, edeg, vedges)
        dS = wm.delta_S_shape(deltas, trans, ESTAR)
        F2, edeg2, _ = verify_exact(F, faces, edeg, vedges, face, d, e,
                                    deltas, trans)
        fu = followups(F2, edeg2)
        n32 = sum(1 for k, v, _, _ in fu if k == "3-2" and v)
        n44v = [x for x in fu if x[0] == "4-4" and x[1]]
        row = dict(
            count=len(sites),
            face_degs=tuple(sorted(edeg[frozenset(p)]
                                   for p in [(sorted(face)[0], sorted(face)[1]),
                                             (sorted(face)[1], sorted(face)[2]),
                                             (sorted(face)[2], sorted(face)[0])])),
            ill=ill_signature(deltas),
            dS=round(dS, 4),
            transitions=describe_transitions(trans),
            followups_32=n32,
            followups_44=[(x[2], x[3]) for x in n44v],
        )
        rows.append(row)
        print(f"x{row['count']:4d}  face degs {row['face_degs']}  "
              f"-> ill {row['ill']}  dS_shape={row['dS']:+.3f}")
        print(f"       transitions: {', '.join(row['transitions'])}")
        print(f"       follow-ups: {n32} valid 3-2 (undo); "
              f"{len(n44v)} valid 4-4 -> "
              + (", ".join(f"ill{s} dS{ds:+.3f}" for s, ds in
                           row['followups_44'][:4]) or "none"))
        print(f"       [verified exactly]")
    if args.json:
        with open(args.json, "w") as f:
            json.dump([{**r, "face_degs": list(r["face_degs"]),
                        "ill": list(r["ill"])} for r in rows], f, indent=1)
        print(f"\nwrote {os.path.abspath(args.json)}")


if __name__ == "__main__":
    main()
