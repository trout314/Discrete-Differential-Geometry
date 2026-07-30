#!/usr/bin/env python3
"""Composability test for the deg-3-catalyzed deg-4 step: WALK a real deg-4
fragment across the crystal by iterating the 2-move catalyzed reaction.

Per step, on the LIVE manifold:
  * anchor = the current head deg-4 edge; enumerate all (2->3|3->2)^2
    composites in its radius-2 core (first move touching the head, second
    the disturbed region);
  * keep candidates that (a) heal/remove the head, (b) preserve the local
    illegal-degree composition exactly (content-neutral), (c) create >= 1
    new deg-4 OUTSIDE the head's hinge-flip octahedron, (d) have dS <= --max-ds;
  * choose the candidate whose escaped landing maximises displacement along
    the current direction (greedy directional persistence), apply it through
    do_bistellar_move, and advance the head to the landing edge.

After the walk, every step is UNDONE in reverse and the facet set must equal
the original exactly (invertibility check).  Positions for displacement
reporting are the snapshot's harmonic chart (static; vertices persist since
no 1->4/4->1 moves are involved).
"""
import argparse
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
import tip_retract_search as T
import defect_viewer as dv

CELL = 1e6


def facs(m):
    return [tuple(sorted(int(x) for x in t)) for t in np.asarray(m.facets())]


def mid(X, P, e, ref):
    d0 = X[e[0]] - ref
    d0 -= np.round(d0 / P) * P
    d1 = X[e[1]] - ref
    d1 -= np.round(d1 / P) * P
    return ref + (d0 + d1) / 2


def as_bistellar(mv):
    if mv[0] == "23":
        _, tri, x, y = mv
        return list(tri), [x, y]
    _, (u, w), a, b, c = mv
    return [u, w], [a, b, c]


def step_candidates(m, head, max_ds):
    """Enumerate content-neutral escaped-landing composites at `head`."""
    F = facs(m)
    Vball = T.ball(F, head, T.R_BALL)
    Vcore = frozenset(T.ball(F, head, T.R_CORE))
    region = [t for t in F if set(t) <= Vball]
    P = T.Patch([frozenset(t) for t in region])
    init_edeg = dict(P.edeg)
    init_comp = Counter(d for d in init_edeg.values() if d not in (5, 6))
    octa = set(head)
    for t in P.tets:
        if set(head) <= t:
            octa |= t

    def disturbed():
        f = set(head)
        for ed in P.touched:
            if P.edeg.get(ed, 0) != init_edeg.get(ed, 0):
                f |= set(ed)
        return f

    def apply_mv(mv):
        if mv[0] == "23":
            _, tri, x, y = mv
            P.apply_23(tri, x, y)
        else:
            _, (u, w), a, b, c = mv
            P.apply_32(u, w, a, b, c)

    def undo_mv(mv):
        if mv[0] == "23":
            _, tri, x, y = mv
            P.undo_23(tri, x, y)
        else:
            _, (u, w), a, b, c = mv
            P.undo_32(u, w, a, b, c)

    out = []
    for m1 in P.moves(Vcore, set(head)):
        apply_mv(m1)
        for m2 in P.moves(Vcore, disturbed()):
            apply_mv(m2)
            hd = P.edeg.get(head, 0)
            if hd in (5, 6) or hd == 0:                    # head healed/gone
                comp = Counter(d for d in P.edeg.values() if d not in (5, 6))
                if comp == init_comp:                      # content-neutral
                    new4 = [ed for ed in P.touched
                            if P.edeg.get(ed, 0) == 4
                            and init_edeg.get(ed, 0) != 4
                            and not set(ed) <= octa]
                    if new4:
                        ds = T.LAM * T.dS_full(init_edeg, P.edeg,
                                               frozenset(P.touched))
                        if ds <= max_ds:
                            out.append((m1, m2, ds, new4))
            undo_mv(m2)
        undo_mv(m1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snap", default=os.path.join(
        _ROOT, "data/mgas/lam35r_snap15000.mfd"))
    ap.add_argument("--head", type=int, nargs=2, default=[3282, 3318],
                    help="starting deg-4 edge (default: the verified tip)")
    ap.add_argument("--steps", type=int, default=8)
    ap.add_argument("--max-ds", type=float, default=1.0)
    args = ap.parse_args()

    m = ddg.Manifold.load(args.snap, 3)
    fac0 = frozenset(facs(m))
    X, P = dv.positions(args.snap, np.asarray(m.facets()))
    head = tuple(sorted(args.head))
    ref = X[head[0]].astype(float)
    pos0 = mid(X, P, head, ref)
    direction = None
    applied = []                    # (center, cocenter) in order, for undo
    total = np.zeros(3)
    print(f"walk from head {head}  (snapshot "
          f"{os.path.basename(args.snap)}; dS cap {args.max_ds})\n")

    for k in range(args.steps):
        cands = step_candidates(m, head, args.max_ds)
        if not cands:
            print(f"step {k + 1}: NO candidate (head {head}) -- walk ends")
            break
        best = None
        for m1, m2, ds, new4 in cands:
            for ed in new4:
                d = mid(X, P, ed, ref) - mid(X, P, head, ref)
                score = (float(np.dot(d, direction))
                         if direction is not None else float(np.linalg.norm(d)))
                if best is None or score > best[0]:
                    best = (score, m1, m2, ds, ed, d)
        score, m1, m2, ds, land, d = best
        for mv in (m1, m2):
            cen, coc = as_bistellar(mv)
            assert m.has_bistellar_move(cen, coc), f"invalid at step {k+1}"
            m.do_bistellar_move(cen, coc)
            applied.append((cen, coc))
        dirn = d / max(np.linalg.norm(d), 1e-12)
        direction = dirn if direction is None else \
            0.5 * (direction + dirn) / max(np.linalg.norm(
                0.5 * (direction + dirn)), 1e-12)
        total += d
        print(f"step {k + 1}: {head} -> {land}  |d|={np.linalg.norm(d)/CELL:.3f} "
              f"cells  dS={ds:+.3f}  ({len(cands)} candidates, "
              f"align={score/max(np.linalg.norm(d),1e-12):+.2f})")
        head = land

    n = len(applied) // 2
    print(f"\nwalked {n} steps, net displacement "
          f"{np.linalg.norm(total)/CELL:.3f} cells "
          f"(path length ~{n * 0.25:.2f}); head now {head}")
    pairs, degs = m.illegal_edges()
    print(f"final illegal composition: "
          f"{dict(Counter(int(x) for x in degs))}")

    # invertibility: undo everything, demand exact restoration
    for cen, coc in reversed(applied):
        m.do_bistellar_move(coc, cen)
    ok = frozenset(facs(m)) == fac0
    print(f"undo all {len(applied)} moves: facet set restored EXACTLY: {ok}")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
