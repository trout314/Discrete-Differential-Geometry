#!/usr/bin/env python3
"""Mine an FK state for FK->FK block moves: extract closed-star regions up to K
interior vertices, bucket by framed boundary, and flag any framed-boundary type
that is filled by >=2 non-isomorphic FK-legal balls (a definite move).

Pruning cascade (cheap -> exact): every region gets a cheap boundary key and a
cheap ball key; only within a boundary-key bucket that already shows >=2 distinct
cheap ball keys do we compute the exact framed_boundary_canon / ball_canon. A
periodic crystal fills each hole identically -> zero moves (the no-false-positive
sanity). A disordered glass fills recurring boundaries several ways -> moves.

Usage: python glass_mine.py STATE.mfd K [alphabet]
  K        = max interior vertices (2 or 3)
  alphabet = 'fk' (strict, vacuum moves) or 'deg4' (FK + one deg-4 edge; the
             minimal disclination unit -> creation/transport moves). Default fk.
Moves are classified by filling species: vac<->vac (vacuum), vac<->deg4
(creation/annihilation), deg4<->deg4 (transport/reconfiguration)."""
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "python"))
sys.path.insert(0, _HERE)
import fk_moves as fk
from discrete_differential_geometry import Manifold


def cheap_boundary_key(a):
    fr = sorted(a["edeg"][e] for e in a["bedges"])
    return (len(a["bverts"]), len(a["bedges"]), len(a["bfaces"]), tuple(fr))


def cheap_ball_key(a):
    ie = sorted(a["edeg"][e] for e in a["iedges"])
    iv = sorted(sum(1 for e in a["edeg"] if v in e and a["edeg"][e] == 6)
                for v in a["iverts"])
    return (len(a["iverts"]), len(a["tets"]), tuple(ie), tuple(iv))


def regions(v2t, adjv, V, K):
    """Yield connected vertex-sets (as frozensets) of size 1..K."""
    seen = set()
    for v in range(V):
        if v in v2t:
            yield frozenset((v,))
    for v in range(V):
        for w in adjv[v]:
            if v < w:
                yield frozenset((v, w))
    if K >= 3:
        for v in range(V):
            nb = sorted(adjv[v])
            for i in range(len(nb)):
                for j in range(i + 1, len(nb)):
                    s = frozenset((nb[i], v, nb[j]))
                    if len(s) == 3 and s not in seen:
                        seen.add(s)
                        yield s


def mine(facets, K, alphabet="fk", verbose=False):
    tets = [tuple(sorted(int(x) for x in t)) for t in facets]
    V = max(max(t) for t in tets) + 1
    v2t = defaultdict(list)
    adjv = defaultdict(set)
    edeg = defaultdict(int)
    for i, t in enumerate(tets):
        for v in t:
            v2t[v].append(i)
        for a, b in combinations(t, 2):
            adjv[a].add(b); adjv[b].add(a)
            edeg[(a, b)] += 1
    # FK-legal vertices in the full state. A closed-star region's interior legality
    # is fixed by these; deg4 alphabet allows S to also hold the (<=2) endpoints of
    # one deg-4 edge. Sound pre-filter: #illegal interior verts must be 0 (or 2 for
    # deg4); exact fk.interior_ok(a, alphabet) confirms.
    legal = set()
    for v in v2t:
        degs = [edeg[tuple(sorted((v, w)))] for w in adjv[v]]
        if all(d in (5, 6) for d in degs) and degs.count(6) in (0, 2, 3, 4):
            legal.add(v)
    allowed_nill = {0} if alphabet == "fk" else {0, 2}

    buckets = defaultdict(list)     # cheap bkey -> [(cheap ball key, tets, species)]
    n_reg = n_fk = 0
    for S in regions(v2t, adjv, V, K):
        n_reg += 1
        if verbose and n_reg % 50000 == 0:
            print(f"    ...{n_reg} regions, {n_fk} accepted, {len(buckets)} btypes",
                  flush=True)
        if len(S - legal) not in allowed_nill:
            continue
        ti = set()
        for v in S:
            ti.update(v2t[v])
        rt = [tets[i] for i in ti]
        a = fk.analyze(rt)
        if not S <= a["iverts"]:
            continue                # accidental extra interior vertex; skip
        acc, species = fk.interior_ok(a, alphabet)
        if not acc:
            continue
        n_fk += 1
        buckets[cheap_boundary_key(a)].append((cheap_ball_key(a), rt, species))

    # exact verification only where cheap ball keys already disagree
    moves = []
    n_candidate_buckets = 0
    for bk, items in buckets.items():
        if len({ck for ck, _, _ in items}) < 2:
            continue
        n_candidate_buckets += 1
        by_fb = defaultdict(dict)   # framed boundary -> {ball_canon: (tets, species)}
        for ck, rt, sp in items:
            a = fk.analyze(rt)
            by_fb[fk.framed_boundary_canon(a)].setdefault(
                fk.ball_canon(a), (rt, sp))
        for fbc, balls in by_fb.items():
            if len(balls) >= 2:
                moves.append(list(balls.values()))   # [(tets, species), ...]

    return dict(n_regions=n_reg, n_fk_legal=n_fk, n_boundary_types=len(buckets),
                n_candidate_buckets=n_candidate_buckets, moves=moves)


if __name__ == "__main__":
    path, K = sys.argv[1], int(sys.argv[2])
    alphabet = sys.argv[3] if len(sys.argv) > 3 else "fk"
    fac = np.asarray(Manifold.load(path, 3).facets())
    r = mine(fac, K, alphabet=alphabet, verbose=True)
    print(f"{os.path.basename(path)}  K={K}  alphabet={alphabet}")
    print(f"  regions scanned      : {r['n_regions']}")
    print(f"  accepted closed stars: {r['n_fk_legal']}")
    print(f"  distinct boundary types: {r['n_boundary_types']}")
    print(f"  candidate buckets    : {r['n_candidate_buckets']}")
    print(f"  MOVES (>=2 fillings of one framed boundary): {len(r['moves'])}")
    kind = {("vac", "vac"): "vacuum", ("vac", "deg4"): "CREATION",
            ("deg4", "deg4"): "TRANSPORT"}
    for i, mv in enumerate(r["moves"][:15]):
        sp = tuple(sorted(s for _, s in mv))
        label = kind.get(sp, f"{sp}")
        sizes = [len(fk.analyze(t)["iverts"]) for t, _ in mv]
        tets = [len(t) for t, _ in mv]
        print(f"    move {i} [{label}]: {len(mv)} fillings, species {sp}, "
              f"interior sizes {sizes}, tets {tets}")
