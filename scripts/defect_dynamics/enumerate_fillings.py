#!/usr/bin/env python3
"""Enumerate FK-legal-alphabet fillings of ONE fixed framed boundary, to find the
minimal FK->FK block move source-independently (the definitive route: the
dynamics never presents the minimal defect cleanly, but enumeration constructs it).

Take a closed-star region of an FK crystal as the seed filling; its boundary
2-sphere + framing is FIXED. Explore the space of fillings by INTERIOR bistellar
flips that keep the boundary fixed (2->3, 3->2, 4-4; vertex-count-preserving),
BFS-bounded. Collect every filling whose interior satisfies the alphabet
(fk/deg3/deg4/e34). >=2 distinct fillings (by fixed-boundary canon) of the same
boundary = a move; classify by species (vac<->defect = creation, etc.).

Usage: enumerate_fillings.py STATE.mfd seed_v1 seed_v2 alphabet [grow] [cap]
  seed_v1 seed_v2 : an FK EDGE; its closed star is the seed (2 interior verts, so
                    the edge between them can host a fully-interior defect).
"""
import os
import sys
from collections import defaultdict, deque
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "python"))
sys.path.insert(0, _HERE)
import fk_moves as fk
from discrete_differential_geometry import Manifold


def has_edge(tets, u, w):
    return any(u in t and w in t for t in tets)


def _cycle(pairs):
    """order 4 link edges (pairs) into a 4-cycle [c0,c1,c2,c3]."""
    adj = defaultdict(list)
    for a, b in pairs:
        adj[a].append(b); adj[b].append(a)
    start = next(iter(adj))
    order = [start]
    prev = None
    while len(order) < 4:
        nxt = [x for x in adj[order[-1]] if x != prev]
        if not nxt:
            return None
        prev = order[-1]
        order.append(nxt[0] if nxt[0] not in order else nxt[-1])
    return order


def neighbors(tets, bfaces):
    """(rem, add) for each boundary-preserving interior bistellar move."""
    fm = defaultdict(list)
    et = defaultdict(list)
    for t in tets:
        for f in combinations(t, 3):
            fm[f].append(t)
        for e in combinations(t, 2):
            et[e].append(t)
    out = []
    # 2->3 : interior triangle shared by 2 tets
    for f, ts in fm.items():
        if len(ts) == 2 and f not in bfaces:
            d = next(iter(set(ts[0]) - set(f)))
            e = next(iter(set(ts[1]) - set(f)))
            if d != e and not has_edge(tets, d, e):
                out.append((set(ts),
                            {tuple(sorted((*p, d, e))) for p in combinations(f, 2)}))
    for e, ts in et.items():
        u, w = e
        on_bdry = any(u in bf and w in bf for bf in bfaces)
        if on_bdry:
            continue
        if len(ts) == 3:                              # 3->2 on interior deg-3 edge
            opp = set().union(*(set(t) - {u, w} for t in ts))
            if len(opp) == 3:
                out.append((set(ts), {tuple(sorted((*opp, u))),
                                      tuple(sorted((*opp, w)))}))
        elif len(ts) == 4:                            # 4-4 on interior deg-4 edge
            cyc = _cycle([tuple(set(t) - {u, w}) for t in ts])
            if cyc:
                c0, c1, c2, c3 = cyc
                if not has_edge(tets, c0, c2):
                    out.append((set(ts),
                                {tuple(sorted((c0, c2, u, c1))),
                                 tuple(sorted((c0, c2, c1, w))),
                                 tuple(sorted((c0, c2, w, c3))),
                                 tuple(sorted((c0, c2, c3, u)))}))
    return out


def _n_illegal(a):
    return sum(1 for e in a["iedges"] if a["edeg"][e] not in (5, 6))


def enumerate_fillings(seed, bfaces, alphabet, max_tets, max_ill=5, cap=40000):
    """Best-first over interior-flip fillings, prioritizing FEW illegal edges (so
    the rare near-vacuum fillings are reached first). Labels are fixed, so two
    fillings are identical iff same tet set -> dedup by frozenset, no canon."""
    import heapq
    seed = frozenset(seed)
    seen = {seed}
    heap = [(0, 0, seed)]                          # (n_illegal, tiebreak, tets)
    found = {}
    n = tie = 0
    while heap and n < cap:
        _, _, cur = heapq.heappop(heap); n += 1
        a = fk.analyze(list(cur))
        ok, sp = fk.interior_ok(a, alphabet)
        if ok:
            found.setdefault(fk.ball_canon(a), (sp, list(cur)))
        for rem, add in neighbors(cur, bfaces):
            nxt = frozenset((cur - rem) | add)
            if len(nxt) > max_tets or nxt in seen:
                continue
            an = fk.analyze(list(nxt))
            ni = _n_illegal(an)
            if ni > max_ill or nxt != frozenset(an["tets"]):
                continue
            seen.add(nxt); tie += 1
            heapq.heappush(heap, (ni, tie, nxt))
    return found, n, len(seen)


if __name__ == "__main__":
    path = sys.argv[1]
    v1, v2 = int(sys.argv[2]), int(sys.argv[3])
    alphabet = sys.argv[4] if len(sys.argv) > 4 else "e34"
    grow = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    cap = int(sys.argv[6]) if len(sys.argv) > 6 else 40000

    fac = [tuple(sorted(int(x) for x in t)) for t in Manifold.load(path, 3).facets()]
    seed = [t for t in fac if v1 in t or v2 in t]     # closed star of edge (v1,v2)
    a0 = fk.analyze(seed)
    bfaces = a0["bfaces"]
    print(f"seed edge ({v1},{v2}): {len(seed)} tets, boundary {len(bfaces)} faces, "
          f"interior {sorted(a0['iverts'])}, seed species "
          f"{fk.interior_ok(a0, alphabet)[1]}")
    found, nexpl, nseen = enumerate_fillings(seed, bfaces, alphabet,
                                             max_tets=len(seed) + grow, cap=cap)
    print(f"explored {nexpl} fillings ({nseen} distinct), "
          f"alphabet-satisfying: {len(found)}")
    species = sorted(sp for sp, _ in found.values())
    print(f"  species present: {species}")
    if len(found) >= 2:
        print(f"  >>> MOVE: {len(found)} distinct fillings of this framed boundary")
        for bc, (sp, tets) in found.items():
            print(f"        filling species={sp}, {len(tets)} tets")
    else:
        print("  no move (only one alphabet filling reachable)")
