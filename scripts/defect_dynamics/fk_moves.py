#!/usr/bin/env python3
"""Core machinery for the exhaustive FK->FK block-move search.

A move-piece is a triangulated 3-ball (a set of tets). We classify:
  * boundary faces  = triangles in exactly ONE tet; boundary verts/edges lie on them.
  * interior verts  = not on any boundary face (coordination complete -> FK-checkable).
  * interior edges  = not on any boundary face (degree complete -> must be in {5,6}).
  * FRAMING: for each boundary edge, its degree WITHIN the piece (interior tets on it).

Two pieces are a MOVE iff they have ISOMORPHIC FRAMED BOUNDARIES (boundary
triangulation + per-edge framing, up to boundary automorphism) but are
NON-ISOMORPHIC as balls. Same framing => regluing to any exterior preserves every
boundary degree => FK-legality preserved with no illegal intermediate.

Detection: exact canonical form of (a) the framed boundary and (b) the whole
ball, via colored-graph individualization-refinement. Group by framed-boundary
canon; a group with >1 distinct ball canon = a move class.

Self-tests (no crystals): iso-invariance; the octahedron's 3 diagonal fillings
(a genuine framed-boundary match -> MOVE); the 2-3 Pachner pair (framing differs
-> NOT a match).
"""
from collections import Counter, defaultdict
from itertools import combinations


# ---------------------------------------------------------------- ball analysis
def analyze(tets):
    tets = [tuple(sorted(t)) for t in tets]
    face = Counter()
    for t in tets:
        for f in combinations(t, 3):
            face[f] += 1
    bfaces = {f for f, c in face.items() if c == 1}
    bverts = {v for f in bfaces for v in f}
    edeg = Counter()
    for t in tets:
        for e in combinations(t, 2):
            edeg[e] += 1
    bedges = {e for f in bfaces for e in combinations(f, 2)}
    allv = {v for t in tets for v in t}
    iverts = allv - bverts
    iedges = set(edeg) - bedges
    return dict(tets=tets, bfaces=bfaces, bverts=bverts, iverts=iverts,
                edeg=edeg, bedges=bedges, iedges=iedges, allv=allv)


def fk_interior_ok(a):
    """All interior edges deg {5,6}; all interior vertices Frank-Kasper."""
    for e in a["iedges"]:
        if a["edeg"][e] not in (5, 6):
            return False
    for v in a["iverts"]:
        d = [a["edeg"][e] for e in a["edeg"] if v in e]
        if any(x not in (5, 6) for x in d):
            return False
        if d.count(6) not in (0, 2, 3, 4):
            return False
    return True


# --------------------------------------------------- generic canonical labeling
def _refine(adj, icolor):
    """Stable color refinement; icolor maps node->int, returns node->int."""
    col = dict(icolor)
    while True:
        sig = {v: (col[v], tuple(sorted(col[u] for u in adj[v]))) for v in col}
        rank = {s: i for i, s in enumerate(sorted(set(sig.values())))}
        new = {v: rank[sig[v]] for v in col}
        if len(set(new.values())) == len(set(col.values())):
            return new
        col = new


def _signature(adj, orig, order):
    idx = {v: i for i, v in enumerate(order)}
    edges = sorted({tuple(sorted((idx[u], idx[v])))
                    for u in adj for v in adj[u]})
    return (tuple(orig[v] for v in order), tuple(edges))


def _canon(adj, orig, icolor):
    col = _refine(adj, icolor)
    classes = defaultdict(list)
    for v, c in col.items():
        classes[c].append(v)
    nonsing = [c for c, vs in classes.items() if len(vs) > 1]
    if not nonsing:
        return _signature(adj, orig, sorted(col, key=lambda v: col[v]))
    target = min(nonsing)                      # smallest color class with a tie
    m = max(col.values())
    best = None
    for v in classes[target]:                  # individualize each, take lex-min
        c2 = dict(col)
        c2[v] = m + 1
        cand = _canon(adj, orig, c2)
        if best is None or cand < best:
            best = cand
    return best


def canonical(adj, color):
    """Canonical hashable form of a colored graph, invariant under
    color-preserving isomorphism (individualization-refinement, lex-min)."""
    o2i = {c: i for i, c in enumerate(sorted(set(color.values())))}
    base = {v: o2i[color[v]] for v in color}
    return _canon(adj, color, base)


# --------------------------------------------------- framed-boundary / ball canon
def _vert_invariant(a, v):
    """Local FK-ish invariant of a boundary vertex (its in-piece decoration)."""
    inc = sorted(a["edeg"][e] for e in a["edeg"] if v in e)
    ntet = sum(1 for t in a["tets"] if v in t)
    nbf = sum(1 for f in a["bfaces"] if v in f)
    return ("bv", tuple(inc), ntet, nbf)


def framed_boundary_canon(a):
    """Canonical form of the framed boundary: boundary verts + boundary faces +
    boundary-edge framing, as a colored incidence graph."""
    adj = defaultdict(set)
    color = {}
    for v in a["bverts"]:
        n = ("V", v)
        color[n] = _vert_invariant(a, v)
    for i, f in enumerate(sorted(a["bfaces"])):
        fn = ("F", i)
        color[fn] = ("face",)
        for v in f:
            adj[fn].add(("V", v)); adj[("V", v)].add(fn)
    for j, e in enumerate(sorted(a["bedges"])):
        en = ("E", j)
        color[en] = ("edge", a["edeg"][e])          # framing on the edge
        for v in e:
            adj[en].add(("V", v)); adj[("V", v)].add(en)
    for n in color:
        adj.setdefault(n, set())
    return canonical(adj, color)


def ball_canon(a):
    """Canonical form of the whole ball (vertex-tet incidence), boundary verts
    colored by their framed decoration so it refines consistently."""
    adj = defaultdict(set)
    color = {}
    for v in a["allv"]:
        n = ("V", v)
        color[n] = _vert_invariant(a, v) if v in a["bverts"] else ("iv",)
    for i, t in enumerate(a["tets"]):
        tn = ("T", i)
        color[tn] = ("tet",)
        for v in t:
            adj[tn].add(("V", v)); adj[("V", v)].add(tn)
    for n in color:
        adj.setdefault(n, set())
    return canonical(adj, color)


def filling_canon(a):
    """Canonical form of a filling REL FIXED BOUNDARY: boundary vertices
    individualized by identity (the exterior pins them), interior refinable.
    Two fillings of the SAME (same-labelled) boundary are the same move-state iff
    equal. This is the notion a physical move uses -- unlike ball_canon, it does
    NOT quotient by boundary symmetry (AB- vs CE-diagonal are distinct here)."""
    adj = defaultdict(set)
    color = {}
    for v in a["allv"]:
        color[("V", v)] = ("bv", v) if v in a["bverts"] else ("iv",)
    for i, t in enumerate(a["tets"]):
        tn = ("T", i)
        color[tn] = ("tet",)
        for v in t:
            adj[tn].add(("V", v)); adj[("V", v)].add(tn)
    for n in color:
        adj.setdefault(n, set())
    return canonical(adj, color)


def find_moves(regions, require_fk=True):
    """Detect moves among fillings that SHARE a boundary labelling (e.g. the
    output of a filling enumerator on one fixed boundary, or same-hole samples).
    Groups by framed-boundary type, then by fixed-boundary filling_canon; a group
    with >1 distinct filling = a move class. NOTE: cross-location mining needs a
    boundary-alignment step (canonical relabelling) before this applies."""
    buckets = defaultdict(list)
    for i, tets in enumerate(regions):
        a = analyze(tets)
        if require_fk and not fk_interior_ok(a):
            continue
        buckets[framed_boundary_canon(a)].append((i, filling_canon(a)))
    moves = []
    for fb, items in buckets.items():
        distinct = defaultdict(list)
        for i, fc in items:
            distinct[fc].append(i)
        if len(distinct) > 1:
            moves.append([grp[0] for grp in distinct.values()])
    return moves


# ------------------------------------------------------------------- self-tests
def _octahedron_fillings():
    # apexes 0,1; equator cycle 2-3-4-5. Fill along each internal diagonal.
    eq = [2, 3, 4, 5]
    ring = [(eq[i], eq[(i + 1) % 4]) for i in range(4)]
    diagAB = [(0, 1, a, b) for a, b in ring]
    diagCE = [(2, 4, x, y) for x, y in [(0, 3), (3, 1), (1, 5), (5, 0)]]
    diagDF = [(3, 5, x, y) for x, y in [(0, 2), (2, 1), (1, 4), (4, 0)]]
    return diagAB, diagCE, diagDF


def _selftest():
    ok = True
    # (1) iso-invariance: a ball and a relabeled copy are the SAME piece
    A = [(0, 1, 2, 3), (0, 1, 2, 4)]
    relabel = {0: 10, 1: 11, 2: 12, 3: 13, 4: 14}
    B = [tuple(relabel[v] for v in t) for t in A]
    same = ball_canon(analyze(A)) == ball_canon(analyze(B))
    print(f"[1] iso-invariance (relabel -> same canon): {same}")
    ok &= same

    # (2) octahedron: 3 diagonal fillings of ONE fixed framed boundary.
    # Abstractly ISOMORPHIC balls (ball_canon collapses to 1) but 3 distinct
    # fixed-boundary fillings (filling_canon) -> a move class of 3.
    fills = _octahedron_fillings()
    one_boundary = len({framed_boundary_canon(analyze(f)) for f in fills}) == 1
    iso_balls = len({ball_canon(analyze(f)) for f in fills}) == 1
    distinct_fills = len({filling_canon(analyze(f)) for f in fills}) == 3
    mv = find_moves(fills, require_fk=False)
    found = bool(mv) and len(mv[0]) == 3
    print(f"[2] octahedron: one framed boundary={one_boundary}, "
          f"abstractly-iso balls={iso_balls}, distinct fixed-boundary fillings="
          f"{distinct_fills}, move-class of 3 found={found}")
    ok &= one_boundary and iso_balls and distinct_fills and found

    # (3) 2-3 Pachner: framing differs -> NOT the same framed boundary
    two = [(0, 1, 2, 3), (0, 1, 2, 4)]
    three = [(0, 1, 3, 4), (1, 2, 3, 4), (0, 2, 3, 4)]
    fb2, fb3 = framed_boundary_canon(analyze(two)), framed_boundary_canon(analyze(three))
    differ = fb2 != fb3
    e01 = analyze(two)["edeg"][(0, 1)], analyze(three)["edeg"][(0, 1)]
    print(f"[3] 2-3 Pachner: framed boundaries differ={differ}  "
          f"(edge 01 framing {e01[0]} vs {e01[1]})")
    ok &= differ

    print(f"\nALL SELF-TESTS {'PASSED' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    _selftest()
