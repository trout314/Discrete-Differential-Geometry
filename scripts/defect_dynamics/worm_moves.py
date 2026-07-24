#!/usr/bin/env python3
"""Exact combinatorics of the V-preserving bistellar moves (2-3, 3-2, 4-4) --
the worm-move alphabet -- as a REUSABLE LIBRARY over bare facet arrays.

Crystal-independent by design: every function takes a plain (N, 4) facet array
(any labels) plus the derived tables from `build_tables`.  Drivers decide what
to enumerate over (a reference crystal, a defected state, or an abstract
exhaustive search); this module owns validity, exact integer degree deltas,
per-vertex (n6, m) junction transitions, and the closed-form full-action dS.

Move arithmetic (derivations 2026-07-24 session):
  2-3 on face {a,b,c} with apexes d,e (valid iff edge de absent):
      creates edge de at degree 3; face edges ab, bc, ca each -1;
      side edges xd, xe (x in {a,b,c}) each +1.   dN3=+1, dE=+1, dV=0.
  3-2 on a degree-3 edge de with link {a,b,c} (valid iff face abc absent):
      exact inverse of the above.                 dN3=-1, dE=-1, dV=0.
  4-4 on a degree-4 edge ef with link cycle a-b-c-d, choosing diagonal ac
      (valid iff edge ac absent):
      removes ef (deg 4), creates ac (deg 4); equator edges ab, bc, cd, da
      each +1; pole edges eb, fb, ed, fd each -1. all counts unchanged.

Vertex counters follow the sampler convention (source/sampler.d):
  n6(v) = # incident edges with degree >= 6
  m(v)  = # incident edges with degree not in {5, 6}   (impurity valence)
Legal classes (m=0): n6=0 Z12, 2 Z14, 3 Z15, 4 Z16.
"""
from collections import defaultdict
from itertools import combinations

import numpy as np

# ---------------------------------------------------------------------------
# derived tables
# ---------------------------------------------------------------------------


def build_tables(F):
    """faces: frozenset(3)->[tet ids]; edeg: frozenset(2)->degree;
    vedges: vertex -> set of incident edges (frozensets)."""
    F = np.asarray(F)
    faces = defaultdict(list)
    edeg = defaultdict(int)
    for t, tet in enumerate(F):
        tet = [int(x) for x in tet]
        for tri in combinations(tet, 3):
            faces[frozenset(tri)].append(t)
        for e in combinations(tet, 2):
            edeg[frozenset(e)] += 1
    vedges = defaultdict(set)
    for e in edeg:
        for v in e:
            vedges[v].add(e)
    return dict(faces), dict(edeg), dict(vedges)


def vertex_counters(v, edeg, vedges):
    """(n6, m) of vertex v under the sampler convention."""
    degs = [edeg[e] for e in vedges[v]]
    return (sum(1 for d in degs if d >= 6),
            sum(1 for d in degs if d not in (5, 6)))


# ---------------------------------------------------------------------------
# site enumeration + validity
# ---------------------------------------------------------------------------


def two_three_sites(F, faces, edeg):
    """Yield (face frozenset, d, e, valid) for every interior face."""
    F = np.asarray(F)
    for face, tets in faces.items():
        if len(tets) != 2:
            continue                              # boundary (not our case)
        t0, t1 = tets
        d = int(next(x for x in F[t0] if int(x) not in face))
        e = int(next(x for x in F[t1] if int(x) not in face))
        yield face, d, e, frozenset((d, e)) not in edeg


def three_two_sites(F, faces, edeg):
    """Yield (edge frozenset, link_verts, valid) for every degree-3 edge."""
    F = np.asarray(F)
    tets_of = defaultdict(list)
    for t, tet in enumerate(F):
        for e in combinations([int(x) for x in tet], 2):
            tets_of[frozenset(e)].append(t)
    for e, d3 in edeg.items():
        if d3 != 3:
            continue
        link = set()
        for t in tets_of[e]:
            link |= {int(x) for x in F[t]} - set(e)
        if len(link) != 3:
            continue                              # pinched: invalid
        yield e, frozenset(link), frozenset(link) not in faces


def four_four_sites(F, faces, edeg):
    """Yield (edge ef, cycle [a,b,c,d], diagonal, valid) -- two diagonal
    choices per degree-4 edge with a proper 4-cycle link."""
    F = np.asarray(F)
    tets_of = defaultdict(list)
    for t, tet in enumerate(F):
        for e in combinations([int(x) for x in tet], 2):
            tets_of[frozenset(e)].append(t)
    for e, d4 in edeg.items():
        if d4 != 4:
            continue
        # build the link cycle from the 4 tets: neighbors in the cycle share a tet
        pairs = [frozenset({int(x) for x in F[t]} - set(e)) for t in tets_of[e]]
        if any(len(p) != 2 for p in pairs):
            continue
        adj = defaultdict(set)
        for p in pairs:
            a, b = sorted(p)
            adj[a].add(b)
            adj[b].add(a)
        verts = sorted(adj)
        if len(verts) != 4 or any(len(adj[v]) != 2 for v in verts):
            continue                              # not a simple 4-cycle
        a = verts[0]
        b = min(adj[a])
        cyc = [a, b]
        while len(cyc) < 4:
            nxt = next(x for x in adj[cyc[-1]] if x != cyc[-2])
            cyc.append(nxt)
        for diag in (frozenset((cyc[0], cyc[2])), frozenset((cyc[1], cyc[3]))):
            yield e, cyc, diag, diag not in edeg


# ---------------------------------------------------------------------------
# fast targeted enumeration (numpy; for worm heads / small illegal sets):
# sites at SPECIFIC candidate edges, without building full python tables
# ---------------------------------------------------------------------------


def edge_degrees_np(F):
    """(eu (E,2) int64 sorted pairs, edeg (E,) counts) -- vectorized."""
    F = np.asarray(F, np.int64)
    pairs = np.vstack([np.sort(F[:, [i, j]], axis=1)
                       for i in range(4) for j in range(i + 1, 4)])
    eu, cnt = np.unique(pairs, axis=0, return_counts=True)
    return eu, cnt


def _tets_with_pair(F, a, b):
    F = np.asarray(F)
    m = ((F == a).any(1)) & ((F == b).any(1))
    return np.nonzero(m)[0]


def _contains_simplex(F, verts):
    F = np.asarray(F)
    m = np.ones(len(F), bool)
    for v in verts:
        m &= (F == v).any(1)
    return bool(m.any())


def three_two_sites_at(F, edges3):
    """As three_two_sites but only at the given deg-3 edges (pairs)."""
    F = np.asarray(F)
    for a, b in edges3:
        ts = _tets_with_pair(F, a, b)
        if len(ts) != 3:
            continue
        link = sorted(set(F[ts].ravel().tolist()) - {int(a), int(b)})
        if len(link) != 3:
            continue
        yield frozenset((int(a), int(b))), frozenset(link), \
            not _contains_simplex(F, link)


def four_four_sites_at(F, edges4):
    """As four_four_sites but only at the given deg-4 edges (pairs)."""
    F = np.asarray(F)
    for a, b in edges4:
        ts = _tets_with_pair(F, a, b)
        if len(ts) != 4:
            continue
        pairs = [frozenset(set(F[t].tolist()) - {int(a), int(b)})
                 for t in ts]
        if any(len(p) != 2 for p in pairs):
            continue
        adj = defaultdict(set)
        for p in pairs:
            x, y = sorted(p)
            adj[x].add(y)
            adj[y].add(x)
        verts = sorted(adj)
        if len(verts) != 4 or any(len(adj[v]) != 2 for v in verts):
            continue
        v0 = verts[0]
        v1 = min(adj[v0])
        cyc = [v0, v1]
        while len(cyc) < 4:
            cyc.append(next(x for x in adj[cyc[-1]] if x != cyc[-2]))
        for diag in (frozenset((cyc[0], cyc[2])), frozenset((cyc[1], cyc[3]))):
            yield frozenset((int(a), int(b))), cyc, diag, \
                not _contains_simplex(F, sorted(diag))


# ---------------------------------------------------------------------------
# exact degree deltas  {edge: (old, new)}; created edges get old=None,
# removed edges new=None
# ---------------------------------------------------------------------------


def delta_two_three(face, d, e, edeg):
    abc = sorted(face)
    out = {}
    for pair in combinations(abc, 2):
        k = frozenset(pair)
        out[k] = (edeg[k], edeg[k] - 1)
    for x in abc:
        for apex in (d, e):
            k = frozenset((x, apex))
            out[k] = (edeg[k], edeg[k] + 1)
    out[frozenset((d, e))] = (None, 3)
    return out


def delta_three_two(edge, link, edeg):
    d, e = sorted(edge)
    abc = sorted(link)
    out = {}
    for pair in combinations(abc, 2):
        k = frozenset(pair)
        out[k] = (edeg[k], edeg[k] + 1)
    for x in abc:
        for apex in (d, e):
            k = frozenset((x, apex))
            out[k] = (edeg[k], edeg[k] - 1)
    out[frozenset(edge)] = (3, None)
    return out


def delta_four_four(edge, cyc, diag, edeg):
    a, c = sorted(diag)
    b, dd = sorted(set(cyc) - diag)
    out = {frozenset(edge): (4, None), frozenset(diag): (None, 4)}
    for i in range(4):
        k = frozenset((cyc[i], cyc[(i + 1) % 4]))
        out[k] = (edeg[k], edeg[k] + 1)
    for pole in sorted(edge):
        for x in (b, dd):
            k = frozenset((pole, x))
            out[k] = (edeg[k], edeg[k] - 1)
    return out


# ---------------------------------------------------------------------------
# per-vertex junction transitions
# ---------------------------------------------------------------------------

_CLASS = {(0, 0): "Z12", (2, 0): "Z14", (3, 0): "Z15", (4, 0): "Z16"}


def class_name(n6, m):
    if m > 0:
        return f"ill(n6={n6},m={m})"
    return _CLASS.get((n6, m), f"n6={n6}")


def vertex_transitions(deltas, edeg, vedges):
    """{vertex: ((n6, m) before, (n6, m) after)} for all support vertices."""
    verts = sorted({v for e in deltas for v in e})
    out = {}
    for v in verts:
        n6b, mb = vertex_counters(v, edeg, vedges)
        dn6 = dm = 0
        for e, (old, new) in deltas.items():
            if v not in e:
                continue
            if old is not None and new is not None:
                dn6 += (new >= 6) - (old >= 6)
                dm += (new not in (5, 6)) - (old not in (5, 6))
            elif old is None:                     # created edge
                dn6 += (new >= 6)
                dm += (new not in (5, 6))
            else:                                 # removed edge
                dn6 -= (old >= 6)
                dm -= (old not in (5, 6))
        out[v] = ((n6b, mb), (n6b + dn6, mb + dm))
    return out


# ---------------------------------------------------------------------------
# closed-form full-action dS (lambda factored out; see CONVENTIONS.md)
#   S_shape = (e*/6) [ sum_e (deg-e*)^2 - x(1-x) E ] + zleg sum_v U + cimp sum_v m^2
#   dS_vol  = c_N (2 (N - N*) + 1) dN3   (quote separately; 0 on-target)
# ---------------------------------------------------------------------------


def U_zleg(n6):
    if n6 == 0 or 2 <= n6 <= 4:
        return 0
    if n6 in (1, 5):
        return 1
    return (n6 - 4) ** 2


def delta_S_shape(deltas, trans, estar=5.105025, zleg=0.6, cimp=1.0):
    x = estar - int(estar)
    dS = 0.0
    dE = 0
    for e, (old, new) in deltas.items():
        if old is not None and new is not None:
            dS += (estar / 6.0) * ((new - estar) ** 2 - (old - estar) ** 2)
        elif old is None:
            dS += (estar / 6.0) * ((new - estar) ** 2 - x * (1 - x))
            dE += 1
        else:
            dS += (estar / 6.0) * (-((old - estar) ** 2) + x * (1 - x))
            dE -= 1
    for v, ((n6b, mb), (n6a, ma)) in trans.items():
        dS += zleg * (U_zleg(n6a) - U_zleg(n6b))
        dS += cimp * (ma ** 2 - mb ** 2)
    return dS


# ---------------------------------------------------------------------------
# apply moves to a facet array (for exact verification / search recursion)
# ---------------------------------------------------------------------------


def _remove_add(F, remove_tets, add_tets):
    F = np.asarray(F)
    keep = np.ones(len(F), bool)
    keep[list(remove_tets)] = False
    return np.vstack([F[keep], np.array(sorted(map(sorted, add_tets)))])


def apply_two_three(F, faces, face, d, e):
    a, b, c = sorted(face)
    return _remove_add(F, faces[face],
                      [(a, b, d, e), (b, c, d, e), (c, a, d, e)])


def apply_three_two(F, faces, edge, link):
    d, e = sorted(edge)
    a, b, c = sorted(link)
    tets = set()
    for pair in ((a, b), (b, c), (c, a)):
        tets |= set(faces[frozenset((pair[0], pair[1], d))]) & \
                set(faces[frozenset((pair[0], pair[1], e))])
    return _remove_add(F, tets, [(a, b, c, d), (a, b, c, e)])


def apply_four_four(F, faces, edge, cyc, diag):
    ef = sorted(edge)
    old = set()
    for i in range(4):
        old |= set(faces[frozenset((cyc[i], cyc[(i + 1) % 4], ef[0]))]) & \
               set(faces[frozenset((cyc[i], cyc[(i + 1) % 4], ef[1]))])
    a, c = sorted(diag)
    b, dd = sorted(set(cyc) - diag)
    new = [(a, c, b, ef[0]), (a, c, b, ef[1]), (a, c, dd, ef[0]),
           (a, c, dd, ef[1])]
    return _remove_add(F, old, new)
