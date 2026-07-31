"""Dressed generators: deterministic constructors for deep composite moves.

A dressed generator is a single-region composite move sequence achieving a
declared net effect -- the building blocks of the bilocal move program and
the unconditional f-vector generators (see notes/bilocal-worm-design.md).
This module provides the GENERAL machinery:

  * Goal        -- what the composite must achieve, expressed as a
                   progress metric + termination predicate;
  * staged_search -- the deterministic staged depth-limited DFS engine
                   (with cage moves and exact rollback) that cracked both
                   the edge-removal and vertex-collapse problems;
  * goals       -- EdgeRemoval, VertexCollapse (to Z=4, ready for 4->1),
                   VertexLegalize (post-1->4 growth: the insertion's
                   second phase);
  * net_effect / reverse -- the composite's f-vector change and canonical
                   net facet key; the exact inverse sequence.

DETERMINISM CONTRACT: given the same state and the same (goal, depth,
nodecap), the search returns the same composite. Every candidate list is
canonically sorted and every tie is broken by a canonical key -- NEVER by
dict/set iteration order, which is insertion-history-dependent in Python
and would silently break the reverse-check reversibility scheme (see the
bilocal design doc's bookkeeping section).

The 1<->4 boundary ops run through the sampler: the capi supports
targeted vertex-changing moves (4->1 returns the removed label to the
sampler's pool; 1->4 consumes its caller-chosen label), and the Live
index mirrors them, so composites price in place with no dup + fresh-
sampler repricing dance.
"""
from __future__ import annotations

import os
import sys
from collections import Counter
from itertools import combinations

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# region + alphabet (canonically ordered)
# ---------------------------------------------------------------------------

def star_verts(L, vv):
    out = set()
    for t in L.v2t.get(vv, ()):
        out.update(t)
    return out


def region_of(L, verts):
    """Union of closed stars of `verts`."""
    out = set()
    for vv in verts:
        out |= star_verts(L, vv)
    return out


def alphabet(L, region):
    """All valid elementary moves supported in `region`, canonically sorted.

    Returns [(center, cocenter)] with 3->2 entries for deg-3 edges inside
    the region and 2->3 entries for interior faces with non-adjacent
    apexes. Sorted so the DFS expansion order is a function of the state
    alone (determinism contract)."""
    out = []
    seen_e = set()
    for vv in sorted(region):
        for t in L.v2t.get(vv, ()):
            ts = tuple(sorted(t))
            for e in combinations(ts, 2):
                if e in seen_e or not (e[0] in region and e[1] in region):
                    continue
                seen_e.add(e)
                if L.edeg.get(e, 0) == 3:
                    tets = [q for q in L.v2t[e[0]] if e[1] in q]
                    if len(tets) != 3:
                        continue
                    lk = sorted({x for q in tets for x in q} - set(e))
                    if len(lk) == 3:
                        out.append((e, tuple(lk)))
    seen_f = set()
    for vv in sorted(region):
        for t in L.v2t.get(vv, ()):
            for f in combinations(tuple(sorted(t)), 3):
                if f in seen_f or not set(f) <= region:
                    continue
                seen_f.add(f)
                tets = [q for q in L.v2t[f[0]] if set(f) <= set(q)]
                if len(tets) != 2:
                    continue
                apex = sorted((set(tets[0]) | set(tets[1])) - set(f))
                if len(apex) != 2 or tuple(apex) in L.edeg:
                    continue
                out.append((f, tuple(apex)))
    out.sort()
    return out


# ---------------------------------------------------------------------------
# goals
# ---------------------------------------------------------------------------

class Goal:
    """Progress-metric interface for the staged search.

    progress(L): comparable; SMALLER is better; stages commit strictly
                 improving sequences.
    done(L):     terminal success.
    region(L):   vertex set the stage alphabet is confined to.
    """

    def progress(self, L):
        raise NotImplementedError

    def done(self, L):
        raise NotImplementedError

    def region(self, L):
        raise NotImplementedError


class EdgeRemoval(Goal):
    """Delete edge e: drive d(e) to 3, consume with the axis 3->2."""

    def __init__(self, e):
        self.e = tuple(sorted(e))

    def progress(self, L):
        return L.edeg.get(self.e, 0)        # 0 = gone = best

    def done(self, L):
        return self.e not in L.edeg

    def region(self, L):
        return region_of(L, self.e)


class VertexCollapse(Goal):
    """Drive Z(v) down to 4 (the 4->1 itself is a boundary op)."""

    def __init__(self, v):
        self.v = v

    def progress(self, L):
        return max(sum(1 for e in L.edeg if self.v in e), 4)

    def done(self, L):
        return sum(1 for e in L.edeg if self.v in e) == 4

    def region(self, L):
        return region_of(L, [self.v])


class VertexLegalize(Goal):
    """Grow a fresh (post-1->4) vertex's star until it carries no illegal
    edge -- the insertion generator's second phase. Progress = (number of
    illegal edges at v, then their total deficit)."""

    def __init__(self, v):
        self.v = v

    def _ill(self, L):
        return [(e, d) for e, d in L.edeg.items()
                if self.v in e and d not in (5, 6)]

    def progress(self, L):
        ill = self._ill(L)
        return (len(ill), sum(abs(6 - d) for _, d in ill))

    def done(self, L):
        return not self._ill(L)

    def region(self, L):
        return region_of(L, [self.v])


# ---------------------------------------------------------------------------
# the deterministic staged search engine
# ---------------------------------------------------------------------------

def staged_search(L, obj, goal, depth=4, nodecap=150000,
                  max_stages=64):
    """Run staged DFS until goal.done or no stage improves progress.

    Deterministic given (state, goal, depth, nodecap): canonical
    alphabet order, ties broken by (progress, round(dS, 9), move-seq key).

    Returns dict: ok, moves, ds_net, barrier, nodes, stages. On failure
    (ok=False) all moves are rolled back; on success they are LEFT
    APPLIED (callers roll back with `undo` if measuring only).
    """
    o0 = obj()
    all_moves = []
    smax = o0
    total_nodes = 0
    stages = 0

    while not goal.done(L) and stages < max_stages:
        p0 = goal.progress(L)
        best = {"key": None, "seq": None}
        nodes = [0]
        o_base = obj()

        def dfs(dep, seq):
            if nodes[0] > nodecap:
                return
            p = goal.progress(L)
            if p < p0:
                # canonical total order: progress, then exact dS, then the
                # move sequence itself -- a pure function of the state
                key = (p, round(obj() - o_base, 9), tuple(seq))
                if best["key"] is None or key < best["key"]:
                    best["key"] = key
                    best["seq"] = list(seq)
            if dep == 0 or goal.done(L):
                return
            for cen, coc in alphabet(L, goal.region(L)):
                nodes[0] += 1
                try:
                    L.do(cen, coc)
                except Exception:
                    continue
                seq.append((cen, coc))
                dfs(dep - 1, seq)
                seq.pop()
                L.do(coc, cen)

        dfs(depth, [])
        total_nodes += nodes[0]
        stages += 1
        if best["seq"] is None:
            break
        for cen, coc in best["seq"]:
            L.do(cen, coc)
            all_moves.append((cen, coc))
            smax = max(smax, obj())

    ok = goal.done(L)
    res = {"ok": ok, "moves": list(all_moves), "ds_net": obj() - o0,
           "barrier": smax - o0, "nodes": total_nodes, "stages": stages}
    if not ok:
        undo(L, all_moves)
    return res


def undo(L, moves):
    for cen, coc in reversed(moves):
        L.do(coc, cen)


def reverse(moves):
    """The exact inverse composite (apply to the end state to restore)."""
    return [(coc, cen) for cen, coc in reversed(moves)]


def net_effect(moves):
    """(d_f3, canonical net facet diff) of a composite of 2<->3 moves.

    d_f3 = +1 per 2->3, -1 per 3->2 (f0 changes only via the 1<->4
    boundary ops, which are outside this alphabet)."""
    net = Counter()
    df3 = 0
    for cen, coc in moves:
        if len(cen) == 3:                     # 2->3
            df3 += 1
            a, b, c = cen
            x, y = coc
            for t in ((a, b, c, x), (a, b, c, y)):
                net[tuple(sorted(t))] -= 1
            for e in ((a, b), (b, c), (a, c)):
                net[tuple(sorted(e + tuple(coc)))] += 1
        else:                                  # 3->2
            df3 -= 1
            u, w = cen
            a, b, c = coc
            for e in ((a, b), (b, c), (a, c)):
                net[tuple(sorted((u, w) + e))] -= 1
            for t in ((a, b, c, u), (a, b, c, w)):
                net[tuple(sorted(t))] += 1
    return df3, frozenset((t, k) for t, k in net.items() if k != 0)


# ---------------------------------------------------------------------------
# closed-form (Tier-1) constructors: octahedron-vertex <-> deg-4-edge
# ---------------------------------------------------------------------------
# The 4-move transmutation between a deg-4 edge (its cavity's diagonal
# triangulation) and a Z=6 vertex with octahedral link (the cavity coned).
# Frame-parameterized, no search; E->V is valid at EVERY deg-4 edge with no
# ambient preconditions; the matching-slot reverse is the exact inverse
# (validated 6/6 with machine-precision dS antisymmetry, transmute_lab.py).
# The 1<->4 boundary ops run THROUGH the sampler (the capi supports
# targeted 1<->4 with label-pool bookkeeping); pass a
# `sampler_factory(manifold) -> (sampler, Live)` to price the copy.

def cycle_of(L, e):
    """Link 4-cycle of a deg-4 edge, in cyclic order (canonical start)."""
    tets = [t for t in L.v2t[e[0]] if e[1] in t]
    assert len(tets) == 4, "not a deg-4 edge"
    pairs = [tuple(sorted(set(t) - set(e))) for t in tets]
    adj = {}
    for a, b in pairs:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    start = min(adj)
    cyc = [start, sorted(adj[start])[0]]
    while len(cyc) < 4:
        nxt = [x for x in adj[cyc[-1]] if x != cyc[-2]]
        cyc.append(nxt[0])
    return cyc


def link3_of(L, e):
    """The 3 link vertices of a deg-3 edge, read live."""
    tets = [t for t in L.v2t[e[0]] if e[1] in t]
    assert len(tets) == 3
    lk = sorted({x for t in tets for x in t} - set(e))
    assert len(lk) == 3
    return lk


def edge_to_vertex(m0, e, label, sampler_factory, slot=0):
    """E->V transmutation on a copy of m0. Returns (sampler, Live, v,
    frame) with v = the new octahedral-link vertex.

    slot in 0..3 picks which adjacent cycle pair survives to the final
    3->2 co-face (the Hastings slot count of the constructor)."""
    import worm_deg4_slide as _W
    L0 = _W.Live(m0)
    u, w = tuple(sorted(e))
    cyc = cycle_of(L0, (u, w))
    a, b = cyc[slot % 4], cyc[(slot + 1) % 4]
    c, d = cyc[(slot + 2) % 4], cyc[(slot + 3) % 4]
    m2 = m0.dup() if hasattr(m0, "dup") else m0.copy()
    s2, L2 = sampler_factory(m2)
    L2.do(sorted([u, w, c, d]), [label])      # 1->4 through the sampler
    L2.do(tuple(sorted((u, w, c))), (label, b))
    L2.do(tuple(sorted((u, w, d))), (label, a))
    L2.do(tuple(sorted((u, w))), (label, a, b))
    return s2, L2, label, {"e": (u, w), "cyc": cyc, "slot": slot}


def vertex_to_edge(view, v, diag_pair, keep_pair, sampler_factory):
    """V->E transmutation on a copy: diagonal = diag_pair (3 choices per
    octahedral vertex), final-face pair = keep_pair (adjacent link pair
    whose two link faces have the diagonal as apexes). Co-centers read
    live. Returns (sampler, Live, new_edge)."""
    m3 = view.dup()
    s3, L3 = sampler_factory(m3)
    u, w = diag_pair
    p2, p3 = keep_pair
    L3.do(tuple(sorted((v, p2, p3))), (u, w))            # link flip
    L3.do(tuple(sorted((v, p2))),
          tuple(link3_of(L3, tuple(sorted((v, p2))))))
    L3.do(tuple(sorted((v, p3))),
          tuple(link3_of(L3, tuple(sorted((v, p3))))))
    nb = sorted({x for t in L3.v2t[v] for x in t} - {v})
    assert len(nb) == 4
    L3.do([v], nb)                            # 4->1 through the sampler
    return s3, L3, (u, w)
