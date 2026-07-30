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

The 1<->4 boundary ops are NOT executable through a sampler driver (the
capi forbids vertex-changing targeted moves); callers finish a vertex
removal/insertion by duplicating at the Manifold level and repricing with
a fresh sampler (see vertex_removal_v2.py for the pattern).
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
