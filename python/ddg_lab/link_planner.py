"""(link, collar) planner: 2D planning for deep vertex composites.

Moves at v act on lk(v) as 2D bistellar moves (bilocal design doc 3.-1).
This module PLANS a link reduction entirely in the decorated 2D state --
link faces + ambient degrees of spokes/link-edges/chords + m-counters --
and only then LIFTS the plan to 3D moves. Every link flip carries a
FLAVOR: '23' (2->3: needs the target chord absent; df3 = +1; expels a
tet) or '32' (3->2: needs the flipped edge at ambient degree exactly 3;
df3 = -1; absorbs the outside tet). Vertex deletion = 3->2 on a spoke at
link-degree 3.

The planner's dS bookkeeping (pins + c_imp*sum m^2, the m^2-only action)
is EXACT: the validation criterion is planned dS == executed dS to
machine precision. v1 plans with flips + deletions only (no insertions).
"""
from __future__ import annotations

import heapq
import os
import sys
from itertools import combinations

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _e(a, b):
    return (a, b) if a < b else (b, a)


class CollarState:
    """Decorated link of v: faces, ambient edge degrees, m-counters."""

    def __init__(self, L, v, cimp, pin_coef, e_target, f30, f10,
                 f3_target):
        self.v = v
        self.cimp = cimp
        self.pin = pin_coef
        self.et = e_target
        self.f3t = f3_target
        self.faces = frozenset(
            tuple(sorted(set(t) - {v})) for t in L.v2t[v])
        verts = sorted({x for f in self.faces for x in f})
        self.deg = {}                      # tracked ambient edge degrees
        for u in verts:
            self.deg[_e(v, u)] = L.edeg[_e(v, u)]
        for x, y in combinations(verts, 2):
            d = L.edeg.get(_e(x, y), 0)
            if d:
                self.deg[_e(x, y)] = d
        self.m = {}                        # m-counter per tracked vertex
        for u in verts + [v]:
            self.m[u] = sum(1 for ed, d in L.edeg.items()
                            if u in ed and d not in (5, 6))
        self.f3 = f30
        self.f1 = f10
        self.dS = 0.0

    # -- pieces of the exact action ---------------------------------------
    def _glob(self, f3, f1):
        return (0.1 * (f3 - self.f3t) ** 2
                + self.pin * (f1 - 6.0 * f3 / self.et) ** 2)

    def _shift(self, edges, delta, created=(), removed=()):
        """Apply degree shifts; return m^2 delta; update m-counters."""
        dm2 = 0.0
        for ed in removed:
            d0 = self.deg.pop(ed)
            if d0 not in (5, 6):
                for u in ed:
                    dm2 += (self.m[u] - 1) ** 2 - self.m[u] ** 2
                    self.m[u] -= 1
        for ed in created:
            self.deg[ed] = 3
            for u in ed:
                dm2 += (self.m[u] + 1) ** 2 - self.m[u] ** 2
                self.m[u] += 1
        for ed in edges:
            d0 = self.deg[ed]
            d1 = d0 + delta[ed]
            self.deg[ed] = d1
            ill0, ill1 = d0 not in (5, 6), d1 not in (5, 6)
            if ill0 != ill1:
                step = 1 if ill1 else -1
                for u in ed:
                    dm2 += (self.m[u] + step) ** 2 - self.m[u] ** 2
                    self.m[u] += step
        return dm2

    def _move(self, sides_minus, sides_plus, created, removed, df3):
        g0 = self._glob(self.f3, self.f1)
        delta = {}
        for ed in sides_minus:
            delta[ed] = delta.get(ed, 0) - 1
        for ed in sides_plus:
            delta[ed] = delta.get(ed, 0) + 1
        dm2 = self._shift(list(delta), delta, created, removed)
        self.f3 += df3
        self.f1 += df3 + len(created) - len(removed) - df3  # net edges
        # edges: +1 per created, -1 per removed
        self.f1 = self.f1  # (f1 alr. adjusted via created/removed count)
        g1 = self._glob(self.f3, self.f1)
        self.dS += self.cimp * dm2 + (g1 - g0)

    # -- 2D operations -----------------------------------------------------
    def link_edges(self):
        out = set()
        for f in self.faces:
            for p in combinations(f, 2):
                out.add(_e(*p))
        return out

    def flip_candidates(self):
        """[(a,b),(x,y)] flippable link edges with their opposite pair."""
        wing = {}
        for f in self.faces:
            for p in combinations(f, 2):
                wing.setdefault(_e(*p), []).append(
                    next(iter(set(f) - set(p))))
        out = []
        for ed, ws in wing.items():
            if len(ws) == 2 and _e(*ws) not in self.link_edges():
                out.append((ed, _e(*ws)))
        return out

    def flip(self, ab, xy, flavor):
        v = self.v
        a, b = ab
        x, y = xy
        f1_, f2_ = tuple(sorted((a, b, x))), tuple(sorted((a, b, y)))
        self.faces = frozenset(
            (self.faces - {f1_, f2_})
            | {tuple(sorted((x, y, a))), tuple(sorted((x, y, b)))})
        if flavor == "23":       # center (v,a,b), coc (x,y)
            self._move(
                sides_minus=[_e(v, a), _e(v, b), _e(a, b)],
                sides_plus=[_e(v, x), _e(v, y), _e(a, x), _e(a, y),
                            _e(b, x), _e(b, y)],
                created=[_e(x, y)], removed=[], df3=+1)
        else:                    # '32': axis (a,b), coc (v,x,y)
            self._move(
                sides_minus=[_e(v, a), _e(v, b), _e(a, x), _e(a, y),
                             _e(b, x), _e(b, y)],
                sides_plus=[_e(v, x), _e(v, y), _e(x, y)],
                created=[], removed=[_e(a, b)], df3=-1)

    def delete(self, u):
        """3->2 on spoke (v,u); u's link nbrs (x,y,z)."""
        v = self.v
        nbrs = sorted({w for f in self.faces if u in f for w in f} - {u})
        x, y, z = nbrs
        self.faces = frozenset(
            (f for f in self.faces if u not in f)) | frozenset(
            [tuple(sorted((x, y, z)))])
        self._move(
            sides_minus=[_e(v, x), _e(v, y), _e(v, z),
                         _e(u, x), _e(u, y), _e(u, z)],
            sides_plus=[_e(x, y), _e(y, z), _e(x, z)],
            created=[], removed=[_e(v, u)], df3=-1)

    def key(self):
        return (self.faces,
                tuple(sorted(self.deg.items())), round(self.dS, 9))

    def clone(self):
        import copy
        c = object.__new__(CollarState)
        c.__dict__.update(self.__dict__)
        c.deg = dict(self.deg)
        c.m = dict(self.m)
        return c


def plan_collapse(L, v, cimp, pin_coef, e_target, f30, f10, f3_target,
                  target_z=4, nodecap=20000, optimize="z"):
    """optimize='z': fast Z-first (first valid plan). optimize='cost':
    cost-first best-first (priority dS, then Z) -- returns the cheapest
    plan found before nodecap; exact-optimal when step costs are
    nonnegative along the frontier (deletion refunds make this a strong
    heuristic rather than a guarantee)."""
    """Best-first plan reducing |lk(v)| to target_z. Returns
    (ops, planned_dS) with ops = [('flip', ab, xy, flavor) |
    ('delete', u)], or (None, None)."""
    s0 = CollarState(L, v, cimp, pin_coef, e_target, f30, f10, f3_target)
    Z0 = len({x for f in s0.faces for x in f})
    cnt = 0
    def prio(z, ds, ln):
        return ((ds, z, ln) if optimize == "cost" else (z, ds, ln))
    heap = [(prio(Z0, 0.0, 0), cnt, s0, [])]
    seen = set()
    while heap and cnt < nodecap:
        (_, _, _), _, st, ops = heapq.heappop(heap)
        Z = len({x for f in st.faces for x in f})
        if Z <= target_z:
            return ops, st.dS
        k = st.key()
        if k in seen:
            continue
        seen.add(k)
        # deletions first
        for f in st.faces:
            pass
        verts = {x for fc in st.faces for x in fc}
        for u in sorted(verts):
            if st.deg.get(_e(st.v, u)) == 3:
                nx = st.clone()
                nx.delete(u)
                cnt += 1
                heapq.heappush(heap, (
                    prio(len({x for fc in nx.faces for x in fc}),
                         round(nx.dS, 9), len(ops) + 1),
                    cnt, nx, ops + [("delete", u)]))
        for ab, xy in st.flip_candidates():
            for flavor in ("23", "32"):
                if flavor == "23" and xy in st.deg:
                    continue                  # chord present: blocked
                if flavor == "32" and st.deg.get(ab) != 3:
                    continue
                nx = st.clone()
                nx.flip(ab, xy, flavor)
                cnt += 1
                heapq.heappush(heap, (
                    prio(len({x for fc in nx.faces for x in fc}),
                         round(nx.dS, 9), len(ops) + 1),
                    cnt, nx, ops + [("flip", ab, xy, flavor)]))
    return None, None


def lift(L, v, ops):
    """Execute a 2D plan as 3D moves through the Live driver. Returns the
    3D move list. Raises on any validity failure (repair = v2)."""
    moves = []
    for op in ops:
        if op[0] == "delete":
            u = op[1]
            tets = [t for t in L.v2t[v] if u in t]
            lk = sorted({x for t in tets for x in t} - {v, u})
            cen, coc = _e(v, u), tuple(lk)
        else:
            _, ab, xy, flavor = op
            if flavor == "23":
                cen = tuple(sorted((v,) + ab))
                coc = xy
            else:
                cen, coc = ab, tuple(sorted((v,) + xy))
        L.do(cen, coc)
        moves.append((cen, coc))
    return moves
