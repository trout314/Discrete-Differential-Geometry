#!/usr/bin/env python3
"""Incremental defect bookkeeping over the accepted-move event stream.

Every defect script in the tree currently rebuilds the same thing from
scratch: marshal all facets -> edges_from_facets -> vertex_classes -> flood
components. That is O(N) per sample (~800 ms at N3 = 58k, measured), which is
only ~3% of runtime at the production 150-sweep cadence but 83% at per-sweep
resolution, and it caps how finely defects can be followed.

This module maintains the same quantities incrementally instead, driven by the
accepted-move stream (event_replay.event_changes). The stream is exact -- the
replay reproduces the live manifold move for move -- and sparse: about 5
accepted moves per sweep out of ~23M attempted, so per-move bookkeeping costs
nothing. Everything is reconstructable from the move record, so no additional
state is needed in the D core.

WHAT COUNTS AS A DEFECT
-----------------------
A defect vertex is one that is either

  * incident to an ILLEGAL edge (degree not in {5, 6}), or
  * of non-FK coordination: its number of incident degree>=6 edges ("six
    edges") is not in {0, 2, 3, 4}.

The second clause is new, and it closes a real blind spot: the six-edges form
the disclination network, and n6 counts the disclination lines meeting at a
vertex -- 0 = Z12, 2 = a line passing through (Z14), 3 = a 3-fold node (Z15),
4 = a 4-fold node (Z16). A vertex with all-legal edges but n6 = 5 is a defect
that the FK potential actively penalises (U = zleg * dist^2(n6, {0,2,3,4}))
yet which contributes NOTHING to an imp>0 census: it appears in no complex,
and mobile_gas's `legalvert = mean(imp == 0)` scores it as legal.

The two anomalies differ in kind and are counted separately:
  n6 == 1  a disclination line terminating in the bulk -- topologically
           forbidden, and indeed never observed in any state checked;
  n6 >= 5  a high-order node, allowed but geometrically disfavoured (Z16 is
           the Frank-Kasper upper bound). Observed, but rare.

Measured on three lam=0.40 snapshots: 0, 2 and 7 such vertices against 46-63
illegal-edge vertices, and every one of them ATTACHED to an existing complex
rather than forming its own -- component counts were unchanged (8 -> 8) and
only sizes grew. So broadening the definition fattens defects slightly; it
does not percolate.

For an all-legal vertex there is an exact identity, n6 = Z - 12, where Z is
the coordination (number of neighbours): the link is a triangulated 2-sphere
on Z vertices, so it has 2Z - 4 triangles and its degrees sum to 6Z - 12,
whence n5 + n6 = Z and 5*n5 + 6*n6 = 6Z - 12. So "all-legal with n6 > 4" is
exactly "coordination Z >= 17". `audit` checks this identity as a free
consistency test on the incremental counters.
"""
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import event_replay as er
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

FK_N6 = frozenset((0, 2, 3, 4))     # Z12 / Z14 / Z15 / Z16


def _ek(u, w):
    return (u, w) if u < w else (w, u)


def _legal(d):
    return d == 5 or d == 6


class Complex:
    """One connected defect: its vertices, its label, and its closed star.

    The label keeps the historical edge-degree signature as the primary
    species name -- so "(3,4,4)" still means exactly what it always has --
    with the coordination anomalies reported alongside rather than folded in.
    A bare knot is sig == (3,4,4) with nodes == (); one decorated by a 5-fold
    node is (3,4,4) with nodes == (5,).
    """
    __slots__ = ("verts", "sig", "nodes")

    def __init__(self, verts, sig, nodes):
        self.verts = verts          # sorted list of defect vertices
        self.sig = sig              # sorted illegal edge degrees within it
        self.nodes = nodes          # sorted non-FK n6 values within it

    @property
    def key(self):
        return (self.sig, self.nodes)

    def __repr__(self):
        return (f"Complex(n={len(self.verts)}, sig={self.sig}, "
                f"nodes={self.nodes})")


class DefectState:
    """Defect bookkeeping maintained move by move.

    Seed from a manifold (O(N) once), then feed it the drained event batches.
    Components and closed stars are computed on demand from the maintained
    state -- O(|defect vertices|), not O(N).
    """

    def __init__(self, mfd, keep_tets=True):
        fac = np.asarray(mfd.facets())
        self.tets = {tuple(sorted(int(x) for x in t)) for t in fac}
        self.edeg = Counter()
        self.v2t = defaultdict(set)
        for t in self.tets:
            for e in combinations(t, 2):
                self.edeg[e] += 1
            for v in t:
                self.v2t[v].add(t)
        self.n6 = defaultdict(int)
        self.imp = defaultdict(int)
        for e, d in self.edeg.items():
            if d >= 6:
                self.n6[e[0]] += 1
                self.n6[e[1]] += 1
            if not _legal(d):
                self.imp[e[0]] += 1
                self.imp[e[1]] += 1
        # The illegal edges, maintained too: a complex's signature is read
        # off THIS set, not by filtering every edge in the manifold (which
        # would be O(E) per component and was the dominant cost when this was
        # first written).
        self.ill_edges = {e for e, d in self.edeg.items()
                          if d > 0 and not _legal(d)}
        self.defect = {v for v in self.v2t if self._is_defect(v)}
        self._keep_tets = keep_tets
        self.n_events = 0

    # -- defect predicate ---------------------------------------------------

    def _is_defect(self, v):
        return self.imp[v] > 0 or self.n6[v] not in FK_N6

    def _refresh(self, v):
        if v in self.v2t and self.v2t[v] and self._is_defect(v):
            self.defect.add(v)
        else:
            self.defect.discard(v)

    # -- incremental update -------------------------------------------------

    def apply(self, ev):
        """Fold one accepted-move event into the state."""
        rem, add = er.event_changes(ev)
        # Edges whose degree may change; snapshot BEFORE so a single edge
        # touched by several tets of the same move is counted once.
        touched = set()
        for t in rem + add:
            for p in combinations(t, 2):
                touched.add(p)
        before = {p: self.edeg.get(p, 0) for p in touched}

        for t in rem:
            if self._keep_tets:
                self.tets.discard(t)
            for v in t:
                self.v2t[v].discard(t)
            for p in combinations(t, 2):
                self.edeg[p] -= 1
        for t in add:
            if self._keep_tets:
                self.tets.add(t)
            for v in t:
                self.v2t[v].add(t)
            for p in combinations(t, 2):
                self.edeg[p] += 1

        verts = set()
        for p in touched:
            d0, d1 = before[p], self.edeg.get(p, 0)
            if d0 == d1:
                continue
            for v in p:
                verts.add(v)
            self.n6[p[0]] += (d1 >= 6) - (d0 >= 6)
            self.n6[p[1]] += (d1 >= 6) - (d0 >= 6)
            ill0 = d0 > 0 and not _legal(d0)
            ill1 = d1 > 0 and not _legal(d1)
            self.imp[p[0]] += ill1 - ill0
            self.imp[p[1]] += ill1 - ill0
            if ill1 and not ill0:
                self.ill_edges.add(p)
            elif ill0 and not ill1:
                self.ill_edges.discard(p)
            if d1 == 0:
                del self.edeg[p]
                self.ill_edges.discard(p)
        for t in rem + add:
            verts.update(t)
        for v in verts:
            if not self.v2t.get(v):          # vertex destroyed (4->1)
                self.v2t.pop(v, None)
                self.n6.pop(v, None)
                self.imp.pop(v, None)
                self.defect.discard(v)
            else:
                self._refresh(v)
        self.n_events += 1

    def apply_all(self, events):
        for ev in events:
            self.apply(ev)

    # -- queries ------------------------------------------------------------

    def neighbours(self, v):
        """Vertices sharing a tet with v -- from the incidence map, no scan."""
        out = set()
        for t in self.v2t.get(v, ()):
            out.update(t)
        out.discard(v)
        return out

    def components(self):
        """Connected components of defect vertices, as Complex objects."""
        seen, out = set(), []
        for s0 in self.defect:
            if s0 in seen:
                continue
            stack, comp = [s0], []
            seen.add(s0)
            while stack:
                u = stack.pop()
                comp.append(u)
                for w in self.neighbours(u):
                    if w in self.defect and w not in seen:
                        seen.add(w)
                        stack.append(w)
            cv = set(comp)
            sig = sorted(self.edeg[e] for e in self.ill_edges
                         if e[0] in cv and e[1] in cv)
            # Only ALL-LEGAL vertices contribute here. A vertex touching an
            # illegal edge is already a defect and already described by `sig`;
            # what `nodes` adds is the coordination-only defects, the ones the
            # historical imp>0 census could not see at all. (Restricting this
            # also keeps the label consistent with the census counters, which
            # are likewise over all-legal vertices.)
            nodes = sorted(self.n6[v] for v in comp
                           if self.imp[v] == 0 and self.n6[v] not in FK_N6)
            out.append(Complex(sorted(comp), tuple(sig), tuple(nodes)))
        return out

    def star(self, verts):
        """Closed star of a vertex set: every tet incident to any of them.

        Faces are implicit (a simplicial complex is determined by its
        facets), so the tet set IS the closed star."""
        out = set()
        for v in verts:
            out |= self.v2t.get(v, set())
        return out

    def star_boundary(self, verts):
        """Triangles of the closed star lying on its boundary -- those in
        exactly one star tet. This is the surface through which disclination
        flux would be measured."""
        cnt = Counter()
        for t in self.star(verts):
            for f in combinations(t, 3):
                cnt[f] += 1
        return {f for f, c in cnt.items() if c == 1}

    # -- audit --------------------------------------------------------------

    def audit(self, mfd):
        """Compare the incremental state against a from-scratch rebuild.

        Returns a list of discrepancy strings (empty when clean). Also checks
        the identity n6 == Z - 12 on every all-legal vertex, which is an
        independent test of the counters that does not go through
        vertex_classes at all."""
        problems = []
        fac = np.asarray(mfd.facets())
        if self._keep_tets:
            live = {tuple(sorted(int(x) for x in t)) for t in fac}
            if live != self.tets:
                problems.append(
                    f"tet set differs by {len(live ^ self.tets)}")
        eu, edg, V = edges_from_facets(fac)
        ref_edeg = {tuple(int(x) for x in e): int(d) for e, d in zip(eu, edg)}
        ref_ill = {e for e, d in ref_edeg.items() if not _legal(d)}
        if ref_ill != self.ill_edges:
            problems.append(
                f"illegal-edge set differs by {len(ref_ill ^ self.ill_edges)}")
        mine = {e: d for e, d in self.edeg.items() if d > 0}
        if ref_edeg != mine:
            diff = set(ref_edeg.items()) ^ set(mine.items())
            problems.append(f"edge degrees differ on {len(diff)} entries")
        n6r, impr, adj = vertex_classes(fac)
        lab = np.unique(fac)
        for i, v in enumerate(lab):
            v = int(v)
            if self.n6.get(v, 0) != int(n6r[i]):
                problems.append(f"n6[{v}] {self.n6.get(v,0)} != {int(n6r[i])}")
            if self.imp.get(v, 0) != int(impr[i]):
                problems.append(f"imp[{v}] {self.imp.get(v,0)} != {int(impr[i])}")
            if len(problems) > 8:
                break
        # independent check: for an all-legal vertex, n6 == Z - 12, with the
        # coordination Z recovered from the star size (2Z - 4 tets at v)
        for v in list(self.v2t)[:2000]:
            if self.imp.get(v, 0):
                continue
            z = (len(self.v2t[v]) + 4) // 2
            if self.n6.get(v, 0) != z - 12:
                problems.append(
                    f"identity n6 == Z-12 broken at {v}: "
                    f"n6={self.n6.get(v,0)} Z={z}")
                break
        return problems


class Worldlines:
    """Per-move identity over defect complexes.

    One Pachner move touches at most 5 vertices (center u coCenter), so every
    complex not meeting that set is untouched -- identical vertex set,
    identical structure -- and its identity is preserved trivially. Only the
    complexes meeting those few vertices can change, and they inherit the
    label held by the majority of their previous members. At one-move
    granularity that is essentially exact; the inference that frame-based
    linkage has to make over a 150-sweep gap does not arise.

    Merges and splits are recorded rather than silently becoming deaths and
    births: a component whose majority label is already claimed by another
    component this step is a split-off and gets a fresh id.
    """

    def __init__(self, state):
        self.state = state
        self.label = {}          # vertex -> track id
        self.born = {}           # track id -> clock
        self.died = {}           # track id -> clock
        self.merges = 0
        self.splits = 0
        self._next = 0

    def step(self, clock):
        """Relabel after the state has advanced. Returns [(Complex, tid)]."""
        comps = self.state.components()
        claimed, newlabel, out = set(), {}, []
        for c in comps:
            prev = Counter(self.label[x] for x in c.verts if x in self.label)
            t = None
            if prev:
                cand = prev.most_common(1)[0][0]
                if cand in claimed:
                    self.splits += 1
                else:
                    t = cand
                    if len(prev) > 1:
                        self.merges += len(prev) - 1
            if t is None:
                t = self._next
                self._next += 1
                self.born[t] = clock
            claimed.add(t)
            for x in c.verts:
                newlabel[x] = t
            out.append((c, t))
        for t in set(self.label.values()) - claimed:
            self.died.setdefault(t, clock)
        self.label = newlabel
        return out


def census(state):
    """Summary counters for a time series -- the replacement for the
    per-chunk `vertex_classes` block copied across the producer scripts."""
    comps = state.components()
    anom = [v for v in state.defect if state.imp[v] == 0]
    return {
        "n_defect": len(state.defect),
        "n_illedge_vertices": sum(1 for v in state.defect if state.imp[v]),
        "n_nonfk_all_legal": len(anom),
        "n_nonfk_n6_1": sum(1 for v in anom if state.n6[v] == 1),
        "n_nonfk_n6_ge5": sum(1 for v in anom if state.n6[v] >= 5),
        "n_components": len(comps),
        "sizes": sorted((len(c.verts) for c in comps), reverse=True),
        "species": Counter(c.key for c in comps),
    }
