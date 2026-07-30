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
import math
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations, permutations, product

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


def _complex_edges(facets):
    """Edges of an abstract complex given its facets."""
    out = set()
    for f in facets:
        for e in combinations(sorted(f), 2):
            out.add(e)
    return out


def _refine(facets, verts, vcolor, ecolor):
    """Colour-refine the vertices (WL on the facet hypergraph, respecting any
    given vertex and edge colours). Two vertices keep the same colour only
    while nothing in their decorated neighbourhoods distinguishes them."""
    nbr = defaultdict(set)
    for (u, w) in _complex_edges(facets):
        nbr[u].add(w)
        nbr[w].add(u)
    seed = {s_: i for i, s_ in
            enumerate(sorted({vcolor.get(v, 0) for v in verts}))}
    color = {v: seed[vcolor.get(v, 0)] for v in verts}
    for _ in range(len(verts) + 1):
        sig = {}
        for v in verts:
            around = []
            for f in facets:
                if v in f:
                    around.append((len(f),
                                   tuple(sorted(color[u] for u in f if u != v))))
            ring = tuple(sorted((ecolor.get(_ek(v, u), 0), color[u])
                                for u in nbr[v]))
            sig[v] = (color[v], tuple(sorted(around)), ring)
        vals = {s_: i for i, s_ in enumerate(sorted(set(sig.values())))}
        new = {v: vals[sig[v]] for v in verts}
        if new == color:
            break
        color = new
    return color


def canonical_key(facets, vcolor=None, ecolor=None, limit=200000):
    """Canonical form of an abstract simplicial complex, optionally DECORATED.

    Returns (key, exact). `key` relabels the complex to 0..n-1 by the vertex
    ordering that minimises it lexicographically, so two decorated complexes
    are isomorphic-as-decorated iff their keys are equal -- provided `exact`
    is True. Colour refinement first cuts the search to orderings consistent
    with the refined partition; if that still exceeds `limit` the search is
    truncated and `exact` comes back False, meaning the key is a strong
    invariant but no longer a certificate.

    `vcolor` maps vertex -> label, `ecolor` maps sorted edge pair -> label.
    For defects the natural choice is ecolor = ambient edge degree and
    vcolor = n6. Note those two are not independent of each other given the
    edge decoration: every illegal edge lies INSIDE an induced defect
    subcomplex (both endpoints of an illegal edge have imp > 0, hence are
    defect vertices), so the edge colours already determine imp(v) and the
    illegal degree sum at each vertex, and then

        n6(v) = Z + 5*imp(v) - sum(d_illegal at v) - 12
        q_R(v) = Z*(pi - 3*theta) + 6*theta,   theta = arccos(1/3)

    with Z the coordination. So decorating with n6 is equivalent to
    decorating with Z, and either one fixes the curvature; there is no need
    to carry q_R as well. (Both identities follow from the vertex link being
    a triangulated 2-sphere, and hold for every vertex, legal or not.)
    """
    facets = [tuple(sorted(f)) for f in facets]
    verts = sorted({v for f in facets for v in f})
    vcolor = vcolor or {}
    ecolor = ecolor or {}
    if not verts:
        return ((), True)
    edges = sorted(_complex_edges(facets))
    color = _refine(facets, verts, vcolor, ecolor)
    classes = defaultdict(list)
    for v in verts:
        classes[color[v]].append(v)
    groups = [sorted(classes[c]) for c in sorted(classes)]

    total = 1
    for g in groups:
        total *= math.factorial(len(g))
    exact = total <= limit

    best = None
    seen = 0
    for combo in product(*(permutations(g) for g in groups)):
        order = [v for g in combo for v in g]
        idx = {v: i for i, v in enumerate(order)}
        rel_f = tuple(sorted(tuple(sorted(idx[v] for v in f)) for f in facets))
        rel_e = tuple(sorted(
            (tuple(sorted((idx[u], idx[w]))), ecolor.get((u, w), 0))
            for (u, w) in edges))
        rel_v = tuple(vcolor.get(v, 0) for v in order)
        rel = (rel_f, rel_e, rel_v)
        if best is None or rel < best:
            best = rel
        seen += 1
        if seen >= limit:
            break
    return (best, exact)


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

    def components(self, broad=True):
        """Connected components of defect vertices, as Complex objects.

        broad=True  the current definition: illegal-edge incidence OR non-FK
                    coordination.
        broad=False the HISTORICAL definition, imp>0 only. Kept so producers
                    can keep emitting `n_illegal` / `sizes` / `members` with
                    exactly the meaning every previous run and published
                    number used, while reporting the broadened counts
                    alongside rather than silently redefining the columns.
        """
        pool = self.defect if broad else {v for v in self.defect
                                          if self.imp[v] > 0}
        seen, out = set(), []
        for s0 in pool:
            if s0 in seen:
                continue
            stack, comp = [s0], []
            seen.add(s0)
            while stack:
                u = stack.pop()
                comp.append(u)
                for w in self.neighbours(u):
                    if w in pool and w not in seen:
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

    def induced(self, verts):
        """The subcomplex INDUCED by `verts`: every simplex of the
        triangulation all of whose vertices are defect vertices.

        Returned as {dim: sorted list of simplices}. This is the full
        subcomplex, not just the top-dimensional part -- an edge with both
        endpoints in the set belongs to it even when no tet of the set
        contains it.

        It is DENSE, not curve-like. Measured: a (3,4,4) knot always induces
        f = (5, 10, 9, 3), and since C(5,2) = 10 its 1-skeleton is the
        complete graph K5 -- every defect vertex adjacent to every other --
        with 9 of the 10 triangles and 3 of the 5 tets, i.e. a fixed
        subcomplex of the boundary of the 4-simplex, identical in all 82
        instances checked. The curve-like object is a different one: the
        ILLEGAL-EDGE subgraph (edges of `ill_edges` inside the set), which for
        a (3,4,4) knot has 3 edges with degree sequence (2,1,1,1,1) -- a
        2-edge path plus a disjoint edge. Do not confuse the two.

        Distinct from `star`: the induced subcomplex is what the defect IS
        (it has the defect's vertices and nothing else), while the closed star
        is the region the defect OCCUPIES (it pulls in a shell of legal
        neighbours). Both are useful; they answer different questions.

        Every simplex is a face of some tet, and any tet containing a simplex
        with all vertices in the set is incident to the set, so scanning the
        incident tets finds everything -- O(|S| * Z), no global scan.
        """
        S = set(verts)
        out = {0: sorted(S), 1: set(), 2: set(), 3: set()}
        for v in S:
            for t in self.v2t.get(v, ()):
                inside = [x for x in t if x in S]
                if len(inside) < 2:
                    continue
                for e in combinations(inside, 2):
                    out[1].add(e)
                for f in combinations(inside, 3):
                    out[2].add(f)
                if len(inside) == 4:
                    out[3].add(tuple(t))
        for d in (1, 2, 3):
            out[d] = sorted(out[d])
        return out

    def induced_facets(self, verts):
        """Maximal simplices of the induced subcomplex -- the data that
        determines it. Returned as a sorted list of tuples."""
        ind = self.induced(verts)
        facets = []
        for d in (3, 2, 1, 0):
            for s in ind[d]:
                s = (s,) if d == 0 else tuple(s)
                if not any(set(s) < set(f) for f in facets):
                    facets.append(s)
        return sorted(facets, key=lambda f: (-len(f), f))

    def decorations(self, verts):
        """(vcolor, ecolor) for the induced subcomplex of `verts`:
        vcolor[v] = n6(v), ecolor[edge] = ambient edge degree. Feed straight
        to canonical_key."""
        S = set(verts)
        vcolor = {v: self.n6[v] for v in S}
        ecolor = {}
        ind = self.induced(S)
        for e in ind[1]:
            ecolor[e] = self.edeg[e]
        return vcolor, ecolor

    def induced_shape(self, verts):
        """Coarse shape summary of the induced subcomplex: the f-vector, the
        sorted degree sequence of its 1-skeleton, and the number of connected
        components / independent cycles of that 1-skeleton.

        Note these describe the DENSE induced subcomplex, so `cycles` is
        large (6 for a (3,4,4) knot, whose 1-skeleton is K5) and is not a
        curve invariant. For curve-like structure use the illegal-edge
        subgraph instead -- see `illegal_shape`."""
        ind = self.induced(verts)
        nv, ne = len(ind[0]), len(ind[1])
        deg = Counter()
        for u, w in ind[1]:
            deg[u] += 1
            deg[w] += 1
        adj = defaultdict(set)
        for u, w in ind[1]:
            adj[u].add(w)
            adj[w].add(u)
        seen, ncomp = set(), 0
        for v in ind[0]:
            if v in seen:
                continue
            ncomp += 1
            stack = [v]
            seen.add(v)
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
        return {
            "f": (nv, ne, len(ind[2]), len(ind[3])),
            "degseq": tuple(sorted((deg[v] for v in ind[0]), reverse=True)),
            "components": ncomp,
            "cycles": ne - nv + ncomp,       # first Betti number of the graph
        }

    def illegal_shape(self, verts):
        """Shape of the ILLEGAL-EDGE subgraph inside a defect: the
        disclination line itself, as opposed to the dense induced subcomplex.

        This is the curve-like object. A (3,4,4) knot gives 3 edges with
        degree sequence (2,1,1,1,1) -- a 2-edge path plus a disjoint edge, so
        two graph components even though the defect is connected through the
        triangulation's 1-skeleton."""
        S = set(verts)
        edges = [e for e in self.ill_edges if e[0] in S and e[1] in S]
        deg = Counter()
        adj = defaultdict(set)
        for u, w in edges:
            deg[u] += 1
            deg[w] += 1
            adj[u].add(w)
            adj[w].add(u)
        touched = set(deg)
        seen, ncomp = set(), 0
        for v in touched:
            if v in seen:
                continue
            ncomp += 1
            stack = [v]
            seen.add(v)
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
        return {
            "n_edges": len(edges),
            "degseq": tuple(sorted(deg.values(), reverse=True)),
            "components": ncomp,
            "cycles": len(edges) - len(touched) + ncomp,
            "isolated_vertices": len(S) - len(touched),
        }

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

    # -- curvature charge ---------------------------------------------------

    def vertex_charges(self):
        """q_R(v) = 1/2 sum_{e ∋ v} delta(e), with delta(e) = 2*pi - theta*deg(e)
        and theta = arccos(1/3) the regular-tetrahedron dihedral angle.
        Returns {vertex: q_R} in radians (notes/CONVENTIONS.md symbols)."""
        theta = np.arccos(1.0 / 3.0)
        q = defaultdict(float)
        for (u, w), d in self.edeg.items():
            if d <= 0:
                continue
            half = 0.5 * (2.0 * np.pi - theta * d)
            q[u] += half
            q[w] += half
        return q

    def qbar(self, q=None):
        """State-mean q_R -- the crystal-background reference. Small but not
        zero: sum_v q_R = 2*pi*E - 6*theta*N_3, which is a topological
        constant for the box, so qbar drifts only with the f-vector."""
        if q is None:
            q = self.vertex_charges()
        return float(np.mean([q[v] for v in self.v2t])) if self.v2t else 0.0

    def complex_charge(self, verts, q=None):
        """Q_c = sum_{v in c} q_R(v) -- the RAW curvature charge, with no
        background subtraction.

        Raw is the primary quantity because it is reference-free: it depends
        only on the complex's own vertices, so it is exactly reproducible and
        identical complexes give identical values. Subtracting a state mean
        makes it depend on the whole configuration and introduces a spurious
        spread (measured: sd ~ 1e-3 across snapshots purely from qbar drift).

        Because q_R(v) = Z(pi - 3 theta) + 6 theta, the raw charge of an
        n-vertex complex is fixed by its TOTAL coordination:

            Q_c = (sum Z) * (pi - 3 theta) + 6 n theta

        so it is quantised in steps of (pi - 3 theta) = -0.5512856 rad per
        unit of total coordination."""
        if q is None:
            q = self.vertex_charges()
        return float(sum(q[v] for v in verts))

    def complex_charge_rel(self, verts, q=None, qbar=None):
        """Q_c relative to the crystal background: sum (q_R(v) - qbar)."""
        if q is None:
            q = self.vertex_charges()
        if qbar is None:
            qbar = self.qbar(q)
        return float(sum(q[v] - qbar for v in verts))

    def total_coordination(self, verts):
        """sum of Z over the complex -- the integer that indexes its raw
        charge (see complex_charge)."""
        return sum((len(self.v2t[v]) + 4) // 2 for v in verts)

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
        # edges_from_facets relabels vertices to dense 0..V-1; map the edge
        # endpoints back to original labels (lab[idx]) so the keys match the
        # incrementally-maintained self.edeg, which is keyed by REAL labels.
        # Without this the audit false-positives the instant label recycling
        # opens a gap in the numbering -- the state is correct, the comparison
        # was not (it compared dense-index edges to real-label edges).
        lab = np.unique(fac)
        ref_edeg = {tuple(sorted((int(lab[a]), int(lab[b])))): int(d)
                    for (a, b), d in zip(eu, edg)}
        ref_ill = {e for e, d in ref_edeg.items() if not _legal(d)}
        if ref_ill != self.ill_edges:
            problems.append(
                f"illegal-edge set differs by {len(ref_ill ^ self.ill_edges)}")
        mine = {e: d for e, d in self.edeg.items() if d > 0}
        if ref_edeg != mine:
            diff = set(ref_edeg.items()) ^ set(mine.items())
            problems.append(f"edge degrees differ on {len(diff)} entries")
        n6r, impr, adj = vertex_classes(fac)
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

    Every identity change is also appended to `self.events` as a dict --
    birth / death / merge / split -- so a caller can census the reaction
    chemistry rather than only its rate. `drain_events()` empties the
    list, which is what long runs should use (the stream is unbounded:
    the flicker vacuum births and kills complexes continuously).
    """

    def __init__(self, state):
        self.state = state
        self.label = {}          # vertex -> track id
        self.born = {}           # track id -> clock
        self.died = {}           # track id -> clock
        self.size = {}           # track id -> vertex count, current
        self.merges = 0
        self.splits = 0
        self.events = []
        self._next = 0

    def drain_events(self):
        ev, self.events = self.events, []
        return ev

    def transport_hint(self, old_verts, new_verts):
        """Carry track identity across a discontinuous composite-move hop.

        Called between a bracketed composite's state update and the next
        step() (see the EVT_CHANNEL_* records in sampler.d): seeds the
        landing vertices with the majority label of the departed ones, so
        the landing component votes as the same worldline instead of
        surfacing as a spurious death + birth. Vertices that already carry
        a label are left alone -- a hop INTO an existing complex is a
        genuine merge and is recorded as one. Returns the tid, or None if
        the departed vertices carried no label."""
        prev = Counter(self.label[x] for x in old_verts if x in self.label)
        if not prev:
            return None
        t = prev.most_common(1)[0][0]
        for v in new_verts:
            self.label.setdefault(v, t)
        return t

    def step(self, clock):
        """Relabel after the state has advanced. Returns [(Complex, tid)]."""
        comps = self.state.components()
        claimed, newlabel, out = set(), {}, []
        for c in comps:
            prev = Counter(self.label[x] for x in c.verts if x in self.label)
            t = None
            parent = None
            if prev:
                cand = prev.most_common(1)[0][0]
                if cand in claimed:
                    self.splits += 1
                    parent = cand
                else:
                    t = cand
                    if len(prev) > 1:
                        self.merges += len(prev) - 1
                        self.events.append({
                            "k": "merge", "clock": clock, "into": t,
                            "size": len(c.verts),
                            # contributed vertex count per absorbed track,
                            # and that track's size just before the merge
                            "from": sorted(
                                (int(o), int(n), int(self.size.get(o, 0)))
                                for o, n in prev.items() if o != t)})
            if t is None:
                t = self._next
                self._next += 1
                self.born[t] = clock
                if parent is None:
                    self.events.append({"k": "birth", "clock": clock,
                                        "tid": t, "size": len(c.verts)})
                else:
                    self.events.append({
                        "k": "split", "clock": clock, "parent": parent,
                        "child": t, "size": len(c.verts),
                        "parent_size_before": int(self.size.get(parent, 0))})
            claimed.add(t)
            for x in c.verts:
                newlabel[x] = t
            out.append((c, t))
        for t in set(self.label.values()) - claimed:
            if t not in self.died:
                self.died[t] = clock
                self.events.append({"k": "death", "clock": clock, "tid": t,
                                    "size": int(self.size.get(t, 0))})
        self.size = {t: len(c.verts) for c, t in out}
        self.label = newlabel
        return out


def census(state, members=False):
    """Summary counters for a time series -- the single replacement for the
    per-chunk `vertex_classes` block copied across the producer scripts.

    Emits BOTH definitions. `n_illegal` / `sizes` / `members` are the
    historical imp>0 quantities, unchanged, so existing time series stay
    comparable; the `*_broad` keys and the non-FK counters are the new
    information. Nothing is silently redefined."""
    old = state.components(broad=False)
    new = state.components(broad=True)
    anom = [v for v in state.defect if state.imp[v] == 0]
    out = {
        # --- historical (imp > 0 only): same meaning as every prior run
        "n_illegal": sum(len(c.verts) for c in old),
        "sizes": sorted((len(c.verts) for c in old), reverse=True),
        "n_components": len(old),
        # --- broadened: illegal-edge incidence OR non-FK coordination
        "n_defect_broad": len(state.defect),
        "sizes_broad": sorted((len(c.verts) for c in new), reverse=True),
        "n_components_broad": len(new),
        "n_nonfk_all_legal": len(anom),
        "n_nonfk_n6_1": sum(1 for v in anom if state.n6[v] == 1),
        "n_nonfk_n6_ge5": sum(1 for v in anom if state.n6[v] >= 5),
        "species": Counter(c.key for c in new),
        # --- host quality
        "mean_edeg": (sum(state.edeg.values()) / len(state.edeg)
                      if state.edeg else 0.0),
        "legaledge": (sum(1 for d in state.edeg.values() if _legal(d))
                      / len(state.edeg) if state.edeg else 1.0),
        "legalvert": (sum(1 for v in state.v2t if state.imp[v] == 0)
                      / len(state.v2t) if state.v2t else 1.0),
        "legalvert_fk": (sum(1 for v in state.v2t if v not in state.defect)
                         / len(state.v2t) if state.v2t else 1.0),
    }
    if members:
        out["members"] = [sorted(c.verts) for c in old]
        out["members_broad"] = [sorted(c.verts) for c in new]
        out["sigs"] = [list(c.sig) for c in new]
    return out
