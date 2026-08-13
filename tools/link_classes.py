#!/usr/bin/env python3
"""Exact classification of pure-{5,6} vertex links, and the FK legality detector.

Why this exists
---------------
In a closed 3-manifold triangulation, lk(v) is a triangulated S^2 and
deg_M({v,u}) = deg_lk(v)(u).  So "every edge at v has degree 5 or 6" says lk(v)
is a 5/6-triangulated 2-sphere.  Euler gives the local sum rule
sum_{e ni v}(6 - deg e) = 12, hence n_5 = 12 always, deg v = 12 + n_6, and the
number of tets at v is 20 + 2*n_6.  Minimum degree 5 makes the link graph
3-connected planar, so by Steinitz/Whitney the graph determines the sphere and
these links are exactly the **duals of the fullerenes C_{20+2 n_6}**.

The consequence the census cares about: **n_6 does NOT determine lk(v)**.  It
does for n_6 <= 3, but at n_6 = 4 there are TWO classes -- the T_d Friauf
polyhedron (the Frank-Kasper Z16) and a D2 isomer that no FK phase admits.
Bucketing vertices by n_6 alone, as `pure56 & (n6 == 4)`, silently counts the
D2 links as Z16.  A state can therefore have zero illegal edges (f_FK = 1) and
still not be Frank-Kasper.

Enumerated here exhaustively (`enumerate_5_6_spheres`, no external tables):

    n_6   V    classes  |Aut|                cofacial six-pairs at v
    ---------------------------------------------------------------------
    0     12   1        120                  0        Z12  (icosahedron)
    1     13   0        --                   --       IMPOSSIBLE (no C22)
    2     14   1        24                   0        Z14  (D6d)
    3     15   1        12                   0        Z15  (D3h)
    4     16   2        24 / 4               0 / 2    Z16 T_d (Friauf) / D2
    5     17   3        20, 4, 4             5, 3, 2
    6     18   6        4, 2, 12, 2, 12, 6   6, 5, 6, 4, 6, 3
    7     19   6        ...                  min 4
    8     20   15       ...                  min 4

Two degree-6 vertices are adjacent in lk(v) iff the corresponding six-edges at
v are **cofacial** -- they span a triangle with v.  Reading the last column:
the four FK links Z12/Z14/Z15/Z16-T_d are precisely the 5/6-links whose
six-edges are pairwise non-cofacial, and from n_6 = 5 up no class achieves it.
So non-cofaciality is a purely local predicate equivalent to FK legality, and
it is what caps n_6 at 4.  That is the detector wired into the census: count,
per vertex, the pairs of degree-6 edges sharing a triangle with v.

Usage
-----
    python tools/link_classes.py --enumerate 6           # regenerate the table
    python tools/link_classes.py data/tcp_reference/T3_R_m2_N1360.mfd
    python tools/link_classes.py --verify <seed.mfd>     # cross-check detector
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (n_6) -> number of isomorphism classes of 5/6-triangulated 2-spheres.
# Verified by enumerate_5_6_spheres() for n_6 <= 8; equals the fullerene
# isomer count of C_{20+2 n_6}.
CLASS_COUNTS = {0: 1, 1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 6, 8: 15}

# Names for the FK-legal links, keyed by n_6 (cofacial pair count == 0).
FK_NAMES = {0: "Z12", 2: "Z14", 3: "Z15", 4: "Z16"}

# Largest n_6 admitting a link with pairwise non-cofacial six-edges.
MAX_FK_N6 = 4


# --------------------------------------------------------------------------
# Fast manifold-level detector
# --------------------------------------------------------------------------
_TRI_OF_TET = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
_EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def relabel(F):
    """(T, V): facets relabeled to 0..V-1, matching fk_skeleton.edges_from_facets."""
    F = np.asarray(F, np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    return inv.reshape(F.shape), len(lab)


def cofacial_six_pairs(F, eu, edeg, V):
    """Per-vertex count of pairs of degree-6 edges at v that span a triangle with v.

    Equivalently: the number of edges in the subgraph of lk(v) induced on its
    degree-6 vertices.  Zero for every Frank-Kasper link (Z12/Z14/Z15/Z16-T_d);
    positive exactly on the non-FK 5/6-links, the smallest being the D2 isomer
    at n_6 = 4.

    `eu, edeg, V` must come from fk_skeleton.edges_from_facets(F) -- the vertex
    relabeling is recomputed here identically, so labels agree.
    """
    T, V2 = relabel(F)
    assert V2 == V, "vertex count disagrees with the supplied edge census"

    # edge -> degree lookup by lexicographic key (eu rows are sorted a<b)
    key = eu[:, 0].astype(np.int64) * V + eu[:, 1]
    order = np.argsort(key, kind="stable")
    key_s, deg_s = key[order], edeg[order]

    def edge_deg(a, b):
        lo = np.minimum(a, b).astype(np.int64)
        hi = np.maximum(a, b).astype(np.int64)
        idx = np.searchsorted(key_s, lo * V + hi)
        return deg_s[idx]

    # unique triangles of the complex
    tri = np.sort(np.vstack([T[:, list(t)] for t in _TRI_OF_TET]), axis=1)
    tri = np.unique(tri, axis=0)

    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    s_ab = edge_deg(a, b) == 6
    s_ac = edge_deg(a, c) == 6
    s_bc = edge_deg(b, c) == 6

    cof = np.zeros(V, dtype=np.int64)
    np.add.at(cof, a, (s_ab & s_ac).astype(np.int64))
    np.add.at(cof, b, (s_ab & s_bc).astype(np.int64))
    np.add.at(cof, c, (s_ac & s_bc).astype(np.int64))
    return cof


def cofacial_six_pairs_by_label(F):
    """{original vertex label: cofacial degree-6 pairs}, for callers that never
    relabel (DefectState & friends key everything by the .mfd's own labels)."""
    F = np.asarray(F, np.int64)
    eu, edeg, V = _edges(F)
    return dict(zip(np.unique(F).tolist(), cofacial_six_pairs(F, eu, edeg, V).tolist()))


def link_class_census(F, eu, edeg, V):
    """Per-vertex FK link classification.

    Returns (labels, info) where `labels` is an int array of LinkClass codes and
    `info` carries the raw per-vertex arrays.  Class codes:

        0 Z12, 1 Z14, 2 Z15, 3 Z16 (T_d, Friauf)   -- Frank-Kasper legal
        4 Z16_D2   pure {5,6}, n_6 = 4, six-edges cofacial  -- NOT FK
        5 Z17plus  pure {5,6}, n_6 >= 5                     -- NOT FK (hub)
        6 impure   some incident edge of degree <=4 or >=7  -- illegal edges
    """
    def incident(mask):
        c = np.zeros(V, dtype=np.int64)
        np.add.at(c, eu[mask, 0], 1)
        np.add.at(c, eu[mask, 1], 1)
        return c

    n6 = incident(edeg == 6)
    pure56 = (incident(edeg <= 4) == 0) & (incident(edeg >= 7) == 0)
    cof = cofacial_six_pairs(F, eu, edeg, V)

    labels = np.full(V, 6, dtype=np.int8)          # impure
    p = pure56
    labels[p & (n6 == 0)] = 0
    labels[p & (n6 == 2)] = 1
    labels[p & (n6 == 3)] = 2
    labels[p & (n6 == 4) & (cof == 0)] = 3
    labels[p & (n6 == 4) & (cof > 0)] = 4
    labels[p & (n6 >= 5)] = 5

    # Theorems, cheap to assert: n_6 = 1 is impossible in the pure sector, and
    # the unique links at n_6 = 0, 2, 3 have no cofacial six-pairs.
    assert not np.any(p & (n6 == 1)), "pure-{5,6} vertex with n_6 = 1 (impossible)"
    assert not np.any(p & (n6 <= 3) & (cof > 0)), \
        "cofacial six-pair at n_6 <= 3 -- detector or complex is inconsistent"
    assert np.all(cof[p & (n6 == 4)] <= 2), \
        "n_6 = 4 link with >2 cofacial six-pairs -- no such class exists"

    info = dict(n6=n6, pure56=pure56, cofacial=cof)
    return labels, info


CLASS_NAMES = ["Z12", "Z14", "Z15", "Z16", "Z16_D2", "Z17plus", "impure"]
FK_CODES = (0, 1, 2, 3)


def summarize(labels):
    """Fractions by class name + f_FK (historical, n_6-bucketed) and f_FK_strict."""
    V = len(labels)
    cnt = np.bincount(labels, minlength=7)
    out = {name: float(cnt[i]) / V for i, name in enumerate(CLASS_NAMES)}
    out["n_Z16_D2"] = int(cnt[4])
    out["f_FK_strict"] = float(cnt[list(FK_CODES)].sum()) / V
    out["f_FK"] = out["f_FK_strict"] + out["Z16_D2"]   # what n_6 bucketing gives
    return out


# --------------------------------------------------------------------------
# Link-level classification (slow, exact -- for cross-checks and spot audits)
# --------------------------------------------------------------------------
def extract_link(F, v):
    """Triangles of lk(v): the opposite face of every tet containing v."""
    F = np.asarray(F, np.int64)
    rows = F[np.any(F == v, axis=1)]
    return [tuple(sorted(int(x) for x in r if x != v)) for r in rows]


def link_degree_map(tris):
    """Vertex -> degree in a triangulated 2-sphere given by its triangle list."""
    adj = defaultdict(set)
    for a, b, c in tris:
        adj[a].update((b, c))
        adj[b].update((a, c))
        adj[c].update((a, b))
    return {v: len(s) for v, s in adj.items()}, adj


def classify_link_exact(tris):
    """(name, n_6, cofacial_pairs) for a single link given as a triangle list."""
    deg, adj = link_degree_map(tris)
    ds = sorted(set(deg.values()))
    if not set(ds) <= {5, 6}:
        return "impure", None, None
    six = {v for v, d in deg.items() if d == 6}
    n6 = len(six)
    cof = sum(1 for v in six for w in adj[v] if w in six and w > v)
    if n6 == 4:
        return ("Z16" if cof == 0 else "Z16_D2"), n6, cof
    if n6 in FK_NAMES:
        return FK_NAMES[n6], n6, cof
    return f"Z{12 + n6}", n6, cof


# --------------------------------------------------------------------------
# Exhaustive enumeration of 5/6-triangulated 2-spheres
# --------------------------------------------------------------------------
class SphereEnum:
    """Backtracking generator for triangulated S^2 with all degrees in {5,6}.

    Grows a triangle set, maintaining for each vertex the link graph (adjacency
    among its neighbours); that graph must stay a disjoint union of paths and
    close into a single cycle exactly when the vertex is finished.  Since every
    degree is 5 or 6 with n_5 = 12 forced, chi = 2 automatically -- a connected
    closed result is a sphere, so no genus filtering is needed.
    """

    def __init__(self, k):
        self.k = k
        self.V = 12 + k
        self.reset()
        self.solutions = []
        self.nodes = 0

    def reset(self):
        self.tris = []
        self.nv = 0
        self.nbrs, self.link, self.closed = [], [], []
        self.ecount = defaultdict(int)
        self.n6 = 0

    def new_vertex(self):
        self.nbrs.append(set())
        self.link.append({})
        self.closed.append(False)
        self.nv += 1
        return self.nv - 1

    def pop_vertex(self):
        self.nbrs.pop()
        self.link.pop()
        self.closed.pop()
        self.nv -= 1

    def _component(self, v, x):
        lk, seen, stack = self.link[v], {x}, [x]
        while stack:
            for w in lk[stack.pop()]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    def can_add(self, a, b, c):
        if len({a, b, c}) != 3:
            return False
        for x, y in ((a, b), (b, c), (a, c)):
            if self.ecount[frozenset((x, y))] >= 2:
                return False
        for v in (a, b, c):
            if self.closed[v]:
                return False
        for v, (x, y) in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
            nb, lk = self.nbrs[v], self.link[v]
            if len(nb) + (x not in nb) + (y not in nb) > 6:
                return False
            if len(lk.get(x, ())) >= 2 or len(lk.get(y, ())) >= 2:
                return False
            if y in lk.get(x, ()):
                return False
            if x in lk and y in lk and y in self._component(v, x):
                # closing the link cycle: legal only if it uses every neighbour
                d = len(nb)
                if d not in (5, 6) or (d == 6 and self.n6 >= self.k):
                    return False
                if len(self._component(v, x)) != d:
                    return False
        return True

    def add(self, a, b, c):
        undo = []
        for x, y in ((a, b), (b, c), (a, c)):
            self.ecount[frozenset((x, y))] += 1
        for v, (x, y) in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
            nb, lk = self.nbrs[v], self.link[v]
            for u in (x, y):
                if u not in nb:
                    nb.add(u)
                    lk[u] = set()
                    undo.append(('n', v, u))
            lk[x].add(y)
            lk[y].add(x)
            if all(len(s) == 2 for s in lk.values()) and len(self._component(v, x)) == len(nb):
                self.closed[v] = True
                undo.append(('c', v))
                if len(nb) == 6:
                    self.n6 += 1
                    undo.append(('6',))
        self.tris.append((a, b, c))
        return undo

    def remove(self, a, b, c, undo):
        self.tris.pop()
        for op in reversed(undo):
            if op[0] == 'n':
                self.nbrs[op[1]].discard(op[2])
                self.link[op[1]].pop(op[2], None)
            elif op[0] == 'c':
                self.closed[op[1]] = False
            else:
                self.n6 -= 1
        for v, (x, y) in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
            lk = self.link[v]
            if x in lk:
                lk[x].discard(y)
            if y in lk:
                lk[y].discard(x)
        for x, y in ((a, b), (b, c), (a, c)):
            self.ecount[frozenset((x, y))] -= 1

    def candidates(self, a, b):
        out = [c for c in range(self.nv)
               if c not in (a, b) and not self.closed[c] and self.can_add(a, b, c)]
        if self.nv < self.V:
            out.append(-1)
        return out

    def search(self):
        self.nodes += 1
        open_e = [e for e, n in self.ecount.items() if n == 1]
        if not open_e:
            if self.nv == self.V and all(self.closed) and self.n6 == self.k:
                self.solutions.append([tuple(sorted(t)) for t in self.tris])
            return
        best, bestc = None, None
        for e in open_e:                       # most-constrained open edge
            a, b = tuple(e)
            cs = self.candidates(a, b)
            if not cs:
                return
            if bestc is None or len(cs) < len(bestc):
                best, bestc = (a, b), cs
                if len(cs) == 1:
                    break
        a, b = best
        for c in bestc:
            fresh = c == -1
            if fresh:
                c = self.new_vertex()
                if not self.can_add(a, b, c):
                    self.pop_vertex()
                    continue
            undo = self.add(a, b, c)
            self.search()
            self.remove(a, b, c, undo)
            if fresh:
                self.pop_vertex()

    def run(self):
        # symmetry breaking: vertex 0 carries its full star, degree d0 in {5,6}
        for d0 in (5, 6):
            if d0 == 6 and self.k == 0:
                continue
            self.reset()
            v0 = self.new_vertex()
            ring = [self.new_vertex() for _ in range(d0)]
            undos = [self.add(v0, ring[i], ring[(i + 1) % d0]) for i in range(d0)]
            self.search()
            for i in reversed(range(d0)):
                self.remove(v0, ring[i], ring[(i + 1) % d0], undos[i])
        return self.solutions


def enumerate_5_6_spheres(k, dedup=True):
    """All triangulated S^2 with degrees in {5,6} and exactly k degree-6 vertices.

    Returns a list of triangle lists, one per isomorphism class when dedup.
    """
    e = SphereEnum(k)
    sols = e.run()
    if not dedup:
        return sols
    import networkx as nx
    reps, graphs = [], []
    for s in sols:
        G = nx.Graph()
        for a, b, c in s:
            G.add_edges_from([(a, b), (b, c), (a, c)])
        if not any(nx.is_isomorphic(G, H) for H in graphs):
            graphs.append(G)
            reps.append(s)
    return reps


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def _load_facets(path):
    sys.path.insert(0, os.path.join(_ROOT, "python"))
    from discrete_differential_geometry import Manifold
    return np.asarray(Manifold.load(path, 3).facets(), np.int64)


def _edges(F):
    sys.path.insert(0, os.path.join(_ROOT, "scripts"))
    from fk_skeleton import edges_from_facets
    return edges_from_facets(F)


def cmd_enumerate(kmax):
    import networkx as nx
    print(f"{'n_6':>4} {'V':>4} {'classes':>8}   |Aut| / cofacial six-pairs")
    for k in range(kmax + 1):
        reps = enumerate_5_6_spheres(k)
        detail = []
        for s in reps:
            G = nx.Graph()
            for a, b, c in s:
                G.add_edges_from([(a, b), (b, c), (a, c)])
            aut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
            _, _, cof = classify_link_exact(s)
            detail.append(f"{aut}/{cof}")
        exp = CLASS_COUNTS.get(k)
        flag = "" if exp is None or exp == len(reps) else f"  !! expected {exp}"
        print(f"{k:>4} {12+k:>4} {len(reps):>8}   {', '.join(detail) or '--'}{flag}")
        sys.stdout.flush()


def cmd_census(paths, verify=False, sample=200):
    for path in paths:
        F = _load_facets(path)
        eu, edeg, V = _edges(F)
        labels, info = link_class_census(F, eu, edeg, V)
        s = summarize(labels)
        print(f"\n{os.path.basename(path)}   N = {V}")
        print("  " + "  ".join(f"{n}={s[n]:.4f}" for n in CLASS_NAMES))
        print(f"  f_FK (n_6 bucketed) = {s['f_FK']:.6f}"
              f"   f_FK_strict = {s['f_FK_strict']:.6f}"
              f"   Z16_D2 sites = {s['n_Z16_D2']}")
        if s["n_Z16_D2"]:
            bad = np.flatnonzero(labels == 4)
            print(f"  !! non-FK Z16 (D2 isomer) at vertices {bad[:10].tolist()}"
                  f"{' ...' if len(bad) > 10 else ''}")
        if verify:
            rng = np.random.default_rng(0)
            vs = rng.choice(V, size=min(sample, V), replace=False)
            lab_to_name = dict(enumerate(CLASS_NAMES))
            bad = 0
            for v in vs:
                name, _, cof = classify_link_exact(extract_link(F, int(v)))
                want = lab_to_name[int(labels[v])]
                got = name if name != "impure" else "impure"
                if got.startswith("Z1") and got not in CLASS_NAMES:
                    got = "Z17plus"
                if got != want or (cof is not None and cof != info["cofacial"][v]):
                    bad += 1
                    print(f"  MISMATCH v={v}: fast={want} slow={name} "
                          f"cof fast={info['cofacial'][v]} slow={cof}")
            print(f"  verify: {len(vs)} links re-derived from scratch, "
                  f"{bad} mismatches")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*", help=".mfd files to classify")
    ap.add_argument("--enumerate", type=int, metavar="KMAX",
                    help="regenerate the class table for n_6 = 0..KMAX")
    ap.add_argument("--verify", action="store_true",
                    help="cross-check the fast detector against per-vertex link extraction")
    ap.add_argument("--sample", type=int, default=200,
                    help="vertices to re-derive under --verify (default 200)")
    a = ap.parse_args()
    if a.enumerate is not None:
        cmd_enumerate(a.enumerate)
    if a.paths:
        cmd_census(a.paths, verify=a.verify, sample=a.sample)
    if a.enumerate is None and not a.paths:
        ap.error("give .mfd paths and/or --enumerate KMAX")
