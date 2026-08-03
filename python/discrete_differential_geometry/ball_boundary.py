"""Exact boundary distance map of the elementary 2->3 ball.

A 2->3 Pachner move replaces the two tets of a bipyramid by three tets on
the SAME five vertices, so the two complexes are identical outside a ball
B and share ``dB`` -- six unit equilateral triangles forming a 2-sphere.
Every metric consequence of the defect is therefore carried by one
compactly supported object, the pair of intrinsic boundary distance maps

    d_B^R, d_B^D : dB x dB -> R,

because for any x, y in the host

    d(x, y) = min( d_out(x, y),
                   min_{p,q in dB} [ d_out(x,p) + d_B(p,q) + d_out(q,y) ] )

with d_out (distance in the closed complement of B) IDENTICAL before and
after. Everything except d_B is common mode and cancels in d_D - d_R.

The map is UNIVERSAL: all tets are congruent regular unit simplices, so
it depends only on the combinatorics of the move, not on the host, the
vertex labels, or where in the crystal the move is applied. Compute it
once and it applies to every 2->3 in every triangulation.

Labels (fixed): A, B, C = the shared face; P, Q = the two apexes.

    reference R : tets ABCP, ABCQ            internal face ABC
    defected  D : tets ABPQ, ACPQ, BCPQ      internal edge PQ (degree 3)
    boundary dB : ABP ACP BCP ABQ ACQ BCQ    (identical in both)

Exactness
---------
Distances are computed by rigid development in the edge-sqrt(2) rational
embedding of ``development.py``, so every developed coordinate is in Q^3
and every squared length is in Q. A candidate path is a straight chord
through a strip of tets; it is ACCEPTED only if it crosses each internal
face of the strip in the interior, in order -- a sequence of exact
rational sign tests. Consecutive crossings then bound a sub-segment whose
endpoints both lie in one (convex) tet, so containment in the strip is
proved, not sampled. Lengths are returned as exact ``Fraction`` squares
in unit-edge units; take a square root only to report.

The minimum is over all simple tet strips, which for these two tiny
configurations is a complete enumeration -- so d_B is exact, not bounded,
provided no minimizer touches the singular edge PQ. That last case is
detected (``missing`` in the returned summary), never assumed away; it
does not occur, for the reason that PQ carries cone angle 3*theta =
211.59 deg < 2*pi and is therefore positively curved.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations

from .development import _SEED, develop_path, _sub, _dot, _cross

def _perms3(xs):
    from itertools import permutations
    return tuple(permutations(xs))


A, B, C, P, Q = 0, 1, 2, 3, 4

TETS_R = ((A, B, C, P), (A, B, C, Q))
TETS_D = ((A, B, P, Q), (A, C, P, Q), (B, C, P, Q))
BOUNDARY = ((A, B, P), (A, C, P), (B, C, P),
            (A, B, Q), (A, C, Q), (B, C, Q))
VERTEX_NAMES = {A: "A", B: "B", C: "C", P: "P", Q: "Q"}

#: the configuration's own symmetry group: S3 on the shared face (A,B,C)
#: times S2 on the apexes (P,Q), order 12. Every element preserves TETS_R,
#: TETS_D and BOUNDARY setwise, so the intrinsic distance maps must be
#: EQUIVARIANT under it: d(gx, gy) == d(x, y) exactly, for both configs. That
#: is a free and rather sharp correctness check on the whole development --
#: an error in a strip placement, a chord test, or the node grid would almost
#: certainly break it -- and `defect_boundary_map.validate` asserts it.
#:
#: These are the same 12 relabellings `move_site_census` canonicalises a 2->3
#: site over. Note the distinction: there the group acts on a site embedded in
#: a host, where it is only part of the story (the host's Aut decides which
#: sites are equivalent); here the ball is the whole world, so the group is
#: the complete symmetry of the object.
CONFIG_SYMMETRY = tuple(
    {A: fp[0], B: fp[1], C: fp[2], P: ap[0], Q: ap[1]}
    for fp in _perms3((A, B, C)) for ap in ((P, Q), (Q, P)))

#: squared edge length of the development embedding (unit edge -> divide by 2)
EMBED_EDGE_SQ = Fraction(2)


def face_name(tri):
    return "".join(VERTEX_NAMES[v] for v in tri)


# --------------------------------------------------------------------------
# boundary nodes: a shared Steiner grid on dB
# --------------------------------------------------------------------------

def boundary_nodes(order: int):
    """Barycentric grid of order `order` on the six boundary triangles.

    Returns a list of nodes, each a dict {vertex: Fraction weight} summing
    to 1 with only NONZERO weights kept -- so a point on a shared edge or
    a vertex of dB is produced once, with the same key, no matter which
    triangle generated it. That sharing is what makes the map splice into
    a global Steiner graph built at the same order.
    """
    seen = {}
    for tri in BOUNDARY:
        for i in range(order + 1):
            for j in range(order + 1 - i):
                k = order - i - j
                w = {v: Fraction(c, order)
                     for v, c in zip(tri, (i, j, k)) if c}
                seen.setdefault(node_key(w), w)
    return [seen[k] for k in sorted(seen)]


def node_key(w):
    return tuple(sorted((v, wt) for v, wt in w.items()))


def node_name(w):
    return "+".join(f"{wt}{VERTEX_NAMES[v]}" for v, wt in sorted(w.items()))


def node_faces(w):
    """The boundary triangles whose closure contains this node."""
    s = set(w)
    return tuple(t for t in BOUNDARY if s <= set(t))


# --------------------------------------------------------------------------
# tet strips and their exact developments
# --------------------------------------------------------------------------

def strips(tets):
    """All simple dual-graph paths (tet strips) in a configuration."""
    adj = {t: [u for u in tets if u != t and len(set(t) & set(u)) == 3]
           for t in tets}
    out = []

    def walk(path):
        out.append(tuple(path))
        for u in adj[path[-1]]:
            if u not in path:
                walk(path + [u])

    for t in tets:
        walk([t])
    return out


def strip_placements(strip):
    """Developed positions, one dict per tet of the strip, in a common
    exact rational frame. Vertex labels revisited by the strip get one
    placement PER TET (development.develop_path keeps them separate), so
    a cone that fails to close cannot silently self-identify."""
    seed = {v: _SEED[i] for i, v in enumerate(strip[0])}
    return develop_path(strip[0], seed, list(strip))


def _pos(w, place):
    return tuple(sum(wt * place[v][k] for v, wt in w.items())
                 for k in range(3))


def _in_triangle(X, q1, q2, q3):
    """Exact: is coplanar X inside triangle (q1,q2,q3), boundary included?"""
    n = _cross(_sub(q2, q1), _sub(q3, q1))
    a1 = _dot(_cross(_sub(q2, X), _sub(q3, X)), n)
    a2 = _dot(_cross(_sub(q3, X), _sub(q1, X)), n)
    a3 = _dot(_cross(_sub(q1, X), _sub(q2, X)), n)
    return a1 >= 0 and a2 >= 0 and a3 >= 0


def strip_length_sq(strip, places, wx, wy):
    """Exact squared length (unit-edge units) of the straight chord from
    node wx in strip[0] to node wy in strip[-1], or None if that chord is
    not contained in the strip.

    Containment proof: the chord crosses each internal face in its
    interior at parameters t_0 <= t_1 <= ... , and each sub-segment then
    has both endpoints in the closure of one convex tet.
    """
    if not set(wx) <= set(strip[0]) or not set(wy) <= set(strip[-1]):
        return None
    p0 = _pos(wx, places[0])
    p1 = _pos(wy, places[-1])
    d = _sub(p1, p0)
    if d == (0, 0, 0):
        return Fraction(0)
    last = Fraction(0)
    for i in range(len(strip) - 1):
        shared = [v for v in strip[i] if v in strip[i + 1]]
        q1, q2, q3 = (places[i][v] for v in shared)
        n = _cross(_sub(q2, q1), _sub(q3, q1))
        den = _dot(d, n)
        if den == 0:                      # chord parallel to the face
            return None
        t = _dot(_sub(q1, p0), n) / den
        if t < last or t > 1:
            return None
        X = tuple(p0[k] + t * d[k] for k in range(3))
        if not _in_triangle(X, q1, q2, q3):
            return None
        last = t
    return _dot(d, d) / EMBED_EDGE_SQ


class BallBoundaryMap:
    """Exact intrinsic distance maps on dB for one configuration.

    ``config`` is "R" (two tets) or "D" (three tets). ``dist_sq(x, y)``
    returns an exact Fraction; ``matrix()`` returns the full float table
    on the order-n boundary grid.
    """

    def __init__(self, config: str, order: int):
        if config not in ("R", "D"):
            raise ValueError("config must be 'R' or 'D'")
        self.config = config
        self.order = order
        self.tets = TETS_R if config == "R" else TETS_D
        self.strips = strips(self.tets)
        self._places = {s: strip_placements(s) for s in self.strips}
        self.nodes = boundary_nodes(order)
        self.keys = [node_key(w) for w in self.nodes]
        self._index = {k: i for i, k in enumerate(self.keys)}

    def index(self, w):
        return self._index[node_key(w)]

    def node_permutation(self, g):
        """Index permutation induced on the boundary grid by a relabelling.

        `g` maps vertex -> vertex (an element of CONFIG_SYMMETRY). The grid is
        defined by barycentric weights, so relabelling a node's weights lands
        on another grid node of the same order -- the permutation is total.
        """
        return [self.index({g[v]: w for v, w in nd.items()})
                for nd in self.nodes]

    def dist_sq(self, wx, wy):
        """Exact squared intrinsic distance through B, or None if every
        strip chord failed (a minimizer touching the singular edge)."""
        best = None
        for s in self.strips:
            L = strip_length_sq(s, self._places[s], wx, wy)
            if L is not None and (best is None or L < best):
                best = L
        return best

    def matrix(self):
        """(float matrix, list of unresolved index pairs)."""
        import numpy as np
        n = len(self.nodes)
        M = np.zeros((n, n))
        missing = []
        for i in range(n):
            for j in range(i + 1, n):
                L = self.dist_sq(self.nodes[i], self.nodes[j])
                if L is None:
                    missing.append((i, j))
                    M[i, j] = M[j, i] = np.nan
                else:
                    M[i, j] = M[j, i] = float(L) ** 0.5
        return M, missing


def steiner_interior(config, order, grid=None):
    """INDEPENDENT upper bound on d_B: a Steiner graph over the interior
    of B at subdivision `grid` (default = `order`).

    Nodes are barycentric grid points on ALL faces of the configuration,
    internal ones included, keyed by weight dict so points shared between
    tets are identified; every pair of nodes on the same tet is joined by
    the exact straight chord in that tet's regular embedding. Any path it
    represents is a real path in B, so its distances are >= d_B -- and
    unlike the strip enumeration it can also represent paths that touch
    the singular edge PQ or hug dB. Agreement at the boundary nodes is
    therefore evidence that no such path beats the strip minimum.

    Returns (matrix over boundary_nodes(order), node list).
    """
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
    grid = grid or order
    tets = TETS_R if config == "R" else TETS_D
    keys, pos_of = {}, {}
    wmin = {}                    # (lo, hi) -> min weight; coo SUMS duplicates
    for tet in tets:
        place = {v: _SEED[k] for k, v in enumerate(tet)}
        local = []
        for tri in combinations(tet, 3):
            for i in range(grid + 1):
                for j in range(grid + 1 - i):
                    k = grid - i - j
                    w = {v: Fraction(c, grid)
                         for v, c in zip(tri, (i, j, k)) if c}
                    kk = node_key(w)
                    if kk not in keys:
                        keys[kk] = len(keys)
                    pos_of[kk] = _pos(w, place)
                    local.append(kk)
        local = sorted(set(local))
        for ka, kb in combinations(local, 2):
            d = _sub(pos_of[ka], pos_of[kb])
            w = float(_dot(d, d) / EMBED_EDGE_SQ) ** 0.5
            e = (keys[ka], keys[kb]) if keys[ka] < keys[kb] else (keys[kb], keys[ka])
            if e not in wmin or w < wmin[e]:
                wmin[e] = w
    n = len(keys)
    ii = [e[0] for e in wmin]; jj = [e[1] for e in wmin]
    ww = [wmin[e] for e in wmin]
    G = coo_matrix((ww + ww, (ii + jj, jj + ii)), shape=(n, n)).tocsr()
    Dm = shortest_path(G, directed=False)
    nodes = boundary_nodes(order)
    sel = [keys[node_key(w)] for w in nodes]
    return Dm[np.ix_(sel, sel)], nodes


def boundary_surface_lengths(order):
    """Distances measured ALONG dB (the six-triangle sphere), as a graph
    upper bound: a Steiner graph on the boundary surface only. Used to
    check that the solid genuinely shortcuts its own boundary."""
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
    nodes = boundary_nodes(order)
    idx = {node_key(w): i for i, w in enumerate(nodes)}
    seedplace = {}
    ii, jj, ww = [], [], []
    for tri in BOUNDARY:
        pos = {v: _SEED[k] for k, v in enumerate(tri)}     # flat triangle
        seedplace[tri] = pos
        loc = [w for w in nodes if set(w) <= set(tri)]
        for wx, wy in combinations(loc, 2):
            a, b = _pos(wx, pos), _pos(wy, pos)
            d = _sub(a, b)
            ii.append(idx[node_key(wx)])
            jj.append(idx[node_key(wy)])
            ww.append(float(_dot(d, d) / EMBED_EDGE_SQ) ** 0.5)
    n = len(nodes)
    G = coo_matrix((ww + ww, (ii + jj, jj + ii)), shape=(n, n)).tocsr()
    return shortest_path(G, directed=False), nodes
