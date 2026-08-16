"""Frank-Kasper / disclination-skeleton census primitives.

Exact combinatorial backbone: for every vertex v the link is a triangulated
S^2, and the degree of edge {v,u} equals the degree of u in link(v); Euler for
the link gives the LOCAL SUM RULE  sum_{e ni v} (6 - deg e) = 12  exactly.
Hence a vertex whose incident edges are all degree 5 has exactly 12 of them
(icosahedral link, Z12), and pure-{5,6} vertices with two/three/four 6-edges
are the FK coordinations Z14/Z15/Z16. A pure-{5,6} vertex with exactly ONE
6-edge cannot exist (no S^2 triangulation with degree census 5^12 6^1), so
"disclination lines cannot end" is a theorem in the pure-{5,6} sector; line
endpoints require a compensating deg<=4 or deg>=7 edge at the same vertex.

Promoted from scripts/fk_skeleton.py (2026-08 cleanup, Phase 1a); the
g-ladder ensemble census CLI stays there. See `link_classes` for the exact
5/6-link classification (Z16 T_d vs D2) these censuses lean on.
"""
import numpy as np

EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def edges_from_facets(F):
    """(eu, edeg, V) from a facet array: unique edges (relabeled 0..V-1), degrees."""
    F = np.asarray(F, np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    T = inv.reshape(F.shape)
    V = len(lab)
    epairs = np.sort(np.vstack([T[:, [i, j]] for i, j in EDGE_PAIRS]), axis=1)
    eu, edeg = np.unique(epairs, axis=0, return_counts=True)
    return eu, edeg, V


def load_edges(path):
    """Load a seed; return (eu, edeg, V) as in edges_from_facets."""
    from ._manifold import Manifold
    return edges_from_facets(Manifold.load(path, 3).facets())


def vertex_classes(facets):
    """Per-vertex (n6, m_impure) and adjacency list.

    (Promoted from scripts/dopant_pairs.py, which re-exports it.)
    """
    eu, edeg, V = edges_from_facets(facets)
    n6 = np.zeros(V, int)
    imp = np.zeros(V, int)
    adj = [[] for _ in range(V)]
    for (a, b), d in zip(eu, edeg):
        adj[a].append(b)
        adj[b].append(a)
        if d >= 6:
            n6[a] += 1
            n6[b] += 1
        if d < 5 or d > 6:
            imp[a] += 1
            imp[b] += 1
    return n6, imp, adj


def vertex_class_census(eu, edeg, V, facets=None):
    """Per-vertex counts of incident edges by degree class + Z-classification.

    n_6 alone does NOT determine lk(v): at n_6 = 4 there are two 5/6-spheres,
    the T_d Friauf polyhedron (the FK Z16) and a D2 isomer no FK phase admits
    (see `link_classes`).  Pass `facets` to split them -- the discriminator
    is whether any two degree-6 edges at v are cofacial.  Keys "Z12".."Z16"
    keep their historical n_6-bucketed meaning either way, so fFK is unchanged;
    with `facets` the census additionally reports "Z16_Td", "Z16_D2" and
    "FK_strict" = fFK - fZ16_D2, which is the honest FK fraction.
    """
    def incident(mask):
        c = np.zeros(V, dtype=np.int64)
        np.add.at(c, eu[mask, 0], 1)
        np.add.at(c, eu[mask, 1], 1)
        return c

    n_le4 = incident(edeg <= 4)
    n5 = incident(edeg == 5)
    n6 = incident(edeg == 6)
    n_ge7 = incident(edeg >= 7)

    # Exact local sum rule (vertex link = S^2): sum(6 - deg) over incident = 12.
    charge = np.zeros(V, dtype=np.int64)
    np.add.at(charge, eu[:, 0], 6 - edeg)
    np.add.at(charge, eu[:, 1], 6 - edeg)
    assert np.all(charge == 12), "link sum rule violated -- not a 3-manifold?"

    pure56 = (n_le4 == 0) & (n_ge7 == 0)
    fz = {}
    for name, k in (("Z12", 0), ("Z14", 2), ("Z15", 3), ("Z16", 4)):
        fz[name] = float(np.mean(pure56 & (n6 == k)))
    n_broken = int(np.sum(pure56 & (n6 == 1)))  # forbidden by S^2 combinatorics
    fz["pure56_other"] = float(np.mean(pure56)) - sum(fz.values())
    fz["impure"] = float(np.mean(~pure56))

    if facets is not None:
        from .link_classes import cofacial_six_pairs
        cof = cofacial_six_pairs(facets, eu, edeg, V)
        assert not np.any(pure56 & (n6 <= 3) & (cof > 0)), \
            "cofacial six-pair at n_6 <= 3 -- no such 5/6-sphere exists"
        z16 = pure56 & (n6 == 4)
        fz["Z16_D2"] = float(np.mean(z16 & (cof > 0)))
        fz["Z16_Td"] = fz["Z16"] - fz["Z16_D2"]
        # keyed "FK_strict" so the callers' f-prefix convention yields fFK_strict
        fz["FK_strict"] = fz["Z12"] + fz["Z14"] + fz["Z15"] + fz["Z16_Td"]
    return fz, n_broken


def skeleton_stats(eu, edeg, V):
    """Structure of the deg>=6 disclination network."""
    sk = eu[edeg >= 6]
    if len(sk) == 0:
        return dict(n_edges=0, frac_val1=np.nan, frac_val2=np.nan,
                    n_comp=0, largest_frac=np.nan, cyclomatic=0)
    val = np.zeros(V, dtype=np.int64)
    np.add.at(val, sk[:, 0], 1)
    np.add.at(val, sk[:, 1], 1)
    touched = val > 0
    vh = np.bincount(val[touched])
    ntouch = int(touched.sum())

    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in sk:
        a, b = int(a), int(b)
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = {}
    for a, b in sk:
        r = find(int(a))
        roots[r] = roots.get(r, 0) + 1
    comp_sizes = sorted(roots.values(), reverse=True)
    return dict(
        n_edges=int(len(sk)),
        frac_val1=float(vh[1] / ntouch) if len(vh) > 1 else 0.0,
        frac_val2=float(vh[2] / ntouch) if len(vh) > 2 else 0.0,
        n_comp=len(comp_sizes),
        largest_frac=float(comp_sizes[0] / len(sk)),
        cyclomatic=int(len(sk) - ntouch + len(comp_sizes)),
    )


def bfs_orders(eu, V, n_centers, rng):
    adj = [[] for _ in range(V)]
    for a, b in eu:
        adj[a].append(b)
        adj[b].append(a)
    mmax = V // 6
    orders = []
    for _ in range(n_centers):
        src = int(rng.integers(V))
        seen = np.zeros(V, bool)
        seen[src] = True
        order = [src]
        frontier = [src]
        while len(order) < mmax and frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if not seen[w]:
                        seen[w] = True
                        nxt.append(w)
            rng.shuffle(nxt)
            order.extend(nxt)
            frontier = nxt
        orders.append(np.array(order[:mmax]))
    return orders, mmax


def window_ratio(q, orders, mgrid, rng, n_shuf=4):
    """Var over centers of window charge, observed / value-shuffled null."""
    def var(field):
        sums = np.array([np.cumsum(field[o])[mgrid - 1] for o in orders])
        return sums.var(axis=0)

    vr = var(q)
    vs = np.mean([var(rng.permutation(q)) for _ in range(n_shuf)], axis=0)
    return vr / vs
