"""Exact combinatorial symmetry of a 3-manifold triangulation.

Computes ``Aut(K)`` -- the group of vertex relabellings carrying the facet set
to itself -- as explicit permutations, together with the orbit maps that most
of this codebase has been approximating with brute force plus invariant
bucketing (Weisfeiler-Leman colours, degree signatures, rounded coordinates).
Those approximations are one-sided: a colour invariant can only ever prove two
objects INEQUIVALENT, so a class count from bucketing is a bound, not the
answer.  The orbits here are exact.

The algorithm
-------------
An automorphism of a CONNECTED 3-manifold triangulation is determined by the
image of a single ordered tetrahedron (a "frame").  Crossing a shared face
forces the image of the new apex: the three face vertices keep their images,
and the neighbouring tet across the image face is unique.  So::

    Aut(K) = { developments of one base frame that close globally }

which is the same forced-propagation argument that ``crystal_grains.py`` uses
to grow crystalline grains -- ``develop_partial`` (grain growth, mismatch =
boundary) and ``develop_total`` (automorphism, mismatch = failure) are the
same traversal under two mismatch policies, and both live here.

Two consequences make this cheap.  First, ``Aut`` acts FREELY on frames (an
automorphism fixing a frame is forced to be the identity), so every frame
orbit has size exactly ``|Aut|`` and ``24 * nT / |Aut|`` is an integer -- a
free self-check.  Second, once some automorphisms are known, every frame in
the base frame's orbit is accounted for and needs no development at all, so
the number of full developments is (generators found) + (failed candidates),
not ``|Aut|``.  The R crystal at m=3 (24462 tets) takes on the order of a
minute; results are cached to a ``.sym.npz`` sidecar.

WHY Aut(K) MAY LEGITIMATELY EXCEED THE PHYSICAL CRYSTAL'S SYMMETRY
------------------------------------------------------------------
Read this before quoting ``|Aut|`` as a space-group order.  ``K`` is an
abstract simplicial complex; the crystal is a decorated point set in R^3.  The
group computed here is guaranteed only to CONTAIN the image of the crystal's
space group acting on the supercell torus, and it can be strictly larger for
three independent reasons:

  1. METRIC BLINDNESS.  ``K`` records which sites are Delaunay-adjacent, not
     where they sit.  A relabelling that preserves every adjacency while
     distorting distances is an automorphism of ``K`` and not an isometry of
     the crystal.  Nothing rules one out.
  2. SUPERCELL PERIODICITY.  We triangulate the m x m x m torus, not R^3.  The
     quotient turns the lattice translations into a FINITE group and can
     manufacture coincidences the infinite crystal does not have (e.g. a BC
     chain that closes only by wrap-around).  Anything derived from |Aut| is
     a property of the torus at that m, not of the crystal.
  3. DELAUNAY TIE-BREAKING.  ``tcp_reference.py`` perturbs every site by ~1e-6
     before tiling, to resolve co-spherical degeneracies consistently across
     periodic images.  That perturbation breaks the point group METRICALLY,
     while the connectivity it produces may retain it -- so ``K`` can be more
     symmetric than the point set it was built from.

For sampler physics ``Aut(K)`` is nevertheless the RIGHT group: the sampler
sees ``K``, so ``Aut(K)`` is exactly the set of relabellings under which the
action, every move rate, and every defect species is invariant.  Just do not
identify it with crystallography without checking the factorisation.

Measured 2026-08-02, for what it is worth: A15 (m=3), C15 (m=3) and R (m=2,3)
all come out EQUAL to (supercell translations) x (point group x centering),
with no accidental automorphisms --

    A15 m3  |Aut| = 1296 = 27 x 48            (Pm-3n)
    C15 m3  |Aut| = 5184 = 27 x 48 x 4        (Fd-3m, fcc centering)
    R   m2  |Aut| =  144 =  8 x  6 x 3        (R-3, rhombohedral centering)
    R   m3  |Aut| =  486 = 27 x  6 x 3

so for these crystals the two groups coincide.  That is an empirical result
at these supercells, not a theorem.

Terminology
-----------
"Orbit" in ``scripts/defect_dynamics/`` almost always means a BC CHAIN (the
closed sliding-window walk of ``worm_helix.bc_orbit``), never a group orbit.
This module says "chain" for the walk and "orbit" only for group orbits;
``chain_orbits`` is the classification of chains under ``Aut``.
"""
from __future__ import annotations

import collections
import hashlib
import itertools
import os
import warnings

import numpy as np

__all__ = [
    "TriView", "CrystalSymmetry",
    "develop_partial", "develop_total",
    "PERMS", "PERM_INDEX", "E6",
]

#: the six edges of a tet, as local vertex-index pairs
E6 = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

#: the 24 orderings of a tet's vertices; a FRAME is (tet id, index into this)
PERMS = tuple(itertools.permutations(range(4)))
PERM_INDEX = {p: i for i, p in enumerate(PERMS)}


# ---------------------------------------------------------------------------
# combinatorial view
# ---------------------------------------------------------------------------


class TriView:
    """Adjacency tables of a 3-manifold triangulation, in internal ids.

    Vertices are relabelled to ``0..V-1`` (``labels[i]`` is the original id,
    ``index[orig]`` the inverse) unless ``relabel=False``, which asserts the
    input is already contiguous from 0.  Tet rows are NOT sorted -- only the
    keys of ``tetid`` and ``face2`` are -- so a caller may keep its own row
    order (``crystal_grains`` does).
    """

    def __init__(self, tets, relabel=True):
        T = np.asarray(tets, dtype=np.int64)
        if T.ndim != 2 or T.shape[1] != 4:
            raise ValueError(f"expected an (nT, 4) tet array, got {T.shape}")
        if relabel:
            labels, inv = np.unique(T, return_inverse=True)
            self.labels = labels.astype(np.int64)
            self.tets = inv.reshape(T.shape).astype(np.int32)
        else:
            self.labels = np.arange(int(T.max()) + 1, dtype=np.int64)
            self.tets = T.astype(np.int32)
        self.V = int(len(self.labels))
        self.nT = int(len(self.tets))
        self.index = {int(v): i for i, v in enumerate(self.labels)}

        face2 = {}
        for t, tv in enumerate(self.tets):
            tv = (int(tv[0]), int(tv[1]), int(tv[2]), int(tv[3]))
            for i in range(4):
                face2.setdefault(tuple(sorted(tv[:i] + tv[i + 1:])),
                                 []).append((t, tv[i]))
        self.face2 = face2
        self.tetid = {tuple(sorted(int(x) for x in tv)): t
                      for t, tv in enumerate(self.tets)}
        self._nbr = None
        self._edeg = None
        self._adj = None
        self._colours = None

    # -- lazily built tables -------------------------------------------------

    @property
    def nbr(self):
        """``(nbr, napex)``: across the face of tet ``t`` opposite local vertex
        ``i``, the neighbouring tet id and its apex vertex."""
        if self._nbr is None:
            nbr = np.full((self.nT, 4), -1, np.int32)
            nap = np.full((self.nT, 4), -1, np.int32)
            for t, tv in enumerate(self.tets):
                tv = (int(tv[0]), int(tv[1]), int(tv[2]), int(tv[3]))
                for i in range(4):
                    pair = self.face2[tuple(sorted(tv[:i] + tv[i + 1:]))]
                    if len(pair) != 2:
                        raise ValueError(
                            f"face of tet {t} has {len(pair)} cofacets; "
                            "this is not a closed 3-manifold")
                    o = pair[0] if pair[0][0] != t else pair[1]
                    nbr[t, i], nap[t, i] = o
            self._nbr = (nbr, nap)
        return self._nbr

    @property
    def edge_degree(self):
        """``{(u, w) sorted: number of tets on that edge}`` (the hinge degree)."""
        if self._edeg is None:
            d = collections.Counter()
            for tv in self.tets:
                tv = (int(tv[0]), int(tv[1]), int(tv[2]), int(tv[3]))
                for i, j in E6:
                    a, b = tv[i], tv[j]
                    d[(a, b) if a < b else (b, a)] += 1
            self._edeg = dict(d)
        return self._edeg

    @property
    def adj(self):
        if self._adj is None:
            a = collections.defaultdict(set)
            for (u, w) in self.edge_degree:
                a[u].add(w)
                a[w].add(u)
            self._adj = {v: sorted(s) for v, s in a.items()}
        return self._adj

    @property
    def edges(self):
        return sorted(self.edge_degree)

    @property
    def faces(self):
        return sorted(self.face2)

    @property
    def tet_keys(self):
        return sorted(self.tetid)

    def colours(self, rounds=20):
        """1-WL refinement of the edge-degree-labelled 1-skeleton.

        Aut-invariant by construction, so safe as a lossless PREFILTER on
        candidate frames -- it can only rule candidates out.  (It is exactly
        the invariant ``move_site_census`` used to bucket with; here it merely
        accelerates the exact computation.)
        """
        if self._colours is not None:
            return self._colours
        deg, adj = self.edge_degree, self.adj
        col = {v: len(ns) for v, ns in adj.items()}
        for _ in range(rounds):
            new = {}
            for v, ns in adj.items():
                new[v] = (col[v], tuple(sorted(
                    (col[u], deg[(v, u) if v < u else (u, v)]) for u in ns)))
            pal = {c: i for i, c in enumerate(sorted(set(new.values())))}
            nxt = {v: pal[new[v]] for v in new}
            if len(pal) == len(set(col.values())) and nxt == col:
                break
            col = nxt
        self._colours = col
        return col

    def is_dual_connected(self):
        """Is the dual (tet-adjacency) graph connected?

        The whole method rests on "an automorphism is determined by the image
        of one frame", which is FALSE across components -- a disconnected
        complex can have automorphisms permuting components that no single
        development can see. Callers must refuse rather than report the group
        they happen to find.
        """
        nbr, _ = self.nbr
        seen = np.zeros(self.nT, bool)
        seen[0] = True
        stack = [0]
        n = 1
        while stack:
            t = stack.pop()
            for i in range(4):
                u = int(nbr[t, i])
                if not seen[u]:
                    seen[u] = True
                    n += 1
                    stack.append(u)
        return n == self.nT

    def frame_window(self, fid):
        """Frame id -> the ordered 4-tuple of internal vertex ids."""
        t, pi = divmod(int(fid), 24)
        tv = self.tets[t]
        return tuple(int(tv[k]) for k in PERMS[pi])

    def frame_id(self, window):
        """Ordered 4-tuple of internal vertex ids -> frame id."""
        w = tuple(int(x) for x in window)
        t = self.tetid[tuple(sorted(w))]
        tv = [int(x) for x in self.tets[t]]
        return t * 24 + PERM_INDEX[tuple(tv.index(x) for x in w)]

    def digest(self):
        """Content hash of the facet set (cache key)."""
        f = np.sort(self.tets, axis=1)
        f = f[np.lexsort(f.T[::-1])]
        return hashlib.sha1(np.ascontiguousarray(f, np.int32).tobytes()).hexdigest()


# ---------------------------------------------------------------------------
# development (the shared workhorse)
# ---------------------------------------------------------------------------


def develop_partial(src, dst, seed_tet, sig0, accept):
    """Grow a single-valued covering map ``src -> dst`` from one seed tet.

    ``sig0`` maps the four vertices of ``src`` tet ``seed_tet`` to vertices of
    a ``dst`` tet.  Crossing a face FORCES the neighbour's image; the neighbour
    is accepted iff ``accept(src_tet_vertices, sig)`` holds, and a mismatch or
    a single-valuedness conflict simply makes that face a BOUNDARY (this is
    the grain-growth policy -- see ``develop_total`` for the strict one).

    Returns ``{src tet id: (dst tet id, sig)}`` over every tet consistently
    reached, where each ``sig`` is the 4-vertex correspondence of that tet.
    """
    seed_rt = tuple(sorted(int(sig0[int(v)]) for v in src.tets[seed_tet]))
    assign = {seed_tet: (dst.tetid[seed_rt], dict(sig0))}
    q = collections.deque([seed_tet])
    while q:
        t = q.popleft()
        rtid, sig = assign[t]
        tv = [int(v) for v in src.tets[t]]
        for apex in tv:
            face = [v for v in tv if v != apex]
            snb = [x for x in src.face2.get(tuple(sorted(face)), ())
                   if x[0] != t]
            if not snb:
                continue
            t2, apex2 = snb[0]
            if t2 in assign:
                continue                                   # already fixed
            rface = tuple(sorted(sig[v] for v in face))
            rnb = [x for x in dst.face2.get(rface, ()) if x[0] != rtid]
            if not rnb:
                continue
            rtid2, rapex2 = rnb[0]
            sig2 = {v: sig[v] for v in face}
            sig2[apex2] = rapex2
            if accept([int(v) for v in src.tets[t2]], sig2):
                assign[t2] = (rtid2, sig2)
                q.append(t2)
    return assign


def develop_total(view, base, img):
    """Extend the frame correspondence ``base -> img`` to a global automorphism.

    ``base`` and ``img`` are ordered 4-tuples of internal vertex ids.  Returns
    the vertex permutation as an ``int32`` array of length ``V``, or ``None``
    if propagation hits a conflict, fails to cover every tet, or is not a
    bijection.  Unlike ``develop_partial`` a mismatch is FATAL, because an
    automorphism must be globally single-valued.
    """
    tets, (nbr, napex) = view.tets, view.nbr
    sig = np.full(view.V, -1, np.int32)
    for i in range(4):
        v, w = int(base[i]), int(img[i])
        if sig[v] != -1 and sig[v] != w:
            return None
        sig[v] = w
    seen = np.full(view.nT, -1, np.int32)
    t0 = view.tetid[tuple(sorted(int(x) for x in base))]
    seen[t0] = view.tetid[tuple(sorted(int(x) for x in img))]
    stack = [t0]
    while stack:
        t = stack.pop()
        r = int(seen[t])
        tv, rv = tets[t], tets[r]
        for i in range(4):
            sx = int(sig[int(tv[i])])
            j = -1
            for k in range(4):
                if int(rv[k]) == sx:
                    j = k
                    break
            if j < 0:
                return None                     # image tet lacks the image vertex
            t2, a2 = int(nbr[t, i]), int(napex[t, i])
            r2, b2 = int(nbr[r, j]), int(napex[r, j])
            if sig[a2] == -1:
                sig[a2] = b2
            elif sig[a2] != b2:
                return None                     # holonomy: not single-valued
            if seen[t2] == -1:
                seen[t2] = r2
                stack.append(t2)
            elif seen[t2] != r2:
                return None
    if (sig < 0).any() or (seen < 0).any():
        return None                             # disconnected: did not cover
    if len(np.unique(sig)) != view.V:
        return None                             # not injective
    return sig


# ---------------------------------------------------------------------------
# the group
# ---------------------------------------------------------------------------


def _compose(g, h):
    """(g o h)(v) = g[h[v]]."""
    return g[h]


def _closure(gens, V):
    """All elements generated by ``gens`` (BFS over left multiplication)."""
    ident = np.arange(V, dtype=np.int32)
    elems = {ident.tobytes(): ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                p = _compose(g, h)
                k = p.tobytes()
                if k not in elems:
                    elems[k] = p
                    nxt.append(p)
        frontier = nxt
    return list(elems.values())


def _frame_signature(view, window, col):
    """Aut-invariant decoration of an ordered tet (lossless candidate filter)."""
    deg = view.edge_degree
    w = window
    return (tuple(col[v] for v in w),
            tuple(deg[(w[i], w[j]) if w[i] < w[j] else (w[j], w[i])]
                  for i, j in E6))


def _find_automorphisms(view, progress=None):
    """Generators of ``Aut(view)``, by development from a fixed base frame."""
    col = view.colours()
    base = tuple(int(x) for x in view.tets[0])
    bsig = _frame_signature(view, base, col)
    base_multiset = (tuple(sorted(bsig[0])), tuple(sorted(bsig[1])))

    gens, orbit = [], {base}
    tried = 0
    for t in range(view.nT):
        tv = [int(x) for x in view.tets[t]]
        # cheap tet-level screen before generating 24 orderings
        cand = _frame_signature(view, tuple(tv), col)
        if (tuple(sorted(cand[0])), tuple(sorted(cand[1]))) != base_multiset:
            continue
        for p in PERMS:
            w = (tv[p[0]], tv[p[1]], tv[p[2]], tv[p[3]])
            if w in orbit:
                continue                        # already an image: free skip
            if _frame_signature(view, w, col) != bsig:
                continue
            tried += 1
            g = develop_total(view, base, w)
            if g is None:
                continue
            gens.append(g)
            elems = _closure(gens, view.V)
            orbit = {tuple(int(e[v]) for v in base) for e in elems}
            if progress:
                progress(f"    generator {len(gens)}: |<gens>| = {len(elems)}")
    order = len(_closure(gens, view.V)) if gens else 1
    return gens, order, tried


class CrystalSymmetry:
    """``Aut(K)`` with orbit maps, stabilizers and BC-chain classes.

    Public methods take and return ORIGINAL vertex labels; internal storage is
    relabelled to ``0..V-1``.  Build with :meth:`compute` or, to reuse a cached
    sidecar, :meth:`for_manifold_path`.
    """

    def __init__(self, view, generators, order):
        self.view = view
        self.generators = [np.asarray(g, np.int32) for g in generators]
        self.order = int(order)
        self._elements = None
        self._orbits = {}
        self._chains = None

    # -- construction --------------------------------------------------------

    @classmethod
    def compute(cls, tets, relabel=True, progress=None):
        view = TriView(tets, relabel=relabel)
        if not view.is_dual_connected():
            raise ValueError(
                "triangulation is disconnected: an automorphism is NOT "
                "determined by the image of one frame across components, so "
                "development would silently report a subgroup (typically the "
                "trivial one). Split into components and handle the "
                "component-permuting part separately.")
        gens, order, tried = _find_automorphisms(view, progress=progress)
        if progress:
            progress(f"    |Aut| = {order} from {len(gens)} generators "
                     f"({tried} developments)")
        sym = cls(view, gens, order)
        sym._check_free_action()
        return sym

    @classmethod
    def for_manifold_path(cls, path, cache=True, progress=None):
        """Symmetry of a ``.mfd`` file, cached to ``<path>.sym.npz``.

        The cache is keyed by a content hash of the facet set, so a regenerated
        or relabelled crystal invalidates it rather than silently reusing the
        wrong group.
        """
        from ._manifold import Manifold
        tets = np.asarray(Manifold.load(path, 3).facets())
        view = TriView(tets)
        side = path + ".sym.npz"
        if cache and os.path.exists(side):
            sym = cls._load_sidecar(side, view)
            if sym is not None:
                return sym
        sym = cls.compute(tets, progress=progress)
        if cache:
            try:
                np.savez_compressed(
                    side, digest=sym.view.digest(), order=sym.order,
                    gens=(np.stack(sym.generators) if sym.generators
                          else np.zeros((0, sym.view.V), np.int32)))
            except OSError as exc:
                warnings.warn(f"could not write symmetry cache {side}: {exc}")
        return sym

    @classmethod
    def _load_sidecar(cls, side, view):
        """Load cached generators, VERIFYING them; ``None`` means recompute.

        A cache of a mathematical object is only worth having if using it is
        indistinguishable from recomputing it, so nothing here is taken on
        trust: the content hash must match, every generator must be checked to
        permute the facet set, the group order is re-derived from the closure
        rather than read, and the free-action invariant is re-asserted. A
        sidecar that fails any of these is discarded WITH A WARNING, never
        silently half-used.
        """
        try:
            z = np.load(side)
            digest, gens, order = str(z["digest"]), z["gens"], int(z["order"])
        except Exception as exc:
            warnings.warn(f"unreadable symmetry cache {side} ({exc}); "
                          "recomputing")
            return None
        if digest != view.digest():
            return None                     # different triangulation: expected
        gens = [np.asarray(gens[i], np.int32) for i in range(len(gens))]
        facets = {tuple(sorted(int(x) for x in tv)) for tv in view.tets}
        for g in gens:
            if len(g) != view.V or len(np.unique(g)) != view.V:
                warnings.warn(f"symmetry cache {side} holds a non-permutation; "
                              "recomputing")
                return None
            if {tuple(sorted(int(g[int(x)]) for x in tv))
                    for tv in view.tets} != facets:
                warnings.warn(f"symmetry cache {side} holds a generator that is "
                              "not an automorphism; recomputing")
                return None
        sym = cls(view, gens, len(_closure(gens, view.V)))
        if sym.order != order:
            warnings.warn(f"symmetry cache {side} claims |Aut| = {order} but "
                          f"its generators close to {sym.order}; recomputing")
            return None
        try:
            sym._check_free_action()
        except AssertionError as exc:
            warnings.warn(f"symmetry cache {side} fails the free-action check "
                          f"({exc}); recomputing")
            return None
        return sym

    def _check_free_action(self):
        """``Aut`` acts freely on frames, so ``24*nT`` must divide by ``|Aut|``.

        A failure means the development argument broke (disconnected complex,
        or a bug), not that the crystal is unusual -- so this is an assert.
        """
        n = 24 * self.view.nT
        if n % self.order:
            raise AssertionError(
                f"free action violated: 24*nT = {n} not divisible by "
                f"|Aut| = {self.order}")

    # -- elements ------------------------------------------------------------

    @property
    def elements(self):
        """All ``|Aut|`` permutations (internal ids).  Built on demand."""
        if self._elements is None:
            self._elements = _closure(self.generators, self.view.V)
        return self._elements

    def permutation(self, g):
        """Group element as ``{original label: original label}``."""
        lab = self.view.labels
        return {int(lab[i]): int(lab[g[i]]) for i in range(self.view.V)}

    def act(self, g, obj):
        """Apply ``g`` to a vertex, or to any tuple/list of vertices.

        Simplices (edges, faces, tets) come back sorted; ordered tuples such as
        frames and chain windows keep their order, so pass a tuple when order
        matters and a sorted tuple when it does not.
        """
        ix, lab = self.view.index, self.view.labels
        if np.isscalar(obj):
            return int(lab[g[ix[int(obj)]]])
        return tuple(int(lab[g[ix[int(v)]]]) for v in obj)

    # -- orbits --------------------------------------------------------------

    def _orbit_table(self, kind):
        if kind in self._orbits:
            return self._orbits[kind]
        v = self.view
        if kind == "vertex":
            items = list(range(v.V))
            act = lambda g, x: int(g[x])
        elif kind == "edge":
            items = v.edges
            act = lambda g, e: (min(int(g[e[0]]), int(g[e[1]])),
                                max(int(g[e[0]]), int(g[e[1]])))
        elif kind == "face":
            items = v.faces
            act = lambda g, f: tuple(sorted(int(g[x]) for x in f))
        elif kind == "tet":
            items = v.tet_keys
            act = lambda g, t: tuple(sorted(int(g[x]) for x in t))
        else:
            raise ValueError(f"unknown kind {kind!r}")
        pos = {x: i for i, x in enumerate(items)}
        par = list(range(len(items)))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        for g in self.generators:
            for x in items:
                a, b = find(pos[x]), find(pos[act(g, x)])
                if a != b:
                    par[a] = b
        groups = collections.defaultdict(list)
        for x in items:
            groups[find(pos[x])].append(x)
        reps = sorted(groups, key=lambda r: -len(groups[r]))
        oid = {}
        members = []
        for i, r in enumerate(reps):
            members.append(groups[r])
            for x in groups[r]:
                oid[x] = i
        self._orbits[kind] = (oid, members, act)
        return self._orbits[kind]

    def _to_internal(self, kind, obj):
        ix = self.view.index
        if kind == "vertex":
            return ix[int(obj)]
        key = tuple(sorted(ix[int(x)] for x in obj))
        return key

    def _to_external(self, kind, obj):
        lab = self.view.labels
        if kind == "vertex":
            return int(lab[obj])
        return tuple(int(lab[x]) for x in obj)

    def orbit_id(self, kind, obj):
        """Orbit index of a vertex / edge / face / tet (original labels)."""
        oid, _, _ = self._orbit_table(kind)
        return oid[self._to_internal(kind, obj)]

    def n_orbits(self, kind):
        return len(self._orbit_table(kind)[1])

    def orbit_members(self, kind, i):
        _, members, _ = self._orbit_table(kind)
        return [self._to_external(kind, x) for x in members[i]]

    def orbit_sizes(self, kind):
        return [len(m) for m in self._orbit_table(kind)[1]]

    def orbit_representatives(self, kind):
        _, members, _ = self._orbit_table(kind)
        return [self._to_external(kind, m[0]) for m in members]

    def orbit_id_map(self, kind):
        """``{object: orbit index}`` over the whole complex, original labels."""
        oid, _, _ = self._orbit_table(kind)
        return {self._to_external(kind, x): i for x, i in oid.items()}

    def stabilizer_order(self, kind, obj):
        """``|Aut| / |orbit|`` -- exact, and free of any element enumeration."""
        _, members, _ = self._orbit_table(kind)
        return self.order // len(members[self.orbit_id(kind, obj)])

    def stabilizer(self, kind, obj):
        """The stabilizer subgroup's elements (materializes ``Aut``)."""
        _, _, act = self._orbit_table(kind)
        x = self._to_internal(kind, obj)
        return [g for g in self.elements if act(g, x) == x]

    # -- BC chains -----------------------------------------------------------

    def _chain_tables(self):
        """Enumerate every BC chain as a cycle of the sliding-window map.

        The map (v0,v1,v2,v3) -> (v1,v2,v3, other apex of face v1v2v3) is a
        BIJECTION on frames (its inverse prepends the other apex of v0v1v2), so
        the ``24*nT`` frames partition into disjoint cycles -- the complete set
        of directed BC chains, obtained in one linear pass.  This is the global
        version of ``worm_helix.bc_orbit``, which walks one chain from one
        hand-picked seed tet.
        """
        if self._chains is not None:
            return self._chains
        v = self.view
        tets, (nbr, napex) = v.tets, v.nbr
        pos = [None] * v.nT
        for t in range(v.nT):
            tv = tets[t]
            pos[t] = {int(tv[0]): 0, int(tv[1]): 1, int(tv[2]): 2, int(tv[3]): 3}
        nxt = np.empty(24 * v.nT, np.int64)
        for t in range(v.nT):
            tv = [int(x) for x in tets[t]]
            for pi, p in enumerate(PERMS):
                w1, w2, w3 = tv[p[1]], tv[p[2]], tv[p[3]]
                t2 = int(nbr[t, p[0]])
                a2 = int(napex[t, p[0]])
                q = pos[t2]
                nxt[t * 24 + pi] = t2 * 24 + PERM_INDEX[
                    (q[w1], q[w2], q[w3], q[a2])]
        chain_of = np.full(24 * v.nT, -1, np.int64)
        chains = []
        for f in range(24 * v.nT):
            if chain_of[f] >= 0:
                continue
            c, x = [], f
            cid = len(chains)
            while chain_of[x] < 0:
                chain_of[x] = cid
                c.append(x)
                x = int(nxt[x])
            chains.append(np.array(c, np.int64))
        self._chains = (nxt, chain_of, chains)
        return self._chains

    @property
    def chains(self):
        """All directed BC chains, each as an array of frame ids."""
        return self._chain_tables()[2]

    def chain_vertices(self, i):
        """Chain ``i`` as its cyclic vertex sequence (original labels), matching
        ``worm_helix.bc_orbit``'s convention of taking window[0] per step."""
        v = self.view
        return [int(v.labels[v.frame_window(f)[0]]) for f in self.chains[i]]

    def chain_of_window(self, window):
        """Which chain an ordered vertex 4-tuple (original labels) lies on."""
        w = tuple(self.view.index[int(x)] for x in window)
        return int(self._chain_tables()[1][self.view.frame_id(w)])

    def chain_orbits(self):
        """``(orbit id per chain, [members], [representative chain ids])``."""
        if "chain" in self._orbits:
            return self._orbits["chain"]
        _, chain_of, chains = self._chain_tables()
        v = self.view
        par = list(range(len(chains)))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        for g in self.generators:
            for i, c in enumerate(chains):
                w = v.frame_window(int(c[0]))
                j = int(chain_of[v.frame_id(tuple(int(g[x]) for x in w))])
                a, b = find(i), find(j)
                if a != b:
                    par[a] = b
        groups = collections.defaultdict(list)
        for i in range(len(chains)):
            groups[find(i)].append(i)
        reps = sorted(groups, key=lambda r: -len(groups[r]))
        oid = np.empty(len(chains), np.int64)
        members = []
        for k, r in enumerate(reps):
            members.append(groups[r])
            for i in groups[r]:
                oid[i] = k
        self._orbits["chain"] = (oid, members, [m[0] for m in members])
        return self._orbits["chain"]

    # -- reporting -----------------------------------------------------------

    def summary(self):
        oid, members, reps = self.chain_orbits()
        lens = collections.Counter(len(c) for c in self.chains)
        return dict(
            V=self.view.V, nT=self.view.nT, order=self.order,
            n_generators=len(self.generators),
            frame_orbits=24 * self.view.nT // self.order,
            orbits={k: self.n_orbits(k)
                    for k in ("vertex", "edge", "face", "tet")},
            n_chains=len(self.chains),
            chain_lengths=dict(sorted(lens.items())),
            n_chain_orbits=len(members),
            chain_orbit_sizes=[len(m) for m in members],
        )
