"""Pythonic wrapper for Manifold."""

from __future__ import annotations

import ctypes
from typing import Sequence

import numpy as np

from . import _dlang
from ._simplicial_complex import SimplicialComplex

_lib = _dlang._lib


class Manifold:
    """A combinatorial manifold of a fixed dimension.

    Parameters
    ----------
    dimension : int
        The dimension of the manifold (2, 3, or 4).
    facets : sequence of sequences of int
        The facets (maximal simplices). Each must have dimension+1 vertices.
    """

    def __init__(self, dimension: int, facets, *, _handle=None):
        if _handle is not None:
            self._handle = _handle
            return

        facets_list = [list(f) for f in facets]
        flat = []
        for f in facets_list:
            flat.extend(f)
        arr = (ctypes.c_int * len(flat))(*flat)
        self._handle = _lib.ddg_manifold_from_facets(dimension, arr, len(facets_list))

    def __del__(self):
        if hasattr(self, "_handle") and self._handle is not None:
            _lib.ddg_manifold_free(self._handle)
            self._handle = None

    @classmethod
    def standard_sphere(cls, dim: int) -> Manifold:
        """Create the standard sphere triangulation of the given dimension."""
        handle = _lib.ddg_manifold_standard_sphere(dim)
        obj = cls.__new__(cls)
        obj._handle = handle
        return obj

    @classmethod
    def load(cls, path: str, dim: int) -> Manifold:
        """Load a manifold from a .mfd file."""
        handle = _lib.ddg_manifold_load(path.encode(), dim)
        obj = cls.__new__(cls)
        obj._handle = handle
        return obj

    def copy(self) -> Manifold:
        """Return a deep copy."""
        obj = Manifold.__new__(Manifold)
        obj._handle = _lib.ddg_manifold_copy(self._handle)
        return obj

    # -- Properties --

    @property
    def dimension(self) -> int:
        return _lib.ddg_manifold_dimension(self._handle)

    @property
    def num_facets(self) -> int:
        return _lib.ddg_manifold_num_facets(self._handle)

    @property
    def f_vector(self) -> np.ndarray:
        buf = (ctypes.c_long * 10)()
        n = _lib.ddg_manifold_f_vector(self._handle, buf, 10)
        return np.array(buf[:n], dtype=np.int64)

    @property
    def euler_characteristic(self) -> int:
        return _lib.ddg_manifold_euler_characteristic(self._handle)

    @property
    def is_orientable(self) -> bool:
        return bool(_lib.ddg_manifold_is_orientable(self._handle))

    # -- Data access --

    def facets(self) -> np.ndarray:
        """Return all facets as ndarray of shape (n, dim+1)."""
        dim = self.dimension
        count = _lib.ddg_manifold_facets(self._handle, None)
        if count == 0:
            return np.empty((0, dim + 1), dtype=np.intc)
        buf = (ctypes.c_int * (count * (dim + 1)))()
        _lib.ddg_manifold_facets(self._handle, buf)
        return np.frombuffer(buf, dtype=np.intc).reshape(count, dim + 1).copy()

    def simplices(self, dim: int) -> np.ndarray:
        """Return all simplices of given dimension as ndarray of shape (n, dim+1)."""
        count = _lib.ddg_manifold_simplices(self._handle, dim, None)
        if count == 0:
            return np.empty((0, dim + 1), dtype=np.intc)
        buf = (ctypes.c_int * (count * (dim + 1)))()
        _lib.ddg_manifold_simplices(self._handle, dim, buf)
        return np.frombuffer(buf, dtype=np.intc).reshape(count, dim + 1).copy()

    def degree(self, simplex: Sequence[int]) -> int:
        """Return the degree of a simplex in the manifold."""
        arr = np.asarray(simplex, dtype=np.intc).ravel()
        ptr = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        return _lib.ddg_manifold_degree(self._handle, ptr, len(arr))

    def degree_histogram(self, dim: int) -> np.ndarray:
        """Return degree histogram: result[i] = count of simplices with degree i+1."""
        n = _lib.ddg_manifold_degree_histogram(self._handle, dim, None)
        if n == 0:
            return np.empty(0, dtype=np.int64)
        buf = (ctypes.c_long * n)()
        _lib.ddg_manifold_degree_histogram(self._handle, dim, buf)
        return np.array(buf[:n], dtype=np.int64)

    def mean_degree(self, dim: int) -> float:
        """Return the mean degree of simplices of given dimension."""
        return _lib.ddg_manifold_mean_degree(self._handle, dim)

    def degree_variance(self, dim: int) -> float:
        """Return the variance of degree for simplices of given dimension."""
        return _lib.ddg_manifold_degree_variance(self._handle, dim)

    def do_bistellar_move(self, center: Sequence[int],
                          cocenter: Sequence[int]) -> None:
        """Apply the SPECIFIED bistellar (Pachner) move: replace star(center)
        (= center joined with the boundary of cocenter) by the star of
        cocenter. dim=3 grammar: 1-4 center=facet cocenter=[new vertex];
        2-3 center=face, cocenter=the two apexes; 3-2 center=a degree-3
        edge, cocenter=its 3 link vertices; 4-1 center=a degree-4 vertex,
        cocenter=its 4 neighbors. All preconditions are validated in the D
        core (raises on an invalid move).

        WARNING: bypasses any ManifoldSampler cocycle tracking on this
        manifold -- analysis/catalog use."""
        c = np.asarray(center, dtype=np.intc).ravel()
        cc = np.asarray(cocenter, dtype=np.intc).ravel()
        _lib.ddg_manifold_do_bistellar_move(
            self._handle, c.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            len(c), cc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(cc))

    def has_bistellar_move(self, center: Sequence[int],
                           cocenter: Sequence[int]) -> bool:
        """True iff do_bistellar_move(center, cocenter) would be valid."""
        c = np.asarray(center, dtype=np.intc).ravel()
        cc = np.asarray(cocenter, dtype=np.intc).ravel()
        return bool(_lib.ddg_manifold_has_bistellar_move(
            self._handle, c.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            len(c), cc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(cc)))

    def do_hinge_move(self, removed_edge: Sequence[int],
                      link_cycle: Sequence[int], diagonal: int) -> None:
        """Apply the specified 4-4 hinge move (dim=3): replace the star of
        degree-4 edge `removed_edge` (link cycle `link_cycle`, in cyclic
        order) with the retriangulation using `diagonal` = 0
        (cycle[0]-cycle[2]) or 1 (cycle[1]-cycle[3]). Validated in the D
        core. Same cocycle-tracking warning as do_bistellar_move."""
        re_ = np.asarray(removed_edge, dtype=np.intc).ravel()
        cy = np.asarray(link_cycle, dtype=np.intc).ravel()
        _lib.ddg_manifold_do_hinge_move(
            self._handle, re_.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            cy.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), int(diagonal))

    def has_hinge_move(self, removed_edge: Sequence[int],
                       link_cycle: Sequence[int], diagonal: int) -> bool:
        """True iff do_hinge_move(...) would be valid."""
        re_ = np.asarray(removed_edge, dtype=np.intc).ravel()
        cy = np.asarray(link_cycle, dtype=np.intc).ravel()
        return bool(_lib.ddg_manifold_has_hinge_move(
            self._handle, re_.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            cy.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), int(diagonal)))

    def freeze_vertices(self, vertices, frozen: bool = True) -> None:
        """Freeze (or unfreeze) vertices. The sampler rejects any move whose
        support contains a frozen vertex; since every facet a move adds or
        removes has all its vertices in the support, this preserves the frozen
        set's entire CLOSED STAR exactly -- the facets and the degrees (hence
        Regge curvature) of every simplex meeting a frozen vertex."""
        arr = np.asarray(list(vertices), dtype=np.intc).ravel()
        ptr = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        _lib.ddg_manifold_freeze_vertices(self._handle, ptr, len(arr),
                                          1 if frozen else 0)

    def freeze_facets(self, facets) -> np.ndarray:
        """Freeze every vertex of the given facets (rows of vertex labels).
        This pins the facets AND the curvature on all their hinges (the
        star-closure needed to 'fix the curvature on a shell'). Returns the
        frozen vertex labels."""
        verts = np.unique(np.asarray(facets, dtype=np.intc).ravel())
        self.freeze_vertices(verts)
        return verts

    def clear_frozen(self) -> None:
        """Unfreeze all vertices."""
        _lib.ddg_manifold_clear_frozen(self._handle)

    def vertex_frozen(self, v: int) -> bool:
        """Is this vertex frozen?"""
        return bool(_lib.ddg_manifold_vertex_frozen(self._handle, int(v)))

    def num_frozen(self) -> int:
        """Number of frozen vertices."""
        return int(_lib.ddg_manifold_num_frozen(self._handle))

    def count_valid_moves(self) -> int:
        """Count valid Pachner moves (including stellar subdivisions)."""
        return _lib.ddg_manifold_count_valid_moves(self._handle)

    def valence_census(self, n6_cap: int = 8, m_cap: int = 8) -> np.ndarray:
        """Joint (n6, m) vertex census, shape (n6_cap+1, m_cap+1); dim=3 only.

        n6 = # incident degree>=6 edges, m = impurity valence; final bins
        clamp. Bin (0, 0) is the FK-Z12 bulk; rows n6 >= 1 are the
        disclination-network vertices.
        """
        from .disclination import valence_census_from_handle
        return valence_census_from_handle(self._handle, n6_cap, m_cap)

    def disclination_census(self, host_classes=None) -> dict:
        """Census of the degree>=6 edge graph (the disclination network);
        dim=3 only. host_classes = native host n6 classes (C15: [0, 4])
        enables the host/dopant graft split. See disclination.py for keys.
        """
        from .disclination import disclination_census_from_handle
        return disclination_census_from_handle(self._handle, host_classes)

    def importance_weight(self) -> float:
        """Return 1/V(x), the importance weight correcting the sampler's
        stationary distribution back to exp(-objective(x)).

        The default sampler (pure Metropolis, no Hastings correction) samples
        from pi(x) ~ exp(-obj(x)) * V(x), where V(x) is the number of valid
        Pachner moves. Multiplying observables by this weight corrects for
        the bias.
        """
        return _lib.ddg_manifold_importance_weight(self._handle)

    # -- Moves --

    def do_move(self) -> None:
        """Perform a random Pachner move."""
        _lib.ddg_manifold_do_pachner_move(self._handle)

    # -- I/O --

    def save(self, path: str, comments: list[str] | None = None) -> None:
        """Save to a .mfd file."""
        if comments:
            arr = (ctypes.c_char_p * len(comments))(*[c.encode() for c in comments])
            _lib.ddg_manifold_save_with_comments(self._handle, path.encode(), arr, len(comments))
        else:
            _lib.ddg_manifold_save(self._handle, path.encode())

    def save_edge_graph(self, path: str) -> None:
        """Save the 1-skeleton as an edge list."""
        _lib.ddg_manifold_save_edge_graph(self._handle, path.encode())

    def save_dual_graph(self, path: str) -> None:
        """Save the dual graph."""
        _lib.ddg_manifold_save_dual_graph(self._handle, path.encode())

    def to_simplicial_complex(self) -> SimplicialComplex:
        """Convert to a SimplicialComplex (copies data)."""
        facets_arr = self.facets()
        return SimplicialComplex(facets_arr.tolist())

    def __repr__(self):
        return f"Manifold(dim={self.dimension}, num_facets={self.num_facets})"
