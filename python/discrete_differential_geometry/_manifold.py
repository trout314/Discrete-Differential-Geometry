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

    def count_valid_moves(self) -> int:
        """Count valid Pachner moves (including stellar subdivisions)."""
        return _lib.ddg_manifold_count_valid_moves(self._handle)

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

    def save(self, path: str) -> None:
        """Save to a .mfd file."""
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
