"""Read-only view of a Manifold, for use by ManifoldSampler."""

from __future__ import annotations

import ctypes
from typing import Sequence

import numpy as np

from . import _dlang
from ._manifold import Manifold

_lib = _dlang._lib


class ManifoldView:
    """Read-only view of a manifold.

    Provides all query methods of Manifold but no mutation (do_move, save, etc.).
    The view borrows the underlying handle — it does NOT own it, so it must not
    outlive the object that created it (e.g. a ManifoldSampler).

    Use ``dup()`` to get a mutable, independently-owned Manifold copy.
    """

    def __init__(self, handle):
        self._handle = handle

    # No __del__ — we don't own the handle.

    def dup(self) -> Manifold:
        """Return a deep, mutable copy of this manifold."""
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

    def mean_degree(self, dim: int) -> float:
        """Return the mean degree of simplices of given dimension."""
        return _lib.ddg_manifold_mean_degree(self._handle, dim)

    def count_valid_moves(self) -> int:
        """Count valid Pachner moves (including stellar subdivisions)."""
        return _lib.ddg_manifold_count_valid_moves(self._handle)

    def importance_weight(self) -> float:
        """Return 1/V, the importance weight correcting the sampler's
        stationary distribution back to exp(-objective).
        """
        return _lib.ddg_manifold_importance_weight(self._handle)

    # -- I/O (read-only: saves current state, does not mutate) --

    def save(self, path: str) -> None:
        """Save to a .mfd file."""
        _lib.ddg_manifold_save(self._handle, path.encode())

    def save_edge_graph(self, path: str) -> None:
        """Save the 1-skeleton as an edge list."""
        _lib.ddg_manifold_save_edge_graph(self._handle, path.encode())

    def save_dual_graph(self, path: str) -> None:
        """Save the dual graph."""
        _lib.ddg_manifold_save_dual_graph(self._handle, path.encode())

    def __repr__(self):
        return f"ManifoldView(dim={self.dimension}, num_facets={self.num_facets})"
