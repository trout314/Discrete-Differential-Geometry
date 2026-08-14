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

    def illegal_edges(self) -> tuple[np.ndarray, np.ndarray]:
        """FK-illegal edge set (degree not in {5, 6}; dim 3 only).

        Returns (edges (n, 2) int array, degrees (n,) int array). One O(E)
        scan in D -- cheap enough for per-chunk instrumentation."""
        n = _lib.ddg_manifold_illegal_edges(self._handle, None, None)
        if n == 0:
            return (np.empty((0, 2), dtype=np.intc),
                    np.empty(0, dtype=np.intc))
        pairs = np.empty((n, 2), dtype=np.intc)
        degs = np.empty(n, dtype=np.intc)
        _lib.ddg_manifold_illegal_edges(
            self._handle,
            pairs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            degs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        return pairs, degs

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
        """Weight converting default-sampler samples back to exp(-objective).

        ``mcmcStep``'s proposal is a REJECTION LOOP, so it proposes a move not
        with the raw draw probability but with the conditioned
        ``q(m|x) = c(m)/C(x)``, where ``c = 1`` for the ``1->4`` facet centre
        (whose keep-coin clips) and ``2`` otherwise, and

            C = 2*count_valid_moves() - f3 + 2*count_valid_hinge_moves()

        The applied Hastings factor uses ``f3(x)/f3(y)`` in place of
        ``C(x)/C(y)``, so the chain is stationary for ``exp(-S) * C/f3`` --
        states with more escape routes are over-weighted. This returns the
        correcting weight ``w = f3/C``. Estimate an ``exp(-S)`` expectation as

            sum(w_i * O_i) / sum(w_i)

        Only ratios matter, so the overall scale is irrelevant.

        The bias is small at production size (``C/f3`` is intensive, so the
        drift is O(1/N): measured ``<d ln(C/f3)>`` per move is a few times
        1e-3, and it carries no vertex fugacity), but this makes it exact.

        VALID ONLY FOR THE PURE BISTELLAR+HINGE CHAIN: plain ``run()`` at
        dim 3 with ``set_bistellar_hastings`` ON (the default). With the
        contract/split channel, worm or non-local slide enabled the chain is a
        MIXTURE of kernels with different targets, whose stationary law has no
        closed form -- no per-state weight corrects it. dim 2/4 raise, because
        their sampler applies no Hastings factor at all; use
        ``run(exact=True)`` there, which is exact in every dimension.
        """
        return _lib.ddg_manifold_importance_weight(self._handle)

    # -- I/O (read-only: saves current state, does not mutate) --

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

    def __repr__(self):
        return f"ManifoldView(dim={self.dimension}, num_facets={self.num_facets})"
