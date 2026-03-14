"""Pythonic wrapper for SimplicialComplex."""

from __future__ import annotations

import ctypes
from typing import Sequence

import numpy as np

from . import _dlang

_lib = _dlang._lib


def _to_c_array(simplex):
    """Convert a Python sequence to a ctypes int array."""
    arr = np.asarray(simplex, dtype=np.intc).ravel()
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(arr)


class SimplicialComplex:
    """A simplicial complex with integer vertices.

    Parameters
    ----------
    facets : sequence of sequences of int, optional
        Maximal simplices. All facets in one call must have the same dimension.
    """

    def __init__(self, facets=None, *, _handle=None):
        if _handle is not None:
            self._handle = _handle
            return

        if facets is None or len(facets) == 0:
            self._handle = _lib.ddg_sc_create()
            return

        facets_list = [list(f) for f in facets]
        verts_per_facet = len(facets_list[0])
        flat = []
        for f in facets_list:
            if len(f) != verts_per_facet:
                raise ValueError(
                    "All facets must have the same number of vertices "
                    "when constructing from a list. Use insert_facet() "
                    "for mixed-dimension facets."
                )
            flat.extend(f)
        arr = (ctypes.c_int * len(flat))(*flat)
        self._handle = _lib.ddg_sc_from_facets(arr, len(facets_list), verts_per_facet)

    def __del__(self):
        if hasattr(self, "_handle") and self._handle is not None:
            _lib.ddg_sc_free(self._handle)
            self._handle = None

    @classmethod
    def load(cls, path: str) -> SimplicialComplex:
        """Load a simplicial complex from a .sc file."""
        handle = _lib.ddg_sc_load(path.encode())
        return cls(_handle=handle)

    def copy(self) -> SimplicialComplex:
        """Return a deep copy."""
        return SimplicialComplex(_handle=_lib.ddg_sc_copy(self._handle))

    # -- Properties --

    @property
    def num_facets(self) -> int:
        return _lib.ddg_sc_num_facets(self._handle)

    @property
    def f_vector(self) -> np.ndarray:
        buf = (ctypes.c_long * 20)()
        n = _lib.ddg_sc_f_vector(self._handle, buf, 20)
        return np.array(buf[:n], dtype=np.int64)

    @property
    def euler_characteristic(self) -> int:
        return _lib.ddg_sc_euler_characteristic(self._handle)

    @property
    def is_connected(self) -> bool:
        return bool(_lib.ddg_sc_is_connected(self._handle))

    @property
    def is_pure(self) -> bool:
        return bool(_lib.ddg_sc_is_pure(self._handle))

    # -- Data access --

    def facets(self, dim: int | None = None) -> np.ndarray:
        """Return facets as an ndarray of shape (n, k).

        Parameters
        ----------
        dim : int, optional
            If given, return only facets of this dimension.
        """
        if dim is not None:
            count = _lib.ddg_sc_facets_of_dim(self._handle, dim, None)
            if count == 0:
                return np.empty((0, dim + 1), dtype=np.intc)
            buf = (ctypes.c_int * (count * (dim + 1)))()
            _lib.ddg_sc_facets_of_dim(self._handle, dim, buf)
            return np.frombuffer(buf, dtype=np.intc).reshape(count, dim + 1).copy()

        # All facets — need sizes since they may differ
        count = _lib.ddg_sc_facets(self._handle, None, None)
        if count == 0:
            return np.empty((0,), dtype=object)
        # Get sizes first
        sizes = (ctypes.c_int * count)()
        # Estimate total vertices
        total_verts = count * 20  # generous upper bound
        data = (ctypes.c_int * total_verts)()
        actual = _lib.ddg_sc_facets(self._handle, data, sizes)
        # Check if all same size
        all_same = all(sizes[i] == sizes[0] for i in range(actual))
        if all_same and actual > 0:
            k = sizes[0]
            total = sum(sizes[i] for i in range(actual))
            return np.frombuffer(data, dtype=np.intc, count=total).reshape(actual, k).copy()
        # Mixed sizes — return list of arrays
        result = []
        idx = 0
        for i in range(actual):
            s = sizes[i]
            result.append(np.array(data[idx : idx + s], dtype=np.intc))
            idx += s
        return np.array(result, dtype=object)

    def simplices(self, dim: int) -> np.ndarray:
        """Return all simplices of given dimension as ndarray of shape (n, dim+1)."""
        count = _lib.ddg_sc_simplices(self._handle, dim, None)
        if count == 0:
            return np.empty((0, dim + 1), dtype=np.intc)
        buf = (ctypes.c_int * (count * (dim + 1)))()
        _lib.ddg_sc_simplices(self._handle, dim, buf)
        return np.frombuffer(buf, dtype=np.intc).reshape(count, dim + 1).copy()

    def contains(self, simplex: Sequence[int]) -> bool:
        """Check if a simplex (of any dimension) is in the complex."""
        ptr, n = _to_c_array(simplex)
        return bool(_lib.ddg_sc_contains(self._handle, ptr, n))

    def contains_facet(self, facet: Sequence[int]) -> bool:
        """Check if a facet is in the complex."""
        ptr, n = _to_c_array(facet)
        return bool(_lib.ddg_sc_contains_facet(self._handle, ptr, n))

    def star(self, simplex: Sequence[int]) -> list[np.ndarray]:
        """Return the star of a simplex (facets containing it)."""
        ptr, n = _to_c_array(simplex)
        count = _lib.ddg_sc_star(self._handle, ptr, n, None, None)
        if count == 0:
            return []
        sizes = (ctypes.c_int * count)()
        total_verts = count * 20
        data = (ctypes.c_int * total_verts)()
        _lib.ddg_sc_star(self._handle, ptr, n, data, sizes)
        result = []
        idx = 0
        for i in range(count):
            s = sizes[i]
            result.append(np.array(data[idx : idx + s], dtype=np.intc))
            idx += s
        return result

    def link(self, simplex: Sequence[int]) -> list[np.ndarray]:
        """Return the link of a simplex."""
        ptr, n = _to_c_array(simplex)
        count = _lib.ddg_sc_link(self._handle, ptr, n, None, None)
        if count == 0:
            return []
        sizes = (ctypes.c_int * count)()
        total_verts = count * 20
        data = (ctypes.c_int * total_verts)()
        _lib.ddg_sc_link(self._handle, ptr, n, data, sizes)
        result = []
        idx = 0
        for i in range(count):
            s = sizes[i]
            result.append(np.array(data[idx : idx + s], dtype=np.intc))
            idx += s
        return result

    # -- Mutation --

    def insert_facet(self, simplex: Sequence[int]) -> None:
        """Insert a facet into the complex."""
        ptr, n = _to_c_array(simplex)
        _lib.ddg_sc_insert_facet(self._handle, ptr, n)

    def remove_facet(self, simplex: Sequence[int]) -> None:
        """Remove a facet from the complex."""
        ptr, n = _to_c_array(simplex)
        _lib.ddg_sc_remove_facet(self._handle, ptr, n)

    # -- I/O --

    def save(self, path: str) -> None:
        """Save to a .sc file."""
        _lib.ddg_sc_save(self._handle, path.encode())

    def save_edge_graph(self, path: str) -> None:
        """Save the 1-skeleton as an edge list."""
        _lib.ddg_sc_save_edge_graph(self._handle, path.encode())

    # -- Topology --

    def connected_components(self) -> list[SimplicialComplex]:
        """Return connected components as a list of SimplicialComplex."""
        count = _lib.ddg_sc_connected_components(self._handle, None)
        if count == 0:
            return []
        handles = (ctypes.c_void_p * count)()
        _lib.ddg_sc_connected_components(self._handle, handles)
        return [SimplicialComplex(_handle=handles[i]) for i in range(count)]

    def is_pure_of_dim(self, d: int) -> bool:
        return bool(_lib.ddg_sc_is_pure_of_dim(self._handle, d))

    def is_circle(self) -> bool:
        return bool(_lib.ddg_sc_is_circle(self._handle))

    def is_2_sphere(self) -> bool:
        return bool(_lib.ddg_sc_is_2_sphere(self._handle))

    def is_2_torus(self) -> bool:
        return bool(_lib.ddg_sc_is_2_torus(self._handle))

    def is_orientable_surface_of_genus(self, g: int) -> bool:
        return bool(_lib.ddg_sc_is_orientable_surface_of_genus(self._handle, g))

    def __repr__(self):
        return f"SimplicialComplex(num_facets={self.num_facets})"


def join(sc1: SimplicialComplex, sc2: SimplicialComplex) -> SimplicialComplex:
    """Return the join of two simplicial complexes."""
    handle = _lib.ddg_sc_join(sc1._handle, sc2._handle)
    return SimplicialComplex(_handle=handle)
