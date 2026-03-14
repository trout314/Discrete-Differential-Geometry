"""Pythonic wrapper for the MCMC manifold sampler."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
from typing import Callable, Optional

from . import _dlang
from ._manifold import Manifold

_lib = _dlang._lib


@dataclass
class SamplerParams:
    """Parameters for the Metropolis-Hastings manifold sampler.

    Attributes
    ----------
    num_facets_target : int
        Target number of top-dimensional simplices.
    hinge_degree_target : float
        Target average degree for codimension-2 faces (hinges).
    num_facets_coef : float
        Coefficient for the volume (facet count) penalty.
    num_hinges_coef : float
        Coefficient for the global curvature (hinge count) penalty.
    hinge_degree_variance_coef : float
        Coefficient for the local curvature (hinge degree variance) penalty.
    codim3_degree_variance_coef : float
        Coefficient for the codimension-3 degree variance penalty.
    """

    num_facets_target: int = 100
    hinge_degree_target: float = 4.5
    num_facets_coef: float = 0.1
    num_hinges_coef: float = 0.05
    hinge_degree_variance_coef: float = 0.2
    codim3_degree_variance_coef: float = 0.1


class ManifoldSampler:
    """MCMC sampler for manifold triangulations using Pachner moves.

    Parameters
    ----------
    manifold : Manifold
        The initial manifold (will be copied).
    params : SamplerParams
        Sampling parameters.
    """

    def __init__(self, manifold: Manifold, params: SamplerParams):
        self._params = params
        self._handle = _lib.ddg_sampler_create(
            manifold._handle,
            params.num_facets_target,
            params.hinge_degree_target,
            params.num_facets_coef,
            params.num_hinges_coef,
            params.hinge_degree_variance_coef,
            params.codim3_degree_variance_coef,
        )
        # Hold a reference to keep the callback alive
        self._callback_ref = None

    def __del__(self):
        if hasattr(self, "_handle") and self._handle is not None:
            _lib.ddg_sampler_free(self._handle)
            self._handle = None

    def run(
        self,
        *,
        sweeps: float | None = None,
        moves: int | None = None,
        exact: bool = False,
        callback: Callable[[int, int], bool] | None = None,
    ) -> int:
        """Run the sampler.

        Exactly one of ``sweeps`` or ``moves`` must be given.

        Parameters
        ----------
        sweeps : float, optional
            Number of sweeps (1 sweep = num_facets attempted moves).
        moves : int, optional
            Raw number of attempted moves.
        exact : bool, optional
            If True, use exact Hastings correction (execute-then-undo with
            countValidBistellarMoves). Slower but samples from the exact
            target distribution exp(-objective). Default False uses the
            approximate nFacets ratio.
        callback : callable, optional
            Called periodically with ``(moves_done, moves_total)``.
            Return ``True`` to stop early.

        Returns
        -------
        int
            Number of moves accepted.
        """
        if (sweeps is None) == (moves is None):
            raise ValueError("Exactly one of 'sweeps' or 'moves' must be given")

        if sweeps is not None:
            mfd_handle = _lib.ddg_sampler_get_manifold(self._handle)
            num_facets = _lib.ddg_manifold_num_facets(mfd_handle)
            moves = int(sweeps * num_facets)

        if callback is not None:

            @_dlang.CALLBACK_TYPE
            def c_callback(done, total, user_data):
                try:
                    return 1 if callback(done, total) else 0
                except Exception:
                    return 1  # stop on exception

            self._callback_ref = c_callback
        else:
            c_callback = _dlang.CALLBACK_TYPE()
            self._callback_ref = c_callback

        run_fn = _lib.ddg_sampler_run_exact if exact else _lib.ddg_sampler_run
        return run_fn(self._handle, moves, c_callback, None)

    @property
    def manifold(self) -> Manifold:
        """Get the current manifold state (copies the manifold)."""
        mfd_handle = _lib.ddg_sampler_get_manifold(self._handle)
        # Create a copy so the Python Manifold owns its handle
        copy_handle = _lib.ddg_manifold_copy(mfd_handle)
        obj = Manifold.__new__(Manifold)
        obj._handle = copy_handle
        return obj

    # -- Direct queries (no copy) --

    @property
    def f_vector(self) -> "np.ndarray":
        """Get the f-vector of the sampler's current manifold (no copy)."""
        import numpy as np
        buf = (ctypes.c_long * 10)()
        n = _lib.ddg_sampler_f_vector(self._handle, buf, 10)
        return np.array(buf[:n], dtype=np.int64)

    def importance_weight(self) -> float:
        """Get the importance weight of the sampler's current manifold (no copy)."""
        return _lib.ddg_sampler_importance_weight(self._handle)

    def simplices(self, dim: int) -> "np.ndarray":
        """Get simplices of given dimension from the sampler's current manifold (no copy)."""
        import numpy as np
        count = _lib.ddg_sampler_simplices(self._handle, dim, None)
        if count == 0:
            return np.empty((0, dim + 1), dtype=np.intc)
        buf = (ctypes.c_int * (count * (dim + 1)))()
        _lib.ddg_sampler_simplices(self._handle, dim, buf)
        return np.frombuffer(buf, dtype=np.intc).reshape(count, dim + 1).copy()

    def degree(self, simplex) -> int:
        """Get the degree of a simplex in the sampler's current manifold (no copy)."""
        import numpy as np
        arr = np.asarray(simplex, dtype=np.intc).ravel()
        ptr = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        return _lib.ddg_sampler_degree(self._handle, ptr, len(arr))
