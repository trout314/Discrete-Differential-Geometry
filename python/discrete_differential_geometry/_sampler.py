"""Pythonic wrapper for the MCMC manifold sampler."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from . import _dlang
from ._manifold import Manifold
from ._manifold_view import ManifoldView

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
    hinge_move_prob : float
        Probability of proposing a hinge (4-4) move per step (dim=3 only).
    """

    num_facets_target: int = 100
    hinge_degree_target: float = 4.5
    num_facets_coef: float = 0.1
    num_hinges_coef: float = 0.05
    hinge_degree_variance_coef: float = 0.2
    codim3_degree_variance_coef: float = 0.1
    hinge_move_prob: float = 0.0


@dataclass
class SamplerStats:
    """Cumulative MCMC statistics from a sampler."""

    total_tried: int
    total_accepted: int
    hinge_tries: int
    hinge_accepts: int
    bistellar_tries: np.ndarray   # per move type, indexed by coCenter.length - 1
    bistellar_accepts: np.ndarray


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
        self._handle = _lib.ddg_sampler_create_ext(
            manifold._handle,
            params.num_facets_target,
            params.hinge_degree_target,
            params.num_facets_coef,
            params.num_hinges_coef,
            params.hinge_degree_variance_coef,
            params.codim3_degree_variance_coef,
            params.hinge_move_prob,
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
            target distribution exp(-objective). Default False uses pure
            Metropolis (no Hastings); use importance_weight() to correct.
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

    # -- Parameter setters --

    def set_num_facets_target(self, target: int) -> None:
        """Update the target number of facets."""
        _lib.ddg_sampler_set_num_facets_target(self._handle, target)
        self._params.num_facets_target = target

    def set_hinge_move_prob(self, prob: float) -> None:
        """Update the hinge move probability."""
        _lib.ddg_sampler_set_hinge_move_prob(self._handle, prob)
        self._params.hinge_move_prob = prob

    def set_num_facets_coef(self, coef: float) -> None:
        """Update the volume penalty coefficient."""
        _lib.ddg_sampler_set_num_facets_coef(self._handle, coef)
        self._params.num_facets_coef = coef

    def set_num_hinges_coef(self, coef: float) -> None:
        """Update the global curvature penalty coefficient."""
        _lib.ddg_sampler_set_num_hinges_coef(self._handle, coef)
        self._params.num_hinges_coef = coef

    def set_hinge_degree_variance_coef(self, coef: float) -> None:
        """Update the local curvature penalty coefficient."""
        _lib.ddg_sampler_set_hinge_degree_variance_coef(self._handle, coef)
        self._params.hinge_degree_variance_coef = coef

    def set_codim3_degree_variance_coef(self, coef: float) -> None:
        """Update the codimension-3 degree variance coefficient."""
        _lib.ddg_sampler_set_codim3_degree_variance_coef(self._handle, coef)
        self._params.codim3_degree_variance_coef = coef

    def set_hinge_degree_target(self, target: float) -> None:
        """Update the target hinge degree."""
        _lib.ddg_sampler_set_hinge_degree_target(self._handle, target)
        self._params.hinge_degree_target = target

    # -- Statistics --

    def get_stats(self) -> SamplerStats:
        """Return cumulative MCMC statistics."""
        dim = self.manifold.dimension
        n = dim + 1
        tt = ctypes.c_long()
        ta = ctypes.c_long()
        ht = ctypes.c_long()
        ha = ctypes.c_long()
        bt = (ctypes.c_long * n)()
        ba = (ctypes.c_long * n)()
        _lib.ddg_sampler_get_stats(
            self._handle,
            ctypes.byref(tt), ctypes.byref(ta),
            ctypes.byref(ht), ctypes.byref(ha),
            bt, ba, n,
        )
        return SamplerStats(
            total_tried=tt.value,
            total_accepted=ta.value,
            hinge_tries=ht.value,
            hinge_accepts=ha.value,
            bistellar_tries=np.array(bt[:n], dtype=np.int64),
            bistellar_accepts=np.array(ba[:n], dtype=np.int64),
        )

    def reset_stats(self) -> None:
        """Reset cumulative statistics counters."""
        _lib.ddg_sampler_reset_stats(self._handle)

    @property
    def current_objective(self) -> float:
        """Return the current objective function value."""
        return _lib.ddg_sampler_current_objective(self._handle)

    # -- Ramped growth --

    def ramped_grow(
        self,
        target_facets: int,
        step_size: int = 500,
        eq_sweeps_per_step: int = 5,
        callback: Callable[[int, int], None] | None = None,
    ) -> None:
        """Grow the manifold to target_facets via ramped growth with equilibration.

        At each step, the sampler's num_facets_target is increased by step_size,
        MCMC is run until the manifold reaches the step target, then equilibration
        sweeps are run.

        Parameters
        ----------
        target_facets : int
            Final target number of facets.
        step_size : int
            Increase num_facets_target by this much each step.
        eq_sweeps_per_step : int
            Equilibration sweeps to run after reaching each step target.
        callback : callable, optional
            Called with ``(current_facets, step_target)`` after each step completes.
        """
        current = self.manifold.num_facets
        step_target = current

        while step_target < target_facets:
            step_target = min(step_target + step_size, target_facets)
            self.set_num_facets_target(step_target)

            # Run until we reach the step target
            while self.manifold.num_facets < step_target:
                self.run(sweeps=1)

            # Equilibrate
            self.run(sweeps=eq_sweeps_per_step)

            if callback is not None:
                callback(self.manifold.num_facets, step_target)

    @property
    def manifold(self) -> ManifoldView:
        """Read-only view of the sampler's current manifold.

        Returns a lightweight ManifoldView that provides all query methods
        (f_vector, degree, simplices, importance_weight, etc.) but no
        mutation. The view is only valid while this sampler is alive.

        Use ``sampler.manifold.dup()`` to get a mutable, independent copy.
        """
        mfd_handle = _lib.ddg_sampler_get_manifold(self._handle)
        return ManifoldView(mfd_handle)
