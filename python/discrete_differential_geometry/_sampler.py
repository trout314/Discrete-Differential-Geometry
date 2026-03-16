"""Pythonic wrapper for the MCMC manifold sampler."""

from __future__ import annotations

import ctypes
import sys
import time
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
    """

    num_facets_target: int = 100
    hinge_degree_target: float = 4.5
    num_facets_coef: float = 0.1
    num_hinges_coef: float = 0.05
    hinge_degree_variance_coef: float = 0.2
    codim3_degree_variance_coef: float = 0.1


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
        callback_interval: int | None = None,
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
            Return ``True`` to stop early. The callback can access
            ``self.manifold``, ``self.current_objective``, etc. to
            inspect the current state of the sampler.
        callback_interval : int, optional
            Number of moves between callback invocations. Default 1000.
            Increase for large triangulations to reduce overhead.

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

        if callback_interval is not None:
            _lib.ddg_sampler_set_callback_interval(self._handle, callback_interval)

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

    def run_with_display(
        self,
        *,
        sweeps: float | None = None,
        moves: int | None = None,
        exact: bool = False,
        update_seconds: float = 0.5,
        callback_interval: int = 10000,
        histogram_dim: int | None = None,
        histogram_width: int = 40,
    ) -> int:
        """Run the sampler with a live terminal display.

        Parameters
        ----------
        sweeps, moves, exact : same as ``run()``.
        update_seconds : float
            Minimum seconds between display updates.
        callback_interval : int
            Moves between D-side callbacks. Default 10000.
        histogram_dim : int, optional
            If given, show a degree histogram for simplices of this dimension.
            Typically 0 (vertex degrees) or dim-2 (hinge degrees).
        histogram_width : int
            Character width of histogram bars.

        Returns
        -------
        int
            Number of moves accepted.
        """
        t_start = time.monotonic()
        last_display = 0.0
        stats_before = self.get_stats()
        dim = self.manifold.dimension

        def _display(done, total):
            nonlocal last_display
            now = time.monotonic()
            if now - last_display < update_seconds:
                return False
            last_display = now

            elapsed = now - t_start
            mfd = self.manifold
            stats = self.get_stats()
            tried = stats.total_tried - stats_before.total_tried
            accepted = stats.total_accepted - stats_before.total_accepted
            accept_pct = 100.0 * accepted / tried if tried > 0 else 0.0

            # Progress line
            pct = 100.0 * done / total if total > 0 else 0.0
            eta = ""
            if done > 0:
                eta_secs = elapsed * (total - done) / done
                if eta_secs < 120:
                    eta = f"  ETA {eta_secs:.0f}s"
                else:
                    eta = f"  ETA {eta_secs / 60:.1f}m"

            lines = []
            lines.append(
                f"  {done:,}/{total:,} moves ({pct:.1f}%)  "
                f"{elapsed:.1f}s elapsed{eta}"
            )

            # State line
            nf = mfd.num_facets
            obj = self.current_objective
            dv = mfd.degree_variance(0)
            lines.append(
                f"  facets={nf:,}  obj={obj:.1f}  "
                f"vtx_deg_var={dv:.1f}  accept={accept_pct:.1f}%"
            )

            # Acceptance breakdown
            bt = stats.bistellar_tries - stats_before.bistellar_tries
            ba = stats.bistellar_accepts - stats_before.bistellar_accepts
            parts = []
            for i in range(dim + 1):
                if bt[i] > 0:
                    r = 100.0 * ba[i] / bt[i]
                    parts.append(f"{i+1}->{dim+2-i-1}:{r:.0f}%")
            move_str = "  " + "  ".join(parts)

            ht = stats.hinge_tries - stats_before.hinge_tries
            ha = stats.hinge_accepts - stats_before.hinge_accepts
            if ht > 0:
                move_str += f"  hinge:{100.0 * ha / ht:.0f}%"
            lines.append(move_str)

            # Optional histogram
            if histogram_dim is not None:
                h = mfd.degree_histogram(histogram_dim)
                if len(h) > 0:
                    max_count = max(h) if max(h) > 0 else 1
                    lines.append(f"  deg histogram (dim={histogram_dim}):")
                    for i, count in enumerate(h):
                        if count == 0:
                            continue
                        deg = i + 1
                        bar_len = int(histogram_width * count / max_count)
                        bar = "\u2588" * bar_len
                        lines.append(f"    {deg:3d}: {bar} {count}")

            # Clear and print
            output = "\n".join(lines)
            # Move cursor up to overwrite previous display
            if hasattr(_display, "_prev_lines"):
                sys.stdout.write(f"\033[{_display._prev_lines}A\033[J")
            sys.stdout.write(output + "\n")
            sys.stdout.flush()
            _display._prev_lines = len(lines)

            return False

        # Print blank lines to reserve space for first update
        print()
        _display._prev_lines = 1

        result = self.run(
            sweeps=sweeps,
            moves=moves,
            exact=exact,
            callback=_display,
            callback_interval=callback_interval,
        )

        # Final display
        _display(1, 1)
        return result

    # -- Parameter setters --

    def set_num_facets_target(self, target: int) -> None:
        """Update the target number of facets."""
        _lib.ddg_sampler_set_num_facets_target(self._handle, target)
        self._params.num_facets_target = target

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
