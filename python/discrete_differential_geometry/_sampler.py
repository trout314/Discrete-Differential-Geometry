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
    hinge_degree_target_coef : float
        Coefficient for the fixed-target hinge penalty
        ``sum_e (deg_e - hinge_degree_target)^2`` (extensive, strictly local).
        0 = off. Unlike the variance penalties, this is a per-hinge coupling
        on an extensive sum, not a coefficient on a mean variance.
    codim3_degree_target_coef : float
        Coefficient for the fixed-target codimension-3 penalty
        ``sum_v (deg_v - codim3_degree_target)^2`` (extensive, strictly
        local). 0 = off.
    codim3_degree_target : float
        Constant target for the codim-3 fixed-target penalty. Choose it
        consistently with the pinned f-vector (see
        :func:`vertex_degree_target` for dimension 3), else it fights the
        edge pin.
    """

    num_facets_target: int = 100
    hinge_degree_target: float = 4.5
    num_facets_coef: float = 0.1
    num_hinges_coef: float = 0.05
    hinge_degree_variance_coef: float = 0.2
    codim3_degree_variance_coef: float = 0.1
    hinge_degree_target_coef: float = 0.0
    codim3_degree_target_coef: float = 0.0
    codim3_degree_target: float = 0.0


def vertex_degree_target(edge_degree_target: float) -> float:
    """Vertex-degree target consistent with an edge-degree target (dim 3).

    For S^3 triangulations the Dehn-Sommerville relations (f2 = 2*f3,
    f0 = f1 - f3) pin the mean vertex degree once the mean edge degree is
    fixed: Dbar = 4/(6/dbar - 1). Using any other codim3_degree_target makes
    the fixed-target vertex penalty fight the edge pin.
    """
    return 4.0 / (6.0 / edge_degree_target - 1.0)


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
        # Fixed-target penalties are post-create setters so the C create
        # signature (and every existing caller) stays unchanged.
        if params.hinge_degree_target_coef:
            _lib.ddg_sampler_set_hinge_degree_target_coef(
                self._handle, params.hinge_degree_target_coef)
        if params.codim3_degree_target:
            _lib.ddg_sampler_set_codim3_degree_target(
                self._handle, params.codim3_degree_target)
        if params.codim3_degree_target_coef:
            _lib.ddg_sampler_set_codim3_degree_target_coef(
                self._handle, params.codim3_degree_target_coef)
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

    def set_hinge_degree_target_coef(self, coef: float) -> None:
        """Update the fixed-target hinge penalty coefficient (0 = off)."""
        _lib.ddg_sampler_set_hinge_degree_target_coef(self._handle, coef)
        self._params.hinge_degree_target_coef = coef

    def set_codim3_degree_target_coef(self, coef: float) -> None:
        """Update the fixed-target codim-3 penalty coefficient (0 = off)."""
        _lib.ddg_sampler_set_codim3_degree_target_coef(self._handle, coef)
        self._params.codim3_degree_target_coef = coef

    def set_codim3_degree_target(self, target: float) -> None:
        """Update the codim-3 fixed-target degree (see vertex_degree_target)."""
        _lib.ddg_sampler_set_codim3_degree_target(self._handle, target)
        self._params.codim3_degree_target = target

    def set_n6_potential(self, zleg_coef: float, imp_coef: float = 0.0,
                         tilt=None) -> None:
        """Configure the vertex 6-valence potential (dim=3 only; 0,0,None = off).

        Per-vertex energy on n6 = #incident edges with degree >= 6 and
        m = #incident edges with degree outside {5, 6}:

            U(n6) = zleg_coef * dist^2(n6, {0,2,3,4}) + tilt[n6]  (tilt: n6 <= 4)
            V(m)  = imp_coef * m^2

        By the link sum rule, zero energy <=> the vertex is exactly a
        Frank-Kasper coordination (Z12/Z14/Z15/Z16). ``tilt`` (length-5
        sequence, default zeros) applies chemical potentials to the legal
        classes: lowering tilt[2] vs tilt[4] favors Z14 (A15-type) over Z16
        (Laves-type) stoichiometry.
        """
        if tilt is None:
            tilt_ptr = None
        else:
            vals = list(tilt)
            if len(vals) != 5:
                raise ValueError("tilt must have length 5")
            tilt_ptr = (ctypes.c_double * 5)(*vals)
        _lib.ddg_sampler_set_n6_potential(
            self._handle, zleg_coef, imp_coef, tilt_ptr)

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

    # -- Per-vertex move-attribution counters (measured combinatorial lapse) --

    def track_move_counts(self, enable: bool = True) -> None:
        """Enable/disable per-vertex move-attribution counters (dim=3 only).

        Off by default (small per-proposal overhead). Every proposed / valid /
        accepted move distributes total weight 1 uniformly over its support
        vertices: the 5 bistellar-ball vertices (= the one 4-simplex the move
        glues on, so 1/5 each) or the 6 support vertices of a 4-4 hinge move.
        Bistellar and hinge accepted ledgers are kept separate so any 4-volume
        convention (e.g. a 4-4 move = 2 stacked 4-simplices) can be applied
        downstream. Enabling does not clear existing counts; see
        :meth:`reset_move_counts`.
        """
        _lib.ddg_sampler_track_move_counts(self._handle, 1 if enable else 0)

    def reset_move_counts(self) -> None:
        """Zero the per-vertex move-attribution counters."""
        _lib.ddg_sampler_reset_move_counts(self._handle)

    def move_counts(self) -> dict[str, np.ndarray]:
        """Return the per-vertex move-attribution counters.

        Returns a dict of equal-length arrays keyed by:
        ``vertex`` (labels, sorted), ``proposed`` (concrete move formed, post
        proposal-thinning, pre validity), ``valid`` (passed the validity check,
        i.e. a counted "try"), ``accepted_bistellar`` and ``accepted_hinge``
        (accepted moves by type). Ledger sums equal event counts. Note a 1-4
        move's created-vertex label is attributed like any other; intersect
        with surviving vertices in analysis.
        """
        n = _lib.ddg_sampler_move_counts(self._handle, None, None, None, None, None)
        if n <= 0:
            empty_i = np.empty(0, dtype=np.int32)
            empty_d = np.empty(0, dtype=np.float64)
            return dict(vertex=empty_i, proposed=empty_d, valid=empty_d,
                        accepted_bistellar=empty_d, accepted_hinge=empty_d)
        labels = (ctypes.c_int * n)()
        prop = (ctypes.c_double * n)()
        valid = (ctypes.c_double * n)()
        acc_b = (ctypes.c_double * n)()
        acc_h = (ctypes.c_double * n)()
        _lib.ddg_sampler_move_counts(self._handle, labels, prop, valid, acc_b, acc_h)
        return dict(
            vertex=np.array(labels[:n], dtype=np.int32),
            proposed=np.array(prop[:n], dtype=np.float64),
            valid=np.array(valid[:n], dtype=np.float64),
            accepted_bistellar=np.array(acc_b[:n], dtype=np.float64),
            accepted_hinge=np.array(acc_h[:n], dtype=np.float64),
        )

    # -- Role-resolved geometry ledger + event log (see move_geometry module) --

    def track_geometry(self, enable: bool = True) -> None:
        """Enable/disable the role-resolved geometry ledger (dim=3 only).

        Records, per vertex and per edge, participation counts in every
        (move type, role) channel of accepted moves, plus tet birth/death
        aggregates and a lifetime histogram. Role taxonomy, degree-change
        tables, and derived fields live in
        :mod:`discrete_differential_geometry.move_geometry`.
        """
        _lib.ddg_sampler_track_geometry(self._handle, 1 if enable else 0)

    def reset_geometry(self) -> None:
        """Zero the geometry ledger (roles, tet aggregates, and clock)."""
        _lib.ddg_sampler_reset_geometry(self._handle)

    def vertex_role_counts(self) -> dict:
        """Per-vertex role counts: dict(vertex[n], counts[n, 11]).

        Columns follow ``move_geometry.VROLE_NAMES``.
        """
        n = _lib.ddg_sampler_vertex_role_counts(self._handle, None, None)
        if n <= 0:
            return dict(vertex=np.empty(0, np.int32), counts=np.empty((0, 11)))
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_double * (n * 11))()
        _lib.ddg_sampler_vertex_role_counts(self._handle, labels, counts)
        return dict(vertex=np.array(labels[:n], np.int32),
                    counts=np.array(counts[:n * 11]).reshape(n, 11))

    def edge_role_counts(self) -> dict:
        """Per-edge role counts: dict(edge[n, 2], counts[n, 15]).

        Columns follow ``move_geometry.EROLE_NAMES``; edges are sorted pairs.
        """
        n = _lib.ddg_sampler_edge_role_counts(self._handle, None, None, None)
        if n <= 0:
            return dict(edge=np.empty((0, 2), np.int32), counts=np.empty((0, 15)))
        la = (ctypes.c_int * n)()
        lb = (ctypes.c_int * n)()
        counts = (ctypes.c_double * (n * 15))()
        _lib.ddg_sampler_edge_role_counts(self._handle, la, lb, counts)
        edge = np.stack([np.array(la[:n], np.int32), np.array(lb[:n], np.int32)], 1)
        return dict(edge=edge, counts=np.array(counts[:n * 15]).reshape(n, 15))

    def tet_stats(self) -> dict:
        """Tet aggregates: created/destroyed by move type, lifetime histogram.

        Returns dict(created[5], destroyed[5], lifetime_hist[64] (log2 bins of
        age in attempted moves), living, censored_deaths, clock). Move type
        codes: 0:1→4, 1:2→3, 2:3→2, 3:4→1, 4:4-4.
        """
        cr = (ctypes.c_long * 5)()
        de = (ctypes.c_long * 5)()
        lh = (ctypes.c_long * 64)()
        living = ctypes.c_long()
        cens = ctypes.c_long()
        clock = ctypes.c_long()
        _lib.ddg_sampler_tet_stats(self._handle, cr, de, lh,
                                   ctypes.byref(living), ctypes.byref(cens),
                                   ctypes.byref(clock))
        return dict(created=np.array(cr[:5], np.int64),
                    destroyed=np.array(de[:5], np.int64),
                    lifetime_hist=np.array(lh[:64], np.int64),
                    living=living.value, censored_deaths=cens.value,
                    clock=clock.value)

    def enable_event_log(self, capacity_mb: float = 16.0) -> None:
        """Enable the accepted-move event log (dim=3 only); 0 disables.

        One fixed-size record per accepted move (see
        ``move_geometry.EVENT_DTYPE``). Drain regularly with
        :meth:`drain_event_log`; if the buffer fills between drains, records
        are DROPPED and :meth:`event_log_overflowed` reports it.
        """
        _lib.ddg_sampler_event_log_enable(self._handle,
                                          int(capacity_mb * 1024 * 1024))

    def drain_event_log(self) -> np.ndarray:
        """Copy out and clear buffered event records (structured array)."""
        from .move_geometry import EVENT_DTYPE
        used = _lib.ddg_sampler_event_log_drain(self._handle, None, 0)
        if used <= 0:
            return np.empty(0, dtype=EVENT_DTYPE)
        buf = (ctypes.c_ubyte * used)()
        got = _lib.ddg_sampler_event_log_drain(self._handle, buf, used)
        return np.frombuffer(bytes(buf[:got]), dtype=EVENT_DTYPE)

    def event_log_overflowed(self) -> bool:
        """True if records were dropped since last check (clears the flag)."""
        return bool(_lib.ddg_sampler_event_log_overflowed(self._handle))

    # -- Six-edge flip log (disclination-network rewiring stream) --

    def enable_six_flip_log(self, capacity_mb: float = 16.0) -> None:
        """Enable the six-edge flip log (dim=3 only); 0 disables.

        One fixed-size record (``disclination.SIX_FLIP_DTYPE``) per edge
        crossing the degree 5<->6 threshold in an accepted move — the
        complete rewiring history of the disclination network, on the same
        clock as the move event log. Rates run a few records per accepted
        move; drain every few sweeps.
        """
        _lib.ddg_sampler_six_flip_log_enable(self._handle,
                                             int(capacity_mb * 1024 * 1024))

    def drain_six_flip_log(self) -> np.ndarray:
        """Copy out and clear buffered flip records (structured array)."""
        from .disclination import SIX_FLIP_DTYPE
        used = _lib.ddg_sampler_six_flip_log_drain(self._handle, None, 0)
        if used <= 0:
            return np.empty(0, dtype=SIX_FLIP_DTYPE)
        buf = (ctypes.c_ubyte * used)()
        got = _lib.ddg_sampler_six_flip_log_drain(self._handle, buf, used)
        return np.frombuffer(bytes(buf[:got]), dtype=SIX_FLIP_DTYPE)

    def six_flip_log_overflowed(self) -> bool:
        """True if flip records were dropped since last check (clears flag)."""
        return bool(_lib.ddg_sampler_six_flip_log_overflowed(self._handle))

    # -- Integer 1-cocycle tracking (T^3 winding forms) --

    def enable_cocycle(self, edges, omega) -> None:
        """Enable integer 1-cocycle tracking (dim=3 only).

        edges: (n, 2) int array of edge labels; omega: (n, 3) int array,
        omega[i] = winding values of edges[i][0] -> edges[i][1]. Must cover
        the manifold's edge set exactly and be closed on every triangle
        (verified in D; raises otherwise). See ``cocycle.py`` for building
        the initial assignment from reference coordinates.
        """
        e = np.ascontiguousarray(np.asarray(edges, dtype=np.intc))
        w = np.ascontiguousarray(np.asarray(omega, dtype=np.intc))
        if e.ndim != 2 or e.shape[1] != 2 or w.shape != (e.shape[0], 3):
            raise ValueError("edges must be (n, 2), omega must be (n, 3)")
        _lib.ddg_sampler_cocycle_enable(
            self._handle,
            e.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            w.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            len(e))

    def disable_cocycle(self) -> None:
        """Disable cocycle tracking and free its state."""
        _lib.ddg_sampler_cocycle_enable(self._handle, None, None, 0)

    def read_cocycle(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (edges (n, 2) sorted u < v, omega (n, 3)) of the tracked
        cocycle (edge order unspecified)."""
        n = _lib.ddg_sampler_cocycle_read(self._handle, None, None, 0)
        e = np.empty((n, 2), dtype=np.intc)
        w = np.empty((n, 3), dtype=np.intc)
        got = _lib.ddg_sampler_cocycle_read(
            self._handle,
            e.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            w.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), n)
        return e[:got].copy(), w[:got].copy()

    def check_cocycle(self) -> None:
        """Audit the cocycle (edge-set match + closedness on every triangle);
        raises RuntimeError on drift. The production integrity check."""
        _lib.ddg_sampler_cocycle_check(self._handle)

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
