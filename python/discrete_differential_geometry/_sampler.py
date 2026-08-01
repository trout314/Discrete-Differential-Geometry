"""Pythonic wrapper for the MCMC manifold sampler."""

from __future__ import annotations

import ctypes
import sys
import time
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from . import _dlang
from ._manifold import Manifold
from ._manifold_view import ManifoldView

_lib = _dlang._lib

#: Slots per chord in the knot-slide proposal: 2 chord orientations x 6
#: ordered ``(c2, c3)`` picks from the chord's sorted 3-vertex link. Constant
#: by construction (an invalid slot is a rejected proposal, not a smaller
#: menu), so the proposal denominator cancels from the Hastings ratio.
SLIDE_SLOTS = 12


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

    def do_bistellar_move(self, center: Sequence[int],
                          cocenter: Sequence[int]) -> None:
        """Apply a SPECIFIED bistellar move through the sampler.

        Unlike ``Manifold.do_bistellar_move`` on a bare manifold, this acts
        on the sampler's internal state and keeps its bookkeeping coherent:
        the forced cocycle update is applied when cocycle tracking is
        enabled, and the tracked objective is invalidated (recomputed lazily,
        along with the n6-potential state). Use this to apply externally
        proposed compound moves (e.g. the knot slide) mid-run.

        Vertex-changing moves are supported: a 4->1 (``center=[v]``,
        ``cocenter=<4 neighbors>``) returns the removed label to the
        sampler's label pool; a 1->4 (``center=<tet>``,
        ``cocenter=[label]``) consumes its caller-chosen label from the
        pool (a fresh label above the pool is simply used). The opt-in
        geometry ledger / event logs do not see these moves. All
        preconditions are validated in the D core (raises on an invalid
        move)."""
        c = np.asarray(center, dtype=np.intc).ravel()
        cc = np.asarray(cocenter, dtype=np.intc).ravel()
        _lib.ddg_sampler_do_bistellar_move(
            self._handle, c.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            len(c), cc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(cc))

    def set_slide_prob(self, prob: float) -> None:
        """Enable the knot-slide move class (dim = 3 only).

        ``prob`` is the probability of proposing a knot slide rather than the
        ordinary 3->2 bistellar move, given that the unified proposal has
        landed on a degree-3 edge. 0 (the default) disables slides.

        Acceptance is plain Metropolis on the exact action change, with no
        Hastings correction (k_f = k_r = 1). Slides run entirely inside the D
        sampler and participate in objective tracking, cocycle updates, the
        geometry ledger and move counters like any other move type. See
        :meth:`set_slide_clean_only` to restrict the class.

        Note this is the probability that a whole MCMC step is a slide
        attempt, not a conditional probability at a degree-3 edge: the slide
        is an independent channel, so a value near 1 starves the thermal
        channel entirely."""
        self._recorded_slide_prob = prob
        _lib.ddg_sampler_set_slide_prob(self._handle, float(prob))

    def set_slide_clean_only(self, clean_only: bool) -> None:
        """Restrict the slide move class to CLEAN (species-preserving) slides.

        Off by default. Cleanliness is not needed for detailed balance -- the
        proposal is symmetric and k_f = k_r = 1 whether or not the slide
        preserves the illegal-degree multiset (verified on 124 transitions,
        66 of them dirty) -- so by default every valid slide is proposed and
        dirty ones are turned away by the Metropolis test on their real
        energy cost rather than excluded by fiat.

        Turning it ON selects a subclass with a special analytic property: a
        clean slide preserves the ENTIRE edge-degree multiset, so under a
        volume pin plus an edge-degree term alone its dS vanishes identically.
        That makes the clean class an exactly zero-energy orbit -- every
        proposal accepted, defect species conserved -- which is useful for
        studying pure transport at fixed species, and degenerate if you want
        the defect population to equilibrate."""
        _lib.ddg_sampler_set_slide_clean_only(
            self._handle, 1 if clean_only else 0)

    def slide_stats(self) -> tuple[int, int]:
        """``(tries, accepts)`` for the knot-slide move class. ``tries``
        counts proposals that formed a legal clean slide (i.e. reached the
        Metropolis test), not raw slot draws."""
        t = ctypes.c_long()
        a = ctypes.c_long()
        _lib.ddg_sampler_slide_stats(
            self._handle, ctypes.byref(t), ctypes.byref(a))
        return int(t.value), int(a.value)

    def set_worm_prob(self, prob: float) -> None:
        """Enable the deg-4 worm channel (dim=3): each step proposes the
        catalysed 2-move transport move with this probability, anchored
        uniformly on the live deg-4 edge set.  Inert while a cocycle is
        attached.  0 disables (default)."""
        self._recorded_worm_prob = prob
        _lib.ddg_sampler_set_worm_prob(self._handle, float(prob))

    def worm_stats(self) -> tuple[int, int, int]:
        """(tries, accepts, no_candidate_rejections) for the worm channel."""
        t = ctypes.c_long()
        a = ctypes.c_long()
        n = ctypes.c_long()
        _lib.ddg_sampler_worm_stats(self._handle, ctypes.byref(t),
                                    ctypes.byref(a), ctypes.byref(n))
        return t.value, a.value, n.value

    def worm_enum(self, a: int, b: int, with_ds: bool = True):
        """Enumerate worm candidates at deg-4 edge (a, b): returns
        (landings[n,2], dS[n]) -- the crossval window into the D move."""
        cap = 512
        la = (ctypes.c_int * cap)()
        lb = (ctypes.c_int * cap)()
        ds = (ctypes.c_double * cap)() if with_ds else None
        n = _lib.ddg_sampler_worm_at(
            self._handle, int(a), int(b), 0, -1, la, lb, ds, cap)
        land = np.array([[la[i], lb[i]] for i in range(n)], dtype=np.int64)
        dsa = (np.array([ds[i] for i in range(n)]) if with_ds
               else np.full(int(n), np.nan))
        return land, dsa

    def worm_at(self, a: int, b: int, cand: int, commit: bool = False):
        """Trial (or commit) worm candidate `cand` at anchor (a, b);
        returns its exact dS, or None if the candidate is invalid."""
        ds = (ctypes.c_double * 1)()
        r = _lib.ddg_sampler_worm_at(
            self._handle, int(a), int(b), 2 if commit else 1, int(cand),
            None, None, ds, 1)
        return float(ds[0]) if r != 0 or not commit else None

    def set_worm_f0(self, utab: dict, ufb, ufbc: float, z0: float,
                    lmax: int = 4096, zeta: float = 0.05,
                    aof: float = 0.5, ph: float = 0.45,
                    pg: float = 0.49, bcf: float = 0.01,
                    bc4: float = 0.05, maxstep: int = 100000,
                    ucap_hi: float = 35.0,
                    ucap_lo: float = -20.0, mu: float = 1.5,
                    f0_ref: "tuple | None" = None) -> None:
        """Configure the f0 worm channel (scheme C, design doc 3.2).

        ``utab`` maps spoke-degree multisets (tuples) to umbrella
        values -- either plain floats (frozen table), or triples
        ``(cum_dS, df1, df3)`` for the f-ADAPTIVE tube: ``cum_dS`` the
        cumulative action change measured at build time, ``df1``/``df3``
        the corridor state's exact (f1, f3) offset from the corridor
        start, and ``f0_ref = (f1, f3)`` the build-time f-vector
        reference (required with triples). The D engine then recompiles
        the frozen scalar table at each episode open from the current
        f-vector -- the global-pin part of every tube value self-adjusts
        across the whole f0 range, so one tube build serves many
        commits. Keys are packed into degree-bucket counts (3..9+, 8
        bits each) -- distinct multisets with any degree above 9 may
        collide, resolved by min value. ``ufb`` is the 6-coefficient
        linear fallback (n3, n4, n5, n6, n7plus, Zdeficit) with constant
        ``ufbc`` and Z reference ``z0``. Mixture weights: ``aof`` =
        open-flag share of openings; ``ph``/``pg``/``bcf``/``bc4`` =
        head / global-repair / close-flag / close-41 shares of open
        steps. ``maxstep`` caps episodes (capped walks are exactly
        undone in the D core)."""
        packed = {}
        adaptive = False
        for ms, u in utab.items():
            k = 0
            for d in ms:
                k += 1 << (8 * min(max(int(d) - 3, 0), 6))
            if isinstance(u, (tuple, list)):
                trip = (float(u[0]), float(u[1]), float(u[2]))
                adaptive = True
            else:
                trip = (float(u), 0.0, 0.0)
            if k not in packed or trip[0] < packed[k][0]:
                packed[k] = trip
        if adaptive and f0_ref is None:
            raise ValueError("f0_ref=(f1, f3) is required with "
                             "(cum, df1, df3) table values")
        tf1, tf3 = (float(f0_ref[0]), float(f0_ref[1])) \
            if f0_ref is not None else (0.0, 0.0)
        ks = sorted(packed)
        keys = np.array(ks, dtype=np.uint64)
        vals = np.array([packed[k][0] for k in ks], dtype=np.float64)
        df1 = np.array([packed[k][1] for k in ks], dtype=np.float64)
        df3 = np.array([packed[k][2] for k in ks], dtype=np.float64)
        as_p = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fb = np.asarray(ufb, dtype=np.float64)
        assert fb.shape == (6,)
        _lib.ddg_sampler_worm_f0_config(
            self._handle,
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            as_p(vals),
            as_p(df1) if adaptive else None,
            as_p(df3) if adaptive else None,
            tf1, tf3, len(keys),
            fb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            float(ufbc), float(z0), int(lmax), float(zeta),
            float(aof), float(ph), float(pg), float(bcf), float(bc4),
            int(maxstep), float(ucap_hi), float(ucap_lo), float(mu))

    def worm_f0_episode(self) -> dict:
        """Run one f0-worm episode; returns a diagnostics dict. The
        tracked objective stays exact across episodes."""
        out = np.zeros(12, dtype=np.float64)
        changed = _lib.ddg_sampler_worm_f0_episode(
            self._handle,
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return {"changed": bool(changed), "opened": int(out[0]),
                "head": int(out[1]), "steps": int(out[2]),
                "closed": {0: None, 1: "cf", 2: "c4", 3: "undone"}[
                    int(out[3])],
                "dS": float(out[4]), "umax": float(out[5]),
                "nH": int(out[6]), "accH": int(out[7]),
                "nG": int(out[8]), "accG": int(out[9]),
                "zmin": int(out[10]), "nZ4": int(out[11]),
                # diagnostics for the BEST close attempt of the
                # episode (the pair has no f census of its own, so
                # res.df carries them): dS41*100, uRawKeep*100,
                # log(Wc/wc)*100, z of the kept ball. "umax" is that
                # attempt's log-alpha, NOT the open's U sum.
                "cdiag": tuple(int(out[12 + k]) for k in range(4))}

    def set_worm_pair(self, zeta2: float = 1.0, bcp: float = 0.05,
                      chain_k: int = 20) -> None:
        """Configure the BILOCAL (two-ball) sector: ``zeta2`` the pair
        fugacity (pass ``float("nan")`` to auto-calibrate it from the
        proposal density -- for the strict chord channel the open makes
        no move, so the balanced value is just minus the log proposal
        density and needs no probe). ``bcp`` is the close-pair share --
        but under auto-zeta2 it is reinterpreted: p_close is derived as
        maxstep/(1/bcp), i.e. the mean episode length becomes
        ``maxstep * bcp`` and abandonment falls to ~e^(-1/bcp). p_close
        cannot be tuned independently once zeta2 is auto, because it
        already sits inside zeta2* = -(log n3 + log nSite + log
        p_close); the two are derived together. The
        umbrella, caps and seed bias come from :meth:`set_worm_f0`,
        which must be called first."""
        _lib.ddg_sampler_worm_pair_config(
            self._handle, float(zeta2), float(bcp), int(chain_k))

    def worm_pair_episode(self) -> dict:
        """Run one bilocal episode: a vertex is created at one ball and
        destroyed at the other, so the net f-change vanishes EXACTLY
        (both closures are f-neutral -- this channel transports a
        vertex, it cannot change f0) and the global pins cost exactly
        zero (measured 0.0e+00 at every separation).

        The umbrella is ROLE-SIGNED: the created ball carries -U and the
        adopted ball +U, so one star grows while the other collapses and
        the U content is conserved open-to-close; the close prices both
        balls under the REVERSE episode's roles. ``zeta2`` must be a
        calibrated finite value (see f0_worm.calib_zeta2) -- auto is
        chord-only -- and ``bcp`` IS p_close, so it sets the mean
        episode length. ``closed`` is "transport" (the adopted vertex
        was deleted, so a vertex moved) or "roundtrip" (the created one
        was deleted, committing only the walk's corridor work)."""
        out = np.zeros(16, dtype=np.float64)
        changed = _lib.ddg_sampler_worm_pair_episode(
            self._handle,
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return {"changed": bool(changed), "opened": int(out[0]),
                "head": int(out[1]), "steps": int(out[2]),
                "closed": {0: None, 3: "undone", 4: "transport",
                           6: "roundtrip"}.get(int(out[3])),
                "dS": float(out[4]), "umax": float(out[5]),
                "nH": int(out[6]), "accH": int(out[7]),
                "nG": int(out[8]), "accG": int(out[9]),
                "zmin": int(out[10]), "nZ4": int(out[11]),
                # diagnostics for the BEST close attempt of the
                # episode (the pair has no f census of its own, so
                # res.df carries them): dS41*100, uRawKeep*100,
                # log(Wc/wc)*100, z of the kept ball. "umax" is that
                # attempt's log-alpha, NOT the open's U sum.
                "cdiag": tuple(int(out[12 + k]) for k in range(4))}

    def worm_chord_strict_episode(self) -> dict:
        """One STRICT-CLOSURE chord episode. Both marks are pure flags,
        so the sector boundary makes no move and the whole f-change
        comes from the walk. The close fires only when one mark is a
        degree-3 chord, the other is ABSENT, and the absent one's region
        carries no degree-3 edge -- the flicker relocated AND its old
        neighbourhood came out clean, i.e. a reaction rather than a bare
        relocation. ``closed == "strict"`` marks a committed one."""
        out = np.zeros(16, dtype=np.float64)
        changed = _lib.ddg_sampler_worm_chord_strict(
            self._handle,
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return {"changed": bool(changed), "opened": int(out[0]),
                "steps": int(out[2]),
                "closed": {0: None, 3: "undone", 7: "strict"}.get(
                    int(out[3])),
                "dS": float(out[4]), "nH": int(out[6]),
                "accH": int(out[7]), "nG": int(out[8]),
                "accG": int(out[9]),
                "df": tuple(int(out[12 + k]) for k in range(4))}

    def set_worm_chord(self, ctab: dict, offpen: float = 0.0) -> None:
        """Upload the CHORD carrier's umbrella, replayed from measured
        catalysed paths. ``ctab`` maps the chord's local signature --
        the sorted degrees of every edge at either endpoint -- to the
        cumulative dS along the path. Empty dict disables it (U = 0)."""
        packed = {}
        for ms, u in ctab.items():
            # ms = (region deg-3 count, *endpoint edge degrees). The
            # count rides in the top byte, which wf0Key leaves free
            # (it uses buckets 0..6, bits 0..55) -- see wf0ChordKey.
            k = (int(ms[0]) & 0xFF) << 56
            for d in ms[1:]:
                k += 1 << (8 * min(max(int(d) - 3, 0), 6))
            packed[k] = min(packed.get(k, float(u)), float(u))
        ks = sorted(packed)
        keys = np.array(ks, dtype=np.uint64)
        vals = np.array([packed[k] for k in ks], dtype=np.float64)
        _lib.ddg_sampler_worm_chord_config(
            self._handle,
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            vals.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            len(keys), float(offpen))

    def worm_chord_episode(self) -> dict:
        """Run one CHORD (2<->3) bilocal episode -- the flicker carrier.

        A chord is created at one ball and annihilated at the other, so
        the net f-change vanishes and the pins cost nothing. Told apart
        by the closure: ``closed == "chord"`` is a committed transport
        (the flicker relocated); an episode that does work and then
        annihilates the chord it made is catalysis. No umbrella is used
        -- a chord's closure condition (degree exactly 3) is what a 2->3
        creates, so there is no barrier to flatten. ``set_worm_pair``'s
        ``zeta2`` acts as the LOG chemical potential for this carrier."""
        out = np.zeros(12, dtype=np.float64)
        changed = _lib.ddg_sampler_worm_chord_episode(
            self._handle,
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return {"changed": bool(changed), "opened": int(out[0]),
                "head": int(out[1]), "steps": int(out[2]),
                "closed": {0: None, 3: "undone", 5: "chord"}.get(
                    int(out[3])),
                "dS": float(out[4]), "nH": int(out[6]),
                "accH": int(out[7]), "nG": int(out[8]),
                "accG": int(out[9]),
                "df": tuple(int(out[12 + k]) for k in range(4))}

    def set_nonlocal_slide_prob(self, prob: float, max_step: int = 8) -> None:
        """Enable the NON-LOCAL slide channel (dim = 3 only).

        Each MCMC step proposes it with probability ``prob``: pick a degree-3
        chord uniformly (1/n_3) from the live defect set, then annihilate it
        (3->2) and re-create it a uniform 1..``max_step`` tets down its BC
        chain (2->3). Direction is the slot orientation, so the step
        distribution is symmetric about 0; acceptance is Metropolis on the
        exact dS with the 1/n_3 Hastings factor n_3/(n_3 + dn_3).

        This is an EQUILIBRIUM-sampling move -- it teleports the excitation
        across same-rung sites, defeating the washboard caging. ``max_step`` =
        4 recovers the local slide's displacement (the physical-kinetics
        limit). 0 disables the channel (the default). While enabled, the D
        sampler keeps a live degree-3 chord set (rebuilt at each run start,
        maintained incrementally); it costs nothing when disabled."""
        self._recorded_nonlocal_slide = (prob, max_step)
        _lib.ddg_sampler_set_nonlocal_slide_prob(
            self._handle, float(prob), int(max_step))

    def nonlocal_slide_stats(self) -> tuple[int, int]:
        """``(tries, accepts)`` for the non-local slide channel. ``tries``
        counts proposals that formed a legal move (reached Metropolis)."""
        t = ctypes.c_long()
        a = ctypes.c_long()
        _lib.ddg_sampler_nonlocal_slide_stats(
            self._handle, ctypes.byref(t), ctypes.byref(a))
        return int(t.value), int(a.value)

    def deg3_count(self) -> int:
        """Size of the live degree-3 chord set (diagnostic). Populated only
        while the non-local channel is enabled and a run has started."""
        return int(_lib.ddg_sampler_deg3_count(self._handle))

    def slide_at(self, a: int, b: int, slot: int,
                 commit: bool = False) -> "float | None":
        """Attempt the knot slide at chord ``(a, b)`` in slot ``slot``.

        Scripted / crossval entry into the same D code path the sampler's
        slide move uses -- NOT a sampling path (it bypasses the Metropolis
        test). Returns the exact action change if the slot yields a legal
        clean slide, or ``None`` if it does not. With ``commit=False`` the
        state is restored exactly; with ``commit=True`` the slide is applied
        unconditionally and all bookkeeping advances as for an accepted move.

        There are :data:`SLIDE_SLOTS` = 12 slots per chord: 2 chord
        orientations x 6 ordered picks of ``(c2, c3)`` from the chord's sorted
        3-vertex link."""
        ds = ctypes.c_double()
        rc = _lib.ddg_sampler_slide_at(
            self._handle, int(a), int(b), int(slot),
            1 if commit else 0, ctypes.byref(ds))
        return float(ds.value) if rc == 1 else None

    def slide_at2(self, a: int, b: int, slot: int, commit: bool = False):
        """Like :meth:`slide_at`, also returning the slot's ARRIVAL chord
        from the frame decode: (dS_or_None, (c4, c8) or None). The arrival
        is reported whenever the frame derives, even for rejected slides --
        drivers use it to identify exact inverses (FPKMC walk-back)."""
        ds = ctypes.c_double()
        c4 = ctypes.c_int(-1)
        c8 = ctypes.c_int(-1)
        rc = _lib.ddg_sampler_slide_at2(
            self._handle, int(a), int(b), int(slot),
            1 if commit else 0, ctypes.byref(ds),
            ctypes.byref(c4), ctypes.byref(c8))
        arr = (int(c4.value), int(c8.value)) if c4.value >= 0 else None
        return (float(ds.value) if rc == 1 else None), arr

    def nonlocal_slide_at(self, a: int, b: int, slot: int, steps: int,
                          commit: bool = False):
        """Non-local slide: annihilate the degree-3 chord ``(a, b)`` with a
        3->2 and re-create it ``steps`` tets down the BC chain with a 2->3.

        ``slot`` (0..11) selects the walk (orientation x link order). Returns
        ``(dS, dn3, arrival_chord)`` if the slot yields a legal move, else
        ``None`` -- ``dS`` the exact action change, ``dn3`` the exact change in
        the degree-3 edge count (for the 1/n_3 Metropolis factor
        n3/(n3+dn3)), ``arrival_chord`` the (sorted) target chord. With
        ``commit=False`` the state is restored exactly; with ``commit=True``
        the move is applied unconditionally. ``steps=4`` reproduces the local
        slide's displacement; larger steps are the non-local sampling move.
        See sampler.tryNonlocalSlide. Equilibrium-sampling move -- for physical
        kinetics keep steps small (=4)."""
        ds = ctypes.c_double()
        dn3 = ctypes.c_long()
        ta = ctypes.c_int(-1)
        tb = ctypes.c_int(-1)
        rc = _lib.ddg_sampler_nonlocal_slide_at(
            self._handle, int(a), int(b), int(slot), int(steps),
            1 if commit else 0, ctypes.byref(ds), ctypes.byref(dn3),
            ctypes.byref(ta), ctypes.byref(tb))
        if rc != 1:
            return None
        return (float(ds.value), int(dn3.value),
                tuple(sorted((int(ta.value), int(tb.value)))))

    def site_survey(self, chain) -> dict:
        """Washboard site survey along a BC chain (notes/FPKMC_DESIGN.md).

        For each window k (chain[k:k+5]) the D core creates the knot through
        the ordinary speculative-delta path, measures every slide slot in
        trial-only mode, undoes the creation, and audits exact restoration
        of manifold, potential state and objective.

        Returns dict of arrays over the ``n = len(chain) - 4`` windows:
          ``dS_create``  (n,)      creation action delta; NaN = invalid site
          ``slot_dS``    (n, 12)   per-slot trial dS; NaN = not a legal slide
          ``slot_dest``  (n, 12, 2) per-slot destination chord (c4, c8);
                                    -1 where the frame derivation failed
          ``slot_clean`` (n, 12)   1.0 species-preserving, 0.0 dirty,
                                    NaN illegal
        Destinations vs chain arithmetic classify each slot as
        chain-forward / chain-backward / off-chain (the slot census).
        MEASURED on R m4: 1 clean chain slot each direction (I1 exact)
        plus 2-3 clean OFF-CHAIN slots per site -- the clean slide network
        is a branching graph (see notes/FPKMC_DESIGN.md R3)."""
        ch = np.ascontiguousarray(np.asarray(chain, dtype=np.intc).ravel())
        if ch.size < 5:
            raise ValueError("chain must have at least 5 vertices")
        stride = 1 + 4 * SLIDE_SLOTS
        nwin = ch.size - 4
        out = np.empty(nwin * stride, dtype=np.float64)
        got = _lib.ddg_sampler_site_survey(
            self._handle, ch.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            int(ch.size), out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        assert got == nwin, f"survey returned {got} windows, expected {nwin}"
        rows = out.reshape(nwin, stride)
        slots = rows[:, 1:].reshape(nwin, SLIDE_SLOTS, 4)
        dest = slots[:, :, 1:3]
        dest = np.where(np.isnan(dest), -1, dest).astype(int)
        return {"dS_create": rows[:, 0].copy(),
                "slot_dS": slots[:, :, 0].copy(),
                "slot_dest": dest,
                "slot_clean": slots[:, :, 3].copy()}

    def slide_graph_scan(self, root_chord, dS_max: float = 20.0,
                         max_depth: int = 6, max_nodes: int = 20000,
                         max_edges: int = 200000, blocked_verts=None) -> dict:
        """Bounded DFS over defect states reachable by LEGAL slides (dirty
        included) from the current state, rooted at the degree-3 chord
        ``root_chord``. Exact per-edge dS; exact rollback (audited in D).
        Nodes are keyed exactly by the edge-degree overlay vs the current
        state. Expansion prunes at cum dS > dS_max; capped scans raise.

        FP mode: ``blocked_verts`` (another defect's vertex labels) adds a
        per-node ``dock`` flag -- 1 iff the node's knot complex touches the
        one-tet-layer neighborhood of the blocked set (computed once at
        the root state). Dock nodes are enumerated but never expanded
        (absorbing). Without blocked_verts the scan is unchanged and
        ``dock`` is all zeros.

        Complete-interior guarantee (FP): a node with depth < max_depth,
        dS <= dS_max and dock == 0 was fully expanded -- every legal slide
        out of it appears in the edge list. All other nodes are frontier.

        Returns dict: nodes (dS, depth, n_chords, sig, chord, dock) and
        edges (src, dst, dS, chord, slot). See notes/FPKMC_DESIGN.md M2."""
        a, b = int(root_chord[0]), int(root_chord[1])
        nd = ctypes.POINTER(ctypes.c_double)
        ni = ctypes.POINTER(ctypes.c_int)
        node_dS = np.empty(max_nodes, np.float64)
        node_depth = np.empty(max_nodes, np.intc)
        node_nch = np.empty(max_nodes, np.intc)
        node_sig = np.empty(max_nodes, np.int64)
        node_ca = np.empty(max_nodes, np.intc)
        node_cb = np.empty(max_nodes, np.intc)
        node_dock = np.zeros(max_nodes, np.intc)
        e_src = np.empty(max_edges, np.intc)
        e_dst = np.empty(max_edges, np.intc)
        e_dS = np.empty(max_edges, np.float64)
        e_ca = np.empty(max_edges, np.intc)
        e_cb = np.empty(max_edges, np.intc)
        e_slot = np.empty(max_edges, np.intc)
        n_edges = ctypes.c_long()
        if blocked_verts is not None and len(blocked_verts):
            bv = np.ascontiguousarray(
                [int(v) for v in blocked_verts], dtype=np.intc)
            bv_ptr, bv_n = bv.ctypes.data_as(ni), len(bv)
        else:
            bv_ptr, bv_n = None, 0
        n = _lib.ddg_sampler_slide_graph_scan(
            self._handle, a, b, float(dS_max), int(max_depth),
            int(max_nodes), int(max_edges),
            node_dS.ctypes.data_as(nd), node_depth.ctypes.data_as(ni),
            node_nch.ctypes.data_as(ni),
            node_sig.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
            node_ca.ctypes.data_as(ni), node_cb.ctypes.data_as(ni),
            e_src.ctypes.data_as(ni), e_dst.ctypes.data_as(ni),
            e_dS.ctypes.data_as(nd), e_ca.ctypes.data_as(ni),
            e_cb.ctypes.data_as(ni), e_slot.ctypes.data_as(ni),
            ctypes.byref(n_edges),
            bv_ptr, bv_n, node_dock.ctypes.data_as(ni))
        if n < 0:
            err = _lib.ddg_last_error()
            raise RuntimeError(err.decode() if err else "graph scan failed")
        err = _lib.ddg_last_error()
        if err and b"CAPPED" in err:
            raise RuntimeError(err.decode())
        m = int(n_edges.value)
        return {"dS": node_dS[:n].copy(), "depth": node_depth[:n].copy(),
                "n_chords": node_nch[:n].copy(), "sig": node_sig[:n].copy(),
                "chord": np.stack([node_ca[:n], node_cb[:n]], 1),
                "dock": node_dock[:n].copy(),
                "edge_src": e_src[:m].copy(), "edge_dst": e_dst[:m].copy(),
                "edge_dS": e_dS[:m].copy(),
                "edge_chord": np.stack([e_ca[:m], e_cb[:m]], 1),
                "edge_slot": e_slot[:m].copy()}

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
        self._recorded_n6_potential = {
            "zleg_coef": float(zleg_coef), "imp_coef": float(imp_coef),
            "tilt": [float(x) for x in (tilt if tilt is not None else [0.0] * 5)]}
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

    def enable_cocycle_positions(self, box) -> None:
        """Enable the incrementally-maintained vertex lift (needs a cocycle).

        Integrates the cocycle over a spanning forest ONCE to seed a position
        per vertex, then keeps it exact move-by-move at O(1). This replaces
        re-deriving positions from the whole cochain on every sample
        (``cocycle.tree_positions``), which is O(V+E) per call and rebuilds
        the spanning tree each time -- so it can reassign a persisting vertex,
        the "gauge glitch" that consumers otherwise detect and discard. Here
        the gauge is fixed for the lifetime of the run, so glitches cannot
        arise and displacement needs no filtering.

        ``box`` is the torus period per axis, in the same units as the cocycle
        values (e.g. ``scale * mcell`` for a cocycle built by
        ``cocycle.build_from_positions``)."""
        b = np.ascontiguousarray(
            np.broadcast_to(np.asarray(box, dtype=np.intc), (3,)))
        _lib.ddg_sampler_cocycle_enable_positions(
            self._handle, b.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))

    def vertex_positions(self, verts) -> np.ndarray:
        """Lift positions of ``verts`` as an (n, 3) int array, in the order
        given. O(1) per vertex -- query only what you need (a chord is two
        lookups) rather than marshalling the whole cochain."""
        v = np.ascontiguousarray(np.asarray(verts, dtype=np.intc).ravel())
        out = np.empty((len(v), 3), dtype=np.intc)
        _lib.ddg_sampler_cocycle_positions(
            self._handle, v.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            len(v), out.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        return out

    def check_cocycle_positions(self) -> None:
        """Audit the lift: ``pos(v) - pos(u) == omega(u->v)`` modulo the torus
        period, on every edge. Raises RuntimeError if it has drifted. The
        position-channel analogue of :meth:`check_cocycle`."""
        _lib.ddg_sampler_cocycle_pos_check(self._handle)

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
