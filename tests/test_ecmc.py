"""Tests for lifted / event-chain transport primitives (ecmc.py).

The load-bearing claim is that a FLIGHT is exactly equivalent to stepping: the
sampled stopping point has the same distribution as running per-step factorized
Metropolis and recording where it first fails. If that holds, skipping the
intermediate positions is event-driven MC rather than an approximation to it.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
for _p in ("scripts", "tools"):
    sys.path.insert(0, os.path.join(_ROOT, _p))

from discrete_differential_geometry.ecmc import (
    face_rung, chain_rungs, uphill_staircase, sample_flight)

C = 0.34034                      # lam * e* / 6, the rung coefficient


def _step_by_step(rungs, start, sigma, rng, beta, c, horizon):
    """Reference implementation: walk one step at a time, accepting each with
    min(1, e^{-beta dS}), and report how many steps survived."""
    L = len(rungs)
    j = 0
    while j < horizon:
        q0 = rungs[(start + sigma * j) % L]
        q1 = rungs[(start + sigma * (j + 1)) % L]
        dS = c * (int(q1) - int(q0))
        if dS > 0 and rng.random() > np.exp(-beta * dS):
            return j
        j += 1
    return horizon


@pytest.mark.parametrize("beta", [0.5, 1.0, 3.0])
def test_flight_matches_step_by_step(beta):
    """Same stopping distribution as stepping -- the whole justification for
    jumping straight to the event."""
    rng = np.random.default_rng(20260802)
    rungs = rng.choice([46, 48, 50, 52], size=60)
    horizon = 40
    for start in (0, 7, 23):
        for sigma in (+1, -1):
            a = np.array([sample_flight(rungs, start, sigma,
                                        np.random.default_rng(s), beta=beta,
                                        c=C, horizon=horizon).steps
                          for s in range(4000)])
            b = np.array([_step_by_step(rungs, start, sigma,
                                        np.random.default_rng(s), beta, C,
                                        horizon)
                          for s in range(4000)])
            # compare means and the full survival curve
            assert abs(a.mean() - b.mean()) < 0.15 * (1 + b.mean())
            for j in range(0, min(horizon, 12)):
                pa, pb = (a >= j).mean(), (b >= j).mean()
                assert abs(pa - pb) < 0.04, (beta, start, sigma, j, pa, pb)


def test_survival_is_exactly_exp_minus_beta_lambda():
    """P(survive j steps) = exp(-beta * Lambda_j) by construction."""
    rng = np.random.default_rng(7)
    rungs = rng.choice([46, 48, 50, 52], size=50)
    beta, start, sigma, horizon = 1.0, 3, +1, 30
    lam = uphill_staircase(rungs, start, sigma, c=C, max_steps=horizon)
    steps = np.array([sample_flight(rungs, start, sigma,
                                    np.random.default_rng(s), beta=beta, c=C,
                                    horizon=horizon).steps
                      for s in range(20000)])
    for j in range(1, 12):
        assert abs((steps >= j).mean() - np.exp(-beta * lam[j])) < 0.02


def test_same_rung_stretches_are_free():
    """A ray that never changes rung accumulates zero uphill, so the flight
    runs to the horizon every time -- the free Q46/Q50 web is ballistic."""
    rungs = np.full(40, 46)
    lam = uphill_staircase(rungs, 0, +1, c=C)
    assert np.all(lam == 0.0)
    for s in range(200):
        r = sample_flight(rungs, 0, +1, np.random.default_rng(s), horizon=25)
        assert r.steps == 25 and r.reason == "horizon" and r.uphill == 0.0


def test_downhill_is_free_too():
    """Only UPHILL steps accumulate; a descending ladder never stops the
    flight (factorized Metropolis accepts every downhill step)."""
    rungs = np.array([52, 50, 48, 46] * 10)
    r = sample_flight(rungs, 0, +1, np.random.default_rng(0), horizon=30)
    assert r.uphill >= 0.0
    lam = uphill_staircase(rungs, 0, +1, c=C, max_steps=30)
    assert lam[-1] > 0          # the ladder wraps, so some uphill exists
    flat = uphill_staircase(np.array([52, 50, 48, 46]), 0, +1, c=C, max_steps=3)
    assert np.all(flat == 0.0)  # strictly descending: no uphill at all


def test_blocker_stops_the_flight_short():
    """The flight must stop BEFORE a blocker even when the washboard would
    have let it run on -- otherwise it tunnels through an obstacle."""
    rungs = np.full(40, 46)                       # free ray: no washboard stop
    r = sample_flight(rungs, 0, +1, np.random.default_rng(0),
                      blockers={5: "other-chord"}, horizon=25)
    assert r.reason == "blocker" and r.steps == 4
    assert r.stop_index == 4 and r.blocker == "other-chord"


def test_blocker_beyond_the_washboard_stop_is_irrelevant():
    rungs = np.array([46, 52] * 20)               # steep: stops almost at once
    hits = [sample_flight(rungs, 0, +1, np.random.default_rng(s),
                          blockers={18: "far"}, horizon=25).reason
            for s in range(200)]
    assert "blocker" not in hits


# ---------------------------------------------------------------------------
# against the real crystal
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def r_chain():
    """A real R chain, built in-memory. The chain comes from CrystalSymmetry's
    own enumeration rather than worm_helix.bc_orbit -- same sequence (asserted
    in test_chain_select.py), no dependency on scripts/."""
    import tcp_reference as tr
    from discrete_differential_geometry import Manifold, CrystalSymmetry
    facets = np.asarray(tr.build_t3_triangulation("r", 2)[0])
    m = Manifold(3, facets.tolist())
    edeg = {tuple(sorted(map(int, e))): m.degree(e) for e in m.simplices(1)}
    sym = CrystalSymmetry.compute(facets)
    return sym.chain_vertices(0), edeg, sym


def test_rung_ladder_on_R_has_exactly_four_values(r_chain):
    """Q in {46,48,50,52} on R -- the measured 4-rung ladder."""
    verts, edeg, _ = r_chain
    q = chain_rungs(verts, edeg)
    assert set(q.tolist()) <= {46, 48, 50, 52}
    assert len(set(q.tolist())) == 4


def test_rung_is_constant_on_aut_face_orbits(r_chain):
    """Q is a function of the FACE, so it must be constant on Aut orbits --
    which is what lets the rung sequence be cached per chain CLASS."""
    import tcp_reference as tr
    facets = np.asarray(tr.build_t3_triangulation("r", 2)[0])
    _, edeg, sym = r_chain
    f2a = {}
    for t in facets:
        t = tuple(sorted(int(x) for x in t))
        for i in range(4):
            f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
    fo = sym.orbit_id_map("face")
    by_orbit = {}
    for f, ap in f2a.items():
        if len(ap) != 2:
            continue
        q = face_rung(f, tuple(ap), edeg)
        by_orbit.setdefault(fo[f], set()).add(q)
    assert all(len(v) == 1 for v in by_orbit.values()), \
        "Q differs within an Aut face orbit -- impossible if Q is a face function"
    assert set().union(*by_orbit.values()) <= {46, 48, 50, 52}


# ---------------------------------------------------------------------------
# pluggable event resolution
# ---------------------------------------------------------------------------


def _flight(reason, blocker=None, stop=3):
    from discrete_differential_geometry.ecmc import FlightResult
    return FlightResult(3, reason, stop, 0.0, blocker)


def test_horizon_always_refreshes_whatever_the_policy():
    """A horizon stop is not an event: nothing to reflect off or hand to, and
    reflecting there would trap a 2-cycle."""
    from discrete_differential_geometry.ecmc import resolve_event, EVENT_POLICIES
    rng = np.random.default_rng(0)
    for pol in EVENT_POLICIES:
        a, s, kind = resolve_event(_flight("horizon"), "chordA", 3, rng,
                                   policy=pol, reverse_slot=lambda s: -s,
                                   refresh=lambda r: ("fresh", 7))
        assert (a, s, kind) == ("fresh", 7, "refresh"), pol


def test_policies_dispatch_distinctly():
    from discrete_differential_geometry.ecmc import resolve_event
    rng = np.random.default_rng(0)
    rev = lambda s: -s
    fl = _flight("blocker", blocker="chordB")
    assert resolve_event(fl, "chordA", 3, rng, "reflect", reverse_slot=rev) \
        == ("chordA", -3, "reflect")
    assert resolve_event(fl, "chordA", 3, rng, "handoff_sigma") \
        == ("chordB", 3, "handoff_sigma")
    assert resolve_event(fl, "chordA", 3, rng, "handoff_chain") \
        == ("chordB", 3, "handoff_chain")
    a, s, k = resolve_event(fl, "chordA", 3, rng, "rotate_slot",
                            reverse_slot=rev,
                            slot_choices=lambda a, i: [1, 2, 3, 4])
    assert a == "chordA" and s in (1, 2, 4) and k == "rotate_slot"


def test_washboard_stop_cannot_hand_off():
    """No blocker exists at a washboard stop, so even a handoff policy must
    reflect rather than invent one."""
    from discrete_differential_geometry.ecmc import resolve_event
    rng = np.random.default_rng(0)
    a, s, kind = resolve_event(_flight("washboard"), "chordA", 3, rng,
                               policy="handoff_sigma", reverse_slot=lambda s: -s)
    assert (a, s, kind) == ("chordA", -3, "reflect")


def test_unknown_policy_is_refused():
    from discrete_differential_geometry.ecmc import resolve_event
    with pytest.raises(ValueError, match="unknown policy"):
        resolve_event(_flight("blocker", "b"), "a", 1,
                      np.random.default_rng(0), policy="steer_to_free")


# ---------------------------------------------------------------------------
# same-rung hops
# ---------------------------------------------------------------------------


def test_hop_is_a_bijection_with_a_deterministic_inverse():
    """T_sigma and T_-sigma are inverse, which is what makes the proposal
    symmetric -- and with dS = 0 that gives acceptance exactly 1."""
    from discrete_differential_geometry.ecmc import hop_target, hop_index
    rng = np.random.default_rng(3)
    rungs = rng.choice([46, 48, 50, 52], size=80)
    for start in range(0, 80, 7):
        for sigma in (+1, -1):
            for k in (1, 2, 3, 5):
                t = hop_target(rungs, start, sigma, k)
                if t is None:
                    continue
                assert rungs[t] == rungs[start]          # same rung
                assert hop_target(rungs, t, -sigma, k) == start % len(rungs)
                assert hop_index(rungs, start, t, sigma) == k


def test_hop_skips_intervening_off_rung_sites():
    """The hop passes OVER barriers rather than through them -- it never
    visits the intermediate sites, so they cost nothing."""
    from discrete_differential_geometry.ecmc import hop_target
    rungs = np.array([46, 52, 52, 52, 52, 46, 50, 50, 46])
    assert hop_target(rungs, 0, +1, 1) == 5      # jumps the four 52s
    assert hop_target(rungs, 0, +1, 2) == 8
    assert hop_target(rungs, 5, -1, 1) == 0


def test_hop_returns_none_when_the_rung_is_unique():
    from discrete_differential_geometry.ecmc import hop_target
    rungs = np.array([46, 50, 50, 50])
    assert hop_target(rungs, 0, +1, 1) is None   # the only 46 on the chain


def test_same_rung_positions_on_a_real_chain(r_chain):
    """On R every rung is well populated, so a hop always has targets."""
    from discrete_differential_geometry.ecmc import (
        same_rung_positions, hop_target)
    verts, edeg, _ = r_chain
    q = chain_rungs(verts, edeg)
    for idx in (0, 17, 101):
        pos = same_rung_positions(q, idx)
        assert len(pos) > 10 and idx % len(q) in set(pos.tolist())
        for sigma in (+1, -1):
            t = hop_target(q, idx, sigma, 1)
            assert t is not None and q[t] == q[idx % len(q)]


def test_hop_never_returns_the_source():
    """The chain is cyclic, so a full lap lands back on the start -- a hop to
    the source is not a move and must never be proposed as one."""
    from discrete_differential_geometry.ecmc import hop_target
    rungs = np.array([46, 50, 50, 50])
    for sigma in (+1, -1):
        for k in range(1, 6):
            assert hop_target(rungs, 0, sigma, k) is None
    rungs2 = np.array([46, 50, 46, 50])           # exactly two 46s
    assert hop_target(rungs2, 0, +1, 1) == 2
    assert hop_target(rungs2, 0, +1, 2) is None   # would be the source again
