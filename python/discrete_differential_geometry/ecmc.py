"""Event-chain (lifted) transport of a degree-3 chord along a BC chain.

The reversible non-local slide picks a fresh slot every proposal, so the chord
random-walks and its displacement is diffusive. LIFTING holds the slot fixed
between events: the chord runs in one direction until it would have been
rejected, and only then turns around. That converts diffusive O(L^2) exploration
into ballistic O(L) -- the best a lift can do (Chen-Lovasz-Pak bound the gain at
quadratic), and the scaling already measured for a single defect
(tau_int 158.4 -> 2.4, growing with L).

Two primitives live here.

RUNG SEQUENCE.  Every triangle has an exact integer 2->3 creation cost

    Q(t) = 2 * (Scross - Sface) + 18

with Sface the three face-edge degrees and Scross the six apex-shell degrees.
A slide preserves N3, hence the edge count and the degree sum, so

    dS_slide(t1 -> t2) = c * (Q(t2) - Q(t1))     EXACT (9e-16 measured)

On R the ladder has only four rungs, Q in {46, 48, 50, 52}, and dQ = 0 slides
are exactly free. Q is a function of the FACE alone, so it is constant on Aut
face orbits -- which is both a correctness check and the reason a chain's rung
sequence can be cached per chain CLASS (14 on R) instead of per chain (756).

FLIGHT SAMPLER.  Under per-step factorized Metropolis the probability of
surviving j steps along the ray is

    P(j) = prod_{i<j} min(1, e^{-beta dS_i}) = exp(-beta * Lambda_j),
    Lambda_j = sum_{i<j} [dS_i]^+

so drawing E ~ Exp(1) and taking the largest j with beta*Lambda_j <= E samples
the stopping point exactly. Lambda is an explicit staircase over a cached
integer sequence, so the whole flight costs one lookup rather than j
accept/reject tests -- and the intermediate positions are never visited, which
is the point of event-driven MC, not an approximation to it.

WHAT THIS MODULE DOES NOT DO.  It samples the flight against the STATIC
washboard only. A real channel must also stop at contacts with other chords and
at windows where the recreating 2->3 is invalid or changes species; omit any of
those and the flight tunnels through an obstacle. See `blockers=` -- the caller
supplies them, and the sampler stops at whichever event comes first.
"""
from __future__ import annotations

import numpy as np

__all__ = ["face_rung", "chain_rungs", "uphill_staircase", "sample_flight",
           "FlightResult", "resolve_event", "EVENT_POLICIES"]


def _d(edeg, a, b):
    return edeg[(a, b)] if a < b else edeg[(b, a)]


def face_rung(face, apexes, edeg):
    """Exact integer 2->3 creation cost Q of a triangle.

    `face` is the three vertices, `apexes` the two apexes of its two tets,
    `edeg` a dict of sorted vertex pairs -> edge degree.
    """
    a, b, c = face
    sface = _d(edeg, a, b) + _d(edeg, b, c) + _d(edeg, a, c)
    scross = sum(_d(edeg, x, y) for x in apexes for y in (a, b, c))
    return 2 * (scross - sface) + 18


def chain_rungs(verts, edeg):
    """Rung sequence Q[k] along a closed BC chain.

    The 2->3 site at chain position k is the face between chain tets k and k+1:
    face (v[k+1], v[k+2], v[k+3]) with apexes (v[k], v[k+4]) -- the
    chain-internal geometry that makes chain-aligned worms possible at all.
    Indices are cyclic, so this is defined at every position of a closed chain.

    Returns an int array of length len(verts).
    """
    v = [int(x) for x in verts]
    L = len(v)
    out = np.empty(L, np.int64)
    for k in range(L):
        f = (v[(k + 1) % L], v[(k + 2) % L], v[(k + 3) % L])
        ap = (v[k], v[(k + 4) % L])
        out[k] = face_rung(f, ap, edeg)
    return out


def uphill_staircase(rungs, start, sigma, c=0.34034, max_steps=None):
    """Cumulative uphill action Lambda_j along the ray, j = 0, 1, 2, ...

    Lambda_0 = 0 and Lambda_{j+1} = Lambda_j + [c*(Q_{j+1} - Q_j)]^+, so
    Lambda is flat exactly where the ray runs along one rung -- the free
    Q46/Q50 web is crossed at zero accumulated cost.

    Returned array has length max_steps+1 (default: once around the chain).
    """
    L = len(rungs)
    n = L if max_steps is None else int(max_steps)
    idx = (int(start) + int(sigma) * np.arange(n + 1)) % L
    q = rungs[idx].astype(np.float64)
    step = c * np.diff(q)
    return np.concatenate([[0.0], np.cumsum(np.maximum(step, 0.0))])


class FlightResult:
    """Outcome of one ballistic flight."""

    __slots__ = ("steps", "reason", "stop_index", "uphill", "blocker")

    def __init__(self, steps, reason, stop_index, uphill, blocker=None):
        self.steps = int(steps)          # successful steps taken
        self.reason = reason             # 'washboard' | 'blocker' | 'horizon'
        self.stop_index = int(stop_index)  # chain position after the flight
        self.uphill = float(uphill)      # action climbed during the flight
        self.blocker = blocker           # whatever the caller passed, if any

    def __repr__(self):
        return (f"FlightResult(steps={self.steps}, reason={self.reason!r}, "
                f"stop_index={self.stop_index}, uphill={self.uphill:.4f})")


def sample_flight(rungs, start, sigma, rng, beta=1.0, c=0.34034,
                  blockers=None, horizon=None):
    """Sample where a ballistic flight stops. Exact, not an approximation.

    Draws E ~ Exp(1) and returns the largest j with beta*Lambda_j <= E, i.e.
    the number of steps that survive per-step factorized Metropolis. The
    equivalence to stepping is exact (see the module docstring), and
    `tests/test_ecmc.py` checks it against a direct step-by-step simulation.

    `blockers` maps chain position -> opaque object (another chord, an invalid
    target, a species change). The flight stops at the FIRST blocker on the ray
    even if the washboard would have let it run further -- that is what makes a
    flight equivalent to a walk rather than a tunnel.

    `horizon` caps the flight length (one lap by default).
    """
    L = len(rungs)
    n = L if horizon is None else int(horizon)
    lam = uphill_staircase(rungs, start, sigma, c=c, max_steps=n)
    E = rng.exponential()

    # first blocker on the ray, if any (j = 1.. so the chord never blocks itself)
    j_block, blk = None, None
    if blockers:
        for j in range(1, n + 1):
            pos = (int(start) + int(sigma) * j) % L
            if pos in blockers:
                j_block, blk = j, blockers[pos]
                break

    # largest j with beta*Lambda_j <= E  (Lambda_0 = 0, so j >= 0 always)
    j_wash = int(np.searchsorted(beta * lam, E, side="right") - 1)
    j_wash = min(j_wash, n)

    if j_block is not None and j_block <= j_wash:
        # the blocker is reached before the washboard would have stopped us;
        # we stop just short of it
        j = j_block - 1
        reason = "blocker"
    elif j_wash >= n:
        j, reason, blk = n, "horizon", None
        return FlightResult(j, reason, (int(start) + int(sigma) * j) % L,
                            float(lam[j]), None)
    else:
        j, reason, blk = j_wash, "washboard", None

    return FlightResult(j, reason, (int(start) + int(sigma) * j) % L,
                        float(lam[j]), blk)


# ---------------------------------------------------------------------------
# event resolution -- deliberately pluggable
# ---------------------------------------------------------------------------

#: What may happen when a flight stops. The choice is a physics question, not
#: a detail, so it is a parameter rather than a hard-coded rule.
#:
#:   reflect        sigma -> -sigma, same chord stays active. Correct and
#:                  simple, but it destroys the transport chain and with it the
#:                  only momentum-like object in the sampler.
#:   handoff_sigma  the blocker becomes active, keeping the direction --
#:                  hard-sphere ECMC, a Newton's-cradle momentum transfer.
#:   handoff_chain  the blocker becomes active on the SAME chain. Differs from
#:                  handoff_sigma exactly when the two chords ride different
#:                  chains through one contact, which is the generic case.
#:   rotate_slot    the same chord stays active but takes a different slot.
#:                  The velocity space here is 12 SLOTS, not +-1, so this is a
#:                  lifting move that reflect cannot express: it changes chain
#:                  without reversing, and it is how a chord reaches the
#:                  crystal-spanning free web (which lives ACROSS chains).
#:   refresh        redraw the lifting variables entirely. Always balance-safe
#:                  (it is just velocity resampling) and what a horizon stop
#:                  should do -- reflecting at a non-event would pair a forward
#:                  flight with an equal backward one and trap a 2-cycle.
EVENT_POLICIES = ("reflect", "handoff_sigma", "handoff_chain", "rotate_slot",
                  "refresh")


def resolve_event(result, active, slot, rng, policy="handoff_sigma",
                  reverse_slot=None, slot_choices=None, refresh=None):
    """Decide the next (active object, slot) after a flight stops.

    Pure bookkeeping over the lifting variables -- the caller supplies the
    problem-specific maps, so this stays testable without a sampler:

      reverse_slot(slot)  -> the slot walking the same chain backwards
      slot_choices(active, stop_index) -> slots available at the stop site
      refresh(rng)        -> a fresh (active, slot)

    A horizon stop always refreshes regardless of `policy`: it is not an
    event, so there is nothing to reflect off or hand to.

    BALANCE WARNING. Whatever policy you pick must pair with the involution
    sigma -> -sigma, or skew detailed balance fails and pi moves. In
    particular a rule that PREFERS cheap continuations (say, "pick a dS-free
    slot if one exists") is not automatically reversible, because freeness of
    the step out of a site is not symmetric under negating the direction. That
    is the same trap as putting selectivity in a closure criterion. Verify any
    new policy with a two-sided test before trusting it.
    """
    if policy not in EVENT_POLICIES:
        raise ValueError(f"unknown policy {policy!r}; choose from "
                         f"{EVENT_POLICIES}")
    if result.reason == "horizon":
        if refresh is None:
            raise ValueError("a horizon stop must refresh; supply refresh=")
        return (*refresh(rng), "refresh")

    if policy == "refresh":
        if refresh is None:
            raise ValueError("policy 'refresh' needs refresh=")
        return (*refresh(rng), "refresh")

    if policy == "reflect" or result.reason == "washboard":
        # a washboard stop has no blocker to hand anything to
        if reverse_slot is None:
            raise ValueError("reflection needs reverse_slot=")
        return active, reverse_slot(slot), "reflect"

    if policy == "rotate_slot":
        if slot_choices is None:
            raise ValueError("policy 'rotate_slot' needs slot_choices=")
        opts = [s for s in slot_choices(active, result.stop_index) if s != slot]
        if not opts:
            return active, reverse_slot(slot), "reflect"
        return active, opts[rng.integers(len(opts))], "rotate_slot"

    # handoff_sigma / handoff_chain: the blocker takes over
    if result.blocker is None:
        raise ValueError("handoff policy but the flight reported no blocker")
    if policy == "handoff_sigma":
        return result.blocker, slot, "handoff_sigma"
    return result.blocker, slot, "handoff_chain"
