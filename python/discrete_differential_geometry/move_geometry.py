"""Role tables, derived fields, and event-log replay for the geometry ledger.

Python mirror of the role taxonomy in ``source/sampler.d`` (GeometryLedger).
Every accepted move partitions its support into orbits ("roles") whose degree
changes are fixed integers, so the role-resolved ledgers are a lossless
generating set for linear geometry-change observables; the named derived
fields below are the standard combinations. The event log records one
fixed-size record per accepted move — the full 4D cobordism, one 4-simplex
(or, for a 4-4 move, one 2-pentachoron stack) per record — and everything
here is reconstructible from it offline via :func:`replay_role_counts`.

Move type codes: 0: 1→4, 1: 2→3, 2: 3→2, 3: 4→1, 4: 4-4 hinge.
Event label layout: bistellar = center then coCenter (sizes implied by type,
unused slots −1); 4-4 = removedEdge (2) then linkCycle (4) rotated so the
added diagonal is (labels[2], labels[4]).
"""
from __future__ import annotations

import numpy as np

THETA_TET = float(np.arccos(1.0 / 3.0))   # regular-tet dihedral angle

MOVE_NAMES = ["1->4", "2->3", "3->2", "4->1", "4-4"]

# Vertex roles (column order of vertex_role_counts) and their degree changes.
VROLE_NAMES = [
    "v23_pole", "v23_equator",
    "v32_pole", "v32_equator",
    "v14_base", "v14_created",
    "v41_base", "v41_destroyed",
    "v44_pole", "v44_diag", "v44_passive",
]
V_DELTA = np.array([+2, 0, -2, 0, +2, 0, -2, 0, -2, +2, 0], dtype=float)
V_BIRTH = np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=float)
V_DEATH = np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype=float)

# Edge roles (column order of edge_role_counts) and their degree changes.
EROLE_NAMES = [
    "e23_equator", "e23_spoke", "e23_created",
    "e32_triangle", "e32_spoke", "e32_destroyed",
    "e14_base", "e14_created",
    "e41_base", "e41_destroyed",
    "e44_destroyed", "e44_created", "e44_polediag", "e44_polepassive",
    "e44_equator",
]
E_DELTA = np.array([-1, +1, 0, +1, -1, 0, +1, 0, -1, 0, 0, 0, 0, -1, +1],
                   dtype=float)
E_BIRTH = np.array([0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0], dtype=float)
E_DEATH = np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0], dtype=float)

# Event record: u64 clock + u32 type + 6 x i32 labels = 36 bytes, packed.
EVENT_DTYPE = np.dtype([("clock", "<u8"), ("type", "<u4"),
                        ("labels", "<i4", (6,))], align=False)
assert EVENT_DTYPE.itemsize == 36

# labels split (nA, nB) per move type: A/B = center/coCenter for bistellar,
# removedEdge/linkCycle for the 4-4 move.
SUPPORT_SPLIT = {0: (4, 1), 1: (3, 2), 2: (2, 3), 3: (1, 4), 4: (2, 4)}


# ---------------------------------------------------------------------------
# Derived fields (standard linear combinations of the ledgers)
# ---------------------------------------------------------------------------

def degree_velocity(counts: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """Net degree change per simplex: counts[n, nroles] @ delta."""
    return counts @ delta


def volume_flux(vcounts: np.ndarray) -> np.ndarray:
    """Signed local volume deposition per vertex (trace-K analog).

    Net participation in volume-changing channels: (1→4) − (4→1). Under the
    volume pin this sums to ≈ 0 globally — the discrete maximal-slicing check.
    """
    i = {n: k for k, n in enumerate(VROLE_NAMES)}
    return (vcounts[:, i["v14_base"]] + vcounts[:, i["v14_created"]]
            - vcounts[:, i["v41_base"]] - vcounts[:, i["v41_destroyed"]])


def deficit_flux(ecounts: np.ndarray) -> np.ndarray:
    """Regge deficit-angle change rate per edge: −θ_tet · Δdeg(e).

    Positive = curvature increasing at that edge. The per-edge curvature
    transport field.
    """
    return -THETA_TET * (ecounts @ E_DELTA)


def lapse(vcounts: np.ndarray) -> np.ndarray:
    """Total accepted participation per vertex (the unnormalized lapse)."""
    return vcounts.sum(axis=1)


# ---------------------------------------------------------------------------
# Event-log replay
# ---------------------------------------------------------------------------

def _support(ev):
    nA, nB = SUPPORT_SPLIT[int(ev["type"])]
    lab = ev["labels"]
    return lab[:nA], lab[nA:nA + nB]


def replay_role_counts(events: np.ndarray):
    """Reconstruct the vertex/edge role ledgers from an event array.

    Returns (vroles, eroles): dicts mapping vertex label / sorted edge tuple
    to length-11 / length-15 count arrays. Exact mirror of the D-side ledger
    (used as a cross-validation of both the log and the ledger).
    """
    vi = {n: k for k, n in enumerate(VROLE_NAMES)}
    ei = {n: k for k, n in enumerate(EROLE_NAMES)}
    vroles: dict[int, np.ndarray] = {}
    eroles: dict[tuple, np.ndarray] = {}

    def bv(v, role):
        a = vroles.get(v)
        if a is None:
            a = vroles[v] = np.zeros(len(VROLE_NAMES))
        a[vi[role]] += 1

    def be(a_, b_, role):
        k = (a_, b_) if a_ < b_ else (b_, a_)
        a = eroles.get(k)
        if a is None:
            a = eroles[k] = np.zeros(len(EROLE_NAMES))
        a[ei[role]] += 1

    for ev in events:
        t = int(ev["type"])
        A, B = _support(ev)
        A = [int(x) for x in A]
        B = [int(x) for x in B]
        if t == 0:                                    # 1->4
            for v in A: bv(v, "v14_base")
            bv(B[0], "v14_created")
            for i in range(4):
                for j in range(i + 1, 4): be(A[i], A[j], "e14_base")
            for v in A: be(v, B[0], "e14_created")
        elif t == 1:                                  # 2->3
            for v in A: bv(v, "v23_equator")
            for v in B: bv(v, "v23_pole")
            for i in range(3):
                for j in range(i + 1, 3): be(A[i], A[j], "e23_equator")
            for c in A:
                for p in B: be(c, p, "e23_spoke")
            be(B[0], B[1], "e23_created")
        elif t == 2:                                  # 3->2
            for v in A: bv(v, "v32_pole")
            for v in B: bv(v, "v32_equator")
            be(A[0], A[1], "e32_destroyed")
            for c in A:
                for q in B: be(c, q, "e32_spoke")
            for i in range(3):
                for j in range(i + 1, 3): be(B[i], B[j], "e32_triangle")
        elif t == 3:                                  # 4->1
            bv(A[0], "v41_destroyed")
            for v in B: bv(v, "v41_base")
            for v in B: be(A[0], v, "e41_destroyed")
            for i in range(4):
                for j in range(i + 1, 4): be(B[i], B[j], "e41_base")
        else:                                         # 4-4; diag = B[0], B[2]
            for v in A: bv(v, "v44_pole")
            bv(B[0], "v44_diag"); bv(B[2], "v44_diag")
            bv(B[1], "v44_passive"); bv(B[3], "v44_passive")
            be(A[0], A[1], "e44_destroyed")
            be(B[0], B[2], "e44_created")
            for p in A:
                be(p, B[0], "e44_polediag"); be(p, B[2], "e44_polediag")
                be(p, B[1], "e44_polepassive"); be(p, B[3], "e44_polepassive")
            for i in range(4):
                be(B[i], B[(i + 1) % 4], "e44_equator")
    return vroles, eroles


def edge_incarnations(events: np.ndarray):
    """Per-incarnation edge records reconstructed from the event log.

    Returns (closed, open_) where closed is a structured array with fields
    (a, b, birth_clock, birth_type, death_clock, death_type) and open_ maps
    still-alive (a, b) -> (birth_clock, birth_type). Edges alive before the
    log started appear only at death and are returned in the third element
    (censored list of (a, b, death_clock, death_type)).
    """
    open_: dict[tuple, tuple] = {}
    closed, censored = [], []

    def born(a, b, clk, t):
        k = (a, b) if a < b else (b, a)
        open_[k] = (clk, t)

    def died(a, b, clk, t):
        k = (a, b) if a < b else (b, a)
        if k in open_:
            bc, bt = open_.pop(k)
            closed.append((k[0], k[1], bc, bt, clk, t))
        else:
            censored.append((k[0], k[1], clk, t))

    for ev in events:
        t = int(ev["type"]); clk = int(ev["clock"])
        A, B = _support(ev)
        A = [int(x) for x in A]; B = [int(x) for x in B]
        if t == 0:
            for v in A: born(v, B[0], clk, t)
        elif t == 1:
            born(B[0], B[1], clk, t)
        elif t == 2:
            died(A[0], A[1], clk, t)
        elif t == 3:
            for v in B: died(A[0], v, clk, t)
        else:
            died(A[0], A[1], clk, t)
            born(B[0], B[2], clk, t)
    rec = np.array(closed, dtype=[("a", "<i4"), ("b", "<i4"),
                                  ("birth_clock", "<u8"), ("birth_type", "<u1"),
                                  ("death_clock", "<u8"), ("death_type", "<u1")])
    return rec, open_, censored
