"""Intrinsic PL geometry via rigid development (CONVENTIONS.md section 6).

The unit-edge triangulation's local geometry is computed by DEVELOPMENT:
unrolling a simply-connected complex of regular tetrahedra into R^3,
placing each tet congruently across a shared face. For regular tets the
step is a single point reflection: the new tet's apex is the mirror
image of the previous tet's dropped vertex through the shared-face
plane.

Everything here is EXACT. We develop in the edge-sqrt(2) rational
embedding -- seed tet (0,0,0),(1,1,0),(1,0,1),(0,1,1) -- where face
planes have rational normals, reflections are rational maps, and all
developed coordinates stay in Q^3. Consequences:

  - every developed tet has all squared edge lengths == 2, checked
    exactly (any deviation is a bug, not roundoff);
  - a BC chain's step map (window k -> k+1) is a rational screw motion
    whose rotation part R satisfies trace(R) == -1/3 exactly
    (twist arccos(-2/3) ~ 131.81 deg -- the ideal Boerdijk-Coxeter
    tetrahelix; the registry's 135 deg is embedding strain);
  - the helix axis (the user-facing definition: axis of the cylinder
    through the stack's vertices) is the rational nullspace of R - I;
  - relative angles between two developed axes are exact rational
    cos^2(theta) values -- intrinsic PL invariants of the local complex
    given the connecting tet path.

Angle comparisons are LOCAL: two axes are comparable only after
developing both stacks in a common frame from a connecting dual-graph
tet path; different paths differ by holonomy (that path-dependence is
curvature, not error -- state the path).

Python + Fraction is deliberate for now: this is offline exact
analysis. A D port is queued for when development enters hot collision
loops (sweep relabeling), per the D-first performance policy.
"""
from __future__ import annotations

from fractions import Fraction
from collections import deque

Vec = tuple  # (Fraction, Fraction, Fraction)

_SEED = ((Fraction(0), Fraction(0), Fraction(0)),
         (Fraction(1), Fraction(1), Fraction(0)),
         (Fraction(1), Fraction(0), Fraction(1)),
         (Fraction(0), Fraction(1), Fraction(1)))
EDGE_SQ = Fraction(2)          # squared edge length of the embedding


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def reflect(p: Vec, q1: Vec, q2: Vec, q3: Vec) -> Vec:
    """Mirror image of p through the plane of (q1, q2, q3). Exact."""
    n = _cross(_sub(q2, q1), _sub(q3, q1))
    t = Fraction(2) * _dot(_sub(p, q1), n) / _dot(n, n)
    return (p[0] - t * n[0], p[1] - t * n[1], p[2] - t * n[2])


def _check_tet(pos, tet):
    for i in range(4):
        for j in range(i + 1, 4):
            d = _sub(pos[tet[i]], pos[tet[j]])
            if _dot(d, d) != EDGE_SQ:
                raise AssertionError(
                    f"development broke congruence on tet {tet}")


def develop_chain(seq):
    """Develop the tet stack of a BC chain vertex sequence (window k =
    seq[k:k+4]). Returns {vertex: position} -- a chain stack is a path
    of tets, so each vertex gets one position. Exact; congruence of
    every placed tet is checked."""
    pos = {seq[i]: _SEED[i] for i in range(4)}
    _check_tet(pos, seq[0:4])
    for k in range(len(seq) - 4):
        drop, face, new = seq[k], seq[k + 1:k + 4], seq[k + 4]
        if new in pos:
            raise AssertionError("chain stack self-overlap (wrap?)")
        pos[new] = reflect(pos[drop], *(pos[v] for v in face))
        _check_tet(pos, seq[k + 1:k + 5])
    return pos


def develop_path(start_tet, start_pos, path):
    """Develop along a dual-graph tet path. start_tet: 4 vertices with
    positions in start_pos ({v: Vec}); path: list of tets (4-tuples),
    consecutive tets sharing a face, path[0] adjacent to (or equal to)
    start_tet. Returns per-tet position dicts [{v: Vec}, ...] -- kept
    per-tet, NOT merged, so revisited vertex labels (holonomy) cannot
    silently conflict."""
    placements = []
    prev_tet, prev_pos = tuple(start_tet), dict(start_pos)
    for tet in path:
        tet = tuple(tet)
        shared = [v for v in tet if v in prev_tet]
        if len(shared) == 4:
            placements.append({v: prev_pos[v] for v in tet})
            prev_tet, prev_pos = tet, placements[-1]
            continue
        if len(shared) != 3:
            raise ValueError(f"path step shares {len(shared)} != 3 verts")
        drop = next(v for v in prev_tet if v not in shared)
        new = next(v for v in tet if v not in shared)
        p = {v: prev_pos[v] for v in shared}
        p[new] = reflect(prev_pos[drop], *(p[v] for v in shared))
        _check_tet(p, tet)
        placements.append(p)
        prev_tet, prev_pos = tet, p
    return placements


def _solve3(M, rhs):
    """Solve 3x3 rational system M x = rhs (columns of M are vectors)."""
    a, b, c = M
    det = _dot(a, _cross(b, c))
    if det == 0:
        raise ValueError("singular step matrix")
    x0 = _dot(rhs, _cross(b, c)) / det
    x1 = _dot(a, _cross(rhs, c)) / det
    x2 = _dot(a, _cross(b, rhs)) / det
    return (x0, x1, x2)


def chain_axis(seq, pos):
    """Exact helix axis of a developed chain stack: the fixed direction
    of the step screw map seq[k+i] -> seq[k+1+i]. Verifies that the
    rotation part R has trace exactly -1/3 (ideal BC twist) and that
    every step gives the same axis. Returns a rational direction vector
    (not normalized)."""
    def step_R(k):
        p = [pos[seq[k + i]] for i in range(4)]
        q = [pos[seq[k + 1 + i]] for i in range(4)]
        M = [_sub(p[i + 1], p[0]) for i in range(3)]
        N = [_sub(q[i + 1], q[0]) for i in range(3)]
        # R M = N  (columns); solve row-wise: R^T columns from M^T x = N_row
        # Build R by solving M^T r_j = N^T_j for each coordinate row j.
        MT = [tuple(M[i][j] for i in range(3)) for j in range(3)]
        R = []
        for j in range(3):
            rhs = tuple(N[i][j] for i in range(3))
            R.append(_solve3(MT, rhs))
        return R  # rows

    def axis_of(R):
        tr = R[0][0] + R[1][1] + R[2][2]
        if tr != Fraction(-1, 3):
            raise AssertionError(f"step trace {tr} != -1/3 (not a BC step)")
        # nullspace of R - I via cross products of two rows
        A = [(R[0][0] - 1, R[0][1], R[0][2]),
             (R[1][0], R[1][1] - 1, R[1][2]),
             (R[2][0], R[2][1], R[2][2] - 1)]
        for i in range(3):
            for j in range(i + 1, 3):
                u = _cross(A[i], A[j])
                if any(x != 0 for x in u):
                    return u
        raise AssertionError("R - I has no rank-2 pair (not a rotation)")

    ks = [0, (len(seq) - 5) // 2, len(seq) - 5]
    axes = []
    for k in ks:
        u = axis_of(step_R(k))
        # canonicalize sign/scale for comparison
        g = None
        for x in u:
            if x != 0:
                g = x
                break
        axes.append(tuple(x / g for x in u))
    if len(set(axes)) != 1:
        raise AssertionError("axis varies along the stack (not rigid?)")
    return axes[0]


def cos2_between(u, v):
    """Exact cos^2 of the (undirected) angle between two rational
    directions in a COMMON developed frame."""
    return Fraction(_dot(u, v)) ** 2 / (_dot(u, u) * _dot(v, v))


def _rotation_between(p, q):
    """Exact rotation R (rows of Fractions) with R (p_i - p_0) =
    (q_i - q_0), for two congruent 4-point placements."""
    M = [_sub(p[i + 1], p[0]) for i in range(3)]
    N = [_sub(q[i + 1], q[0]) for i in range(3)]
    MT = [tuple(M[i][j] for i in range(3)) for j in range(3)]
    R = []
    for j in range(3):
        rhs = tuple(N[i][j] for i in range(3))
        R.append(_solve3(MT, rhs))
    return R


# ---------------------------------------------------------------------------
# Exact rotations in the spinor representation: integer quaternions
# ---------------------------------------------------------------------------

from math import gcd as _gcd


class Quat:
    """Primitive integer quaternion (a, b, c, d) -- the exact spinor
    avatar of a rational rotation (Euler-Rodrigues: R entries are
    quadratic forms over the norm N = a^2+b^2+c^2+d^2). Composition is
    the integer Hamilton product; R(q1 * q2) = R(q1) R(q2).

    The SIGN is the SU(2) double-cover datum: +q and -q give the same
    rotation but opposite spin lifts. Composition preserves signs, so a
    loop composed from per-step lifts carries a well-defined Z2 spin
    holonomy (per-step sign conventions are a Z2 gauge choice; use
    canonical() only for gauge-fixed rotation LABELS, never mid-path).
    """
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a, b, c, d):
        g = _gcd(_gcd(abs(a), abs(b)), _gcd(abs(c), abs(d)))
        if g == 0:
            raise ValueError("zero quaternion")
        if g > 1:
            a, b, c, d = a // g, b // g, c // g, d // g
        self.a, self.b, self.c, self.d = a, b, c, d

    def __repr__(self):
        return f"Quat({self.a},{self.b},{self.c},{self.d})"

    def __eq__(self, other):
        return (self.a, self.b, self.c, self.d) == \
            (other.a, other.b, other.c, other.d)

    def __hash__(self):
        return hash((self.a, self.b, self.c, self.d))

    @property
    def norm(self):
        return self.a ** 2 + self.b ** 2 + self.c ** 2 + self.d ** 2

    def canonical(self):
        """Gauge-fixed sign: first nonzero of (a, b, c, d) positive."""
        for x in (self.a, self.b, self.c, self.d):
            if x != 0:
                return self if x > 0 else Quat(-self.a, -self.b,
                                               -self.c, -self.d)
        raise AssertionError

    def same_rotation(self, other):
        return self.canonical() == other.canonical()

    def __mul__(self, o):
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = o.a, o.b, o.c, o.d
        return Quat(a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                    a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                    a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                    a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)

    def conj(self):
        return Quat(self.a, -self.b, -self.c, -self.d)

    def __pow__(self, k):
        if k < 0:
            return self.conj() ** (-k)
        r = Quat(1, 0, 0, 0)
        base = self
        while k:
            if k & 1:
                r = r * base
            base = base * base
            k >>= 1
        return r

    def cos_phi(self):
        """Exact cos of the rotation angle."""
        return Fraction(2 * self.a ** 2, self.norm) - 1

    def axis(self):
        """Integer rotation axis (b, c, d), gcd-reduced; None if
        identity."""
        b, c, d = self.b, self.c, self.d
        g = _gcd(_gcd(abs(b), abs(c)), abs(d))
        if g == 0:
            return None
        return (b // g, c // g, d // g)

    def to_matrix(self):
        """Rows of Fractions; exact, orthogonal, det +1."""
        a, b, c, d = self.a, self.b, self.c, self.d
        N = self.norm
        return [(Fraction(a * a + b * b - c * c - d * d, N),
                 Fraction(2 * (b * c - a * d), N),
                 Fraction(2 * (b * d + a * c), N)),
                (Fraction(2 * (b * c + a * d), N),
                 Fraction(a * a - b * b + c * c - d * d, N),
                 Fraction(2 * (c * d - a * b), N)),
                (Fraction(2 * (b * d - a * c), N),
                 Fraction(2 * (c * d + a * b), N),
                 Fraction(a * a - b * b - c * c + d * d, N))]

    @staticmethod
    def lift(R):
        """Exact lift of a rational rotation matrix (rows of Fractions)
        to a primitive integer quaternion, canonical sign. The
        roundtrip to_matrix(lift(R)) == R is asserted."""
        tr = R[0][0] + R[1][1] + R[2][2]
        cands = [
            (1 + tr, R[2][1] - R[1][2], R[0][2] - R[2][0],
             R[1][0] - R[0][1]),
            (R[2][1] - R[1][2], 1 + 2 * R[0][0] - tr,
             R[0][1] + R[1][0], R[0][2] + R[2][0]),
            (R[0][2] - R[2][0], R[0][1] + R[1][0],
             1 + 2 * R[1][1] - tr, R[1][2] + R[2][1]),
            (R[1][0] - R[0][1], R[0][2] + R[2][0],
             R[1][2] + R[2][1], 1 + 2 * R[2][2] - tr)]
        for k, v in enumerate(cands):
            if v[k] != 0:
                den = 1
                for x in v:
                    den = den * Fraction(x).denominator // _gcd(
                        den, Fraction(x).denominator)
                q = Quat(*(int(x * den) for x in v)).canonical()
                assert q.to_matrix() == [tuple(r) for r in R], \
                    "quaternion lift roundtrip failed (R not in SO(3,Q)?)"
                return q
        raise AssertionError("all lift candidates vanish")


# ---------------------------------------------------------------------------
# Parallel transport: canonical frames + Wilson lines
# ---------------------------------------------------------------------------

class TransportContext:
    """Exact parallel transport over one host triangulation.

    Every tet gets a CANONICAL frame: its vertices in sorted order --
    with the last two swapped where the global orientation demands it
    (orientation parity is propagated once over the dual graph, so
    frames are consistent manifold-wide and Wilson lines from separate
    calls compose). A Wilson line along a face-adjacent tet path is the
    rotation from the start tet's canonical frame to the end tet's,
    computed through the rigid development; returned as a Quat (spinor
    lift with the canonical per-line sign).

    facets: iterable of 4-tuples/lists of vertex labels.
    """

    def __init__(self, facets):
        self.tets = [tuple(sorted(int(v) for v in f)) for f in facets]
        tset = {frozenset(t) for t in self.tets}
        self.dual = {}
        f2t = {}
        for t in tset:
            for v in t:
                f2t.setdefault(t - {v}, []).append(t)
        for fc, ts in f2t.items():
            for t in ts:
                for u in ts:
                    if u != t:
                        self.dual.setdefault(t, []).append(u)
        self.dual = {t: sorted(ns, key=sorted)
                     for t, ns in self.dual.items()}
        self._orient(tset)

    def _orient(self, tset):
        """Propagate a global orientation: parity[tet] = 1 iff the
        canonical order is sorted-with-last-two-swapped. Two adjacent
        tets are consistent iff they induce OPPOSITE orientations on
        the shared face; removing one vertex from a sorted tuple leaves
        a sorted (identity-permutation) face order, so the induced sign
        is just (-1)^(index of the omitted vertex)."""
        from collections import deque
        self.parity = {}
        for start in sorted(tset, key=sorted):
            if start in self.parity:
                continue
            self.parity[start] = 0
            q = deque([start])
            while q:
                t = q.popleft()
                st = sorted(t)
                for u in self.dual.get(t, ()):
                    it = st.index(next(v for v in t if v not in u))
                    su = sorted(u)
                    iu = su.index(next(v for v in u if v not in t))
                    want = self.parity[t] ^ ((it ^ iu ^ 1) & 1)
                    if u in self.parity:
                        if self.parity[u] != want:
                            raise ValueError("non-orientable complex")
                    else:
                        self.parity[u] = want
                        q.append(u)

    def canon_order(self, tet):
        t = tuple(sorted(tet))
        if self.parity[frozenset(t)]:
            return (t[0], t[1], t[3], t[2])
        return t

    def canon_placement(self, tet):
        order = self.canon_order(tet)
        return {v: _SEED[i] for i, v in enumerate(order)}

    def wilson_line(self, path):
        """Exact transport along a face-adjacent tet path, start
        canonical frame -> end canonical frame. Returns a Quat."""
        path = [frozenset(t) for t in path]
        start = tuple(self.canon_order(path[0]))
        placements = develop_path(start, self.canon_placement(path[0]),
                                  [tuple(sorted(t)) for t in path])
        end_pos = placements[-1]
        order = self.canon_order(path[-1])
        placed = [end_pos[v] for v in order]
        canonical = [_SEED[i] for i in range(4)]
        # g_end: canonical -> placed; transport = g_end^{-1} (g_0 = id)
        Rg = _rotation_between(canonical, placed)
        RgT = [tuple(Rg[i][j] for i in range(3)) for j in range(3)]
        return Quat.lift(RgT)

    def holonomy(self, loop):
        """Wilson line of a closed path (first == last tet): the exact
        loop rotation in the base tet's canonical frame. Returns
        (Quat, cos_phi, axis)."""
        if frozenset(loop[0]) != frozenset(loop[-1]):
            raise ValueError("loop must return to its base tet")
        q = self.wilson_line(loop)
        return q, q.cos_phi(), q.axis()


def chain_step_quat(seq, pos):
    """Spinor lift of the BC-chain step screw's rotation part. Verifies
    the norm-6 signature (cos phi = -2/3, axis norm^2 = 5*a^2)."""
    p = [pos[seq[i]] for i in range(4)]
    q4 = [pos[seq[1 + i]] for i in range(4)]
    R = _rotation_between(p, q4)
    RT = [tuple(R[i][j] for i in range(3)) for j in range(3)]
    q = Quat.lift(RT)
    assert q.cos_phi() == Fraction(-2, 3), "not a BC chain step"
    return q


def chain_transport(seq, pos, k):
    """Exact k-step transport along a developed chain: q_step ** k."""
    return chain_step_quat(seq, pos) ** k
