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
