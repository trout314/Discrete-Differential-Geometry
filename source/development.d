/*******************************************************************************
Intrinsic PL geometry via rigid development (notes/CONVENTIONS.md section 6).

D port of python/discrete_differential_geometry/development.py -- the
hot-loop primitive for intrinsic contact geometry (M5 sweep relabeling,
in-loop FP direction data). A simply-connected complex of regular
tetrahedra is unrolled into R^3 by placing each tet congruently across a
shared face; for unit regular tets each step is ONE point reflection.

Everything is exact: we develop in the edge-sqrt(2) rational embedding
(seed tet (0,0,0),(1,1,0),(1,0,1),(0,1,1)), so all coordinates stay in
Q^3 (Rational!BigInt). Exactness gates, enforced with asserts/exceptions:

  - every placed tet has all squared edge lengths == 2 exactly;
  - a BC-chain step map's rotation part has trace == -1/3 exactly
    (twist arccos(-2/3), the ideal Boerdijk-Coxeter tetrahelix);
  - the helix axis (nullspace of R - I; the axis of the cylinder through
    the stack's vertices) is identical at every sampled step.

Relative angles between two axes developed in a COMMON frame are exact
rational cos^2 values -- PL invariants of the local complex + connecting
path (path-dependence is holonomy, i.e. curvature: state the path).
*******************************************************************************/
module development;

import std.bigint : BigInt;
import std.conv : to;
import std.exception : enforce;

import rational_number : Rational, rational;

alias Q = Rational!BigInt;
alias QVec = Q[3];

Q qint(long n)
{
    return rational(BigInt(n));
}

QVec qvec(long a, long b, long c)
{
    return [qint(a), qint(b), qint(c)];
}

/// Canonical seed tet of the edge-sqrt(2) rational embedding.
immutable long[3][4] seedCoords = [[0, 0, 0], [1, 1, 0], [1, 0, 1],
                                   [0, 1, 1]];

QVec seedVertex(size_t i)
{
    assert(i < 4);
    return qvec(seedCoords[i][0], seedCoords[i][1], seedCoords[i][2]);
}

QVec sub(const QVec a, const QVec b)
{
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

Q dot(const QVec a, const QVec b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

QVec cross(const QVec a, const QVec b)
{
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]];
}

QVec scale(const Q t, const QVec v)
{
    return [t * v[0], t * v[1], t * v[2]];
}

/// Mirror image of p through the plane of (q1, q2, q3). Exact.
QVec reflectPoint(const QVec p, const QVec q1, const QVec q2, const QVec q3)
{
    immutable QVec n = cross(sub(q2, q1), sub(q3, q1));
    immutable Q t = qint(2) * dot(sub(p, q1), n) / dot(n, n);
    return sub(p, scale(t, n));
}

private void checkTet(const QVec[4] p, string what)
{
    foreach (i; 0 .. 4)
        foreach (j; i + 1 .. 4)
        {
            immutable QVec d = sub(p[i], p[j]);
            enforce(dot(d, d) == qint(2),
                    "development broke congruence: " ~ what);
        }
}

/*******************************************************************************
Develop the tet stack of a BC-chain vertex sequence (window k =
seq[k .. k+4]). Returns one exact position per vertex (a chain stack is
a path of tets). Congruence of every placed tet is verified.
*******************************************************************************/
QVec[int] developChain(const int[] seq)
{
    enforce(seq.length >= 5, "need at least one step");
    QVec[int] pos;
    foreach (i; 0 .. 4)
        pos[seq[i]] = seedVertex(i);
    QVec[4] first;
    foreach (i; 0 .. 4) first[i] = pos[seq[i]];
    checkTet(first, "seed");
    foreach (k; 0 .. seq.length - 4)
    {
        immutable int drop = seq[k];
        immutable int nv = seq[k + 4];
        enforce(nv !in pos, "chain stack self-overlap (wrap?)");
        pos[nv] = reflectPoint(pos[drop], pos[seq[k + 1]],
                               pos[seq[k + 2]], pos[seq[k + 3]]);
        QVec[4] t;
        foreach (i; 0 .. 4) t[i] = pos[seq[k + 1 + i]];
        checkTet(t, "step");
    }
    return pos;
}

/// Solve the 3x3 rational system with columns (a, b, c): x s.t.
/// x0*a + x1*b + x2*c = rhs (Cramer).
private QVec solve3(const QVec a, const QVec b, const QVec c,
                    const QVec rhs)
{
    immutable Q det = dot(a, cross(b, c));
    enforce(det != qint(0), "singular step matrix");
    return [dot(rhs, cross(b, c)) / det,
            dot(a, cross(rhs, c)) / det,
            dot(a, cross(b, rhs)) / det];
}

/*******************************************************************************
Exact helix axis of a developed chain stack: the fixed direction of the
step screw map seq[k+i] -> seq[k+1+i]. Verifies trace(R) == -1/3 at each
sampled step and that all sampled steps yield the same axis. Returns a
rational direction (not normalized; sign-canonicalized).
*******************************************************************************/
QVec chainAxis(const int[] seq, const QVec[int] pos)
{
    QVec axisAt(size_t k)
    {
        QVec[4] p, q;
        foreach (i; 0 .. 4)
        {
            p[i] = pos[seq[k + i]];
            q[i] = pos[seq[k + 1 + i]];
        }
        QVec[3] M, N;
        foreach (i; 0 .. 3)
        {
            M[i] = sub(p[i + 1], p[0]);
            N[i] = sub(q[i + 1], q[0]);
        }
        // rows of R: solve M^T r_j = (N^T)_j
        QVec[3] MT, R;
        foreach (j; 0 .. 3)
            MT[j] = [M[0][j], M[1][j], M[2][j]];
        foreach (j; 0 .. 3)
        {
            immutable QVec rhs = [N[0][j], N[1][j], N[2][j]];
            R[j] = solve3(MT[0], MT[1], MT[2], rhs);
        }
        enforce(R[0][0] + R[1][1] + R[2][2] == rational(BigInt(-1), BigInt(3)),
                "step trace != -1/3 (not a BC step)");
        QVec[3] A;
        foreach (j; 0 .. 3)
        {
            A[j] = R[j].dup[0 .. 3];
            A[j][j] = A[j][j] - qint(1);
        }
        foreach (i; 0 .. 3)
            foreach (j; i + 1 .. 3)
            {
                QVec u = cross(A[i], A[j]);
                if (u[0] != qint(0) || u[1] != qint(0) || u[2] != qint(0))
                {
                    // canonicalize: divide by first nonzero component
                    foreach (x; u)
                        if (x != qint(0))
                            return [u[0] / x, u[1] / x, u[2] / x];
                }
            }
        throw new Exception("R - I has no rank-2 row pair (not a rotation)");
    }

    immutable size_t last = seq.length - 5;
    immutable QVec a0 = axisAt(0);
    immutable QVec a1 = axisAt(last / 2);
    immutable QVec a2 = axisAt(last);
    enforce(a0 == a1 && a1 == a2, "axis varies along the stack");
    return a0.dup[0 .. 3];
}

/// Exact cos^2 of the undirected angle between two rational directions
/// in a COMMON developed frame.
Q cos2Between(const QVec u, const QVec v)
{
    immutable Q d = dot(u, v);
    return d * d / (dot(u, u) * dot(v, v));
}

///
@system unittest
{
    // reflection is an involution
    immutable QVec p = qvec(3, -2, 5);
    immutable QVec q1 = qvec(0, 0, 0), q2 = qvec(1, 1, 0),
        q3 = qvec(1, 0, 1);
    assert(reflectPoint(reflectPoint(p, q1, q2, q3), q1, q2, q3) == p);
}

@system unittest
{
    import std.range : iota, array;

    // an abstract 20-vertex chain develops to the ideal tetrahelix:
    // congruence, trace -1/3 and axis constancy are all enforced inside
    int[] seq = iota(20).array.to!(int[]);
    auto pos = developChain(seq);
    auto ax = chainAxis(seq, pos);
    assert(cos2Between(ax, ax) == qint(1));

    // two disjoint sub-windows of the same stack: identical axis
    auto axA = chainAxis(seq[0 .. 9], pos);
    auto axB = chainAxis(seq[11 .. 20], pos);
    assert(cos2Between(axA, axB) == qint(1));

    // a genuinely different direction gives cos^2 < 1
    assert(cos2Between(ax, qvec(1, 0, 0)) != qint(1)
           || cos2Between(ax, qvec(0, 1, 0)) != qint(1));
}
