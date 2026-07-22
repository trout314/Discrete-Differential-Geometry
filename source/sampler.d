/// Shared MCMC sampler core: objective function, move selection, and helpers.
/// Used by manifold_sampler.d, bench_sampler.d, and the C API (ddg_capi.d).
module sampler;

import std.algorithm, std.array, std.conv, std.format, std.math, std.range, std.typecons;
import std.random : uniform, uniform01, rndChoice = choice;
import manifold, manifold_moves, simplicial_complex, utility;

alias isIRof = isInputRangeOf;

// ---------------------------------------------------------------------------
// Penalty computation
// ---------------------------------------------------------------------------

struct Penalty
{
    real volumePenalty;
    real globalCurvPenalty;
    real localCurvPenalty;
    real localSolidAngleCurvPenalty;
    real hingeDegTargetPenalty;
    real coDim3DegTargetPenalty;
}

/// Compute penalties from raw values (no manifold needed).
Penalty penaltiesFromValues(int dim_, P)(
    long nFacets, long nHinges, ulong hingeTotSqDeg,
    long nCoDim3, ulong coDim3TotSqDeg, P params)
{
    enum hingesPerFacet = dim_ * (dim_ + 1) / 2;
    enum coDim3PerFacet = binomial(dim_ + 1, dim_ - 2);

    Penalty penalty;
    penalty.volumePenalty = (nFacets - params.numFacetsTarget) ^^ 2;

    immutable nHingesTarget = hingesPerFacet * nFacets / params.hingeDegreeTarget;
    penalty.globalCurvPenalty = (nHinges - nHingesTarget) ^^ 2;

    immutable degTarget = hingesPerFacet * nFacets / cast(real) nHinges;
    real _;
    real x = modf(degTarget, _);
    real minPenalty = (x - x ^^ 2) * nHinges;

    penalty.localCurvPenalty = (
        degTarget ^^ 2 * nHinges - 2 * degTarget * hingesPerFacet * nFacets + hingeTotSqDeg) - minPenalty;
    // Intensive: divide by count to get mean variance per hinge
    if (nHinges > 0) penalty.localCurvPenalty /= nHinges;

    // Fixed-target hinge penalty: sum_e (deg_e - t)^2 with t a CONSTANT target
    // (params.hingeDegreeTarget), minus the integer-lattice floor x(1-x) per
    // hinge. Extensive and linear in (nFacets, nHinges, totSqDeg), so the
    // action term is strictly local: a move's delta depends only on degrees
    // inside its bistellar ball. Identity vs. the variance form:
    //   sum(d - t)^2 = sum(d - dbar)^2 + nHinges*(dbar - t)^2.
    immutable tHinge = cast(real) params.hingeDegreeTarget;
    x = modf(tHinge, _);
    penalty.hingeDegTargetPenalty =
        tHinge ^^ 2 * nHinges - 2 * tHinge * hingesPerFacet * nFacets + hingeTotSqDeg
        - (x - x ^^ 2) * nHinges;

    static if (dim_ > 2)
    {
        immutable coDim3DegTarget = coDim3PerFacet * nFacets / cast(real) nCoDim3;
        // Codim-3 face degrees are always EVEN: a codim-3 face's link is a
        // 2-sphere, and a triangulated 2-sphere always has an even number of
        // triangles — which equals the face's facet-degree. So the minimum
        // achievable variance splits across the nearest *even* integers, giving
        // floor 4 y(1-y) with y = frac(degTarget/2), not the integer-lattice
        // x(1-x). (The hinge term above keeps x(1-x): edge links are cycles of
        // any length, so edge degrees take any integer.)
        x = modf(coDim3DegTarget / 2.0, _);
        minPenalty = 4.0 * (x - x ^^ 2) * nCoDim3;

        penalty.localSolidAngleCurvPenalty = (
            coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + coDim3TotSqDeg) - minPenalty;
        // Intensive: divide by count to get mean variance per codim-3 face
        if (nCoDim3 > 0) penalty.localSolidAngleCurvPenalty /= nCoDim3;

        // Fixed-target codim-3 penalty: sum_v (deg_v - t)^2 with t a CONSTANT
        // target (params.coDim3DegreeTarget), minus the even-lattice floor
        // 4y(1-y) per face (codim-3 degrees are even; see comment above).
        // Extensive and linear in the counters, hence strictly local.
        immutable tCoDim3 = cast(real) params.coDim3DegreeTarget;
        x = modf(tCoDim3 / 2.0, _);
        penalty.coDim3DegTargetPenalty =
            tCoDim3 ^^ 2 * nCoDim3 - 2 * tCoDim3 * coDim3PerFacet * nFacets + coDim3TotSqDeg
            - 4.0 * (x - x ^^ 2) * nCoDim3;
    }
    else
    {
        penalty.localSolidAngleCurvPenalty = 0;
        penalty.coDim3DegTargetPenalty = 0;
    }

    return penalty;
}

Penalty penalties(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    immutable nFacets = mfd.fVector[dim];
    immutable nHinges = mfd.fVector[dim - 2];
    immutable totSqDeg = mfd.totalSquareDegree(dim - 2);

    static if (dim > 2)
    {
        immutable nCoDim3 = mfd.fVector[dim - 3];
        immutable totSAsqDeg = mfd.totalSquareDegree(dim - 3);
    }
    else
    {
        enum nCoDim3 = 0;
        enum totSAsqDeg = 0;
    }

    return penaltiesFromValues!dim(
        nFacets, nHinges, totSqDeg, nCoDim3, totSAsqDeg, params);
}

real objectiveFromPenalty(P)(Penalty pen, P params)
{
    return params.numFacetsCoef * pen.volumePenalty
        + params.numHingesCoef * pen.globalCurvPenalty
        + params.hingeDegreeVarianceCoef * pen.localCurvPenalty
        + params.coDim3DegreeVarianceCoef * pen.localSolidAngleCurvPenalty
        + params.hingeDegreeTargetCoef * pen.hingeDegTargetPenalty
        + params.coDim3DegreeTargetCoef * pen.coDim3DegTargetPenalty;
}

real objective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    auto pen = mfd.penalties(params);
    return objectiveFromPenalty(pen, params);
}

/// Fixed-target vs variance-about-mean identity:
///   sum(d - t)^2 = sum(d - dbar)^2 + n*(dbar - t)^2
/// checked through both floor conventions on a live manifold.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.7;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    alias BM = BistellarMove!3;
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    auto params = TestParams();
    auto pen = mfd.penalties(params);

    immutable nH = cast(real) mfd.fVector[1];
    immutable nV = cast(real) mfd.fVector[0];

    // Reconstruct raw (floor-free) sums of squared deviations from each form
    real _;
    immutable dbarH = 6.0L * mfd.fVector[3] / nH;
    real xm = modf(dbarH, _);
    real xt = modf(cast(real) params.hingeDegreeTarget, _);
    immutable rawVarH = pen.localCurvPenalty * nH + (xm - xm ^^ 2) * nH;
    immutable rawTgtH = pen.hingeDegTargetPenalty + (xt - xt ^^ 2) * nH;
    assert(isClose(rawTgtH,
        rawVarH + nH * (dbarH - params.hingeDegreeTarget) ^^ 2, 1e-9L));

    immutable dbarV = 4.0L * mfd.fVector[3] / nV;
    real ym = modf(dbarV / 2, _);
    real yt = modf(cast(real) params.coDim3DegreeTarget / 2, _);
    immutable rawVarV = pen.localSolidAngleCurvPenalty * nV + 4 * (ym - ym ^^ 2) * nV;
    immutable rawTgtV = pen.coDim3DegTargetPenalty + 4 * (yt - yt ^^ 2) * nV;
    assert(isClose(rawTgtV,
        rawVarV + nV * (dbarV - params.coDim3DegreeTarget) ^^ 2, 1e-9L));
}

// ---------------------------------------------------------------------------
// Speculative delta
// ---------------------------------------------------------------------------

/******************************************************************************
Compute the objective delta for a bistellar move without executing it.
Enumerates affected sub-simplices and looks up their degrees to compute
speculative totSqDeg deltas.
*/
real speculativeBistellarDelta(int dim, Vertex, P)(
    const ref Manifold!(dim, Vertex) mfd,
    const ref BistellarMove!(dim, Vertex) move,
    real currentObjective,
    P params)
{
    auto center = move.center;
    auto coCenter = move.coCenter;
    immutable cenLen = cast(int) center.length;
    immutable coCenLen = cast(int) coCenter.length;

    // Combined vertex set, sorted (needed for degreeMap lookup via subsetsOfSize)
    Vertex[dim + 2] allVertsBuf;
    allVertsBuf[0 .. cenLen] = center[];
    allVertsBuf[cenLen .. cenLen + coCenLen] = coCenter[];
    auto allVerts = allVertsBuf[0 .. cenLen + coCenLen];
    allVerts.sort();

    // Compute speculative f-vector
    uint[dim + 1] newFVector = mfd.fVector[0 .. dim + 1];
    newFVector[].modifyFVector(move);

    // Compute speculative totSqDeg for dimensions 0 through dim-2.
    // For each sub-simplex s of dimension d in the combined vertex set:
    //   delta(s) = (|C| - |s∩C|) - (|CC| - |s∩CC|)
    //   ΔtotSqDeg[d] += 2*deg(s)*delta + delta²
    long[dim - 1] newTotSqDeg;
    foreach (d; 0 .. dim - 1)
        newTotSqDeg[d] = cast(long) mfd.totalSquareDegree(d);

    // Enumerate subsets of each relevant dimension
    static foreach (d; 0 .. dim - 1)
    {{
        foreach (subset; allVerts[].subsetsOfSize(d + 1))
        {
            // Count how many vertices in this subset are from the center
            int s_C = 0;
            foreach (v; subset)
            {
                if (center.canFind(v)) s_C++;
            }
            int s_CC = d + 1 - s_C;
            int delta = (cenLen - s_C) - (coCenLen - s_CC);

            if (delta == 0) continue;

            long deg = cast(long) mfd.degreeOrZero!d(subset);
            newTotSqDeg[d] += 2 * deg * delta + cast(long) delta * delta;
        }
    }}

    // Compute new objective from speculative values
    static if (dim > 2)
    {
        auto newPen = penaltiesFromValues!dim(
            cast(long) newFVector[dim], cast(long) newFVector[dim - 2],
            cast(ulong) newTotSqDeg[dim - 2],
            cast(long) newFVector[dim - 3],
            cast(ulong) newTotSqDeg[dim - 3],
            params);
    }
    else
    {
        auto newPen = penaltiesFromValues!dim(
            cast(long) newFVector[dim], cast(long) newFVector[dim - 2],
            cast(ulong) newTotSqDeg[dim - 2],
            0, 0, params);
    }

    return objectiveFromPenalty(newPen, params) - currentObjective;
}

///
unittest
{
    import std.random : Mt19937;

    // Test that speculative delta matches actual delta for all move types
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    import manifold_examples : standardSphere;

    // Start from a sphere and do some 1→4 moves to get a nontrivial triangulation
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    // Grow the manifold a bit to create diverse degree distributions
    alias BM = BistellarMove!3;
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Now test speculative delta on all available bistellar moves
    foreach (move; mfd.allBistellarMoves)
    {
        auto currentObj = mfd.objective(params);
        auto specDelta = mfd.speculativeBistellarDelta(move, currentObj, params);

        // Actually execute and compute
        auto mfdCopy = mfd;
        mfdCopy.doMove(move);
        auto actualNewObj = mfdCopy.objective(params);
        auto actualDelta = actualNewObj - currentObj;

        assert(isClose(specDelta, actualDelta, 1e-6),
            "speculative delta mismatch: spec=%s actual=%s for move %s"
            .format(specDelta, actualDelta, move));
    }
}

// ---------------------------------------------------------------------------
// Speculative delta for hinge moves
// ---------------------------------------------------------------------------

/******************************************************************************
Compute the objective delta for a 4-4 hinge move without executing it.
A hinge move preserves the f-vector; only totSqDeg changes for dims 0..dim-2.
Enumerates all sub-simplices of the 6 involved vertices to compute degree deltas.
*/
real speculativeHingeDelta(Vertex, P)(
    const ref Manifold!(3, Vertex) mfd,
    const ref HingeMove!Vertex move,
    real currentObjective,
    P params)
{
    enum dim = 3;

    // Collect the 6 involved vertices (sorted for subset enumeration)
    Vertex[6] allVertsBuf;
    allVertsBuf[0] = move.removedEdge[0];
    allVertsBuf[1] = move.removedEdge[1];
    allVertsBuf[2 .. 6] = move.linkCycle;
    allVertsBuf[].sort();
    auto allVerts = allVertsBuf[];

    // Get old and new facets
    auto oldFacets = move.oldFacets;
    auto newFacets = move.newFacets;

    // f-vector is unchanged by a 4-4 move
    long[dim - 1] newTotSqDeg;
    foreach (d; 0 .. dim - 1)
        newTotSqDeg[d] = cast(long) mfd.totalSquareDegree(d);

    // For each sub-simplex dimension d (0 and 1 for dim=3),
    // enumerate all subsets of the 6 vertices and compute degree deltas.
    static foreach (d; 0 .. dim - 1)
    {{
        foreach (subset; allVerts[].subsetsOfSize(d + 1))
        {
            // Count how many old facets contain this subset
            int oldCount = 0;
            foreach (ref f; oldFacets)
                if (subset.isSubsetOf(f[]))
                    oldCount++;

            // Count how many new facets contain this subset
            int newCount = 0;
            foreach (ref f; newFacets)
                if (subset.isSubsetOf(f[]))
                    newCount++;

            int delta = newCount - oldCount;
            if (delta == 0) continue;

            long deg = cast(long) mfd.degreeOrZero!d(subset);
            newTotSqDeg[d] += 2 * deg * delta + cast(long) delta * delta;
        }
    }}

    // Compute new objective from unchanged f-vector and updated totSqDeg
    auto newPen = penaltiesFromValues!dim(
        cast(long) mfd.fVector[dim], cast(long) mfd.fVector[dim - 2],
        cast(ulong) newTotSqDeg[dim - 2],
        cast(long) mfd.fVector[dim - 3],
        cast(ulong) newTotSqDeg[dim - 3],
        params);

    return objectiveFromPenalty(newPen, params) - currentObjective;
}

///
unittest
{
    // Test that speculativeHingeDelta matches actual delta
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    alias BM = BistellarMove!3;

    // Build a 3-sphere and do 1-4 moves to create degree-4 edges
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Find all degree-4 edges and test speculative delta on valid hinge moves
    int tested = 0;
    foreach (edge; mfd.simplices(1))
    {
        if (mfd.degree(edge) != 4) continue;

        int[2] edgeArr = [edge[0], edge[1]];
        // Pick a start vertex from the link
        // Use any facet containing this edge to find a link vertex
        foreach (facet; mfd.facets)
        {
            if (!edgeArr[].isSubsetOf(facet[]))
                continue;

            // Find the two vertices not in the edge
            int[2] others;
            int oi = 0;
            foreach (v; facet)
                if (v != edgeArr[0] && v != edgeArr[1])
                    others[oi++] = v;

            foreach (diag; 0 .. 2)
            {
                auto hm = mfd.hingeMove(edgeArr, others[0], diag);
                if (!mfd.hasValidHingeMove(hm)) continue;

                auto currentObj = mfd.objective(params);
                auto specDelta = mfd.speculativeHingeDelta(hm, currentObj, params);

                // Actually execute and compute
                auto mfdCopy = mfd;
                mfdCopy.doHingeMove(hm);
                auto actualNewObj = mfdCopy.objective(params);
                auto actualDelta = actualNewObj - currentObj;

                assert(isClose(specDelta, actualDelta, 1e-6),
                    "speculative hinge delta mismatch: spec=%s actual=%s for move %s"
                    .format(specDelta, actualDelta, hm));
                tested++;
            }
            break; // only need one facet per edge
        }
    }
    assert(tested > 0, "should have tested at least one hinge move");
}

// ---------------------------------------------------------------------------
// Vertex 6-valence potential (Z-legality + chemical tilts + impurity valence)
// ---------------------------------------------------------------------------
//
// A strictly-local vertex potential on two per-vertex counters (dim = 3):
//   n6(v) = # incident edges with degree >= 6
//   m(v)  = # incident edges with degree not in {5, 6} ("impurity valence")
// Energy per vertex:
//   U(n6) = zlegCoef * dist^2(n6, {0,2,3,4}) + tilt[n6] (tilt only for n6 <= 4)
//   V(m)  = impCoef * m^2
// Physics: by the link sum rule (sum over incident edges of (6 - deg) = 12,
// link = S^2), a vertex with m = 0 and n6 in {0,2,3,4} is EXACTLY a
// Frank-Kasper coordination (Z12/Z14/Z15/Z16); n6 = 1 is combinatorially
// impossible at m = 0. So U penalizes illegal 6-valences (with quadratic
// analytic tail — n6 is unbounded, a clamped table would leave hubs with no
// restoring force), the tilts are chemical potentials selecting AMONG the
// legal classes (phase selection: A15-type vs Laves-type stoichiometry), and
// V(m) = m^2 penalizes defect CLUSTERING at a vertex (nonlinear on purpose: a
// linear-in-m term is just an edge-level penalty in disguise since
// sum_v m(v) = 2 * #impure-edges; the quadratic makes an octahedral vertex,
// m = 6, cost 36 rather than 6). Zero coefficients = disabled, zero overhead.

struct VertexPot
{
    real zlegCoef = 0;
    real impCoef = 0;
    real[5] tilt = [0, 0, 0, 0, 0];

    bool enabled() const pure nothrow @nogc @safe
    {
        if (zlegCoef != 0 || impCoef != 0) return true;
        foreach (t; tilt) if (t != 0) return true;
        return false;
    }

    real U(long k) const pure nothrow @nogc @safe
    {
        real d2;
        if (k == 0 || (k >= 2 && k <= 4)) d2 = 0;
        else if (k == 1 || k == 5) d2 = 1;
        else d2 = cast(real)(k - 4) ^^ 2;   // analytic quadratic tail
        immutable t = (k >= 0 && k <= 4) ? tilt[cast(size_t) k] : 0.0L;
        return zlegCoef * d2 + t;
    }

    real V(long m) const pure nothrow @nogc @safe
    {
        return impCoef * cast(real)(m * m);
    }
}

/// Per-vertex counter state for the vertex potential. Only nonzero counters
/// are stored; a vertex absent from both maps contributes U(0) + V(0)
/// (= tilt[0]) so `total` must account for the vertex count (see recompute).
struct VertexPotState(Vertex)
{
    int[Vertex] n6;
    int[Vertex] mImp;
    real total = 0;
}

private bool ind6(long d) pure nothrow @nogc @safe { return d >= 6; }
private bool indImp(long d) pure nothrow @nogc @safe
{
    return d > 0 && (d < 5 || d > 6);
}

/// Rebuild the counter state from scratch (init, coefficient changes, tests).
void recomputeVertexPotState(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot)
{
    st.n6.clear;
    st.mImp.clear;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d))
        {
            st.n6.require(e[0], 0)++;
            st.n6.require(e[1], 0)++;
        }
        if (indImp(d))
        {
            st.mImp.require(e[0], 0)++;
            st.mImp.require(e[1], 0)++;
        }
    }
    // total = f0 * U(0)  +  corrections for vertices with nonzero counters.
    st.total = cast(real) mfd.fVector[0] * pot.U(0);
    foreach (v, k; st.n6) st.total += pot.U(k) - pot.U(0);
    foreach (v, m; st.mImp) st.total += pot.V(m);
}

/// Potential delta for a bistellar move (dim 3). Must be called BEFORE the
/// move is applied (reads pre-move degrees). With commit = true, updates the
/// counter state and running total (still pre-move: all reads are speculative).
real potentialBistellarDelta(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    const ref BistellarMove!(3, Vertex) move,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot,
    bool commit)
{
    enum dim = 3;
    auto center = move.center;
    auto coCenter = move.coCenter;
    immutable cenLen = cast(int) center.length;
    immutable coCenLen = cast(int) coCenter.length;

    Vertex[dim + 2] allVertsBuf;
    allVertsBuf[0 .. cenLen] = center[];
    allVertsBuf[cenLen .. cenLen + coCenLen] = coCenter[];
    auto allVerts = allVertsBuf[0 .. cenLen + coCenLen];
    allVerts.sort();
    immutable nv = cast(int) allVerts.length;

    int[dim + 2] dn6;
    int[dim + 2] dm;
    foreach (subset; allVerts.subsetsOfSize(2))
    {
        int s_C = 0;
        foreach (v; subset)
            if (center.canFind(v)) s_C++;
        immutable s_CC = 2 - s_C;
        immutable delta = (cenLen - s_C) - (coCenLen - s_CC);
        if (delta == 0) continue;

        immutable oldDeg = cast(long) mfd.degreeOrZero!1(subset);
        immutable newDeg = oldDeg + delta;
        immutable d6 = (ind6(newDeg) ? 1 : 0) - (ind6(oldDeg) ? 1 : 0);
        immutable dI = (indImp(newDeg) ? 1 : 0) - (indImp(oldDeg) ? 1 : 0);
        if (d6 == 0 && dI == 0) continue;

        foreach (v; subset)
        {
            foreach (i; 0 .. nv)
            {
                if (allVerts[i] == v)
                {
                    dn6[i] += d6;
                    dm[i] += dI;
                    break;
                }
            }
        }
    }

    // 1->4 creates coCenter[0]; 4->1 destroys center[0].
    immutable hasCreated = (coCenLen == 1);
    immutable hasRemoved = (cenLen == 1);

    real dS = 0;
    foreach (i; 0 .. nv)
    {
        immutable v = allVerts[i];
        immutable isCreated = hasCreated && v == coCenter[0];
        immutable isRemoved = hasRemoved && v == center[0];
        if (dn6[i] == 0 && dm[i] == 0 && !isCreated && !isRemoved) continue;

        immutable long oldN6 = st.n6.get(v, 0);
        immutable long oldM = st.mImp.get(v, 0);
        immutable real oldE = isCreated ? 0 : pot.U(oldN6) + pot.V(oldM);
        immutable long newN6 = oldN6 + dn6[i];
        immutable long newM = oldM + dm[i];
        immutable real newE = isRemoved ? 0 : pot.U(newN6) + pot.V(newM);
        dS += newE - oldE;

        if (commit)
        {
            if (isRemoved)
            {
                st.n6.remove(v);
                st.mImp.remove(v);
            }
            else
            {
                if (newN6 != 0) st.n6[v] = cast(int) newN6;
                else st.n6.remove(v);
                if (newM != 0) st.mImp[v] = cast(int) newM;
                else st.mImp.remove(v);
            }
        }
    }
    if (commit) st.total += dS;
    return dS;
}

/// Potential delta for a 4-4 hinge move (dim 3). Same contract as the
/// bistellar version (pre-move reads; commit updates state).
real potentialHingeDelta(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    const ref HingeMove!Vertex move,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot,
    bool commit)
{
    Vertex[6] allVertsBuf;
    allVertsBuf[0] = move.removedEdge[0];
    allVertsBuf[1] = move.removedEdge[1];
    allVertsBuf[2 .. 6] = move.linkCycle;
    allVertsBuf[].sort();
    auto allVerts = allVertsBuf[];

    auto oldFacets = move.oldFacets;
    auto newFacets = move.newFacets;

    int[6] dn6;
    int[6] dm;
    foreach (subset; allVerts.subsetsOfSize(2))
    {
        int oldCount = 0;
        foreach (ref f; oldFacets)
            if (subset.isSubsetOf(f[])) oldCount++;
        int newCount = 0;
        foreach (ref f; newFacets)
            if (subset.isSubsetOf(f[])) newCount++;
        immutable delta = newCount - oldCount;
        if (delta == 0) continue;

        immutable oldDeg = cast(long) mfd.degreeOrZero!1(subset);
        immutable newDeg = oldDeg + delta;
        immutable d6 = (ind6(newDeg) ? 1 : 0) - (ind6(oldDeg) ? 1 : 0);
        immutable dI = (indImp(newDeg) ? 1 : 0) - (indImp(oldDeg) ? 1 : 0);
        if (d6 == 0 && dI == 0) continue;

        foreach (v; subset)
        {
            foreach (i; 0 .. 6)
            {
                if (allVerts[i] == v)
                {
                    dn6[i] += d6;
                    dm[i] += dI;
                    break;
                }
            }
        }
    }

    real dS = 0;
    foreach (i; 0 .. 6)
    {
        if (dn6[i] == 0 && dm[i] == 0) continue;
        immutable v = allVerts[i];
        immutable long oldN6 = st.n6.get(v, 0);
        immutable long oldM = st.mImp.get(v, 0);
        immutable long newN6 = oldN6 + dn6[i];
        immutable long newM = oldM + dm[i];
        dS += (pot.U(newN6) + pot.V(newM)) - (pot.U(oldN6) + pot.V(oldM));

        if (commit)
        {
            if (newN6 != 0) st.n6[v] = cast(int) newN6;
            else st.n6.remove(v);
            if (newM != 0) st.mImp[v] = cast(int) newM;
            else st.mImp.remove(v);
        }
    }
    if (commit) st.total += dS;
    return dS;
}

// ---------------------------------------------------------------------------
// Disclination-network observables (dim = 3)
// ---------------------------------------------------------------------------
//
// Near the flat mean edge degree the curvature-carrying structure is the
// DISCLINATION NETWORK: the graph of degree>=6 edges ("six-edges"). By the
// link sum rule a legal (m = 0) Z14 vertex carries exactly 2 six-edge ends
// (a disclination line passing through), Z15 exactly 3 (a branch node), Z16
// exactly 4 — so the six-edge graph IS the Frank-Kasper skeleton, and its
// components / chains / loops / grafts are the physical line objects. This
// block provides snapshot-cadence censuses of that graph plus the joint
// (n6, m) vertex census. Everything is O(#edges) over a const manifold; no
// sampler state is required, so censuses work on loaded .mfd files too.

struct DisclinationCensus
{
    long nNetVerts;      // vertices with >= 1 incident six-edge
    long nSixEdges;      // edges of degree >= 6
    long nComponents;    // connected components of the six-edge graph
    long giantSize;      // largest component (in vertices)
    long secondSize;     // second-largest component
    long giantDiameter;  // double-BFS pseudo-diameter of the giant component
    long cycleRank;      // first Betti number E - V + C of the network
    long nSegments;      // maximal chains between non-degree-2 network vertices
    long sumSegLen;      // six-edges lying in such chains
    long nPureLoops;     // components that are cycles of degree-2 vertices
    long sumLoopLen;     // six-edges lying in pure loops
    long nEndpoints;     // network vertices with exactly 1 six-edge
    long nFrayVerts;     // network vertices with impurity valence m > 0
    long nImpEndEdges;   // six-edges with >= 1 impure endpoint
    long nZ14, nZ15, nZ16;  // legal (m == 0) network-vertex class populations
    // hostMask-resolved edge categories (bit k of hostMask = "n6-class k is a
    // native host class"). An "FK" endpoint is legal with n6 in {2, 3, 4}.
    long eDopDop;        // both endpoints dopant-FK (legal, class not in mask)
    long eDopHost;       // dopant-FK -- host-FK (the graft edges)
    long eHostHost;      // both endpoints host-FK (the host's own skeleton)
    long eImpAny;        // >= 1 endpoint impure or non-FK (n6 >= 5)
    long[8] netDegCensus;   // network degree 1..7; index 7 clamps >= 8
    long[64] segLenHist;    // chain length; index clamps at 63
    long[32] compSizeHist;  // log2-binned component size (index = bsr(size))
}

/// Slot count of the flattened C-API census layout (see flattenCensus).
enum disclinationCensusSlots = 24 + 8 + 64 + 32;

/// Flatten to the C-API layout: slots 0..20 = the 21 scalar fields in
/// declaration order, 21..23 reserved (zero), 24..31 netDegCensus,
/// 32..95 segLenHist, 96..127 compSizeHist.
void flattenCensus(const ref DisclinationCensus c, long[] outBuf)
in (outBuf.length >= disclinationCensusSlots)
{
    import std.traits : Unqual;
    outBuf[0 .. disclinationCensusSlots] = 0;
    size_t i = 0;
    foreach (field; c.tupleof)
        static if (is(Unqual!(typeof(field)) == long))
            outBuf[i++] = field;
    static assert(DisclinationCensus.tupleof.length == 24);
    assert(i == 21);
    outBuf[24 .. 32] = c.netDegCensus[];
    outBuf[32 .. 96] = c.segLenHist[];
    outBuf[96 .. 128] = c.compSizeHist[];
}

/// Joint (n6, m) vertex census over ALL vertices: outBuf is row-major
/// [min(n6, n6Cap)][min(m, mCap)] with (n6Cap+1)*(mCap+1) slots. Bin (0, 0)
/// is the FK-Z12 bulk; rows n6 >= 1 sum to the disclination-network vertex
/// count. Uncapped, sum_k k * row_k = 2 * #six-edges.
void valenceCensus(Vertex)(const ref Manifold!(3, Vertex) mfd,
    long[] outBuf, int n6Cap, int mCap)
in (outBuf.length >= (n6Cap + 1) * (mCap + 1))
{
    outBuf[0 .. (n6Cap + 1) * (mCap + 1)] = 0;
    int[Vertex] n6, mImp;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d)) { n6.require(e[0], 0)++; n6.require(e[1], 0)++; }
        if (indImp(d)) { mImp.require(e[0], 0)++; mImp.require(e[1], 0)++; }
    }
    long touched = 0;
    foreach (v, k; n6)
    {
        outBuf[min(k, n6Cap) * (mCap + 1) + min(mImp.get(v, 0), mCap)]++;
        ++touched;
    }
    foreach (v, m; mImp)
    {
        if (v in n6) continue;
        outBuf[min(m, mCap)]++;   // n6 == 0 row
        ++touched;
    }
    outBuf[0] += cast(long) mfd.fVector[0] - touched;
}

/// Full disclination-network census (see DisclinationCensus). hostMask marks
/// the host's native n6 classes (e.g. C15: bit 0 | bit 4); pass 0 for no
/// host/dopant split (all legal-legal edges then count as eDopDop).
DisclinationCensus disclinationCensus(Vertex)(
    const ref Manifold!(3, Vertex) mfd, int hostMask = 0)
{
    import core.bitop : bsr;

    DisclinationCensus c;

    // Pass 1: six-edges and per-vertex (n6, m) counters.
    Vertex[2][] six;
    int[Vertex] n6, mImp;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d))
        {
            Vertex[2] pr = [e[0], e[1]];
            six ~= pr;
            n6.require(e[0], 0)++;
            n6.require(e[1], 0)++;
        }
        if (indImp(d))
        {
            mImp.require(e[0], 0)++;
            mImp.require(e[1], 0)++;
        }
    }
    immutable nv = cast(int) n6.length;
    c.nNetVerts = nv;
    c.nSixEdges = cast(long) six.length;
    if (nv == 0) return c;

    // Dense reindex + CSR adjacency.
    int[Vertex] idx;
    auto vlab = new Vertex[nv];
    {
        int k = 0;
        foreach (v, _; n6) { idx[v] = k; vlab[k] = v; ++k; }
    }
    auto ndeg = new int[nv];
    foreach (ref e; six) { ndeg[idx[e[0]]]++; ndeg[idx[e[1]]]++; }
    auto off = new int[nv + 1];
    foreach (i; 0 .. nv) off[i + 1] = off[i] + ndeg[i];
    auto adj = new int[2 * six.length];
    {
        auto fill = off[0 .. nv].dup;
        foreach (ref e; six)
        {
            immutable a = idx[e[0]], b = idx[e[1]];
            adj[fill[a]++] = b;
            adj[fill[b]++] = a;
        }
    }

    // Vertex classification. Note ndeg[i] == n6(vlab[i]) by construction.
    auto imp = new bool[nv];
    foreach (i; 0 .. nv) imp[i] = (vlab[i] in mImp) !is null;
    bool isFK(int i) { return !imp[i] && ndeg[i] >= 2 && ndeg[i] <= 4; }
    bool isHost(int i) { return isFK(i) && ((hostMask >> ndeg[i]) & 1); }

    foreach (i; 0 .. nv)
    {
        c.netDegCensus[min(ndeg[i], 8) - 1]++;
        if (ndeg[i] == 1) c.nEndpoints++;
        if (imp[i]) c.nFrayVerts++;
        else if (ndeg[i] == 2) c.nZ14++;
        else if (ndeg[i] == 3) c.nZ15++;
        else if (ndeg[i] == 4) c.nZ16++;
    }
    foreach (ref e; six)
    {
        immutable a = idx[e[0]], b = idx[e[1]];
        if (imp[a] || imp[b]) c.nImpEndEdges++;
        if (!isFK(a) || !isFK(b)) c.eImpAny++;
        else
        {
            immutable nHost = (isHost(a) ? 1 : 0) + (isHost(b) ? 1 : 0);
            if (nHost == 2) c.eHostHost++;
            else if (nHost == 1) c.eDopHost++;
            else c.eDopDop++;
        }
    }

    // Connected components (iterative DFS with a preallocated stack).
    auto comp = new int[nv];
    comp[] = -1;
    auto work = new int[nv];
    long[] compSizes;
    int giantId = -1;
    foreach (s; 0 .. nv)
    {
        if (comp[s] >= 0) continue;
        immutable cid = cast(int) compSizes.length;
        size_t sp = 0;
        work[sp++] = s;
        comp[s] = cid;
        long size = 0;
        while (sp)
        {
            immutable u = work[--sp];
            ++size;
            foreach (w; adj[off[u] .. off[u + 1]])
                if (comp[w] < 0) { comp[w] = cid; work[sp++] = w; }
        }
        compSizes ~= size;
        if (size > c.giantSize)
        {
            c.secondSize = c.giantSize;
            c.giantSize = size;
            giantId = cid;
        }
        else if (size > c.secondSize)
            c.secondSize = size;
    }
    c.nComponents = cast(long) compSizes.length;
    c.cycleRank = c.nSixEdges - c.nNetVerts + c.nComponents;
    foreach (sz; compSizes)
        c.compSizeHist[min(bsr(cast(ulong) sz), 31)]++;

    // Pseudo-diameter of the giant component: double BFS.
    auto dist = new int[nv];
    int bfsFar(int src)
    {
        dist[] = -1;
        size_t head = 0, tail = 0;
        work[tail++] = src;
        dist[src] = 0;
        int far = src;
        while (head < tail)
        {
            immutable u = work[head++];
            if (dist[u] > dist[far]) far = u;
            foreach (w; adj[off[u] .. off[u + 1]])
                if (dist[w] < 0) { dist[w] = dist[u] + 1; work[tail++] = w; }
        }
        return far;
    }
    foreach (i; 0 .. nv)
        if (comp[i] == giantId)
        {
            immutable a = bfsFar(i);
            immutable b = bfsFar(a);
            c.giantDiameter = dist[b];
            break;
        }

    // Segment decomposition: maximal chains of degree-2 vertices between
    // nodes (network degree != 2), then leftover pure loops (all-degree-2
    // cycles = free disclination loops). The six-graph is simple (edges of a
    // simplicial complex are unique), so chain walking needs no multi-edge
    // guards. Every six-edge lands in exactly one segment or one pure loop.
    bool[ulong] seen;
    ulong ekey(int a, int b)
    {
        return a < b ? (cast(ulong) a << 32) | cast(uint) b
                     : (cast(ulong) b << 32) | cast(uint) a;
    }
    int chainNext(int prev, int cur)
    {
        foreach (w; adj[off[cur] .. off[cur + 1]])
            if (w != prev) return w;
        assert(0, "degree-2 vertex with no forward neighbor");
    }
    foreach (i; 0 .. nv)
    {
        if (ndeg[i] == 2) continue;
        foreach (w0; adj[off[i] .. off[i + 1]])
        {
            if (ekey(i, w0) in seen) continue;
            seen[ekey(i, w0)] = true;
            int prev = i, cur = w0;
            long len = 1;
            while (ndeg[cur] == 2)
            {
                immutable nxt = chainNext(prev, cur);
                seen[ekey(cur, nxt)] = true;
                prev = cur;
                cur = nxt;
                ++len;
            }
            c.nSegments++;
            c.sumSegLen += len;
            c.segLenHist[min(len, 63)]++;
        }
    }
    foreach (i; 0 .. nv)
    {
        if (ndeg[i] != 2) continue;
        foreach (w0; adj[off[i] .. off[i + 1]])
        {
            if (ekey(i, w0) in seen) continue;
            seen[ekey(i, w0)] = true;
            int prev = i, cur = w0;
            long len = 1;
            while (cur != i)
            {
                immutable nxt = chainNext(prev, cur);
                seen[ekey(cur, nxt)] = true;
                prev = cur;
                cur = nxt;
                ++len;
            }
            c.nPureLoops++;
            c.sumLoopLen += len;
        }
    }
    return c;
}

///
unittest
{
    // Boundary of the 4-simplex: every edge has degree 3 — no six-edges, and
    // every vertex has impurity valence 4 (all 4 incident edges are degree 3).
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto c = mfd.disclinationCensus(0);
    c.nNetVerts.shouldEqual(0);
    c.nSixEdges.shouldEqual(0);
    c.nComponents.shouldEqual(0);

    long[5 * 5] vc;
    mfd.valenceCensus(vc[], 4, 4);
    vc[0 * 5 + 4].shouldEqual(5);   // all 5 vertices in bin (n6=0, m=4)
    long tot = 0;
    foreach (x; vc) tot += x;
    tot.shouldEqual(5);
}

// ---------------------------------------------------------------------------
// Integer 1-cocycle tracking (T^3 winding forms / emergent coordinates)
// ---------------------------------------------------------------------------
//
// Maintains representatives of a basis of H^1(M; Z^3) as integer 1-cochains
// omega: oriented edges -> Z^3, CLOSED on every triangle (signed sum = 0).
// For M = T^3 these are the three winding forms; the pairing sum_gamma omega
// over a closed edge loop gamma is its winding vector — gauge-invariant, so
// ANY representative serves for topology (winding, spanning, class checks).
// Geometry (edge direction vectors, Fourier space, nematic order) comes from
// the harmonic representative, computed OFFLINE in Python (graph-Laplacian
// solve); the D side only maintains an exact integer representative.
//
// The update under Pachner moves is FORCED by closedness (the move ball is
// simply connected, so path sums within it are path-independent):
//   1->4  new vertex w in tet {a,b,c,d}: gauge choice places w at a
//         (omega(a->w) = 0, omega(w->v) = omega(a->v)) — pure gauge.
//   2->3  new pole-pole edge: omega(p->q) = omega(p->a) + omega(a->q) for
//         any equator vertex a (choices agree by closedness of the old
//         cocycle on the two old tets).
//   3->2, 4->1  deletions only.
//   4-4   new diagonal: omega(x->z) = omega(x->p) + omega(p->z) via a pole p;
//         then the removed edge is deleted.
// Values are Z^3 so the class is preserved EXACTLY forever; no float drift.

struct CocycleState(Vertex)
{
    int[3][Vertex[2]] omega;   // key: sorted edge (u < v); value = omega(u->v)
    bool enabled;
}

/// omega(a->b) with sign handling for the sorted-key convention.
private int[3] cocGet(Vertex)(const ref CocycleState!Vertex st,
    const Vertex a, const Vertex b)
{
    auto p = mkEdge(a, b) in st.omega;
    assert(p !is null, "cocycle missing edge");
    if (a < b) return *p;
    int[3] r = [-(*p)[0], -(*p)[1], -(*p)[2]];
    return r;
}

/// Store omega(a->b) = val under the sorted-key convention.
private void cocSet(Vertex)(ref CocycleState!Vertex st,
    const Vertex a, const Vertex b, const int[3] val)
{
    if (a < b)
        st.omega[mkEdge(a, b)] = val;
    else
    {
        int[3] r = [-val[0], -val[1], -val[2]];
        st.omega[mkEdge(a, b)] = r;
    }
}

/// Apply the forced cocycle update for an ACCEPTED bistellar move. Reads only
/// edges that persist through the move, so it may run before or after doMove.
void cocycleBistellar(Vertex)(ref CocycleState!Vertex st,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    final switch (cast(int) coCenter.length - 1)
    {
    case 0: // 1->4: place the new vertex at center[0] (gauge choice)
        immutable w = coCenter[0];
        immutable a = center[0];
        immutable int[3] zero = [0, 0, 0];
        cocSet(st, a, w, zero);
        foreach (v; center[1 .. $])
            cocSet(st, w, v, cocGet(st, a, v));
        break;
    case 1: // 2->3: new pole-pole edge, value forced via an equator vertex
        immutable p = coCenter[0], q = coCenter[1], a = center[0];
        int[3] val = cocGet(st, p, a);
        val[] += cocGet(st, a, q)[];
        cocSet(st, p, q, val);
        break;
    case 2: // 3->2: the center edge is destroyed
        st.omega.remove(mkEdge(center[0], center[1]));
        break;
    case 3: // 4->1: the 4 spokes at the removed vertex are destroyed
        immutable w = center[0];
        foreach (v; coCenter)
            st.omega.remove(mkEdge(w, v));
        break;
    }
}

/// Hinge (4-4) analog: diagonal forced via a pole, removed edge deleted.
void cocycleHinge(Vertex)(ref CocycleState!Vertex st,
    const ref HingeMove!Vertex hm)
{
    immutable p = hm.removedEdge[0];
    int[3] val = cocGet(st, hm.addedEdge[0], p);
    val[] += cocGet(st, p, hm.addedEdge[1])[];
    cocSet(st, hm.addedEdge[0], hm.addedEdge[1], val);
    st.omega.remove(mkEdge(hm.removedEdge[0], hm.removedEdge[1]));
}

/// Full audit: key set must equal the manifold's edge set and the cochain
/// must be closed on every triangle. Returns null if clean, else a message.
/// This is the drift check — cheap enough for test/production cadence.
string cocycleProblems(Vertex)(const ref Manifold!(3, Vertex) mfd,
    const ref CocycleState!Vertex st)
{
    long ne = 0;
    foreach (e; mfd.simplices(1))
    {
        ++ne;
        if (mkEdge(e[0], e[1]) !in st.omega)
            return format("cocycle missing edge [%s, %s]", e[0], e[1]);
    }
    if (ne != st.omega.length)
        return format("cocycle has %s edges, manifold has %s",
                      st.omega.length, ne);
    foreach (t; mfd.simplices(2))
    {
        // closedness: omega(a->b) + omega(b->c) - omega(a->c) == 0
        immutable ab = cocGet(st, t[0], t[1]);
        immutable bc = cocGet(st, t[1], t[2]);
        immutable ac = cocGet(st, t[0], t[2]);
        foreach (i; 0 .. 3)
            if (ab[i] + bc[i] - ac[i] != 0)
                return format("cocycle not closed on triangle [%s, %s, %s]",
                              t[0], t[1], t[2]);
    }
    return null;
}

///
unittest
{
    // Closedness is preserved under arbitrary churn. Start from an exact
    // cocycle omega = delta(phi) with random integer phi (closed by
    // construction; on S^3 every closed cochain is exact, so this loses no
    // generality), run mixed MCMC with the cocycle attached, and audit.
    struct TestParams
    {
        int numFacetsTarget = 48;
        real hingeDegreeTarget = 5.1;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    int[int] phi;
    foreach (v; 0 .. 5) phi[v] = uniform(-50, 50);

    CocycleState!int coc;
    coc.enabled = true;
    foreach (e; mfd.simplices(1))
    {
        int[3] val = [phi[e[1]] - phi[e[0]], 2 * (phi[e[1]] - phi[e[0]]), 0];
        coc.omega[mkEdge(e[0], e[1])] = val;
    }
    assert(cocycleProblems(mfd, coc) is null);

    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    foreach (step; 0 .. 3000)
        mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
            hingeTries, hingeAccepts, bTries, bAccepts, null, null, null, null,
            &coc);

    auto prob = cocycleProblems(mfd, coc);
    assert(prob is null, "cocycle drift after churn: " ~ prob);
    assert(hingeAccepts + bAccepts[].sum > 100, "churn too weak to test");
}

/******************************************************************************
Try to propose a hinge move on a 3-manifold. Picks a random facet, a random
edge from that facet, checks for degree 4, picks a random diagonal, and
checks validity. Returns null if no valid move found in this single attempt.

Proposal is symmetric: the forward and reverse proposal probabilities are
equal (same number of containing facets, same edge count per facet, same
f-vector), so no Hastings correction is needed beyond the objective delta.
*/
Nullable!(HingeMove!Vertex) tryProposeHingeMove(Vertex)(
    ref Manifold!(3, Vertex) mfd)
{
    // Pick a random facet
    auto facet = mfd.randomFacetOfDim(3);

    // Pick a random edge from the facet (one of C(4,2)=6 edges)
    static immutable int[2][6] edgePairs = [
        [0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
    auto pair = edgePairs[uniform(0, 6)];

    Vertex[2] edge = [facet[pair[0]], facet[pair[1]]];
    edge[].sort();

    // Check degree 4
    if (mfd.degree(edge[]) != 4)
        return typeof(return).init;

    // Find a start vertex (any facet vertex not in the edge)
    Vertex startVertex = void;
    foreach (v; facet)
        if (v != edge[0] && v != edge[1]) { startVertex = v; break; }

    // Pick a random diagonal and construct the move
    auto hm = mfd.hingeMove(edge, startVertex, uniform(0, 2));

    if (!mfd.hasValidHingeMove(hm))
        return typeof(return).init;

    return nullable(hm);
}

///
unittest
{
    alias BM = BistellarMove!3;

    // Build a manifold with some degree-4 edges
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Try many proposals; at least some should succeed
    int successes = 0;
    foreach (_; 0 .. 200)
    {
        auto result = mfd.tryProposeHingeMove();
        if (!result.isNull) successes++;
    }
    assert(successes > 0, "should have found at least one valid hinge move");
}

/******************************************************************************
Per-vertex move-attribution counters (the "measured combinatorial lapse").

Attribution rule ("Option A"): every event distributes TOTAL WEIGHT 1 uniformly
over the move's support vertices. For a bistellar move the support is the 5
vertices of its bistellar ball (center ∪ coCenter) — the vertex set of the one
4-simplex the move glues onto the triangulation, hence 1/5 per vertex. For a
4-4 hinge move the support is its 6 vertices (removedEdge ∪ linkCycle); its
pentachoron-stack 4-volume is 2, applied in ANALYSIS — the bistellar and hinge
ledgers are kept separate so any weighting convention is a linear combination
downstream. For a 1-4 move the coCenter label is the vertex the move CREATES
(it is attributed like the others; intersect with surviving vertices in
analysis).

Ladder: proposed (concrete move formed, post proposal-thinning, pre validity)
→ valid (passed hasValid*, i.e. counted as a "try") → accepted. Per vertex,
valid/proposed is the kinematic-availability field and accepted/valid the
energetic (Metropolis) field.
*/
struct MoveCounters(Vertex)
{
    double[Vertex] proposed;
    double[Vertex] valid;
    double[Vertex] acceptedBistellar;
    double[Vertex] acceptedHinge;

    void clear()()
    {
        proposed = null;
        valid = null;
        acceptedBistellar = null;
        acceptedHinge = null;
    }
}

/// Distribute total weight 1 uniformly over `support` in `ledger`.
private void addSupport(Vertex)(ref double[Vertex] ledger,
    scope const(Vertex)[] support)
{
    immutable w = 1.0 / support.length;
    foreach (v; support)
        ledger[v] = ledger.get(v, 0.0) + w;
}

///
unittest
{
    double[int] ledger;
    ledger.addSupport([0, 1, 2, 3, 4]);       // one bistellar-ball event
    ledger.addSupport([0, 1, 2, 3, 4]);
    ledger.addSupport([3, 4, 5, 6, 7, 8]);    // one hinge-support event
    ledger[0].shouldEqual(0.4);
    ledger[3].isClose(0.4 + 1.0 / 6).shouldEqual(true);
    // each event contributes total weight 1
    double total = 0;
    foreach (v; ledger.byValue) total += v;
    total.isClose(3.0).shouldEqual(true);
}

/******************************************************************************
Role-resolved geometry ledger (the "maximalist" move-participation record).

Every accepted move partitions its support simplices into orbits ("roles")
under the move's symmetry, and the degree change of a simplex is a FIXED
integer determined by (move type, role) — so the role-resolved ledger is a
lossless generating set for every linear geometry-change observable (volume
flux / trace-K, degree velocity, per-edge deficit flux, the lapse, channel
decompositions). Move type codes follow coCenter.length-1 for bistellar moves,
with 4 for the 4-4 hinge move: 0:1→4, 1:2→3, 2:3→2, 3:4→1, 4:4-4.

Vertex roles and their degree changes (degree = # incident facets):
  2→3: pole(+2)×2 = coCenter, equator(0)×3 = center
  3→2: pole(−2)×2 = center,   equator(0)×3 = coCenter
  1→4: base(+2)×4 = center,   created ×1 = coCenter (born at degree 4)
  4→1: base(−2)×4 = coCenter, destroyed ×1 = center (dies at degree 4)
  4-4: pole(−2)×2 = removedEdge, diag(+2)×2 = addedEdge, passive(0)×2

Edge roles (degree = # incident facets = deficit-angle carrier):
  2→3: equator(−1)×3, spoke(+1)×6, created (pole–pole, born at 3)
  3→2: triangle(+1)×3, spoke(−1)×6, destroyed (center edge, dies at 3)
  1→4: base(+1)×6, created spokes ×4 (born at 3)
  4→1: base(−1)×6, destroyed spokes ×4 (die at 3)
  4-4: destroyed (removedEdge, dies at 4), created (diagonal, born at 4),
       pole–diag(0)×4, pole–passive(−1)×4, equator(+1)×4

Tets have no partial roles (wholesale birth/death only) and are tracked as
AGGREGATES: created/destroyed counts by move type + a log2 lifetime histogram
(age in attempted moves; tets alive when tracking started count as censored).

The optional EVENT LOG appends one fixed-size record per accepted move
(clock: u64, type: u32, labels: 6×i32 = 36 bytes packed): bistellar labels are
center-then-coCenter (support size implied by type, unused slots = -1); 4-4
labels are removedEdge then linkCycle ROTATED so the added diagonal is
(labels[2], labels[4]). The log is the full 4D cobordism, one 4-simplex per
record; roles/ledgers/incarnations are reconstructible offline by replay.
*/
enum VRole
{
    v23Pole, v23Equator,
    v32Pole, v32Equator,
    v14Base, v14Created,
    v41Base, v41Destroyed,
    v44Pole, v44Diag, v44Passive
}

enum ERole
{
    e23Equator, e23Spoke, e23Created,
    e32Triangle, e32Spoke, e32Destroyed,
    e14Base, e14Created,
    e41Base, e41Destroyed,
    e44Destroyed, e44Created, e44PoleDiag, e44PolePassive, e44Equator
}

/// Degree change per role (birth/death roles listed as 0; they are separate
/// channels, not degree increments on a persisting simplex).
immutable int[VRole.max + 1] vRoleDegreeDelta =
    [+2, 0, -2, 0, +2, 0, -2, 0, -2, +2, 0];
immutable int[ERole.max + 1] eRoleDegreeDelta =
    [-1, +1, 0, +1, -1, 0, +1, 0, -1, 0, 0, 0, 0, -1, +1];

/// True for roles that create/destroy the simplex itself.
immutable bool[VRole.max + 1] vRoleIsBirthDeath =
    [false, false, false, false, false, true, false, true, false, false, false];
immutable bool[ERole.max + 1] eRoleIsBirthDeath =
    [false, false, true, false, false, true, false, true, false, true,
     true, true, false, false, false];

enum eventRecordBytes = 36;   // u64 clock + u32 type + 6 x i32 labels, packed
enum sixRecordBytes = 20;     // u64 clock + i32 u + i32 v + i32 dir, packed

struct GeometryLedger(Vertex)
{
    bool trackRoles;    // role-resolved AAs + tet aggregates
    bool logEvents;     // fixed-size event records
    bool logSixFlips;   // six-edge (degree 5<->6 crossing) flip records

    double[Vertex][VRole.max + 1] vertexRoles;
    double[Vertex[2]][ERole.max + 1] edgeRoles;

    // Tets: aggregates only (identities churn too fast to ledger usefully).
    ulong[5] tetsCreated;        // by move type code
    ulong[5] tetsDestroyed;
    ulong[Vertex[4]] tetBirth;   // living tets -> birth clock
    ulong[64] tetLifetimeHist;   // log2-binned age (in attempted moves)
    ulong tetCensoredDeaths;     // destroyed tets born before tracking began

    ulong clock;                 // attempted moves since tracking enabled

    ubyte[] eventBuf;
    size_t eventUsed;
    bool eventOverflow;

    ubyte[] sixBuf;
    size_t sixUsed;
    bool sixOverflow;

    void clearRoles()()
    {
        foreach (ref aa; vertexRoles) aa = null;
        foreach (ref aa; edgeRoles) aa = null;
        tetsCreated[] = 0; tetsDestroyed[] = 0;
        tetBirth = null; tetLifetimeHist[] = 0;
        tetCensoredDeaths = 0;
        clock = 0;
    }
}

private void bump(K)(ref double[K] aa, const K key)
{
    aa[key] = aa.get(key, 0.0) + 1.0;
}

private Vertex[2] mkEdge(Vertex)(const Vertex a, const Vertex b)
{
    Vertex[2] e = a < b ? [a, b] : [b, a];
    return e;
}

private void tetCreate(Vertex)(ref GeometryLedger!Vertex g, Vertex[4] key)
{
    key[].sort();
    g.tetBirth[key] = g.clock;
}

private void tetDestroy(Vertex)(ref GeometryLedger!Vertex g, Vertex[4] key)
{
    import core.bitop : bsr;
    key[].sort();
    if (auto p = key in g.tetBirth)
    {
        immutable age = g.clock - *p;
        g.tetLifetimeHist[age == 0 ? 0 : bsr(age) + 1]++;
        g.tetBirth.remove(key);
    }
    else
        g.tetCensoredDeaths++;   // born before tracking started
}

/// Record an accepted bistellar move. center/coCenter as in BistellarMove.
void recordBistellar(Vertex)(ref GeometryLedger!Vertex g,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    immutable typeCode = cast(int) coCenter.length - 1;
    final switch (typeCode)
    {
    case 0: // 1->4: center = tet, coCenter = created vertex
        foreach (v; center) bump(g.vertexRoles[VRole.v14Base], v);
        bump(g.vertexRoles[VRole.v14Created], coCenter[0]);
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                bump(g.edgeRoles[ERole.e14Base], mkEdge(center[i], center[j]));
        foreach (v; center)
            bump(g.edgeRoles[ERole.e14Created], mkEdge(v, coCenter[0]));
        g.tetsDestroyed[0]++;
        {
            Vertex[4] t = void;
            foreach (i; 0 .. 4) t[i] = center[i];
            tetDestroy(g, t);
        }
        g.tetsCreated[0] += 4;
        foreach (skip; 0 .. 4)
        {
            Vertex[4] t = void; size_t n = 0;
            foreach (i, v; center) if (i != skip) t[n++] = v;
            t[3] = coCenter[0];
            tetCreate(g, t);
        }
        break;
    case 1: // 2->3: center = triangle (equator), coCenter = poles
        foreach (v; center) bump(g.vertexRoles[VRole.v23Equator], v);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v23Pole], v);
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                bump(g.edgeRoles[ERole.e23Equator], mkEdge(center[i], center[j]));
        foreach (c; center)
            foreach (p; coCenter)
                bump(g.edgeRoles[ERole.e23Spoke], mkEdge(c, p));
        bump(g.edgeRoles[ERole.e23Created], mkEdge(coCenter[0], coCenter[1]));
        g.tetsDestroyed[1] += 2;
        foreach (p; coCenter)
        {
            Vertex[4] t = [center[0], center[1], center[2], p];
            tetDestroy(g, t);
        }
        g.tetsCreated[1] += 3;
        foreach (skip; 0 .. 3)
        {
            Vertex[4] t = void; size_t n = 0;
            foreach (i, v; center) if (i != skip) t[n++] = v;
            t[2] = coCenter[0]; t[3] = coCenter[1];
            tetCreate(g, t);
        }
        break;
    case 2: // 3->2: center = edge (poles), coCenter = triangle (equator)
        foreach (v; center) bump(g.vertexRoles[VRole.v32Pole], v);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v32Equator], v);
        bump(g.edgeRoles[ERole.e32Destroyed], mkEdge(center[0], center[1]));
        foreach (c; center)
            foreach (q; coCenter)
                bump(g.edgeRoles[ERole.e32Spoke], mkEdge(c, q));
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                bump(g.edgeRoles[ERole.e32Triangle], mkEdge(coCenter[i], coCenter[j]));
        g.tetsDestroyed[2] += 3;
        foreach (skip; 0 .. 3)
        {
            Vertex[4] t = void; size_t n = 0;
            foreach (i, v; coCenter) if (i != skip) t[n++] = v;
            t[2] = center[0]; t[3] = center[1];
            tetDestroy(g, t);
        }
        g.tetsCreated[2] += 2;
        foreach (p; center)
        {
            Vertex[4] t = [coCenter[0], coCenter[1], coCenter[2], p];
            tetCreate(g, t);
        }
        break;
    case 3: // 4->1: center = destroyed vertex, coCenter = base tet
        bump(g.vertexRoles[VRole.v41Destroyed], center[0]);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v41Base], v);
        foreach (v; coCenter)
            bump(g.edgeRoles[ERole.e41Destroyed], mkEdge(center[0], v));
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                bump(g.edgeRoles[ERole.e41Base], mkEdge(coCenter[i], coCenter[j]));
        g.tetsDestroyed[3] += 4;
        foreach (skip; 0 .. 4)
        {
            Vertex[4] t = void; size_t n = 0;
            foreach (i, v; coCenter) if (i != skip) t[n++] = v;
            t[3] = center[0];
            tetDestroy(g, t);
        }
        g.tetsCreated[3]++;
        {
            Vertex[4] t = void;
            foreach (i; 0 .. 4) t[i] = coCenter[i];
            tetCreate(g, t);
        }
        break;
    }
}

/// Record an accepted 4-4 hinge move.
void recordHinge(Vertex)(ref GeometryLedger!Vertex g,
    const Vertex[2] removedEdge, const Vertex[2] addedEdge,
    const Vertex[4] linkCycleIn)
{
    // Rotate the cycle so the added diagonal is (cycle[0], cycle[2]).
    Vertex[4] lc = linkCycleIn;
    immutable d0 = mkEdge(lc[0], lc[2]);
    if (!(d0[0] == min(addedEdge[0], addedEdge[1])
          && d0[1] == max(addedEdge[0], addedEdge[1])))
        lc = [linkCycleIn[1], linkCycleIn[2], linkCycleIn[3], linkCycleIn[0]];

    foreach (v; removedEdge) bump(g.vertexRoles[VRole.v44Pole], v);
    bump(g.vertexRoles[VRole.v44Diag], lc[0]);
    bump(g.vertexRoles[VRole.v44Diag], lc[2]);
    bump(g.vertexRoles[VRole.v44Passive], lc[1]);
    bump(g.vertexRoles[VRole.v44Passive], lc[3]);

    bump(g.edgeRoles[ERole.e44Destroyed], mkEdge(removedEdge[0], removedEdge[1]));
    bump(g.edgeRoles[ERole.e44Created], mkEdge(lc[0], lc[2]));
    foreach (p; removedEdge)
    {
        bump(g.edgeRoles[ERole.e44PoleDiag], mkEdge(p, lc[0]));
        bump(g.edgeRoles[ERole.e44PoleDiag], mkEdge(p, lc[2]));
        bump(g.edgeRoles[ERole.e44PolePassive], mkEdge(p, lc[1]));
        bump(g.edgeRoles[ERole.e44PolePassive], mkEdge(p, lc[3]));
    }
    foreach (i; 0 .. 4)
        bump(g.edgeRoles[ERole.e44Equator], mkEdge(lc[i], lc[(i + 1) % 4]));

    g.tetsDestroyed[4] += 4;
    foreach (i; 0 .. 4)
    {
        Vertex[4] t = [removedEdge[0], removedEdge[1], linkCycleIn[i],
                       linkCycleIn[(i + 1) % 4]];
        tetDestroy(g, t);
    }
    g.tetsCreated[4] += 4;
    foreach (p; removedEdge)
    {
        Vertex[4] t1 = [p, lc[0], lc[1], lc[2]];
        Vertex[4] t2 = [p, lc[0], lc[2], lc[3]];
        tetCreate(g, t1);
        tetCreate(g, t2);
    }
}

/// Append one fixed-size event record (see eventRecordBytes layout).
void logEvent(Vertex)(ref GeometryLedger!Vertex g, int typeCode,
    scope const(Vertex)[] labelsA, scope const(Vertex)[] labelsB)
{
    if (g.eventUsed + eventRecordBytes > g.eventBuf.length)
    {
        g.eventOverflow = true;
        return;
    }
    import core.stdc.string : memcpy;
    auto p = g.eventBuf.ptr + g.eventUsed;
    immutable ulong clk = g.clock;
    immutable uint tc = cast(uint) typeCode;
    memcpy(p, &clk, 8);
    memcpy(p + 8, &tc, 4);
    int[6] lab = -1;
    size_t n = 0;
    foreach (v; labelsA) lab[n++] = cast(int) v;
    foreach (v; labelsB) lab[n++] = cast(int) v;
    memcpy(p + 12, lab.ptr, 24);
    g.eventUsed += eventRecordBytes;
}

/// Append one six-edge flip record: (clock, u < v, dir). dir = +1 when edge
/// (u, v) crosses degree 5 -> 6 (a disclination-line edge is born), -1 for
/// 6 -> 5 (one dies). The stream is the complete rewiring history of the
/// disclination network: E6(t) = E6(0) + sum of dir over records.
private void logSixFlip(Vertex)(ref GeometryLedger!Vertex g,
    const Vertex a, const Vertex b, int dir)
{
    if (g.sixUsed + sixRecordBytes > g.sixBuf.length)
    {
        g.sixOverflow = true;
        return;
    }
    import core.stdc.string : memcpy;
    auto p = g.sixBuf.ptr + g.sixUsed;
    immutable ulong clk = g.clock;
    immutable int u = cast(int) min(a, b);
    immutable int v = cast(int) max(a, b);
    memcpy(p, &clk, 8);
    memcpy(p + 8, &u, 4);
    memcpy(p + 12, &v, 4);
    memcpy(p + 16, &dir, 4);
    g.sixUsed += sixRecordBytes;
}

/// Emit six-edge flip records for an ACCEPTED bistellar move. Must be called
/// BEFORE the move is applied (reads pre-move degrees). Only persisting edges
/// can cross the 5/6 threshold: per the role tables, created edges are born
/// at degree 3 and destroyed edges die at degree 3, so the six-edge graph
/// changes exactly at the +-1 degree crossings enumerated here.
void sixFlipsBistellar(Vertex)(ref GeometryLedger!Vertex g,
    const ref Manifold!(3, Vertex) mfd,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    void crossing(const Vertex a, const Vertex b, int delta)
    {
        Vertex[2] e = a < b ? [a, b] : [b, a];
        immutable d = cast(long) mfd.degree(e[]);
        if (delta > 0 && d == 5) logSixFlip(g, a, b, +1);
        else if (delta < 0 && d == 6) logSixFlip(g, a, b, -1);
    }
    final switch (cast(int) coCenter.length - 1)
    {
    case 0: // 1->4: the 6 base edges gain a facet
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                crossing(center[i], center[j], +1);
        break;
    case 1: // 2->3: equator edges -1, spokes +1
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                crossing(center[i], center[j], -1);
        foreach (cv; center)
            foreach (p; coCenter)
                crossing(cv, p, +1);
        break;
    case 2: // 3->2: coCenter triangle edges +1, spokes -1
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                crossing(coCenter[i], coCenter[j], +1);
        foreach (cv; center)
            foreach (q; coCenter)
                crossing(cv, q, -1);
        break;
    case 3: // 4->1: the 6 base edges lose a facet
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                crossing(coCenter[i], coCenter[j], -1);
        break;
    }
}

/// Hinge (4-4) analog: equator edges +1, pole-passive edges -1. The removed
/// edge dies at degree 4 and the diagonal is born at 4 — no crossings there.
void sixFlipsHinge(Vertex)(ref GeometryLedger!Vertex g,
    const ref Manifold!(3, Vertex) mfd, const ref HingeMove!Vertex hm)
{
    void crossing(const Vertex a, const Vertex b, int delta)
    {
        Vertex[2] e = a < b ? [a, b] : [b, a];
        immutable d = cast(long) mfd.degree(e[]);
        if (delta > 0 && d == 5) logSixFlip(g, a, b, +1);
        else if (delta < 0 && d == 6) logSixFlip(g, a, b, -1);
    }
    foreach (i; 0 .. 4)
        crossing(hm.linkCycle[i], hm.linkCycle[(i + 1) % 4], +1);
    foreach (p; hm.removedEdge)
        foreach (v; hm.linkCycle)
            if (v != hm.addedEdge[0] && v != hm.addedEdge[1])
                crossing(p, v, -1);
}

///
unittest
{
    // Single controlled moves against a live manifold: measured degree changes
    // must reproduce the (type, role) tables exactly.
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    GeometryLedger!int g;

    long deg(int[] s) { return mfd.degree(s); }

    // --- 1->4 on facet [0,1,2,3], new vertex 5 ---
    int[] c14 = [0,1,2,3];
    auto degBefore = [deg([0]), deg([1]), deg([2]), deg([3])];
    auto edgeBefore = deg([0,1]);
    recordBistellar(g, c14, [5]);
    mfd.doMove(BistellarMove!3([0,1,2,3], [5]));
    foreach (i, v; [0,1,2,3])
        assert(deg([v]) - degBefore[i] == vRoleDegreeDelta[VRole.v14Base]);
    assert(deg([5]) == 4);                                   // born at 4
    assert(deg([0,1]) - edgeBefore == eRoleDegreeDelta[ERole.e14Base]);
    assert(deg([0,5]) == 3);                                 // spoke born at 3
    assert(g.vertexRoles[VRole.v14Created][5] == 1.0);
    assert(g.tetsCreated[0] == 4 && g.tetsDestroyed[0] == 1);

    // --- 2->3 on triangle [0,1,2] (in tets [0,1,2,4],[0,1,2,5]? use valid) ---
    // triangle [0,1,2] now has degree 2 (tets [0,1,2,4] and [0,1,2,5]).
    assert(deg([0,1,2]) == 2);
    auto dp4 = deg([4]); auto dp5 = deg([5]); auto de01 = deg([0,1]);
    recordBistellar(g, [0,1,2], [4,5]);
    mfd.doMove(BistellarMove!3([0,1,2], [4,5]));
    assert(deg([4]) - dp4 == vRoleDegreeDelta[VRole.v23Pole]);
    assert(deg([5]) - dp5 == vRoleDegreeDelta[VRole.v23Pole]);
    assert(deg([0,1]) - de01 == eRoleDegreeDelta[ERole.e23Equator]);
    assert(deg([4,5]) == 3);                                 // pole-pole born at 3
    assert(g.edgeRoles[ERole.e23Created][mkEdge(4,5)] == 1.0);

    // --- role totals: each move contributes its exact multiplicities ---
    double tot(double[int] aa) { double s=0; foreach(v; aa.byValue) s+=v; return s; }
    assert(tot(g.vertexRoles[VRole.v14Base]) == 4.0);
    assert(tot(g.vertexRoles[VRole.v23Pole]) == 2.0);
    assert(tot(g.vertexRoles[VRole.v23Equator]) == 3.0);
    double etot(double[int[2]] aa) { double s=0; foreach(v; aa.byValue) s+=v; return s; }
    assert(etot(g.edgeRoles[ERole.e23Spoke]) == 6.0);
    assert(etot(g.edgeRoles[ERole.e14Created]) == 4.0);

    assert(mfd.findProblems.length == 0);
}

/******************************************************************************
Run one MCMC step using a unified proposal that naturally includes both
bistellar (Pachner) moves and 4-4 hinge moves.

Proposal: pick a random facet, pick a random sub-simplex (center). Check the
degree of the center. For edges (dim-1 center) with degree 4, propose a hinge
move instead of the (invalid) 3→2 bistellar move. All other cases follow the
standard Pachner move logic.

This avoids a separate hinge move probability parameter — hinge moves are
proposed at the natural rate determined by how many degree-4 edges exist.
*/
bool mcmcStep(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    ref Vertex[] unusedVertices,
    P params,
    real hingeMoveProb,  // ignored (kept for API compatibility)
    ref ulong hingeTries,
    ref ulong hingeAccepts,
    ref ulong[4] bistellarTries,
    ref ulong[4] bistellarAccepts,
    MoveCounters!Vertex* counters = null,
    GeometryLedger!Vertex* ledger = null,
    VertexPotState!Vertex* potState = null,
    const(VertexPot)* pot = null,
    CocycleState!Vertex* cocycle = null)
{
    enum dim = 3;
    enum nVerts = dim + 1;
    enum maxMask = (1 << nVerts) - 1;
    alias BM = BistellarMove!(dim, Vertex);

    if (ledger !is null)
        ledger.clock++;               // one tick per attempted move

    // Unified proposal loop
    while (true)
    {
        auto facet = mfd.randomFacetOfDim(dim);

        auto mask = uniform(1, maxMask + 1);
        Vertex[nVerts] centerBuf;
        int centerLen = 0;
        foreach (i; 0 .. nVerts)
        {
            if (mask & (1 << i))
                centerBuf[centerLen++] = facet[i];
        }
        auto center = centerBuf[0 .. centerLen];
        center.sort();

        auto centerDim = centerLen - 1;
        auto centerDeg = mfd.degree(center);

        // --- Edge of degree 4: propose hinge move ---
        if (centerDim == 1 && centerDeg == 4)
        {
            Vertex[2] edge = [center[0], center[1]];

            // Find a start vertex (any facet vertex not in the edge)
            Vertex startVertex = void;
            foreach (v; facet)
                if (v != edge[0] && v != edge[1]) { startVertex = v; break; }

            auto hm = mfd.hingeMove(edge, startVertex, uniform(0, 2));

            // Support = the 6 vertices whose stars the move touches.
            Vertex[6] hingeSupport = void;
            hingeSupport[0 .. 2] = hm.removedEdge[];
            hingeSupport[2 .. 6] = hm.linkCycle[];
            if (counters !is null)
                addSupport(counters.proposed, hingeSupport[]);

            // Frozen-region rejection: every facet this move adds/removes has
            // all its vertices in the support, so rejecting here preserves the
            // frozen set's closed star exactly (facets + hinge degrees).
            if (mfd.anyFrozen(hingeSupport[]))
                continue;

            if (!mfd.hasValidHingeMove(hm))
                continue;

            hingeTries++;
            if (counters !is null)
                addSupport(counters.valid, hingeSupport[]);
            // currentObjective includes the potential total (when enabled); the
            // base speculative delta needs the base-only objective as baseline.
            real baseObj = currentObjective
                - (potState !is null ? potState.total : 0.0L);
            real deltaObj = mfd.speculativeHingeDelta(hm, baseObj, params);
            if (potState !is null)
                deltaObj += mfd.potentialHingeDelta(hm, *potState, *pot, false);
            real logAlpha = -deltaObj;

            if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
            {
                if (potState !is null)
                    mfd.potentialHingeDelta(hm, *potState, *pot, true);
                if (ledger !is null && ledger.logSixFlips)
                    sixFlipsHinge(*ledger, mfd, hm);
                if (cocycle !is null)
                    cocycleHinge(*cocycle, hm);
                mfd.doHingeMove(hm);
                currentObjective += deltaObj;
                hingeAccepts++;
                if (counters !is null)
                    addSupport(counters.acceptedHinge, hingeSupport[]);
                if (ledger !is null)
                {
                    if (ledger.trackRoles)
                        recordHinge(*ledger, hm.removedEdge, hm.addedEdge,
                                    hm.linkCycle);
                    if (ledger.logEvents)
                    {
                        // Rotate cycle so the diagonal is (labels[2], labels[4]).
                        Vertex[4] lc = hm.linkCycle;
                        immutable d0 = mkEdge(lc[0], lc[2]);
                        if (!(d0[0] == min(hm.addedEdge[0], hm.addedEdge[1])
                              && d0[1] == max(hm.addedEdge[0], hm.addedEdge[1])))
                            lc = [hm.linkCycle[1], hm.linkCycle[2],
                                  hm.linkCycle[3], hm.linkCycle[0]];
                        logEvent(*ledger, 4, hm.removedEdge[], lc[]);
                    }
                }
                return true;
            }
            return false;
        }

        // --- Standard bistellar move ---
        if (centerDeg + centerDim != dim + 1)
            continue;

        BM bm;
        if (centerDim == dim)
        {
            if (unusedVertices.empty)
                unusedVertices ~= mfd.fVector[0].to!Vertex;
            bm = BM(center, unusedVertices.back.only);
        }
        else
        {
            auto coCenter = mfd.coCenter(center, facet);
            bm = BM(center, coCenter[]);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        // Support = the 5 vertices of the bistellar ball (one glued 4-simplex).
        Vertex[nVerts + 1] ballBuf = void;
        Vertex[] ball;
        {
            size_t nb = 0;
            foreach (v; bm.center) ballBuf[nb++] = v;
            foreach (v; bm.coCenter) ballBuf[nb++] = v;
            ball = ballBuf[0 .. nb];
        }
        if (counters !is null)
            addSupport(counters.proposed, ball);

        // Frozen-region rejection (see hinge branch). For a 1->4 move the
        // coCenter vertex is new/unused, hence never frozen.
        if (mfd.anyFrozen(ball))
            continue;

        if (!mfd.hasValidMove(bm))
            continue;

        bistellarTries[bm.coCenter.length - 1]++;
        if (counters !is null)
            addSupport(counters.valid, ball);
        real baseObj = currentObjective
            - (potState !is null ? potState.total : 0.0L);
        real deltaObj = mfd.speculativeBistellarDelta(bm, baseObj, params);
        if (potState !is null)
            deltaObj += mfd.potentialBistellarDelta(bm, *potState, *pot, false);
        real logAlpha = -deltaObj;

        if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
        {
            if (potState !is null)
                mfd.potentialBistellarDelta(bm, *potState, *pot, true);
            if (ledger !is null && ledger.logSixFlips)
                sixFlipsBistellar(*ledger, mfd, bm.center, bm.coCenter);
            if (cocycle !is null)
                cocycleBistellar(*cocycle, bm.center, bm.coCenter);
            mfd.doMove(bm);
            if (counters !is null)
                addSupport(counters.acceptedBistellar, ball);
            if (ledger !is null)
            {
                if (ledger.trackRoles)
                    recordBistellar(*ledger, bm.center, bm.coCenter);
                if (ledger.logEvents)
                    logEvent(*ledger, cast(int) bm.coCenter.length - 1,
                             bm.center, bm.coCenter);
            }
            if (bm.coCenter.length == 1)
            {
                // Shrink, then tell the runtime we own the freed slot so the
                // next `~=` reuses the buffer instead of reallocating a fresh
                // int[] every accepted move (that churn was false-pinned by the
                // conservative GC, leaking ~0.14 MB/sweep). Mirrors ddg_capi.d.
                unusedVertices.popBack;
                unusedVertices.assumeSafeAppend;
            }
            if (bm.center.length == 1) unusedVertices ~= bm.center;
            currentObjective += deltaObj;
            bistellarAccepts[bm.coCenter.length - 1]++;
            return true;
        }
        return false;
    }
}

///
unittest
{
    // Integration test: run mixed MCMC for a number of steps, verify manifold integrity
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    // Grow a bit first
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));

    auto currentObj = mfd.objective(params);
    int[] unusedVertices = [8];
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    MoveCounters!int mc;

    int accepted = 0;
    foreach (_; 0 .. 500)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts, &mc))
            accepted++;
    }

    assert(accepted > 0, "should have accepted some moves");
    assert(hingeTries > 0, "should have attempted some hinge moves");

    // Move-counter invariants: every event distributes total weight 1 over its
    // support, so ledger totals must reproduce the scalar tallies exactly.
    double total(double[int] aa)
    {
        double s = 0;
        foreach (v; aa.byValue) s += v;
        return s;
    }
    assert(isClose(total(mc.valid), cast(double)(hingeTries + bTries[].sum), 0.0, 1e-6),
        "valid ledger total != tries");
    assert(isClose(total(mc.acceptedBistellar), cast(double) bAccepts[].sum, 0.0, 1e-6),
        "acceptedBistellar ledger total != bistellar accepts");
    assert(isClose(total(mc.acceptedHinge), cast(double) hingeAccepts, 0.0, 1e-6),
        "acceptedHinge ledger total != hinge accepts");
    assert(total(mc.proposed) >= total(mc.valid) - 1e-9,
        "proposed must dominate valid");

    // Verify manifold integrity after many mixed moves
    assert(mfd.findProblems.length == 0,
        "manifold has problems after mixed MCMC: " ~ mfd.findProblems.to!string);

    // Verify objective tracking is consistent
    auto actualObj = mfd.objective(params);
    assert(isClose(currentObj, actualObj, 1e-4),
        "objective drift: tracked=%s actual=%s".format(currentObj, actualObj));
}

/// Six-flip stream + disclination census: the flip records must exactly
/// account for the six-edge count change, and the census must satisfy its
/// internal identities on a churned manifold.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 64;
        real hingeDegreeTarget = 5.1;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;

    long countSixEdges()
    {
        long n = 0;
        foreach (e; mfd.simplices(1))
            if (mfd.degree(e) >= 6) ++n;
        return n;
    }

    GeometryLedger!int g;
    g.logSixFlips = true;
    g.sixBuf = new ubyte[1 << 20];

    immutable e6Before = countSixEdges();
    foreach (_; 0 .. 4000)
        mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
            hingeTries, hingeAccepts, bTries, bAccepts, null, &g);
    immutable e6After = countSixEdges();

    assert(!g.sixOverflow, "six-flip buffer overflowed in test");
    long net = 0;
    ulong lastClock = 0;
    for (size_t p = 0; p + sixRecordBytes <= g.sixUsed; p += sixRecordBytes)
    {
        import core.stdc.string : memcpy;
        ulong clk; int u, v, dir;
        memcpy(&clk, g.sixBuf.ptr + p, 8);
        memcpy(&u, g.sixBuf.ptr + p + 8, 4);
        memcpy(&v, g.sixBuf.ptr + p + 12, 4);
        memcpy(&dir, g.sixBuf.ptr + p + 16, 4);
        assert(dir == 1 || dir == -1, "bad flip direction");
        assert(u < v, "flip labels not sorted");
        assert(clk >= lastClock, "clock not monotone");
        lastClock = clk;
        net += dir;
    }
    assert(e6After == e6Before + net,
        "six-flip stream does not account for six-edge count change: "
        ~ "%s -> %s but net flips %s".format(e6Before, e6After, net));
    assert(e6After > 0, "test manifold never grew six-edges (weak test)");

    // Census identities on the churned manifold.
    auto c = mfd.disclinationCensus(0);
    c.nSixEdges.shouldEqual(e6After);
    long degTot = 0, vertTot = 0;
    foreach (k; 0 .. 8) { degTot += (k + 1) * c.netDegCensus[k]; vertTot += c.netDegCensus[k]; }
    vertTot.shouldEqual(c.nNetVerts);
    // netDegCensus is clamped at 8+, but a churned 64-facet manifold should
    // not produce network degree >= 8; if it does, skip the handshake check.
    if (c.netDegCensus[7] == 0)
        degTot.shouldEqual(2 * c.nSixEdges);
    (c.sumSegLen + c.sumLoopLen).shouldEqual(c.nSixEdges);
    c.cycleRank.shouldEqual(c.nSixEdges - c.nNetVerts + c.nComponents);
    assert(c.cycleRank >= 0);
    assert(c.giantSize >= c.secondSize);
    assert(c.giantDiameter <= c.giantSize);
    long compTot = 0;
    foreach (x; c.compSizeHist) compTot += x;
    compTot.shouldEqual(c.nComponents);

    // Valence census cross-check: rows n6 >= 1 must reproduce nNetVerts and
    // the six-edge handshake, with caps high enough to avoid clamping.
    enum CAP = 31;
    long[(CAP + 1) * (CAP + 1)] vc;
    mfd.valenceCensus(vc[], CAP, CAP);
    long f0Tot = 0, netTot = 0, handshake = 0;
    foreach (k; 0 .. CAP + 1)
        foreach (m; 0 .. CAP + 1)
        {
            immutable cnt = vc[k * (CAP + 1) + m];
            f0Tot += cnt;
            if (k >= 1) { netTot += cnt; handshake += k * cnt; }
        }
    f0Tot.shouldEqual(mfd.fVector[0]);
    netTot.shouldEqual(c.nNetVerts);
    handshake.shouldEqual(2 * c.nSixEdges);

    // Flattened layout must reproduce the struct (the C API path).
    long[disclinationCensusSlots] flat;
    flattenCensus(c, flat[]);
    flat[0].shouldEqual(c.nNetVerts);
    flat[1].shouldEqual(c.nSixEdges);
    flat[6].shouldEqual(c.cycleRank);
    flat[20].shouldEqual(c.eImpAny);
    flat[21].shouldEqual(0);   // reserved
    long flatDegTot = 0;
    foreach (k; 0 .. 8) flatDegTot += flat[24 + k];
    flatDegTot.shouldEqual(c.nNetVerts);
}

/// Frozen vertices: the sampler must preserve the frozen set's closed star
/// EXACTLY (facets and the degrees of every simplex meeting a frozen vertex)
/// while the rest of the manifold churns freely.
unittest
{
    import std.algorithm : canFind;

    struct TestParams
    {
        int numFacetsTarget = 48;
        real hingeDegreeTarget = 4.8;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    // Grow so there is bulk to churn away from the frozen star.
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));
    mfd.doMove(BM([0,1,2,5], [8]));

    // Freeze vertex 7 and snapshot its closed star + incident edge degrees.
    mfd.setVertexFrozen(7, true);
    assert(mfd.vertexFrozen(7));
    assert(!mfd.vertexFrozen(3));
    assert(mfd.numFrozenVertices == 1);

    int[4][] frozenStar;
    foreach (f; mfd.facets)
        if (f.canFind(7))
        {
            int[4] ff;
            ff[] = f[0 .. 4];
            frozenStar ~= ff;
        }
    size_t[int[2]] edgeDegs;      // degrees of edges at the frozen vertex
    foreach (e; mfd.simplices(1))
        if (e.canFind(7))
        {
            int[2] ee;
            ee[] = e[0 .. 2];
            edgeDegs[ee] = mfd.degree(e);
        }
    immutable frozenVtxDeg = mfd.degree([7]);

    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    int accepted = 0;
    foreach (step; 0 .. 4000)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts))
            accepted++;
    }
    assert(accepted > 50, "sampler froze entirely; test not exercising moves");

    // The frozen closed star is exactly intact.
    foreach (ff; frozenStar)
        assert(mfd.contains(ff[]), "frozen-star facet destroyed");
    size_t starCount = 0;
    foreach (f; mfd.facets)
        if (f.canFind(7)) starCount++;
    starCount.shouldEqual(frozenStar.length);   // no facet added at 7 either
    mfd.degree([7]).shouldEqual(frozenVtxDeg);
    foreach (ee, d; edgeDegs)
        mfd.degree(ee[]).shouldEqual(d);        // hinge curvature at 7 pinned

    // Unfreezing restores full dynamics (flag actually clears).
    mfd.clearFrozenVertices();
    assert(!mfd.vertexFrozen(7));
    assert(mfd.numFrozenVertices == 0);
    assert(mfd.findProblems.length == 0);
}

/// Vertex 6-valence potential: incremental counter state and running total
/// stay exact through mixed MCMC (all bistellar types + hinge moves), and the
/// tracked objective (base + potential) matches a from-scratch recompute.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 24;
        real hingeDegreeTarget = 4.6;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 0.0;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));

    VertexPot pot;
    pot.zlegCoef = 0.3;
    pot.impCoef = 0.05;
    pot.tilt = [0.02L, 0, -0.04L, 0.01L, -0.03L];  // exercise every tilt slot
    assert(pot.enabled);

    // U spot checks: legal set flat at tilt, n6=1/5 pay 1, quadratic tail.
    assert(isClose(pot.U(0), 0.02L));
    assert(isClose(pot.U(1), 0.3L + 0.0L));
    assert(isClose(pot.U(2), -0.04L));
    assert(isClose(pot.U(4), -0.03L));
    assert(isClose(pot.U(5), 0.3L));
    assert(isClose(pot.U(7), 0.3L * 9));
    assert(isClose(pot.V(6), 0.05L * 36));

    VertexPotState!int st;
    mfd.recomputeVertexPotState(st, pot);

    // Counter invariant: sum_v n6 = 2 * #edges(deg >= 6), same for impurity.
    void checkInvariants()
    {
        long s6 = 0, sI = 0;
        foreach (v, k; st.n6) s6 += k;
        foreach (v, m; st.mImp) sI += m;
        long e6 = 0, eI = 0;
        foreach (e; mfd.simplices(1))
        {
            immutable d = cast(long) mfd.degree(e);
            if (d >= 6) e6++;
            if (d < 5 || d > 6) eI++;
        }
        assert(s6 == 2 * e6, "n6 counter sum mismatch");
        assert(sI == 2 * eI, "impurity counter sum mismatch");
    }
    checkInvariants();

    auto currentObj = mfd.objective(params) + st.total;
    int[] unusedVertices = [8];
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;

    int accepted = 0;
    foreach (step; 0 .. 600)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts, null, null,
                &st, &pot))
            accepted++;

        if (step % 100 == 99)
        {
            // Incremental state == from-scratch recompute
            VertexPotState!int fresh;
            mfd.recomputeVertexPotState(fresh, pot);
            assert(fresh.n6 == st.n6,
                "n6 drift at step %s".format(step));
            assert(fresh.mImp == st.mImp,
                "impurity drift at step %s".format(step));
            assert(isClose(fresh.total, st.total, 1e-9L, 1e-9L),
                "potential total drift: tracked=%s fresh=%s"
                .format(st.total, fresh.total));
            checkInvariants();

            // Tracked objective == base + potential, recomputed
            assert(isClose(currentObj, mfd.objective(params) + st.total, 1e-6L),
                "objective drift with potential at step %s".format(step));
        }
    }
    assert(accepted > 0, "should have accepted some moves");
    assert(hingeAccepts + bAccepts[].sum > 0);
    assert(mfd.findProblems.length == 0);
}

// ---------------------------------------------------------------------------
// Move selection
// ---------------------------------------------------------------------------

BistellarMove!(dim, Vertex) chooseRandomMove(int dim, Vertex, P)(
    ref Manifold!(dim, Vertex) manifold, Vertex newVertex, P parameters)
{
    alias BM = BistellarMove!(dim, Vertex);
    enum nVerts = dim + 1;
    enum maxMask = (1 << nVerts) - 1; // 2^(dim+1) - 1

    while(true)
    {
        auto facet = manifold.randomFacetOfDim(dim);

        // Pick a random non-empty subset via bitmask (avoids materializing all subsets)
        auto mask = uniform(1, maxMask + 1);
        Vertex[nVerts] centerBuf;
        int centerLen = 0;
        foreach (i; 0 .. nVerts)
        {
            if (mask & (1 << i))
                centerBuf[centerLen++] = facet[i];
        }
        auto center = centerBuf[0 .. centerLen];
        center.sort();

        auto centerDim = centerLen - 1;
        auto centerDeg = manifold.degree(center);

        if (centerDeg + centerDim != dim + 1)
            continue;

        BM bm;
        if (centerDim == dim)
        {
            bm = BM(center, newVertex.only);
        }
        else
        {
            auto coCenter = manifold.coCenter(center, facet);
            bm = BM(center, coCenter[]);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        if (!manifold.hasValidMove(bm))
            continue;

        return bm;
    }
}
///
@safe unittest
{
    // Smoke test: chooseRandomMove should return without hanging
    auto rp3 = Manifold!3([[1, 2, 3, 7], [1, 2, 3, 11], [1, 2, 6, 9], [1,
            2, 6, 11], [1, 2, 7, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 7,
            10], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 8, 9], [1, 4, 8, 10], [1,
            5, 6, 8], [1, 5, 6, 11], [1, 5, 8, 10], [1, 6, 8, 9], [2, 3, 4, 8],
            [2, 3, 4, 11], [2, 3, 7, 8], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4,
            8, 10], [2, 5, 7, 8], [2, 5, 7, 9], [2, 5, 8, 10], [2, 5, 9, 10],
            [2, 6, 9, 10], [3, 4, 5, 9], [3, 4, 5, 11], [3, 4, 8, 9], [3, 5, 9,
            10], [3, 6, 7, 8], [3, 6, 7, 10], [3, 6, 8, 9], [3, 6, 9, 10], [4,
            5, 6, 7], [4, 5, 6, 11], [4, 5, 7, 9], [4, 6, 7, 10], [5, 6, 7, 8]]);
}

// ---------------------------------------------------------------------------
// Unused vertex management
// ---------------------------------------------------------------------------

Vertex[] getUnusedVertices(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd, Vertex[] initialVertices)
{
    Vertex[] unusedVertices;
    // all gaps in list of vertices should be unused vertices
    if (initialVertices.front != 0)
    {
        unusedVertices ~= initialVertices.front.iota.array;
    }
    foreach (i; 0 .. initialVertices.length - 1)
    {
        if (initialVertices[i] + 1 != initialVertices[i + 1])
        {
            unusedVertices ~= iota(initialVertices[i] + 1, initialVertices[i + 1]).array;
        }
    }
    assert(unusedVertices.all!(v => !mfd.contains(v.only)));
    return unusedVertices;
}

/// Convenience overload: compute initial vertices from the manifold.
Vertex[] getUnusedVertices(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    auto verts = mfd.simplices(0).joiner.array.dup.sort.array;
    if (verts.length == 0) return [];
    return getUnusedVertices(mfd, verts);
}
