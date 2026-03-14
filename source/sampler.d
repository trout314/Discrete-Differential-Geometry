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

    static if (dim_ > 2)
    {
        immutable coDim3DegTarget = coDim3PerFacet * nFacets / cast(real) nCoDim3;
        x = modf(coDim3DegTarget, _);
        minPenalty = (x - x ^^ 2) * nCoDim3;

        penalty.localSolidAngleCurvPenalty = (
            coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + coDim3TotSqDeg) - minPenalty;
        // Intensive: divide by count to get mean variance per codim-3 face
        if (nCoDim3 > 0) penalty.localSolidAngleCurvPenalty /= nCoDim3;
    }
    else
    {
        penalty.localSolidAngleCurvPenalty = 0;
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
        + params.coDim3DegreeVarianceCoef * pen.localSolidAngleCurvPenalty;
}

real objective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    auto pen = mfd.penalties(params);
    return objectiveFromPenalty(pen, params);
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
// Move selection
// ---------------------------------------------------------------------------

BistellarMove!(dim, Vertex) chooseRandomMove(int dim, Vertex, P)(
    Manifold!(dim, Vertex) manifold, Vertex newVertex, P parameters)
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
            bm = BM(center, coCenter);
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
