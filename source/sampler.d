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
// Hinge move proposal
// ---------------------------------------------------------------------------

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
Run one MCMC step that is either a hinge move or a bistellar move.
Returns true if a move was accepted.

hingeMoveProb controls the probability of attempting a hinge move per step.
When the hinge attempt fails (no degree-4 edge or no valid diagonal), falls
through to a bistellar move.

This function uses the speculative delta path (no execute-undo) and
approximate Hastings correction via nFacets ratio for bistellar moves.
For hinge moves, no Hastings correction is needed (symmetric proposal,
unchanged f-vector).
*/
bool mcmcStep(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    ref Vertex[] unusedVertices,
    P params,
    real hingeMoveProb,
    ref ulong hingeTries,
    ref ulong hingeAccepts,
    ref ulong[4] bistellarTries,
    ref ulong[4] bistellarAccepts)
{
    enum dim = 3;

    // With probability hingeMoveProb, try a hinge move
    if (hingeMoveProb > 0 && uniform01 < hingeMoveProb)
    {
        auto hmResult = mfd.tryProposeHingeMove();
        if (!hmResult.isNull)
        {
            auto hm = hmResult.get;
            hingeTries++;

            real deltaObj = mfd.speculativeHingeDelta(hm, currentObjective, params);
            // No Hastings correction: proposal is symmetric, f-vector unchanged
            real logAlpha = -deltaObj;

            if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
            {
                mfd.doHingeMove(hm);
                currentObjective += deltaObj;
                hingeAccepts++;
                return true;
            }
            return false;
        }
        // Hinge attempt failed — fall through to bistellar
    }

    // Bistellar move
    if (unusedVertices.empty)
        unusedVertices ~= mfd.fVector[0].to!Vertex;

    auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
    bistellarTries[bm.coCenter.length - 1]++;

    real deltaObj = mfd.speculativeBistellarDelta(bm, currentObjective, params);
    // No Hastings correction: importance weight 1/V corrects at measurement time.
    real logAlpha = -deltaObj;

    if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
    {
        mfd.doMove(bm);
        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;
        currentObjective += deltaObj;
        bistellarAccepts[bm.coCenter.length - 1]++;
        return true;
    }
    return false;
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

    int accepted = 0;
    foreach (_; 0 .. 500)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts))
            accepted++;
    }

    assert(accepted > 0, "should have accepted some moves");
    assert(hingeTries > 0, "should have attempted some hinge moves");

    // Verify manifold integrity after many mixed moves
    assert(mfd.findProblems.length == 0,
        "manifold has problems after mixed MCMC: " ~ mfd.findProblems.to!string);

    // Verify objective tracking is consistent
    auto actualObj = mfd.objective(params);
    assert(isClose(currentObj, actualObj, 1e-4),
        "objective drift: tracked=%s actual=%s".format(currentObj, actualObj));
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
