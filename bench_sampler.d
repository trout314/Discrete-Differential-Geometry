/// Benchmark for the MCMC sampler hot path.
/// Compares baseline (execute-compute-undo) vs speculative delta with full objective.
module bench_sampler;

version (unittest) {} else {

import std.algorithm, std.array, std.conv, std.datetime.stopwatch, std.math,
    std.range, std.stdio, std.string, std.sumtype;
import std.random : uniform, uniform01, rndChoice = choice;
import manifold, manifold_examples, manifold_moves, simplicial_complex, utility;

enum dim = 3;

struct BenchParams
{
    int numFacetsTarget;
    real hingeDegreeTarget = 4.5;
    real numFacetsCoef = 0.1;
    real numHingesCoef = 0.05;
    real hingeDegreeVarianceCoef = 0.2;
    real coDim3DegreeVarianceCoef = 0.1;
    bool useHingeMoves = false;
}

// --- Penalty computation (same formulas as manifold_sampler.d) ---

private struct Penalty
{
    real volumePenalty;
    real globalCurvPenalty;
    real localCurvPenalty;
    real localSolidAngleCurvPenalty;
}

private Penalty penaltiesFromValues(
    long nFacets, long nHinges, ulong hingeTotSqDeg,
    long nCoDim3, ulong coDim3TotSqDeg, BenchParams params)
{
    enum hingesPerFacet = dim * (dim + 1) / 2;
    enum coDim3PerFacet = binomial(dim + 1, dim - 2);

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

    immutable coDim3DegTarget = coDim3PerFacet * nFacets / cast(real) nCoDim3;
    x = modf(coDim3DegTarget, _);
    minPenalty = (x - x ^^ 2) * nCoDim3;

    penalty.localSolidAngleCurvPenalty = (
        coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + coDim3TotSqDeg) - minPenalty;

    return penalty;
}

real computeObjective(ref const Manifold!dim mfd, BenchParams params)
{
    auto pen = penaltiesFromValues(
        cast(long) mfd.fVector[dim], cast(long) mfd.fVector[dim - 2],
        mfd.totalSquareDegree(dim - 2),
        cast(long) mfd.fVector[dim - 3], mfd.totalSquareDegree(dim - 3),
        params);
    return params.numFacetsCoef * pen.volumePenalty
        + params.numHingesCoef * pen.globalCurvPenalty
        + params.hingeDegreeVarianceCoef * pen.localCurvPenalty
        + params.coDim3DegreeVarianceCoef * pen.localSolidAngleCurvPenalty;
}

real speculativeDelta(
    ref const Manifold!dim mfd,
    ref const BistellarMove!(dim, int) move,
    real currentObjective,
    BenchParams params)
{
    auto center = move.center;
    auto coCenter = move.coCenter;
    immutable cenLen = cast(int) center.length;
    immutable coCenLen = cast(int) coCenter.length;

    int[dim + 2] allVertsBuf;
    allVertsBuf[0 .. cenLen] = center[];
    allVertsBuf[cenLen .. cenLen + coCenLen] = coCenter[];
    auto allVerts = allVertsBuf[0 .. cenLen + coCenLen];
    allVerts.sort();

    size_t[dim + 1] newFVector = mfd.fVector[0 .. dim + 1];
    newFVector[].modifyFVector(move);

    long[dim - 1] newTotSqDeg;
    foreach (d; 0 .. dim - 1)
        newTotSqDeg[d] = cast(long) mfd.totalSquareDegree(d);

    static foreach (d; 0 .. dim - 1)
    {{
        foreach (subset; allVerts[].subsetsOfSize(d + 1))
        {
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

    auto newPen = penaltiesFromValues(
        cast(long) newFVector[dim], cast(long) newFVector[dim - 2],
        cast(ulong) newTotSqDeg[dim - 2],
        cast(long) newFVector[dim - 3],
        cast(ulong) newTotSqDeg[dim - 3],
        params);

    real newObj = params.numFacetsCoef * newPen.volumePenalty
        + params.numHingesCoef * newPen.globalCurvPenalty
        + params.hingeDegreeVarianceCoef * newPen.localCurvPenalty
        + params.coDim3DegreeVarianceCoef * newPen.localSolidAngleCurvPenalty;

    return newObj - currentObjective;
}

// --- Main ---

void main()
{
    writeln("=== MCMC Sampler Benchmark ===");
    writeln("Objective: volume + global curvature + hinge variance + codim-3 variance");
    writeln();
    stdout.flush();

    auto seed = loadManifold!dim("standard_triangulations/dim_3_sphere.mfd");

    foreach (targetTets; [100, 500, 2000])
    {
        auto mfd = grow(seed, targetTets);

        auto mfd1 = mfd;
        writef("  baseline (exec/undo):  ");
        benchmarkBaseline(mfd1, targetTets);

        auto mfd2 = mfd;
        writef("  speculative (full obj):");
        benchmarkSpeculative!false(mfd2, targetTets);

        auto mfd3 = mfd;
        writef("  spec+bitmask proposal: ");
        benchmarkSpeculative!true(mfd3, targetTets);

        writeln();
    }
}

Manifold!dim grow(ref Manifold!dim seed, int targetTets)
{
    auto mfd = seed;
    int nextVertex = 6;

    writef("Growing to %d tets...", targetTets);
    stdout.flush();

    while (mfd.fVector[dim] < targetTets)
    {
        auto facets_ = mfd.asSimplicialComplex.facets(dim);
        auto nFacets = cast(int) facets_.walkLength;
        auto facet = facets_[uniform(0, nFacets)].array;
        auto newV = nextVertex;
        auto move = BistellarMove!dim(facet, [newV]);
        if (!mfd.contains(move.coCenter))
        {
            mfd.doMove(move);
            nextVertex++;
        }
    }

    writefln(" done (%d tets, %d edges, %d vertices)",
        mfd.fVector[dim], mfd.fVector[dim-2], mfd.fVector[0]);
    stdout.flush();
    return mfd;
}

// --- Proposal (original: materializes all subsets) ---
struct ProposedMove
{
    BistellarMove!(dim, int) move;
}

ProposedMove proposeMove(ref Manifold!dim mfd, int newVertex)
{
    alias BM = BistellarMove!(dim, int);

    while (true)
    {
        auto facet = mfd.asSimplicialComplex.randomFacetOfDim(dim);
        auto center = facet.subsets.array.rndChoice.array;

        auto centerLen = cast(int) center.walkLength;
        auto centerDeg = mfd.degree(center);

        if (centerDeg + centerLen - 1 != dim + 1)
            continue;

        BM bm;
        if (centerLen - 1 == dim)
            bm = BM(center, [newVertex]);
        else
        {
            auto coCenter_ = mfd.coCenter(center, facet);
            bm = BM(center, coCenter_);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        alias HM = HingeMove!(dim, int);
        SumType!(BM, HM) wrapped = bm;
        if (!mfd.hasValidMove(wrapped))
            continue;

        return ProposedMove(bm);
    }
}

// --- Optimized proposal: random bitmask instead of materializing subsets ---
ProposedMove proposeMoveOptimized(ref Manifold!dim mfd, int newVertex)
{
    alias BM = BistellarMove!(dim, int);
    enum nVerts = dim + 1;
    enum maxMask = (1 << nVerts) - 1;

    while (true)
    {
        auto facet = mfd.asSimplicialComplex.randomFacetOfDim(dim);

        auto mask = uniform(1, maxMask + 1);
        int[nVerts] centerBuf;
        int centerLen = 0;
        foreach (i; 0 .. nVerts)
        {
            if (mask & (1 << i))
                centerBuf[centerLen++] = facet[i];
        }
        auto center = centerBuf[0 .. centerLen];
        center.sort();

        auto centerDeg = mfd.degree(center);
        if (centerDeg + centerLen - 1 != dim + 1)
            continue;

        BM bm;
        if (centerLen - 1 == dim)
            bm = BM(center, [newVertex]);
        else
        {
            auto coCenter_ = mfd.coCenter(center, facet);
            bm = BM(center, coCenter_);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        alias HM = HingeMove!(dim, int);
        SumType!(BM, HM) wrapped = bm;
        if (!mfd.hasValidMove(wrapped))
            continue;

        return ProposedMove(bm);
    }
}

// --- BASELINE: execute, compute full objective, undo if rejected ---
void benchmarkBaseline(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);
    alias HM = HingeMove!(dim, int);

    int totalAccepted = 0;
    int totalTried = 0;
    int targetAccepted = min(5000, targetTets * 2);

    auto params = BenchParams(targetTets);

    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!int;

        auto proposed = proposeMove(mfd, unusedVertices.back);
        auto bm = proposed.move;
        SumType!(BM, HM) chosenMove = bm;

        mfd.doMove(chosenMove);

        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;

        totalTried++;

        real newObjective = computeObjective(mfd, params);
        real deltaObj = newObjective - currentObjective;

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            mfd.undoMove(chosenMove);
            if (bm.coCenter.length == 1) unusedVertices ~= bm.coCenter;
            if (bm.center.length == 1)
            {
                assert(bm.center.front == unusedVertices.back);
                unusedVertices.popBack;
            }
        }
        else
        {
            totalAccepted++;
            currentObjective = newObjective;
        }
    }

    timer.stop();
    printResult(timer.peek(), totalAccepted, totalTried, targetTets);
}

// --- SPECULATIVE: compute ΔE via speculative delta (full 4-term objective) ---
void benchmarkSpeculative(bool optimizedProposal)(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);
    alias HM = HingeMove!(dim, int);

    int totalAccepted = 0;
    int totalTried = 0;
    int targetAccepted = min(5000, targetTets * 2);

    auto params = BenchParams(targetTets);

    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!int;

        static if (optimizedProposal)
            auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
        else
            auto proposed = proposeMove(mfd, unusedVertices.back);
        auto bm = proposed.move;

        totalTried++;

        real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            // Rejected — nothing to do
        }
        else
        {
            // Accepted — execute the move
            mfd.doMove(bm);

            if (bm.coCenter.length == 1) unusedVertices.popBack;
            if (bm.center.length == 1) unusedVertices ~= bm.center;

            currentObjective += deltaObj;
            totalAccepted++;
        }
    }

    timer.stop();
    printResult(timer.peek(), totalAccepted, totalTried, targetTets);
}

void printResult(T)(T elapsed, int totalAccepted, int totalTried, int targetTets)
{
    auto usPerAccepted = elapsed.total!"usecs" / totalAccepted;
    auto usPerTried = elapsed.total!"usecs" / totalTried;
    auto acceptRate = 100.0 * totalAccepted / totalTried;

    writefln(" %5d tets: %6d acc / %6d try (%.1f%%)  %s  %d μs/acc  %d μs/try",
        targetTets, totalAccepted, totalTried, acceptRate,
        elapsed, usPerAccepted, usPerTried);
    stdout.flush();
}

}
