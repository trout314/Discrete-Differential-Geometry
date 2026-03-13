/// Benchmark for the MCMC sampler hot path.
/// Compares baseline (execute-compute-undo) vs speculative delta with full objective.
module bench_sampler;

version (unittest) {} else {

import std.algorithm, std.array, std.conv, std.datetime.stopwatch, std.math,
    std.range, std.stdio, std.string;
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

    int[] benchSizes = [5000];

    // Grow once with equilibration, benchmarking at each target size
    auto mfd = growWithBenchmarks(seed, benchSizes);
}

Manifold!dim growWithBenchmarks(ref Manifold!dim seed, int[] benchSizes)
{
    auto mfd = seed;
    int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices;

    int finalTarget = benchSizes[$ - 1];

    // Build milestone list: double from 50 up to finalTarget, including bench sizes
    int[] milestones;
    {
        int sz = 50;
        while (sz < finalTarget)
        {
            milestones ~= sz;
            sz = cast(int)(sz * 2);
        }
        // Merge in bench sizes
        foreach (bs; benchSizes)
            if (!milestones.canFind(bs))
                milestones ~= bs;
        milestones.sort();
    }

    // Which bench sizes have we reported?
    int benchIdx = 0;

    foreach (milestone; milestones)
    {
        if (mfd.fVector[dim] >= milestone)
            continue;

        // Growth phase: add vertices via stellar subdivision until we hit milestone
        writef("  growing %d -> %d...", mfd.fVector[dim], milestone);
        stdout.flush();
        while (mfd.fVector[dim] < milestone)
        {
            auto facet = mfd.randomFacetOfDim(dim).array;
            auto newV = nextVertex;
            auto move = BistellarMove!dim(facet, [newV]);
            if (!mfd.contains(move.coCenter))
            {
                mfd.doMove(move);
                nextVertex++;
            }
        }
        writefln(" %d tets", mfd.fVector[dim]);

        // Equilibration phase
        equilibrate(mfd, milestone, unusedVertices, nextVertex);

        // Update nextVertex to be above the max vertex label in the manifold
        nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;

        // Benchmark if this milestone matches (or exceeds) a bench size
        while (benchIdx < benchSizes.length && milestone >= benchSizes[benchIdx])
        {
            int actualTets = cast(int) mfd.fVector[dim];
            writeln();
            writefln("=== Benchmark at ~%d tets (%d actual) ===", benchSizes[benchIdx], actualTets);
            stdout.flush();
            printDegreeStats(mfd);

            auto mfd1 = mfd;
            writef("  speculative (profiled):");
            benchmarkSpeculativeProfiled(mfd1, actualTets);

            version (TrackValidMoves)
            {
                auto mfd2 = mfd;
                writef("  exact incr (V tracked):");
                benchmarkExactHastingsIncremental(mfd2, actualTets);
            }

            writeln();
            stdout.flush();
            benchIdx++;
        }
    }

    return mfd;
}

void equilibrate(ref Manifold!dim mfd, int targetSize,
    ref int[] unusedVertices, ref int nextVertex)
{
    auto eqParams = BenchParams(targetSize);
    auto currentObjective = computeObjective(mfd, eqParams);
    unusedVertices.length = 0;
    unusedVertices ~= nextVertex;

    int currentSize = cast(int) mfd.fVector[dim];
    int eqTarget = currentSize * 10;
    int eqAccepted = 0;

    // Track degree variance over windows to check convergence
    int windowSize = max(1000, currentSize);
    int windowAccepted = 0;
    real[] windowTotSqDegs;

    StopWatch eqTimer;
    eqTimer.start();

    int eqTrials = 0;
    while (eqAccepted < eqTarget)
    {
        if (unusedVertices.empty)
        {
            nextVertex++;
            unusedVertices ~= nextVertex;
        }

        auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
        auto bm = proposed.move;
        eqTrials++;

        // Periodic progress
        if (eqTrials % 100000 == 0)
        {
            writef("\r    eq %d/%d accepted (%d trials, %.1f%%), %d tets   ",
                eqAccepted, eqTarget, eqTrials,
                100.0 * eqAccepted / eqTrials, mfd.fVector[dim]);
            stdout.flush();
        }

        real deltaObj = speculativeDelta(mfd, bm, currentObjective, eqParams);

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            // Rejected
        }
        else
        {
            mfd.doMove(bm);
            if (bm.coCenter.length == 1) unusedVertices.popBack;
            if (bm.center.length == 1) unusedVertices ~= bm.center;
            currentObjective += deltaObj;
            eqAccepted++;
            windowAccepted++;

            // Update nextVertex if a new vertex was consumed
            if (bm.coCenter.length == 1 && bm.coCenter[0] >= nextVertex)
                nextVertex = bm.coCenter[0] + 1;

            // Record degree variance at window boundaries
            if (windowAccepted >= windowSize)
            {
                windowTotSqDegs ~= cast(real) mfd.totalSquareDegree(0);
                windowAccepted = 0;
            }
        }
    }

    eqTimer.stop();

    // Report equilibration convergence
    real finalTotSqDeg = cast(real) mfd.totalSquareDegree(0);
    real avgDeg = cast(real)(dim + 1) * mfd.fVector[dim] / mfd.fVector[0];
    writef("\r    eq %d acc/%d try (%.1f%%, %.1fs), %d tets, deg_var=%.1f",
        eqAccepted, eqTrials, 100.0 * eqAccepted / eqTrials,
        eqTimer.peek.total!"msecs" / 1000.0,
        mfd.fVector[dim],
        finalTotSqDeg / mfd.fVector[0] - avgDeg * avgDeg);

    // Check convergence: compare first and last half of windows
    if (windowTotSqDegs.length >= 4)
    {
        auto nw = windowTotSqDegs.length;
        real firstHalf = 0, secondHalf = 0;
        auto half = nw / 2;
        foreach (i; 0 .. half)
            firstHalf += windowTotSqDegs[i];
        foreach (i; half .. nw)
            secondHalf += windowTotSqDegs[i];
        firstHalf /= half;
        secondHalf /= (nw - half);
        real pctChange = 100.0 * (secondHalf - firstHalf) / firstHalf;
        writefln(", drift=%.2f%%", pctChange);
    }
    else
        writeln();

    stdout.flush();
}

void printDegreeStats(ref const Manifold!dim mfd)
{
    // Vertex degree = number of top-dimensional facets containing the vertex
    auto verts = mfd.simplices(0);
    int nVerts = cast(int) mfd.fVector[0];
    int maxDeg = 0, minDeg = int.max;
    long sumDeg = 0;
    long sumSqDeg = 0;
    int[int] degHist;
    foreach (v; verts)
    {
        int deg = cast(int) mfd.degree(v);
        sumDeg += deg;
        sumSqDeg += deg * deg;
        if (deg > maxDeg) maxDeg = deg;
        if (deg < minDeg) minDeg = deg;
        degHist[deg]++;
    }
    real avgDeg = cast(real) sumDeg / nVerts;
    real variance = cast(real) sumSqDeg / nVerts - avgDeg * avgDeg;

    writefln("  Vertex degree: avg=%.1f, min=%d, max=%d, stddev=%.1f",
        avgDeg, minDeg, maxDeg, sqrt(variance));

    // Show histogram of most common degrees
    int[][] sorted;
    foreach (k, v; degHist)
        sorted ~= [k, v];
    sorted.sort!((a, b) => a[0] < b[0]);
    writef("  Degree histogram:");
    foreach (kv; sorted)
        writef(" %d:%d", kv[0], kv[1]);
    writeln();
    stdout.flush();
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
        auto facet = mfd.randomFacetOfDim(dim);
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

        if (!mfd.hasValidMove(bm))
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
        auto facet = mfd.randomFacetOfDim(dim);

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

        if (!mfd.hasValidMove(bm))
            continue;

        return ProposedMove(bm);
    }
}

// --- BASELINE: execute, compute full objective, undo if rejected ---
void benchmarkBaseline(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);

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
        mfd.doMove(bm);

        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;

        totalTried++;

        real newObjective = computeObjective(mfd, params);
        real deltaObj = newObjective - currentObjective;

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            mfd.undoMove(bm);
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

// --- PROFILED SPECULATIVE: breaks down time by operation ---
void benchmarkSpeculativeProfiled(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);

    int totalAccepted = 0;
    int totalTried = 0;
    int targetAccepted = min(5000, targetTets * 2);

    auto params = BenchParams(targetTets);

    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);

    Duration tPropose, tDelta, tDoMove;
    long nProposeCalls = 0;

    StopWatch timer, sw;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!int;

        sw.reset(); sw.start();
        auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
        sw.stop(); tPropose += sw.peek();
        nProposeCalls++;

        auto bm = proposed.move;
        totalTried++;

        sw.reset(); sw.start();
        real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);
        sw.stop(); tDelta += sw.peek();

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            // Rejected
        }
        else
        {
            sw.reset(); sw.start();
            mfd.doMove(bm);
            sw.stop(); tDoMove += sw.peek();

            if (bm.coCenter.length == 1) unusedVertices.popBack;
            if (bm.center.length == 1) unusedVertices ~= bm.center;

            currentObjective += deltaObj;
            totalAccepted++;
        }
    }

    timer.stop();
    printResult(timer.peek(), totalAccepted, totalTried, targetTets);

    auto total = timer.peek().total!"usecs";
    writefln("    propose: %7d μs (%4.1f%%)  [%d calls, %d μs/call]",
        tPropose.total!"usecs", 100.0 * tPropose.total!"usecs" / total,
        nProposeCalls, tPropose.total!"usecs" / nProposeCalls);
    writefln("    delta:   %7d μs (%4.1f%%)", tDelta.total!"usecs", 100.0 * tDelta.total!"usecs" / total);
    writefln("    doMove:  %7d μs (%4.1f%%)  [%d calls, %d μs/call]",
        tDoMove.total!"usecs", 100.0 * tDoMove.total!"usecs" / total,
        totalAccepted, tDoMove.total!"usecs" / totalAccepted);
    writefln("    other:   %7d μs (%4.1f%%)",
        total - tPropose.total!"usecs" - tDelta.total!"usecs" - tDoMove.total!"usecs",
        100.0 * (total - tPropose.total!"usecs" - tDelta.total!"usecs" - tDoMove.total!"usecs") / total);

    // Report GC stats
    import core.memory : GC;
    auto stats = GC.stats;
    auto profStats = GC.profileStats;
    writefln("    GC: used=%d KB, collections=%d, total_pause=%d μs",
        stats.usedSize / 1024, profStats.numCollections,
        profStats.totalPauseTime.total!"usecs");
    stdout.flush();
}

// --- EXACT HASTINGS: execute, compute objective + valid move count, undo if rejected ---
void benchmarkExactHastings(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);

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

        auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
        auto bm = proposed.move;
        immutable vBefore = cast(real) mfd.countValidBistellarMoves;

        mfd.doMove(bm);

        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;

        totalTried++;

        real newObjective = computeObjective(mfd, params);
        real deltaObj = newObjective - currentObjective;
        immutable vAfter = cast(real) mfd.countValidBistellarMoves;

        real logAlpha = -deltaObj + log(vBefore) - log(vAfter);

        if ((logAlpha < 0) && (uniform01 > exp(logAlpha)))
        {
            mfd.undoMove(bm);
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

// --- EXACT HASTINGS (INCREMENTAL): uses tracked validMoveCount ---
version (TrackValidMoves)
{
void benchmarkExactHastingsIncremental(ref Manifold!dim mfd, int targetTets)
{
    alias BM = BistellarMove!(dim, int);

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

        auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
        auto bm = proposed.move;
        immutable vBefore = cast(real) mfd.validMoveCount;

        mfd.doMove(bm);

        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;

        totalTried++;

        real newObjective = computeObjective(mfd, params);
        real deltaObj = newObjective - currentObjective;
        immutable vAfter = cast(real) mfd.validMoveCount;

        real logAlpha = -deltaObj + log(vBefore) - log(vAfter);

        if ((logAlpha < 0) && (uniform01 > exp(logAlpha)))
        {
            mfd.undoMove(bm);
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
}

// --- EXACT HASTINGS + REJECTION-FREE PROPOSAL: O(1) uniform sampling from valid move array ---
version (TrackValidMoves)
{
void benchmarkRejectionFree(ref Manifold!dim mfd, int targetTets)
{
    import std.random : Mt19937;
    alias BM = BistellarMove!(dim, int);

    int totalAccepted = 0;
    int totalTried = 0;
    int targetAccepted = min(5000, targetTets * 2);

    auto params = BenchParams(targetTets);

    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);
    auto rng = Mt19937(42);

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!int;

        immutable vBefore = cast(real) mfd.validMoveCount;

        // O(1) rejection-free proposal
        auto bm = mfd.sampleValidMove(rng, unusedVertices.back);
        mfd.doMove(bm);

        if (bm.coCenter.length == 1) unusedVertices.popBack;
        if (bm.center.length == 1) unusedVertices ~= bm.center;

        totalTried++;

        real newObjective = computeObjective(mfd, params);
        real deltaObj = newObjective - currentObjective;
        immutable vAfter = cast(real) mfd.validMoveCount;

        real logAlpha = -deltaObj + log(vBefore) - log(vAfter);

        if ((logAlpha < 0) && (uniform01 > exp(logAlpha)))
        {
            mfd.undoMove(bm);
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
