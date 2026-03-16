/// Benchmark for the MCMC sampler hot path.
/// Compares baseline (execute-compute-undo) vs speculative delta with full objective.
module bench_sampler;

version (unittest) {} else {

import std.algorithm, std.array, std.conv, std.datetime.stopwatch, std.math,
    std.range, std.stdio, std.string;
import std.file : exists;
import std.random : uniform, uniform01, rndChoice = choice;
import manifold, manifold_examples, manifold_moves, sampler, simplicial_complex, utility;

enum dim = 3;

struct BenchParams
{
    int numFacetsTarget;
    real hingeDegreeTarget = 4.5;
    real numFacetsCoef = 0.1;
    real numHingesCoef = 0.05;
    real hingeDegreeVarianceCoef = 0.0;
    real coDim3DegreeVarianceCoef = 0.1;
}

// Use shared sampler functions via UFCS (objective, penalties, etc.)
// Aliases for backward compatibility with existing call sites
alias computeObjective = sampler.objective;
alias speculativeDelta = sampler.speculativeBistellarDelta;

// --- Main ---

void main()
{
    auto seed = loadManifold!dim("standard_triangulations/dim_3_sphere.mfd");

    writeln("=== Entropic force measurement ===");
    stdout.flush();

    measureEntropicForce(seed);

    writeln("\ndone");
    stdout.flush();
}

string equilibratedPath(int targetSize)
{
    return format("standard_triangulations/equilibrated_%d.mfd", targetSize);
}

Manifold!dim loadOrBuildEquilibrated(ref Manifold!dim seed, int targetSize)
{
    auto path = equilibratedPath(targetSize);
    if (exists(path))
    {
        writefln("Loading equilibrated manifold from %s", path);
        stdout.flush();
        return loadManifold!dim(path);
    }

    writefln("Building and equilibrating manifold at %d tets (ramped growth)...", targetSize);
    stdout.flush();

    auto mfd = seed;
    int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextVertex];

    // Ramped growth: increase target in steps, equilibrating at each step
    // Scale step size to keep ~100 ramp steps max
    int stepSize = max(500, targetSize / 100);
    enum eqSweepsPerStep = 5;

    for (int currentTarget = stepSize; currentTarget <= targetSize; currentTarget += stepSize)
    {
        // Grow via stellar subdivision to currentTarget
        while (mfd.fVector[dim] < currentTarget)
        {
            auto facet = mfd.randomFacetOfDim(dim).array;
            auto move = BistellarMove!dim(facet, [nextVertex]);
            if (!mfd.contains(move.coCenter))
            {
                mfd.doMove(move);
                nextVertex++;
            }
        }

        // Equilibrate at this step size using shared mcmcStep
        {
            auto params = BenchParams(currentTarget);
            auto currentObjective = computeObjective(mfd, params);
            int sweepSize = max(currentTarget, 1000);
            ulong hT, hA;
            ulong[dim + 1] bT, bA;

            foreach (_; 0 .. eqSweepsPerStep)
            {
                int accepted = 0;
                while (accepted < sweepSize)
                    if (mfd.mcmcStep(currentObjective, unusedVertices, params, 0.0,
                            hT, hA, bT, bA))
                        accepted++;
            }
        }

        if (currentTarget % max(5000, stepSize * 10) == 0 || currentTarget == targetSize)
        {
            writefln("  ramped to %d tets (%d actual), deg_var=%.1f",
                currentTarget, mfd.fVector[dim], mfd.degreeVariance(0));
            stdout.flush();
        }
    }

    // Final equilibration at full target size
    equilibrate(mfd, targetSize, unusedVertices, nextVertex);

    // Save for future runs
    writefln("Saving equilibrated manifold to %s", path);
    stdout.flush();
    mfd.saveTo(path);

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

    ulong hT, hA;
    ulong[dim + 1] bT, bA;

    int eqTrials = 0;
    while (eqAccepted < eqTarget)
    {
        eqTrials++;

        // Periodic progress
        if (eqTrials % 100000 == 0)
        {
            writef("\r    eq %d/%d accepted (%d trials, %.1f%%), %d tets   ",
                eqAccepted, eqTarget, eqTrials,
                100.0 * eqAccepted / eqTrials, mfd.fVector[dim]);
            stdout.flush();
        }

        if (mfd.mcmcStep(currentObjective, unusedVertices, eqParams, 0.0,
                hT, hA, bT, bA))
        {
            eqAccepted++;
            windowAccepted++;

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

        auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
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
void benchmarkSpeculative(ref Manifold!dim mfd, int targetTets)
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

        auto bm = mfd.chooseRandomMove(unusedVertices.back, params);

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
    int[dim + 1] moveTypeCounts; // index by center length - 1

    import core.memory : GC;
    auto gcBefore = GC.profileStats;

    StopWatch timer, sw;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!int;

        sw.reset(); sw.start();
        auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
        sw.stop(); tPropose += sw.peek();
        nProposeCalls++;
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
            moveTypeCounts[bm.center.length - 1]++;

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

    // Report GC stats (delta from benchmark start)
    auto gcAfter = GC.profileStats;
    auto gcCollections = gcAfter.numCollections - gcBefore.numCollections;
    auto gcPause = (gcAfter.totalPauseTime - gcBefore.totalPauseTime).total!"usecs";
    writefln("    GC: collections=%d, pause=%d μs (%.1f%%)",
        gcCollections, gcPause, 100.0 * gcPause / total);
    writef("    Move types:");
    static foreach (i; 0 .. dim + 1)
        writef(" %d→%d:%d", i + 1, dim + 2 - (i + 1), moveTypeCounts[i]);
    writeln();
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

        auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
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

        auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
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

// --- Mixing comparison: bistellar-only vs bistellar+hinge ---

struct MixingResult
{
    int targetTets;
    real hingeMoveProb;
    int totalAccepted;
    int totalTried;
    ulong hingeTries;
    ulong hingeAccepts;
    real finalDegVariance;
    real[] degVarianceTrajectory;
    long elapsedUsecs;
}

MixingResult runMixingSampler(ref Manifold!dim mfd, int targetTets,
    int targetAccepted, real hingeMoveProb)
{
    auto params = BenchParams(targetTets);
    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);

    ulong hingeTries, hingeAccepts;
    ulong[dim + 1] bTries, bAccepts;

    int totalAccepted = 0;
    int totalTried = 0;

    // Record degree variance trajectory at regular intervals
    int windowSize = max(100, targetTets / 10);
    real[] trajectory;

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        totalTried++;
        if (mfd.mcmcStep(currentObjective, unusedVertices, params, hingeMoveProb,
                hingeTries, hingeAccepts, bTries, bAccepts))
        {
            totalAccepted++;

            // Update nextVertex tracking
            auto verts = unusedVertices;
            if (verts.length > 0 && verts[$ - 1] >= nextUnused)
                nextUnused = verts[$ - 1] + 1;

            // Record trajectory point
            if (totalAccepted % windowSize == 0)
            {
                auto totSqDeg = cast(real) mfd.totalSquareDegree(0);
                auto nVerts = cast(real) mfd.fVector[0];
                auto avgDeg = cast(real)(dim + 1) * mfd.fVector[dim] / nVerts;
                trajectory ~= totSqDeg / nVerts - avgDeg * avgDeg;
            }
        }
    }

    timer.stop();

    auto totSqDeg = cast(real) mfd.totalSquareDegree(0);
    auto nVerts = cast(real) mfd.fVector[0];
    auto avgDeg = cast(real)(dim + 1) * mfd.fVector[dim] / nVerts;

    MixingResult result;
    result.targetTets = targetTets;
    result.hingeMoveProb = hingeMoveProb;
    result.totalAccepted = totalAccepted;
    result.totalTried = totalTried;
    result.hingeTries = hingeTries;
    result.hingeAccepts = hingeAccepts;
    result.finalDegVariance = totSqDeg / nVerts - avgDeg * avgDeg;
    result.degVarianceTrajectory = trajectory;
    result.elapsedUsecs = timer.peek.total!"usecs";
    return result;
}

/// Compute lag-1 autocorrelation of a trajectory.
real lag1Autocorrelation(real[] xs)
{
    if (xs.length < 3) return real.nan;
    real mean = 0;
    foreach (x; xs) mean += x;
    mean /= xs.length;

    real c0 = 0, c1 = 0;
    foreach (i; 0 .. xs.length)
        c0 += (xs[i] - mean) ^^ 2;
    foreach (i; 0 .. xs.length - 1)
        c1 += (xs[i] - mean) * (xs[i + 1] - mean);

    if (c0 == 0) return 0;
    return c1 / c0;
}

/// Grow a manifold by 10% via stellar subdivisions (1-4 moves), creating a
/// non-equilibrated state. Returns the grown manifold and the next unused vertex.
Manifold!dim growBy10Percent(ref Manifold!dim mfd, ref int nextVertex)
{
    auto result = mfd;
    int targetTets = cast(int)(mfd.fVector[dim] * 1.1);

    while (result.fVector[dim] < targetTets)
    {
        auto facet = result.randomFacetOfDim(dim).array;
        auto move = BistellarMove!dim(facet, [nextVertex]);
        if (!result.contains(move.coCenter))
        {
            result.doMove(move);
            nextVertex++;
        }
    }
    return result;
}

/// Compute vertex degree variance.
real degVariance(ref const Manifold!dim mfd)
{
    auto totSqDeg = cast(real) mfd.totalSquareDegree(0);
    auto nVerts = cast(real) mfd.fVector[0];
    auto avgDeg = cast(real)(dim + 1) * mfd.fVector[dim] / nVerts;
    return totSqDeg / nVerts - avgDeg * avgDeg;
}

struct EquilibrationResult
{
    real hingeMoveProb;
    int totalAccepted;
    int totalTried;
    ulong hingeTries;
    ulong hingeAccepts;
    real startDegVariance;
    real finalDegVariance;
    real[] degVarianceTrajectory;
    long elapsedUsecs;
}

/// Run equilibration: MCMC until degree variance stabilizes or we hit a sweep limit.
/// Records degree variance trajectory for analysis.
EquilibrationResult runEquilibration(ref Manifold!dim mfd, int targetTets,
    int eqSweeps, real hingeMoveProb, BenchParams paramsOverride = BenchParams.init)
{
    auto params = (paramsOverride.numFacetsTarget > 0)
        ? paramsOverride : BenchParams(targetTets);
    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);

    ulong hingeTries, hingeAccepts;
    ulong[dim + 1] bTries, bAccepts;

    int totalAccepted = 0;
    int totalTried = 0;
    int targetAccepted = eqSweeps * targetTets;

    // Record trajectory at regular intervals (every 0.1 sweeps)
    int windowSize = max(100, targetTets / 10);
    real[] trajectory;
    real startDV = degVariance(mfd);

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        totalTried++;
        if (mfd.mcmcStep(currentObjective, unusedVertices, params, hingeMoveProb,
                hingeTries, hingeAccepts, bTries, bAccepts))
        {
            totalAccepted++;

            if (totalAccepted % windowSize == 0)
                trajectory ~= degVariance(mfd);
        }
    }

    timer.stop();

    EquilibrationResult result;
    result.hingeMoveProb = hingeMoveProb;
    result.totalAccepted = totalAccepted;
    result.totalTried = totalTried;
    result.hingeTries = hingeTries;
    result.hingeAccepts = hingeAccepts;
    result.startDegVariance = startDV;
    result.finalDegVariance = degVariance(mfd);
    result.degVarianceTrajectory = trajectory;
    result.elapsedUsecs = timer.peek.total!"usecs";
    return result;
}

/// Take a seed equilibrated with default BenchParams (hingeDegreeVarianceCoef=0),
/// switch to a new energy function with hingeDegreeVarianceCoef turned on, and
/// measure how fast the manifold relaxes to the new equilibrium.
void benchmarkRelaxation(ref Manifold!dim seedMfd, int seedTets)
{
    enum relaxSweeps = 20;

    // The seed was equilibrated with default BenchParams (coDim3DegreeVarianceCoef=0.1).
    // New target: 10x stronger vertex degree variance penalty.
    auto targetParams = BenchParams(seedTets);
    targetParams.coDim3DegreeVarianceCoef = 1.0;

    writefln("\n=== Relaxation benchmark at %d tets ===", seedTets);
    writefln("  Seed equilibrated with coDim3DegreeVarianceCoef=0.1");
    writefln("  Target energy: coDim3DegreeVarianceCoef=%.1f", targetParams.coDim3DegreeVarianceCoef);
    writefln("  Relaxation: %d sweeps", relaxSweeps);

    writefln("  Starting vertex deg variance: %.2f", degVariance(seedMfd));
    stdout.flush();

    // Run bistellar-only relaxation
    auto mfd1 = seedMfd;
    auto r1 = runEquilibration(mfd1, seedTets, relaxSweeps, 0.0, targetParams);

    // Run with 30% hinge moves
    auto mfd2 = seedMfd;
    auto r2 = runEquilibration(mfd2, seedTets, relaxSweeps, 0.3, targetParams);

    // Report
    writefln("  Bistellar only:   %d acc / %d try (%.1f%%), %.1f ms",
        r1.totalAccepted, r1.totalTried,
        100.0 * r1.totalAccepted / r1.totalTried,
        r1.elapsedUsecs / 1000.0);
    writefln("    deg_var: %.2f -> %.2f", r1.startDegVariance, r1.finalDegVariance);

    writefln("  With 30%% hinge:  %d acc / %d try (%.1f%%), %.1f ms",
        r2.totalAccepted, r2.totalTried,
        100.0 * r2.totalAccepted / r2.totalTried,
        r2.elapsedUsecs / 1000.0);
    writefln("    deg_var: %.2f -> %.2f", r2.startDegVariance, r2.finalDegVariance);
    writefln("    Hinge: %d tried, %d accepted (%.1f%%)",
        r2.hingeTries, r2.hingeAccepts,
        r2.hingeTries > 0 ? 100.0 * r2.hingeAccepts / r2.hingeTries : 0.0);

    if (r1.elapsedUsecs > 0)
        writefln("    Speedup: %.2fx", cast(real) r1.elapsedUsecs / r2.elapsedUsecs);

    auto ac1 = lag1Autocorrelation(r1.degVarianceTrajectory);
    auto ac2 = lag1Autocorrelation(r2.degVarianceTrajectory);
    writefln("    Lag-1 autocorrelation: bistellar=%.4f, with-hinge=%.4f", ac1, ac2);

    // Print deg_var trajectory (first 20 points)
    void printTrajectory(string label, real[] trajectory)
    {
        if (trajectory.length == 0) return;
        writef("    %s trajectory:", label);
        foreach (i, dv; trajectory)
        {
            if (i >= 20) { writef(" ..."); break; }
            if (i % 5 == 0 && i > 0) writef("\n                          ");
            writef(" %.1f", dv);
        }
        writeln();
    }

    printTrajectory("Bistellar", r1.degVarianceTrajectory);
    printTrajectory("With-hinge", r2.degVarianceTrajectory);

    stdout.flush();
}

// ---------------------------------------------------------------------------
// Comprehensive relaxation benchmark
// ---------------------------------------------------------------------------

/// Flat-geometry target mean edge degree for 3-manifolds.
enum real flatEdgeDegree = 2.0 * PI / acos(1.0 / 3.0); // ≈ 5.1044

/// Seed parameters: physically motivated defaults for 3-manifolds.
BenchParams seedParams(int targetTets)
{
    BenchParams p;
    p.numFacetsTarget = targetTets;
    p.hingeDegreeTarget = flatEdgeDegree;
    p.numFacetsCoef = 0.1;
    p.numHingesCoef = 0.1;
    p.hingeDegreeVarianceCoef = 0.0;
    p.coDim3DegreeVarianceCoef = 0.2;
    return p;
}

/// Build an equilibrated seed with physically correct parameters and hinge moves.
Manifold!dim loadOrBuildRelaxSeed(ref Manifold!dim sphere, int targetSize)
{
    auto path = format("standard_triangulations/relax_seed_%d.mfd", targetSize);
    if (exists(path))
    {
        writefln("  Loading seed from %s", path);
        stdout.flush();
        return loadManifold!dim(path);
    }

    writefln("  Building seed at %d tets (ramped growth with hinge moves)...", targetSize);
    stdout.flush();

    auto mfd = sphere;
    int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextVertex];

    int stepSize = max(500, targetSize / 100);
    enum eqSweepsPerStep = 5;
    enum hingeMoveProb = 0.3;

    for (int currentTarget = stepSize; currentTarget <= targetSize; currentTarget += stepSize)
    {
        // Grow via stellar subdivision
        while (mfd.fVector[dim] < currentTarget)
        {
            auto facet = mfd.randomFacetOfDim(dim).array;
            auto move = BistellarMove!dim(facet, [nextVertex]);
            if (!mfd.contains(move.coCenter))
            {
                mfd.doMove(move);
                nextVertex++;
            }
        }

        // Equilibrate with hinge moves
        auto params = seedParams(currentTarget);
        auto currentObjective = computeObjective(mfd, params);
        int sweepSize = max(currentTarget, 1000);
        ulong hT, hA;
        ulong[dim + 1] bT, bA;

        foreach (_; 0 .. eqSweepsPerStep)
        {
            int accepted = 0;
            while (accepted < sweepSize)
            {
                if (unusedVertices.empty)
                    unusedVertices ~= (++nextVertex);
                if (mfd.mcmcStep(currentObjective, unusedVertices, params,
                        hingeMoveProb, hT, hA, bT, bA))
                    accepted++;
            }
        }

        if (currentTarget % max(5000, stepSize * 10) == 0 || currentTarget == targetSize)
        {
            writefln("    ramped to %d tets (%d actual), deg_var=%.1f, mean_edge_deg=%.3f",
                currentTarget, mfd.fVector[dim], degVariance(mfd), meanEdgeDegree(mfd));
            stdout.flush();
        }
    }

    // Final equilibration: 10 sweeps at target size
    {
        auto params = seedParams(targetSize);
        auto currentObjective = computeObjective(mfd, params);
        ulong hT, hA;
        ulong[dim + 1] bT, bA;
        int targetAcc = 10 * targetSize;
        int acc = 0;
        while (acc < targetAcc)
        {
            if (unusedVertices.empty)
                unusedVertices ~= (++nextVertex);
            if (mfd.mcmcStep(currentObjective, unusedVertices, params,
                    hingeMoveProb, hT, hA, bT, bA))
                acc++;
        }
    }

    writefln("  Saving seed to %s", path);
    stdout.flush();
    mfd.saveTo(path);
    return mfd;
}

real meanEdgeDegree(ref const Manifold!dim mfd)
{
    return 6.0 * mfd.fVector[dim] / cast(real) mfd.fVector[1];
}

/// A relaxation scenario: label + perturbed params.
struct RelaxScenario
{
    string label;
    BenchParams params;
}

/// Generate relaxation scenarios from a base seed params.
RelaxScenario[] makeScenarios(int targetTets)
{
    auto base = seedParams(targetTets);
    RelaxScenario[] scenarios;

    // Edge degree perturbations
    {
        auto p = base; p.hingeDegreeTarget = base.hingeDegreeTarget + 0.1;
        scenarios ~= RelaxScenario("edge_deg +0.1", p);
    }
    {
        auto p = base; p.hingeDegreeTarget = base.hingeDegreeTarget - 0.1;
        scenarios ~= RelaxScenario("edge_deg -0.1", p);
    }

    // Volume perturbations
    {
        auto p = base; p.numFacetsTarget = cast(int)(targetTets * 1.05);
        scenarios ~= RelaxScenario("volume +5%", p);
    }
    {
        auto p = base; p.numFacetsTarget = cast(int)(targetTets * 0.95);
        scenarios ~= RelaxScenario("volume -5%", p);
    }

    // Vertex degree variance beta perturbations
    {
        auto p = base; p.coDim3DegreeVarianceCoef = base.coDim3DegreeVarianceCoef * 2.0;
        scenarios ~= RelaxScenario("vtx_var_beta 2x", p);
    }
    {
        auto p = base; p.coDim3DegreeVarianceCoef = base.coDim3DegreeVarianceCoef * 0.5;
        scenarios ~= RelaxScenario("vtx_var_beta 0.5x", p);
    }

    return scenarios;
}

struct RelaxResult
{
    string label;
    real hingeMoveProb;
    int totalAccepted, totalTried;
    ulong hingeTries, hingeAccepts;
    real startObj, finalObj;
    real startDegVar, finalDegVar;
    real startMeanEdgeDeg, finalMeanEdgeDeg;
    int startVolume, finalVolume;
    real[] objTrajectory;
    real[] degVarTrajectory;
    long elapsedUsecs;
}

RelaxResult runRelaxScenario(ref Manifold!dim mfd, BenchParams params,
    int relaxSweeps, real hingeMoveProb)
{
    int targetTets = params.numFacetsTarget;
    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];

    auto currentObjective = computeObjective(mfd, params);
    ulong hingeTries, hingeAccepts;
    ulong[dim + 1] bTries, bAccepts;

    int totalAccepted = 0, totalTried = 0;
    int targetAccepted = relaxSweeps * targetTets;
    int windowSize = max(100, targetTets / 10);

    RelaxResult r;
    r.hingeMoveProb = hingeMoveProb;
    r.startObj = currentObjective;
    r.startDegVar = degVariance(mfd);
    r.startMeanEdgeDeg = meanEdgeDegree(mfd);
    r.startVolume = cast(int) mfd.fVector[dim];

    StopWatch timer;
    timer.start();

    while (totalAccepted < targetAccepted)
    {
        totalTried++;
        if (mfd.mcmcStep(currentObjective, unusedVertices, params, hingeMoveProb,
                hingeTries, hingeAccepts, bTries, bAccepts))
        {
            totalAccepted++;
            if (totalAccepted % windowSize == 0)
            {
                r.objTrajectory ~= currentObjective;
                r.degVarTrajectory ~= degVariance(mfd);
            }
        }
    }

    timer.stop();
    r.totalAccepted = totalAccepted;
    r.totalTried = totalTried;
    r.hingeTries = hingeTries;
    r.hingeAccepts = hingeAccepts;
    r.finalObj = currentObjective;
    r.finalDegVar = degVariance(mfd);
    r.finalMeanEdgeDeg = meanEdgeDegree(mfd);
    r.finalVolume = cast(int) mfd.fVector[dim];
    r.elapsedUsecs = timer.peek.total!"usecs";
    return r;
}

void benchmarkComprehensiveRelaxation(ref Manifold!dim sphere)
{
    int[] sizes = [500, 1000, 2000, 5000, 10_000];
    enum relaxSweeps = 20;

    // Summary table accumulator
    struct SummaryRow
    {
        int size;
        string label;
        real bistellarFinalObj, hingeFinalObj;
        real bistellarFinalDV, hingeFinalDV;
        real bistellarAC, hingeAC;
        long bistellarUsecs, hingeUsecs;
        real hingeAcceptRate;
    }
    SummaryRow[] summaryRows;

    foreach (targetSize; sizes)
    {
        auto seedMfd = loadOrBuildRelaxSeed(sphere, targetSize);
        int actualTets = cast(int) seedMfd.fVector[dim];

        writefln("\n--- Seed: %d tets, mean_edge_deg=%.3f, vtx_deg_var=%.1f ---",
            actualTets, meanEdgeDegree(seedMfd), degVariance(seedMfd));
        stdout.flush();

        auto scenarios = makeScenarios(actualTets);

        foreach (ref scenario; scenarios)
        {
            writefln("\n  Scenario: %s (target_tets=%d, edge_deg=%.3f, vtx_var_beta=%.2f)",
                scenario.label, scenario.params.numFacetsTarget,
                scenario.params.hingeDegreeTarget,
                scenario.params.coDim3DegreeVarianceCoef);
            stdout.flush();

            // Bistellar only
            auto mfd1 = seedMfd;
            auto r1 = runRelaxScenario(mfd1, scenario.params, relaxSweeps, 0.0);

            // With 30% hinge moves
            auto mfd2 = seedMfd;
            auto r2 = runRelaxScenario(mfd2, scenario.params, relaxSweeps, 0.3);

            // Report
            writefln("    Bistellar: obj %.1f->%.1f, dv %.1f->%.1f, "
                ~ "edge_deg %.3f->%.3f, vol %d->%d, %.0fms",
                r1.startObj, r1.finalObj, r1.startDegVar, r1.finalDegVar,
                r1.startMeanEdgeDeg, r1.finalMeanEdgeDeg,
                r1.startVolume, r1.finalVolume,
                r1.elapsedUsecs / 1000.0);
            writefln("    +Hinge:   obj %.1f->%.1f, dv %.1f->%.1f, "
                ~ "edge_deg %.3f->%.3f, vol %d->%d, %.0fms",
                r2.startObj, r2.finalObj, r2.startDegVar, r2.finalDegVar,
                r2.startMeanEdgeDeg, r2.finalMeanEdgeDeg,
                r2.startVolume, r2.finalVolume,
                r2.elapsedUsecs / 1000.0);

            if (r2.hingeTries > 0)
                writefln("    Hinge: %d/%d accepted (%.1f%%)",
                    r2.hingeAccepts, r2.hingeTries,
                    100.0 * r2.hingeAccepts / r2.hingeTries);

            auto ac1 = lag1Autocorrelation(r1.objTrajectory);
            auto ac2 = lag1Autocorrelation(r2.objTrajectory);
            writefln("    Obj lag-1 AC: bistellar=%.4f, +hinge=%.4f", ac1, ac2);

            summaryRows ~= SummaryRow(
                actualTets, scenario.label,
                r1.finalObj, r2.finalObj,
                r1.finalDegVar, r2.finalDegVar,
                ac1, ac2,
                r1.elapsedUsecs, r2.elapsedUsecs,
                r2.hingeTries > 0 ? 100.0 * r2.hingeAccepts / r2.hingeTries : 0);

            stdout.flush();
        }
    }

    // Summary table
    writeln("\n=== Summary ===");
    writefln("%-6s %-18s %10s %10s %10s %10s %8s %8s %6s",
        "Size", "Perturbation", "Bist.Obj", "Hinge.Obj",
        "Bist.DV", "Hinge.DV", "Bist.AC", "Hinge.AC", "H.Acc%");
    foreach (ref row; summaryRows)
    {
        writefln("%6d %-18s %10.1f %10.1f %10.1f %10.1f %8.4f %8.4f %5.1f%%",
            row.size, row.label,
            row.bistellarFinalObj, row.hingeFinalObj,
            row.bistellarFinalDV, row.hingeFinalDV,
            row.bistellarAC, row.hingeAC,
            row.hingeAcceptRate);
    }
    stdout.flush();
}

// ---------------------------------------------------------------------------
// Entropic force measurement
// ---------------------------------------------------------------------------

/// Volume-only energy function. Field names must match what mcmcStep/penalties
/// access via templates.
struct EntropicParams
{
    int numFacetsTarget;
    real hingeDegreeTarget = flatEdgeDegree;
    real numFacetsCoef;
    real numHingesCoef = 0.0;
    real hingeDegreeVarianceCoef = 0.0;
    real coDim3DegreeVarianceCoef = 0.0;
}

struct EntropicForceResult
{
    int nTarget;
    real betaV;
    real meanN, varN;
    real meanInvN;        // <1/N> for Hastings correction
    real entropicForce;   // 2*betaV*(meanN - nTarget) + meanInvN
    real forceErr;        // from batch means
    real finalDegVar;
    int eqSweeps;         // equilibration sweeps used
    bool eqConverged;
    real[] degVarEqTrajectory;   // during equilibration
    real[] degVarProdTrajectory; // during production
    int totalAccepted, totalTried;
    ulong hingeTries, hingeAccepts;
    ulong[dim + 1] bistellarTries, bistellarAccepts;
    long elapsedUsecs;
}

/// Compute standard error via batch means.
real batchMeansError(real[] series, int nBatches)
{
    if (series.length < nBatches) return real.nan;
    auto batchSize = series.length / nBatches;
    real[] batchMeans;
    foreach (b; 0 .. nBatches)
    {
        real sum = 0;
        foreach (i; b * batchSize .. (b + 1) * batchSize)
            sum += series[i];
        batchMeans ~= sum / batchSize;
    }
    // Standard error of the mean of batch means
    real mean = 0;
    foreach (m; batchMeans) mean += m;
    mean /= nBatches;
    real var = 0;
    foreach (m; batchMeans) var += (m - mean) ^^ 2;
    var /= (nBatches - 1);
    return sqrt(var / nBatches);
}

EntropicForceResult runEntropicMeasurement(
    ref Manifold!dim mfd,
    int nTarget,
    real betaV,
    int minEqSweeps,
    int maxEqSweeps,
    int prodSweeps,
    real hingeMoveProb)
{
    auto params = EntropicParams(nTarget, flatEdgeDegree, betaV);

    int nextUnused = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextUnused];
    auto currentObjective = sampler.objective(mfd, params);

    ulong hingeTries, hingeAccepts;
    ulong[dim + 1] bTries, bAccepts;
    int totalAccepted = 0, totalTried = 0;

    // --- Equilibration phase: dynamic convergence of deg_var ---
    // Run in windows of `windowSweeps` sweeps. After each window, record
    // the mean deg_var. Stop when the last `nCheck` window means are within
    // `tolerance` of each other (relative), with at least `minEqSweeps` total.
    enum windowSweeps = 10;
    enum nCheck = 6;       // require 6 consecutive stable windows (= 60 sweeps stable)
    enum tolerance = 0.03; // 3% relative

    int windowSize = windowSweeps * nTarget;  // accepted moves per window
    real[] windowMeans;
    real[] degVarEq;
    int dvSample = max(100, nTarget / 10);
    int eqSweepsDone = 0;
    bool converged = false;

    while (eqSweepsDone < maxEqSweeps)
    {
        // Run one window
        int windowAccepted = 0;
        real dvSum = 0;
        int dvCount = 0;
        while (windowAccepted < windowSize)
        {
            totalTried++;
            if (mfd.mcmcStep(currentObjective, unusedVertices, params, hingeMoveProb,
                    hingeTries, hingeAccepts, bTries, bAccepts))
            {
                totalAccepted++;
                windowAccepted++;
                if (windowAccepted % dvSample == 0)
                {
                    auto dv = degVariance(mfd);
                    degVarEq ~= dv;
                    dvSum += dv;
                    dvCount++;
                }
            }
        }
        eqSweepsDone += windowSweeps;
        if (dvCount > 0)
            windowMeans ~= dvSum / dvCount;

        // Check convergence: last nCheck windows within tolerance
        if (eqSweepsDone >= minEqSweeps && windowMeans.length >= nCheck)
        {
            auto recent = windowMeans[$ - nCheck .. $];
            real rMin = real.max, rMax = -real.max;
            foreach (m; recent) { if (m < rMin) rMin = m; if (m > rMax) rMax = m; }
            real mid = (rMin + rMax) / 2;
            if (mid > 0 && (rMax - rMin) / mid < tolerance)
            {
                converged = true;
                break;
            }
        }
    }

    if (converged)
        writefln("    Equilibrated in %d sweeps (deg_var=%.1f)", eqSweepsDone, degVariance(mfd));
    else
        writefln("    WARNING: not converged after %d sweeps (deg_var=%.1f)", maxEqSweeps, degVariance(mfd));

    // --- Production phase ---
    // Record (nFacets, importance weight) at regular intervals.
    // The importance weight w = 1/countValidMoves is expensive (O(n) recount),
    // so we sample it every `sampleInterval` accepted moves, not every move.
    int prodTarget = prodSweeps * nTarget;
    int sampleInterval = max(1, nTarget);  // ~1 sample per sweep
    int nSamples = prodTarget / sampleInterval;

    auto nFacetsSamples = new real[](nSamples);
    auto weightSamples = new real[](nSamples);
    real[] degVarProd;
    int prodAccepted = 0;
    int sampleIdx = 0;

    StopWatch timer;
    timer.start();

    while (prodAccepted < prodTarget)
    {
        totalTried++;
        if (mfd.mcmcStep(currentObjective, unusedVertices, params, hingeMoveProb,
                hingeTries, hingeAccepts, bTries, bAccepts))
        {
            prodAccepted++;
            totalAccepted++;

            if (prodAccepted % sampleInterval == 0 && sampleIdx < nSamples)
            {
                nFacetsSamples[sampleIdx] = cast(real) mfd.fVector[dim];
                weightSamples[sampleIdx] = 1.0 / cast(real) mfd.countValidMoves;
                sampleIdx++;
            }
            if (prodAccepted % dvSample == 0)
                degVarProd ~= degVariance(mfd);
        }
    }
    nSamples = sampleIdx;  // actual count (may be slightly less due to rounding)

    timer.stop();

    // --- Analysis ---
    // No Hastings correction in the sampler. The importance weight
    // w_i = 1/V(x_i) where V(x) = countValidMoves(x) corrects the chain's
    // biased stationary distribution back to the target exp(-E) * Omega.
    //
    // Importance-weighted mean:
    //   <N>_target = sum(N_i * w_i) / sum(w_i)
    //
    // Entropic force: F(<N>_target) = 2 * betaV * (<N>_target - N_target)

    real sumW = 0, sumNW = 0, sumN2W = 0;
    foreach (i; 0 .. nSamples)
    {
        real n = nFacetsSamples[i];
        real w = weightSamples[i];
        sumW += w;
        sumNW += n * w;
        sumN2W += n * n * w;
    }
    real meanN = sumNW / sumW;
    real varN = sumN2W / sumW - meanN * meanN;

    real force = 2.0 * betaV * (meanN - nTarget);

    // Batch means error: compute importance-weighted force per batch
    int nBatches = 20;
    auto batchSize = nSamples / nBatches;
    real[] batchForces;
    foreach (b; 0 .. nBatches)
    {
        real bSumW = 0, bSumNW = 0;
        foreach (i; b * batchSize .. (b + 1) * batchSize)
        {
            bSumW += weightSamples[i];
            bSumNW += nFacetsSamples[i] * weightSamples[i];
        }
        real bMeanN = bSumNW / bSumW;
        batchForces ~= 2.0 * betaV * (bMeanN - nTarget);
    }
    real forceErr = batchMeansError(batchForces, nBatches);

    EntropicForceResult r;
    r.nTarget = nTarget;
    r.betaV = betaV;
    r.meanN = meanN;
    r.varN = varN;
    r.meanInvN = sumW / nSamples;  // mean importance weight (1/V)
    r.entropicForce = force;
    r.forceErr = forceErr;
    r.finalDegVar = degVariance(mfd);
    r.eqSweeps = eqSweepsDone;
    r.eqConverged = converged;
    r.degVarEqTrajectory = degVarEq;
    r.degVarProdTrajectory = degVarProd;
    r.totalAccepted = totalAccepted;
    r.totalTried = totalTried;
    r.hingeTries = hingeTries;
    r.hingeAccepts = hingeAccepts;
    r.bistellarTries = bTries;
    r.bistellarAccepts = bAccepts;
    r.elapsedUsecs = timer.peek.total!"usecs";
    return r;
}

/// Build a seed equilibrated under a given energy function P.
/// Uses ramped growth with heavy equilibration at each step so the degree
/// distribution adapts naturally to the unconstrained ensemble.
Manifold!dim loadOrBuildSeed(P)(ref Manifold!dim sphere, int targetSize,
    P params, string pathPrefix, int eqSweepsPerStep = 20, real hingeMoveProb = 0.3)
{
    auto path = format("standard_triangulations/%s_%d.mfd", pathPrefix, targetSize);
    if (exists(path))
    {
        writefln("  Loading seed from %s", path);
        stdout.flush();
        return loadManifold!dim(path);
    }

    writefln("  Building %s seed at %d tets...", pathPrefix, targetSize);
    stdout.flush();

    auto mfd = sphere;
    int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
    int[] unusedVertices = [nextVertex];

    // Smaller step size for more gradual adaptation
    int stepSize = max(200, targetSize / 50);

    for (int currentTarget = stepSize; currentTarget <= targetSize; currentTarget += stepSize)
    {
        // Grow via stellar subdivision
        while (mfd.fVector[dim] < currentTarget)
        {
            auto facet = mfd.randomFacetOfDim(dim).array;
            auto move = BistellarMove!dim(facet, [nextVertex]);
            if (!mfd.contains(move.coCenter))
            {
                mfd.doMove(move);
                nextVertex++;
            }
        }

        // Equilibrate with the target energy function
        auto stepParams = params;
        stepParams.numFacetsTarget = currentTarget;
        auto currentObjective = sampler.objective(mfd, stepParams);
        int sweepSize = max(currentTarget, 500);
        ulong hT, hA;
        ulong[dim + 1] bT, bA;

        foreach (_; 0 .. eqSweepsPerStep)
        {
            int accepted = 0;
            while (accepted < sweepSize)
                if (mfd.mcmcStep(currentObjective, unusedVertices, stepParams,
                        hingeMoveProb, hT, hA, bT, bA))
                    accepted++;
        }

        if (currentTarget % max(1000, stepSize * 5) == 0 || currentTarget == targetSize)
        {
            writefln("    %d tets (%d actual), deg_var=%.1f",
                currentTarget, mfd.fVector[dim], degVariance(mfd));
            stdout.flush();
        }
    }

    // Final equilibration: 50 sweeps at target size
    {
        auto finalParams = params;
        finalParams.numFacetsTarget = targetSize;
        auto currentObjective = sampler.objective(mfd, finalParams);
        ulong hT, hA;
        ulong[dim + 1] bT, bA;
        int targetAcc = 50 * targetSize;
        int acc = 0;
        int tried = 0;
        while (acc < targetAcc)
        {
            tried++;
            if (mfd.mcmcStep(currentObjective, unusedVertices, finalParams,
                    hingeMoveProb, hT, hA, bT, bA))
                acc++;
            if (tried % 500000 == 0)
            {
                writef("\r    final eq: %d/%d accepted, deg_var=%.1f   ",
                    acc, targetAcc, degVariance(mfd));
                stdout.flush();
            }
        }
        writefln("\r    final eq done: %d tets, deg_var=%.1f                 ",
            mfd.fVector[dim], degVariance(mfd));
    }

    writefln("  Saving seed to %s", path);
    stdout.flush();
    mfd.saveTo(path);
    return mfd;
}

void measureEntropicForce(ref Manifold!dim sphere)
{
    int[] sizes = [500];
    real[] betaVs = [0.05, 0.1, 0.2];
    enum prodSweeps = 100;
    enum hingeMoveProb = 0.3;


    struct SummaryRow
    {
        int nTarget;
        real betaV;
        real meanN, varN, force, forceErr, degVar;
        int eqSweeps;
        bool eqConverged;
    }
    SummaryRow[] summary;

    foreach (targetSize; sizes)
    {
        // Build volume-only seed at the middle beta_V value
        auto volOnlyParams = EntropicParams(targetSize, flatEdgeDegree, 0.1);
        auto seedMfd = loadOrBuildSeed(sphere, targetSize, volOnlyParams, "vol_only");
        int actualTets = cast(int) seedMfd.fVector[dim];

        writefln("\n--- Target size: %d tets (seed has %d, deg_var=%.1f) ---",
            targetSize, actualTets, degVariance(seedMfd));
        stdout.flush();

        foreach (betaV; betaVs)
        {
            writefln("  beta_V = %.2f:", betaV);
            stdout.flush();

            auto mfd = seedMfd;
            auto r = runEntropicMeasurement(mfd, targetSize, betaV,
                60, 2000, prodSweeps, hingeMoveProb);

            writefln("    <N>_iw = %.2f (target %d, offset %.2f)", r.meanN, r.nTarget,
                r.meanN - r.nTarget);
            writefln("    Var(N)_iw = %.1f (1/(2*beta) = %.1f)", r.varN, 1.0 / (2 * betaV));
            writefln("    F_entropic = %.4f +/- %.4f", r.entropicForce, r.forceErr);
            writefln("    Eq: %d sweeps (%s), deg_var=%.1f, %.0f ms prod",
                r.eqSweeps, r.eqConverged ? "converged" : "NOT converged",
                r.finalDegVar, r.elapsedUsecs / 1000.0);

            // Move acceptance breakdown
            writefln("    Moves: %d tried, %d accepted (%.1f%%)",
                r.totalTried, r.totalAccepted,
                100.0 * r.totalAccepted / r.totalTried);
            writef("      Bistellar:");
            static foreach (i; 0 .. dim + 1)
                writef(" %d->%d: %d/%d", i + 1, dim + 2 - (i + 1),
                    r.bistellarAccepts[i], r.bistellarTries[i]);
            writeln();
            if (r.hingeTries > 0)
                writefln("      Hinge: %d/%d (%.1f%%)",
                    r.hingeAccepts, r.hingeTries,
                    100.0 * r.hingeAccepts / r.hingeTries);

            // Print deg_var trajectory diagnostics
            if (r.degVarEqTrajectory.length > 0)
            {
                writef("    Eq deg_var (last 10):");
                auto eqT = r.degVarEqTrajectory;
                auto start = eqT.length > 10 ? eqT.length - 10 : 0;
                foreach (i; start .. eqT.length) writef(" %.0f", eqT[i]);
                writeln();
            }
            if (r.degVarProdTrajectory.length > 0)
            {
                writef("    Prod deg_var (first 10):");
                foreach (i; 0 .. min(10, r.degVarProdTrajectory.length))
                    writef(" %.0f", r.degVarProdTrajectory[i]);
                writeln();
            }

            summary ~= SummaryRow(r.nTarget, r.betaV, r.meanN, r.varN,
                r.entropicForce, r.forceErr, r.finalDegVar,
                r.eqSweeps, r.eqConverged);
            stdout.flush();
        }
    }

    // Summary table
    writeln("\n=== Entropic force summary ===");
    writefln("%-8s %6s %10s %8s %10s %8s %10s %8s %6s %s",
        "N_target", "beta_V", "<N>", "offset", "Var(N)", "1/2beta", "F_entrop", "F_err", "EqSw", "Conv");
    foreach (ref row; summary)
    {
        writefln("%8d %6.2f %10.2f %8.2f %10.1f %8.1f %10.4f %8.4f %6d %s",
            row.nTarget, row.betaV, row.meanN, row.meanN - row.nTarget,
            row.varN, 1.0 / (2 * row.betaV), row.force, row.forceErr,
            row.eqSweeps, row.eqConverged ? "Y" : "N");
    }

    // Consistency check: for each size, F should be roughly independent of beta_V
    writeln("\n=== Consistency check: F vs beta_V ===");
    foreach (targetSize; sizes)
    {
        writef("  N=%5d:", targetSize);
        foreach (ref row; summary)
            if (row.nTarget == targetSize)
                writef("  beta=%.2f -> F=%.4f+/-%.4f", row.betaV, row.force, row.forceErr);
        writeln();
    }
    stdout.flush();
}

}
