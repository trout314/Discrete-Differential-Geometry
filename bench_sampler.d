/// Benchmark for the MCMC sampler hot path.
/// Compares baseline (execute-compute-undo) vs speculative delta with full objective.
module bench_sampler;

version (unittest) {} else {

import std.algorithm, std.array, std.conv, std.datetime.stopwatch, std.math,
    std.range, std.stdio, std.string;
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
    writeln("=== Growth strategy comparison (3 strategies) ===");
    writeln("strategy,sweep,acc_rate,objective,codim3_var,deg_var_vtx,nTets,elapsed_s");
    stdout.flush();

    auto seed = loadManifold!dim("standard_triangulations/dim_3_sphere.mfd");
    enum targetSize = 10_000;
    enum totalSweeps = 200;

    // --- Strategy 1: Stellar subdivision then equilibrate ---
    {
        writeln("# Strategy 1: stellar growth then equilibrate");
        stdout.flush();

        auto mfd = seed;
        int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
        while (mfd.fVector[dim] < targetSize)
        {
            auto facet = mfd.randomFacetOfDim(dim).array;
            auto move = BistellarMove!dim(facet, [nextVertex]);
            if (!mfd.contains(move.coCenter))
            {
                mfd.doMove(move);
                nextVertex++;
            }
        }

        runEquilibration(mfd, "stellar", targetSize, totalSweeps, nextVertex);
    }

    // --- Strategy 2: Penalty-driven growth (fixed target from start) ---
    {
        writeln("# Strategy 2: penalty-driven growth (fixed target)");
        stdout.flush();

        auto mfd = seed;
        int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
        int[] unusedVertices = [nextVertex];

        auto params = BenchParams(targetSize);
        auto currentObjective = computeObjective(mfd, params);

        StopWatch timer;
        timer.start();

        ulong totalAccepted = 0;
        ulong totalTried = 0;

        bool reachedTarget = false;
        int sweepsSinceTarget = 0;
        int movesInSweep = 0;

        while (sweepsSinceTarget <= totalSweeps)
        {
            if (unusedVertices.empty)
            {
                nextVertex++;
                unusedVertices ~= nextVertex;
            }

            auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
            auto bm = proposed.move;
            totalTried++;

            real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);

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
                totalAccepted++;

                if (bm.coCenter.length == 1 && bm.coCenter[0] >= nextVertex)
                    nextVertex = bm.coCenter[0] + 1;

                if (reachedTarget) movesInSweep++;
            }

            if (!reachedTarget && mfd.fVector[dim] >= targetSize)
            {
                reachedTarget = true;
                movesInSweep = 0;
                logSweep(mfd, "penalty", 0, totalAccepted, totalTried, params, timer);
            }

            if (reachedTarget && movesInSweep >= targetSize)
            {
                sweepsSinceTarget++;
                movesInSweep = 0;
                logSweep(mfd, "penalty", sweepsSinceTarget, totalAccepted, totalTried, params, timer);
            }
        }
    }

    // --- Strategy 3: Ramped target volume ---
    {
        writeln("# Strategy 3: ramped target volume (+500 per step, 5 eq sweeps each)");
        stdout.flush();

        auto mfd = seed;
        int nextVertex = mfd.simplices(0).joiner.array.dup.sort.array.back + 1;
        int[] unusedVertices = [nextVertex];

        enum stepSize = 500;
        enum eqSweepsPerStep = 5;

        StopWatch timer;
        timer.start();

        ulong totalAccepted = 0;
        ulong totalTried = 0;

        // Ramp up from seed size to targetSize in increments
        for (int currentTarget = stepSize; currentTarget <= targetSize; currentTarget += stepSize)
        {
            auto params = BenchParams(currentTarget);
            auto currentObjective = computeObjective(mfd, params);

            // Run until we reach currentTarget, then do eqSweepsPerStep sweeps
            bool reachedStep = (mfd.fVector[dim] >= currentTarget);
            int eqSweepsDone = 0;
            int movesInSweep = 0;
            int sweepSize = max(currentTarget, 1000);

            while (eqSweepsDone < eqSweepsPerStep)
            {
                if (unusedVertices.empty)
                {
                    nextVertex++;
                    unusedVertices ~= nextVertex;
                }

                auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
                auto bm = proposed.move;
                totalTried++;

                real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);

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
                    totalAccepted++;

                    if (bm.coCenter.length == 1 && bm.coCenter[0] >= nextVertex)
                        nextVertex = bm.coCenter[0] + 1;

                    if (reachedStep) movesInSweep++;
                }

                if (!reachedStep && mfd.fVector[dim] >= currentTarget)
                {
                    reachedStep = true;
                    movesInSweep = 0;
                }

                if (reachedStep && movesInSweep >= sweepSize)
                {
                    eqSweepsDone++;
                    movesInSweep = 0;
                }
            }
        }

        // Log sweep 0 at the moment we finish ramping
        auto finalParams = BenchParams(targetSize);
        auto currentObjective = computeObjective(mfd, finalParams);
        logSweep(mfd, "ramped", 0, totalAccepted, totalTried, finalParams, timer);

        // Now run totalSweeps of equilibration at final target
        runEquilibrationFromState(mfd, "ramped", targetSize, totalSweeps,
            nextVertex, unusedVertices, totalAccepted, totalTried, timer);
    }

    writeln("done");
    stdout.flush();
}

void logSweep(ref const Manifold!dim mfd, string label, int sweep,
    ulong totalAccepted, ulong totalTried, BenchParams params, ref StopWatch timer)
{
    auto pen = mfd.penalties(params);
    auto obj = objectiveFromPenalty(pen, params);
    auto elapsed = timer.peek.total!"msecs" / 1000.0;
    writefln("%s,%d,%.4f,%.2f,%.2f,%.3f,%d,%.1f",
        label, sweep,
        totalTried > 0 ? cast(real) totalAccepted / totalTried : 0.0,
        obj, pen.localSolidAngleCurvPenalty,
        mfd.degreeVariance(0), mfd.fVector[dim], elapsed);
    stdout.flush();
}

void runEquilibration(ref Manifold!dim mfd, string label, int targetSize,
    int totalSweeps, int nextVertex)
{
    auto params = BenchParams(targetSize);
    auto currentObjective = computeObjective(mfd, params);
    int[] unusedVertices = [nextVertex];

    ulong totalAccepted = 0;
    ulong totalTried = 0;

    StopWatch timer;
    timer.start();

    foreach (sweep; 0 .. totalSweeps)
    {
        logSweep(mfd, label, sweep, totalAccepted, totalTried, params, timer);

        int accepted = 0;
        while (accepted < targetSize)
        {
            if (unusedVertices.empty)
            {
                nextVertex++;
                unusedVertices ~= nextVertex;
            }

            auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
            auto bm = proposed.move;
            totalTried++;

            real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);

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
                accepted++;
                totalAccepted++;

                if (bm.coCenter.length == 1 && bm.coCenter[0] >= nextVertex)
                    nextVertex = bm.coCenter[0] + 1;
            }
        }
    }
}

void runEquilibrationFromState(ref Manifold!dim mfd, string label, int targetSize,
    int totalSweeps, ref int nextVertex, ref int[] unusedVertices,
    ref ulong totalAccepted, ref ulong totalTried, ref StopWatch timer)
{
    auto params = BenchParams(targetSize);
    auto currentObjective = computeObjective(mfd, params);

    foreach (sweep; 1 .. totalSweeps + 1)
    {
        logSweep(mfd, label, sweep, totalAccepted, totalTried, params, timer);

        int accepted = 0;
        while (accepted < targetSize)
        {
            if (unusedVertices.empty)
            {
                nextVertex++;
                unusedVertices ~= nextVertex;
            }

            auto proposed = proposeMoveOptimized(mfd, unusedVertices.back);
            auto bm = proposed.move;
            totalTried++;

            real deltaObj = speculativeDelta(mfd, bm, currentObjective, params);

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
                accepted++;
                totalAccepted++;

                if (bm.coCenter.length == 1 && bm.coCenter[0] >= nextVertex)
                    nextVertex = bm.coCenter[0] + 1;
            }
        }
    }
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
            writef("  speculative:");
            benchmarkSpeculative!true(mfd1, actualTets);

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
