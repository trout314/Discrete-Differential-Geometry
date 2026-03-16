/// Samples from manifolds (TO DO: better description)
module applications.manifold_sampler;
import std.algorithm, std.array, std.conv, std.datetime, std.datetime.stopwatch,
    std.format, std.getopt, std.math, std.range, std.stdio, std.string,
    std.stdio;
import std.random : rndChoice = choice;
import std.random : uniform, uniform01;
import core.memory : GC;
import algorithms, manifold, manifold_examples, manifold_moves,
    sampler, simplicial_complex, utility;

enum dim = 3;

version (unittest) {} else {
int main(string[] args)
{
    string paramFileName;
    auto helpInformation = getopt(args, std.getopt.config.required,
        "parameterFile|p", &paramFileName);

    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        return 0; // exit with return code zero (success)
    }

    static immutable parametersUsed = [
        ["bool", "disableGC"],
        ["bool", "checkForProblems"],
        ["bool", "saveEdgeGraph"],
        ["bool", "saveDualGraph"],
        ["bool", "doStdoutReport"],
        ["bool", "writeCSVdataFile"],
        ["int", "startingSampleNumber"],
        ["int", "numFacetsTarget"],
        ["int", "movesTriedPerProblemCheck"],
        ["int", "movesTriedPerStdoutReport"],
        ["int", "movesTriedPerGCcollect"],
        ["int", "coDim2DegreeHistogramBins"],
        ["int", "coDim3DegreeHistogramBins"],
        ["int", "coDim2DegreeHistogramBars"],
        ["int", "coDim3DegreeHistogramBars"],
        ["real", "sweepsPerSample"],
        ["real", "sweepsPerCSVline"],
        ["real", "hingeDegreeTarget"],
        ["real", "numFacetsCoef"],
        ["real", "numHingesCoef"],
        ["real", "hingeDegreeVarianceCoef"],
        ["real", "coDim3DegreeVarianceCoef"],
        ["real", "maxSweeps"],
        ["int", "growthStepSize"],
        ["int", "eqSweepsPerStep"],
        ["string", "sampleFilesPrefix"],
        ["string", "initialManifoldFile"]
    ];

    auto params = parseParameterFile!parametersUsed(paramFileName);
    Manifold!dim mfd = loadManifold!dim(params.initialManifoldFile);

    auto initialVertices = mfd.simplices(0).joiner.array.dup.sort.array;
    assert(initialVertices.length > 0, "initial manifold is empty");
    int[] unusedVertices = getUnusedVertices(mfd, initialVertices);

    // bistellarTries[j] counts the (j + 1) -> (dim + 1 - j) moves attempted
    ulong[dim + 1] bistellarTries;

    // bistellarAccepts[j] counts the (j + 1) -> (dim + 1 - j) moves accepted
    ulong[dim + 1] bistellarAccepts;

    ulong hingeTries, hingeAccepts;
    enum hingeMoveProb = 0.3;

    auto sampleThreshold = params.numFacetsTarget * params.sweepsPerSample;
    auto csvLineThreshold = params.numFacetsTarget * params.sweepsPerCSVline;
    
    ulong sampleNumber = params.startingSampleNumber;
    ulong columnReportNumber = 0;

    if (params.disableGC)
    {
        GC.disable;
    }

    auto startTime = Clock.currTime;
    StopWatch timer;
    timer.start;    
    
    //----------------------- RAMPED GROWTH PHASE --------------------------
    // Grow from seed size to numFacetsTarget in steps, equilibrating at each.
    if (mfd.fVector[dim] < params.numFacetsTarget)
    {
        // Start from the first step boundary above current size
        int startTarget = params.growthStepSize;
        while (startTarget <= cast(int) mfd.fVector[dim])
            startTarget += params.growthStepSize;

        for (int stepTarget = startTarget; stepTarget <= params.numFacetsTarget;
             stepTarget += params.growthStepSize)
        {
            // Override numFacetsTarget for this step
            auto stepParams = params;
            stepParams.numFacetsTarget = stepTarget;
            auto stepObjective = mfd.objective(stepParams);

            // Phase 1: MCMC until we reach stepTarget
            while (mfd.fVector[dim] < stepTarget)
            {
                mfd.mcmcStep(stepObjective, unusedVertices, stepParams,
                    hingeMoveProb, hingeTries, hingeAccepts,
                    bistellarTries, bistellarAccepts);
            }

            // Phase 2: equilibration sweeps at this step size
            ulong eqAccepted = 0;
            immutable eqTarget = cast(ulong) params.eqSweepsPerStep * stepTarget;
            while (eqAccepted < eqTarget)
            {
                if (mfd.mcmcStep(stepObjective, unusedVertices, stepParams,
                        hingeMoveProb, hingeTries, hingeAccepts,
                        bistellarTries, bistellarAccepts))
                    ++eqAccepted;
            }

            if (params.doStdoutReport)
            {
                "Ramping: step %d/%d, %d tets".writefln(
                    stepTarget, params.numFacetsTarget, mfd.fVector[dim]);
                stdout.flush();
            }
        }

        // Reset counters so sampling phase starts fresh
        bistellarTries[] = 0;
        bistellarAccepts[] = 0;
    }

    auto currentObjective = mfd.objective(params);
    auto doneSampling = false;

    while (!doneSampling)
    {
        ulong numMovesTried = bistellarTries[].sum + hingeTries;
        ulong numMovesAccepted = bistellarAccepts[].sum + hingeAccepts;
        doneSampling = numMovesAccepted >= params.maxSweeps * params.numFacetsTarget;

        if (unusedVertices.empty)
        {
            /* If there are no unused vertices, then the next unused vertex
            label is just the number of vertices. */
            // TO DO: Make sure this invariant is maintained on load / save, etc.
            unusedVertices ~= mfd.fVector[0].to!int;
        }
        assert(unusedVertices.all!(v => !mfd.contains(v.only)));

        //-------------------------- ATTEMPT MOVE ----------------------------

        bool accepted = false;
        version (TrackValidMoves)
        {
            // Exact Hastings: execute, compute, undo if rejected.
            // Uses incrementally tracked V_before/V_after for correction.
            auto bm = mfd.chooseRandomMove(unusedVertices.back, params);
            ++bistellarTries[bm.coCenter.length - 1];

            immutable vBefore = cast(real) mfd.validMoveCount;

            mfd.doMove(bm);
            if (bm.coCenter.length == 1) unusedVertices.popBack;
            if (bm.center.length == 1) unusedVertices ~= bm.center;

            real newObjective = mfd.objective(params);
            real deltaObj = newObjective - currentObjective;
            immutable vAfter = cast(real) mfd.validMoveCount;

            real logAlpha = -deltaObj + log(vBefore) - log(vAfter);

            if ((logAlpha < 0) && (uniform01 > exp(logAlpha)))
            {
                // Rejected — undo
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
                currentObjective = newObjective;
                ++bistellarAccepts[bm.coCenter.length - 1];
                accepted = true;
            }
        }
        else
        {
            // Speculative path via shared mcmcStep (no Hastings correction;
            // importance weight 1/V corrects at measurement time).
            accepted = mfd.mcmcStep(currentObjective, unusedVertices, params,
                hingeMoveProb, hingeTries, hingeAccepts,
                bistellarTries, bistellarAccepts);
        }

        //--------------------------- MAKE REPORT ----------------------------

        if (params.doStdoutReport && 
            ((numMovesTried % params.movesTriedPerStdoutReport == 0) || doneSampling))
        {
            stdout.write("\033c");  // Clear the screen
            stdout.writeReports(mfd, startTime, timer, bistellarTries[],
                bistellarAccepts[], params);
            stdout.flush();
        }

        //---------------------- WRITE TO CSV DATA FILE ----------------------
        if ((numMovesAccepted == csvLineThreshold) || doneSampling)
        {
            mfd.writeCSVreport(columnReportNumber, startTime, timer, bistellarTries[],
                bistellarAccepts[], params);
            ++columnReportNumber;
            csvLineThreshold += params.numFacetsTarget * params.sweepsPerCSVline;
        }

        //----------------------- SAVE CURRENT MANIFOLD ----------------------
        if ((numMovesAccepted == sampleThreshold) || doneSampling)
        {
            mfd.saveSample(sampleNumber, startTime, timer,
                bistellarTries[], bistellarAccepts[], params);
            ++sampleNumber;
            sampleThreshold += params.numFacetsTarget * params.sweepsPerSample;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if (params.disableGC && (numMovesTried % params.movesTriedPerProblemCheck == 0))
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        if (params.checkForProblems && (numMovesAccepted % params.movesTriedPerProblemCheck == 0))
        {
            auto problems = mfd.findProblems;
            if (!problems.empty)
            {
                problems.each!writeln;
                assert(0);
            }
        }
    }

    "Finished! Time elapsed: %s".writefln(Clock.currTime.to!DateTime - startTime.to!DateTime);
    return 0; // exit with return code zero (success)
}}

void writeReports(M, S, T, W, P)(
    W sink,
    M mfd,
    S startTime,
    T timer,
    ulong[] bistellarTries,
    ulong[] bistellarAccepts,
    P params)
{
    auto numMovesTried = bistellarTries.sum;
    auto numMovesAccepted = bistellarAccepts.sum;

    typeof(timer.peek()) timePerMove;
    if (numMovesTried > 0)
    {
        timePerMove = timer.peek / numMovesTried;
    }
    auto acceptFrac = double(numMovesAccepted) / numMovesTried;

    sink.writeTimingAndTargetsReport(mfd, numMovesAccepted, startTime, timePerMove, acceptFrac, params);
    sink.writeSimplexReport(mfd);
    sink.writeObjectiveReport(mfd, params);
    sink.writeMoveReport(bistellarTries, bistellarAccepts, params);
    sink.writeHistogramReport(mfd, params);
}

void writeCSVreport(M, S, T, P)(M manifold, ulong reportNumber, S startTime, T timer, ulong[] bistellarTries,
    ulong[] bistellarAccepts, P params)
{
    auto file = File(params.sampleFilesPrefix ~ ".dat", "a");
    auto dim = manifold.dimension;

    string[] columnLabels;
    string[] values;

    columnLabels ~= "total_moves_accepted";
    values ~= bistellarAccepts[].sum.to!string;

    columnLabels ~= "objective";
    values ~= manifold.objective(params).to!string;

    columnLabels ~= "volume_penalty";
    values ~= manifold.penalties(params).volumePenalty.to!string;

    columnLabels ~= "global_curvature_penalty";
    values ~= manifold.penalties(params).globalCurvPenalty.to!string;

    columnLabels ~= "local_curvature_variance_penalty";
    values ~= manifold.penalties(params).localCurvPenalty.to!string;

    columnLabels ~= "local_solid_angle_curvature_variance_penalty";
    values ~= manifold.penalties(params).localSolidAngleCurvPenalty.to!string;

    columnLabels ~= "total_moves_tried";
    values ~= bistellarTries[].sum.to!string;

    foreach(n; 1 .. dim+2)
    {
        columnLabels ~= "num_%s_%s_bistellar_accepted".format(n, dim + 2 -n);
        values ~= bistellarAccepts[n-1].to!string;

        columnLabels ~= "num_%s_%s_bistellar_tried".format(n, dim + 2 -n);
        values ~= bistellarTries[n-1].to!string;
    }

    foreach(d; 0 .. dim+1)
    {
        columnLabels ~= "num_%s_simplices".format(d);
        values ~= manifold.fVector[d].to!string;

        columnLabels ~= "deg_mean_%s_simplices".format(d);
        values ~= manifold.meanDegree(d).to!string;

        columnLabels ~= "deg_var_%s_simplices".format(d);
        values ~= manifold.degreeVariance(d).to!string;
    }

    auto maxDeg2 = 2 + params.coDim2DegreeHistogramBins;
    auto hist2 = manifold.degreeHistogram(dim - 2);
    foreach(deg; 3 .. 3 + params.coDim2DegreeHistogramBins)
    {
        columnLabels ~= "codim2_simps_of_deg_%s".format(deg);
        if (deg-1 < hist2.length)
        {
            values ~= hist2[deg-1].to!string;
        }
        else
        {
            values ~= "0";
        }
    }
    auto tail2 = (hist2.length >= maxDeg2) ? hist2[maxDeg2 .. $].sum : 0;
    columnLabels ~= "codim2_deg_tail";
    values ~= tail2.to!string;

    auto maxDeg3 = 2 + 2 * params.coDim3DegreeHistogramBins;
    auto hist3 = manifold.degreeHistogram(dim - 3);
    foreach(deg; iota(4, 4 + 2*params.coDim3DegreeHistogramBins, 2))
    {
        columnLabels ~= "codim3_simps_of_deg_%s".format(deg);
        if (deg - 1 < hist3.length)
        {
            values ~= hist3[deg-1].to!string;
        }
        else
        {
            values ~= "0";
        }
    }
    auto tail3 = (hist3.length >= maxDeg3) ? hist3[maxDeg3 .. $].sum : 0;    
    columnLabels ~= "codim3_deg_tail";
    values ~= tail3.to!string;    
    

    assert(values.length == columnLabels.length);
    if (reportNumber == 0)
    {
        file.writeln(columnLabels.joiner(","));
    }

    file.writeln(values.joiner(","));
}

void writeTimingAndTargetsReport(M, S, T, W, P)(W sink, M mfd, ulong total_moves_accepted, S startTime, T timePerMove, double acceptFrac, P params)
{
    typeof(timePerMove) timePerSweep;
    if (acceptFrac > 0.0)
    {
        timePerSweep = timePerMove * (mfd.fVector[dim] / acceptFrac).to!ulong;
    }

    sink.writeln("-".repeat(80).joiner);
    string timeInfo = "%.1f / %s".format(total_moves_accepted / double(params.numFacetsTarget), params.maxSweeps);
    auto prettyStartTime = startTime.to!string.split('.').front;

    sink.writef("sweeps         : %23-s |", timeInfo);
    if (params.numFacetsCoef > 0.0)
    {
        sink.writef(" target # facets     : %s", params.numFacetsTarget);
    }
    sink.writeln;

    sink.writef("started at     : %23-s |", prettyStartTime);
    if (params.numHingesCoef > 0.0)
    {
        sink.writef(" target hinge degree : %.5f", params.hingeDegreeTarget);
    }
    sink.writeln;

    sink.writefln("Δt/move        : %23-s |", timePerMove.getFirstTwoParts);
    sink.writefln("Δt/sweep       : %23-s |",timePerSweep.getFirstTwoParts);
    sink.writefln("move accept %%  : %23.3-f |", acceptFrac);
}

void writeSimplexReport(int dim, Vertex, W)(W sink, const ref Manifold!(dim, Vertex) mfd)
{
    sink.writeln("-".repeat(80).joiner);
    sink.write("dimension      : ");
    (dim + 1).iota.each!(d => sink.writef("%12s ", d));
    sink.writeln;
    sink.write("# simplices    : ");
    (dim + 1).iota.each!(d => sink.writef("%12s ", mfd.fVector[d]));
    sink.writeln;
    sink.write("degree mean    : ");
    (dim + 1).iota.each!(d => sink.writef("%12.4f ", mfd.meanDegree(d)));
    sink.writeln;
    sink.write("degree std dev : ");
    (dim + 1).iota.each!(d => sink.writef("%12.4f ", mfd.degreeVariance(d).sqrt));
    sink.writeln;
}

void writeObjectiveReport(int dim, Vertex, W, P)(W sink, const ref Manifold!(dim, Vertex) mfd, P params)
{
    sink.writeln("-".repeat(80).joiner);
    sink.writeln("    Penalty", ' '.repeat(10), "Raw",
        ' '.repeat(11), "Coef", ' '.repeat(9), "Value");

    if (params.numFacetsCoef > 0.0)
    {
        auto vp = mfd.penalties(params).volumePenalty;
        sink.writefln("# facets       : %.6e  *  %.4f  =  %.6e",
            vp, params.numFacetsCoef, params.numFacetsCoef * vp);
    }

    if (params.numHingesCoef > 0.0)
    {
        auto gcp = mfd.penalties(params).globalCurvPenalty;
        sink.writefln("# hinges       : %.6e  *  %.4f  =  %.6e",
            gcp, params.numHingesCoef, params.numHingesCoef * gcp);
    }

    static if (dim > 2)
    {
        if (params.hingeDegreeVarianceCoef > 0.0)
        {
            auto lcp = mfd.penalties(params).localCurvPenalty;
            sink.writefln("hinge deg var  : %.6e  *  %.4f  =  %.6e",
                lcp, params.hingeDegreeVarianceCoef, params.hingeDegreeVarianceCoef * lcp);
        }

        if (params.coDim3DegreeVarianceCoef > 0.0)
        {
            auto lsacp = mfd.penalties(params).localSolidAngleCurvPenalty;
            sink.writefln("codim-3 deg var: %.6e  *  %.4f  =  %.6e",
                lsacp, params.coDim3DegreeVarianceCoef, params.coDim3DegreeVarianceCoef * lsacp);
        }
    }
    sink.writefln("total penalty  :                             %.6e",
        mfd.objective(params));
}

void writeMoveReport(S, W, P)(W sink, S bistellarTries, S bistellarAccepts, P params)
{
    sink.writeln("-".repeat(80).joiner);
    sink.writeln("Bistellar Moves      # Accepted          # Tried    %    ");
    foreach (i; 0 .. dim + 1)
    {
        auto accepts = "%10,d".format(bistellarAccepts[i]);
        auto tries = "%10,d".format(bistellarTries[i]);
        sink.writefln("    %2s → %2-s    : %14s / %14s (%5.3f)",
            i + 1, dim + 1 - i, accepts, tries,
            double(bistellarAccepts[i]) / bistellarTries[i]);
    }
    auto totAccepts = bistellarAccepts[].sum;
    auto totTries = bistellarTries[].sum;
    auto accepts = "%10,d".format(totAccepts);
    auto tries = "%10,d".format(totTries);
    sink.writefln("Total Moves    : %14s / %14s (%5.3f)", accepts, tries,
        double(totAccepts) / (totTries));

}

void writeHistogramReport(int dim, Vertex, W, P)(W sink, const ref Manifold!(dim, Vertex) mfd, P params)
{
    sink.writeln("-".repeat(80).joiner);
    auto maxDeg2 = 2 + params.coDim2DegreeHistogramBins;
    auto maxDeg3 = 2 + 2 * params.coDim3DegreeHistogramBins;

    auto hist2 = mfd.degreeHistogram(mfd.dimension - 2);
    auto tail2 = (hist2.length >= maxDeg2) ? hist2[maxDeg2 .. $].sum : 0;
    auto maxDegBin2 = max(hist2.maxElement, tail2);
    auto normedHist2 = hist2.map!(freq => real(freq) / maxDegBin2);

    auto hist3 = mfd.degreeHistogram(mfd.dimension - 3);
    auto tail3 = (hist3.length >= maxDeg3) ? hist3[maxDeg3 .. $].sum : 0;

    static if (dim > 2)
    {
        auto maxDegBin3 = max(hist3.maxElement, tail3);
        auto normedHist3 = hist3.map!(freq => real(freq) / maxDegBin3);
    }

    sink.writeln(' '.repeat(27), "Codimension-2 Degree Histogram");

    static immutable bars = [
        '▏', '▎', '▍', '▌', '▋', '▊', '▉', '█'
    ];

    foreach (bin; iota(0, params.coDim2DegreeHistogramBins))
    {
        sink.writef("%2s ", bin + 3);
        if (bin + 2 < normedHist2.length)
        {
            auto nEighths = (params.coDim2DegreeHistogramBars * normedHist2[bin + 2] * 8.0).to!int;
            auto bar = bars.back.repeat(nEighths / 8).array;
            if (nEighths % 8 > 0)
            {
                bar ~= bars[nEighths % 8];
            }
            sink.writef("%-*s", params.coDim2DegreeHistogramBars, bar);
        }
        sink.writeln;
    }

    auto tailFreq2 = real(tail2) / maxDegBin2;    
    auto nEighths2 = (params.coDim2DegreeHistogramBars * tailFreq2 * 8.0).to!int;
    auto bar2 = bars.back.repeat(nEighths2 / 8).array;
    if (nEighths2 % 8 > 0)
    {
        bar2 ~= bars[nEighths2 % 8];
    }
    sink.writefln(" > %-*s", params.coDim2DegreeHistogramBars, bar2);
    

    static if (dim > 2)
    {
        sink.writeln("-".repeat(80).joiner);
        sink.writeln(' '.repeat(27), "Codimension-3 Degree Histogram");
        foreach (bin; iota(0, params.coDim3DegreeHistogramBins))
        {
            sink.writef("%2s ", 4 + bin * 2);
            if ((3 + bin * 2) < normedHist3.length)
            {
                auto nEighths = (params.coDim3DegreeHistogramBars * normedHist3[3 + bin * 2] * 8.0).to!int;
                auto bar = bars.back.repeat(nEighths / 8).array;
                if (nEighths % 8 > 0)
                {
                    bar ~= bars[nEighths % 8];
                }
                sink.writef("%-*s", params.coDim3DegreeHistogramBars, bar);
            }
            sink.writeln;
        }

        auto tailFreq3 = real(tail3) / maxDegBin3;
        sink.write(" > ");
        auto nEighths3 = (params.coDim3DegreeHistogramBars * tailFreq3 * 8.0).to!int;
        auto bar3 = bars.back.repeat(nEighths3 / 8).array;
        if (nEighths3 % 8 > 0)
        {
            bar3 ~= bars[nEighths3 % 8];
        }
        sink.writefln("%-*s", params.coDim3DegreeHistogramBars, bar3);
    }
}

auto getFirstTwoParts(T)(T time)
{
    return time.to!string.split(",").take(2).join(",").findSplit("and")[0]
        .replace("minute", "min");

}


void saveSample(M, S, T, P)(M manifold, ulong sampleNumber,
                S startTime, T timer, ulong[] bistellarTries, ulong[] bistellarAccepts,
                P parameters)
{
    string prefix = parameters.sampleFilesPrefix ~ "_sample_"
        ~ sampleNumber.to!string;
    
    manifold.saveTo(prefix ~ ".mfd");
        
    auto saveFile = File(prefix ~ ".mfd", "a");       
    saveFile.writeln;
    saveFile.writeReports(manifold, startTime, timer,
        bistellarTries[], bistellarAccepts[], parameters);
    saveFile.writeln(toPrettyString(parameters));

    // To make sure we know history of a sample, we copy the
    // ending info from the initial manifold file 
    File(parameters.initialManifoldFile, "r")
        .byLineCopy
        .find!(line => line.startsWith("---"))
        .each!(line => saveFile.writeln(line));

    if (parameters.saveEdgeGraph)
    {
        manifold.saveEdgeGraphTo(prefix ~ ".edge_graph");
    }

    if (parameters.saveDualGraph)
    {
        manifold.saveDualGraphTo(prefix ~ ".dual_graph");
    }
}