/// Samples from manifolds (TO DO: better description)
module applications.manifold_sampler;
import std.algorithm, std.array, std.conv, std.datetime, std.datetime.stopwatch,
    std.format, std.getopt, std.math, std.range, std.stdio, std.string,
    std.sumtype, std.stdio;
import std.random : rndChoice = choice;
import std.random : uniform, uniform01;
import core.memory : GC;
import algorithms, manifold, manifold_examples, manifold_moves, polygons,
    simplicial_complex, utility;
import unit_threaded;

enum dim = 3;
enum maxHingeMoveDeg = 4;

static assert(maxHingeMoveDeg >= 4, "Hinge move degrees must be at least 4");

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
        ["bool", "useHingeMoves"],
        ["bool", "checkForProblems"],
        ["bool", "saveEdgeGraph"],
        ["bool", "saveDualGraph"],
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
        ["string", "savedFilesPrefix"],
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

    // hingeTries[j] counts the hinge-moves with hinge degree j + 4 tried
    ulong[maxHingeMoveDeg - 3] hingeTries;
 
    // hingeAccepts[j] counts the hinge-moves with hinge degree j + 4 accepted
    ulong[maxHingeMoveDeg - 3] hingeAccepts;

    auto sampleThreshold = params.numFacetsTarget * params.sweepsPerSample;
    auto csvLineThreshold = params.numFacetsTarget * params.sweepsPerCSVline;
    
    ulong sampleNumber = 0;
    ulong columnReportNumber = 0;

    if (params.disableGC)
    {
        GC.disable;
    }

    auto startTime = Clock.currTime;
    StopWatch timer;
    timer.start;    
    
    auto currentObjective = mfd.objective(params);
    auto doneSampling = false;

    while (!doneSampling)
    {
        ulong numMovesTried = bistellarTries[].sum + hingeTries[].sum;
        ulong numMovesAccepted = bistellarAccepts[].sum + hingeAccepts[].sum;
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
        auto chosenMove = mfd.chooseRandomMove(unusedVertices.back, params);
        mfd.doMove(chosenMove);        
        updateUnusedVertices(unusedVertices, chosenMove);
        incrementMoveCounts(bistellarTries, hingeTries, chosenMove);
        incrementMoveCounts(bistellarAccepts, hingeAccepts, chosenMove);

        real newObjective = mfd.objective(params);
        real deltaObj = newObjective - currentObjective;
 
        // Reject move and undo it, if appropriate
        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        {
            mfd.undoMove(chosenMove);
            undoUpdateUnusedVertices(unusedVertices, chosenMove);
            decrementMoveCounts(bistellarAccepts, hingeAccepts, chosenMove);
        }
        else
        {
            currentObjective = newObjective;
        }

        //--------------------------- MAKE REPORT ----------------------------
        if ((numMovesTried % params.movesTriedPerStdoutReport == 0) || doneSampling)
        {
            stdout.write("\033c");  // Clear the screen
            stdout.writeReports(mfd, startTime, timer, bistellarTries[],
                bistellarAccepts[], hingeTries[], hingeAccepts[], params);
            stdout.flush();
        }

        //---------------------- WRITE TO CSV DATA FILE ----------------------
        if ((numMovesAccepted == csvLineThreshold) || doneSampling)
        {
            mfd.writeCSVreport(columnReportNumber, startTime, timer, bistellarTries[],
                bistellarAccepts[], hingeTries[], hingeAccepts[], params);
            ++columnReportNumber;
            csvLineThreshold += params.numFacetsTarget * params.sweepsPerCSVline;
        }

        //----------------------- SAVE CURRENT MANIFOLD ----------------------
        if ((numMovesAccepted == sampleThreshold) || doneSampling)
        {
            mfd.saveSample(sampleNumber, startTime, timer,
                bistellarTries[], bistellarAccepts[], hingeTries[], hingeAccepts[], params);
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
    ulong[] hingeTries,
    ulong[] hingeAccepts,
    P params)
{
    auto numMovesTried = bistellarTries.sum + hingeTries.sum;
    auto numMovesAccepted = bistellarAccepts.sum + hingeAccepts.sum;
    
    typeof(timer.peek()) timePerMove;
    if (numMovesTried > 0)
    {
        timePerMove = timer.peek / numMovesTried;
    }
    auto acceptFrac = double(numMovesAccepted) / numMovesTried;

    sink.writeTimingAndTargetsReport(mfd, numMovesAccepted, startTime, timePerMove, acceptFrac, params);
    sink.writeSimplexReport(mfd);
    sink.writeObjectiveReport(mfd, params);
    sink.writeMoveReport(bistellarTries, bistellarAccepts, hingeTries, hingeAccepts, params);
    sink.writeHistogramReport(mfd, params);
}


Vertex[] getUnusedVertices(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd, Vertex[] initialVertices)
{
    int[] unusedVertices;
    // all gaps in list of vertices should be unusued vertices
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

private struct Penalty
{
    real volumePenalty;
    real globalCurvPenalty;
    real localCurvPenalty;
    real localSolidAngleCurvPenalty;
}

Penalty penalties(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    immutable hingesPerFacet = dim * (dim + 1) / 2;
    immutable coDim3PerFacet = binomial(dim + 1, dim - 2);

    immutable nFacets = mfd.fVector[dim];
    immutable nHinges = mfd.fVector[dim - 2];
    immutable totSqDeg = mfd.totalSquareDegree(dim - 2);

    immutable nHingesTarget = hingesPerFacet * nFacets / params.hingeDegreeTarget;
    immutable degTarget = hingesPerFacet * nFacets / nHinges.to!real;

    static if (dim > 2)
    {
        immutable nCoDim3 = mfd.fVector[dim - 3];
        immutable totSAsqDeg = mfd.totalSquareDegree(dim - 3);
        immutable coDim3DegTarget = coDim3PerFacet * nFacets / nCoDim3.to!real;        
    }

    Penalty penalty;
    penalty.volumePenalty = (nFacets - params.numFacetsTarget) ^^ 2;
    penalty.globalCurvPenalty = (nHinges - nHingesTarget) ^^ 2;

    real _; // dummy for integer part
    real x = modf(degTarget, _); // fractional part
    real minPenalty = (x - x ^^ 2) * nHinges;

    // TO DO: Refer to external paper for this!
    penalty.localCurvPenalty = (
        degTarget ^^ 2 * nHinges - 2 * degTarget * hingesPerFacet * nFacets + totSqDeg) - minPenalty;
    
    static if (dim > 2)
    {
        x = modf(coDim3DegTarget, _); // fractional part
        minPenalty = (x - x ^^ 2) * nCoDim3;

        // TO DO: Refer to arxiv paper for this!
        penalty.localSolidAngleCurvPenalty = (
            coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + totSAsqDeg) - minPenalty;        
    }
    else
    {
        penalty.localSolidAngleCurvPenalty = 0;
    }

    return penalty;
}

real objective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    auto pen = mfd.penalties(params);
    return params.numFacetsCoef * pen.volumePenalty
        + params.numHingesCoef * pen.globalCurvPenalty
        + params.hingeDegreeVarianceCoef * pen.localCurvPenalty
        + params.coDim3DegreeVarianceCoef * pen.localSolidAngleCurvPenalty;
}

void writeCSVreport(M, S, T, P)(M manifold, ulong reportNumber, S startTime, T timer, ulong[] bistellarTries,
    ulong[] bistellarAccepts, ulong[] hingeTries, ulong[] hingeAccepts, P params)
{
    auto file = File(params.savedFilesPrefix ~ ".dat", "a");
    auto dim = manifold.dimension;

    string[] columnLabels;
    string[] values;

    columnLabels ~= "total_moves_accepted";
    values ~= (bistellarAccepts[].sum + hingeAccepts[].sum).to!string;

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
    values ~= (bistellarTries[].sum + hingeTries[].sum).to!string;

    foreach(n; 1 .. dim+2)
    {
        columnLabels ~= "num_%s_%s_bistellar_accepted".format(n, dim + 2 -n);
        values ~= bistellarAccepts[n-1].to!string;

        columnLabels ~= "num_%s_%s_bistellar_tried".format(n, dim + 2 -n);
        values ~= bistellarTries[n-1].to!string;
    }

    foreach(n; 0 .. hingeAccepts.length)
    {
        columnLabels ~= "num_%s_%s_hinge_accepted".format(n+4, (n+2)*(dim-1));
        values ~= hingeAccepts[n].to!string;

        columnLabels ~= "num_%s_%s_hinge_tried".format(n+4, (n+2)*(dim-1));
        values ~= hingeTries[n].to!string;
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

    sink.writefln("Δt/move        : %23-s |", timePerMove.getFirstPart);
    sink.writefln("Δt/sweep       : %23-s |",timePerSweep.getFirstPart);
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

void writeMoveReport(S, T, W, P)(W sink, S bistellarTries, S bistellarAccepts, T hingeTries, T hingeAccepts, P params)
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
    if (params.useHingeMoves)
    {
        sink.writeln("Hinge Moves");
        foreach (i; 0 .. hingeAccepts.length)
        {
            auto accepts = "%10,d".format(hingeAccepts[i]);
            auto tries = "%10,d".format(hingeTries[i]);
            sink.writefln("    %2s → %2-s    : %14s / %14s (%5.3f)", i + 4,
                (i + 2) * (dim - 1), accepts, tries,
                double(hingeAccepts[i]) / hingeTries[i]);
        }
    }
    auto totAccepts = hingeAccepts.sum + bistellarAccepts[].sum;
    auto totTries = hingeTries.sum + bistellarTries[].sum;
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

auto getFirstPart(T)(T time)
{
    return time.to!string.findSplit(",")[0].findSplit("and")[0]
        .replace("minutes", "mins");
}

auto chooseRandomMove(int dim, Vertex, P)(Manifold!(dim, Vertex) manifold, Vertex newVertex, P parameters)
{
    SumType!(BistellarMove!(dim, Vertex), HingeMove!(dim, Vertex)) chosenMove;

    auto done = false;
    while(!done)
    {
        auto facet = manifold.asSimplicialComplex.randomFacetOfDim(dim);
        auto center = facet.subsets.array.rndChoice.array;

        auto centerDim = center.walkLength - 1;
        auto centerDeg = manifold.degree(center);
        auto centerIsHinge = (centerDim == dim - 2);

        // Check if the chosen center simplex has an appropriate degree
        // to be a bistellar move or hinge move
        if ((parameters.useHingeMoves) && (centerDim == dim - 2))
        {
            if (centerDeg > maxHingeMoveDeg)
            {
                continue;
            }
        }
        else if (centerDeg + centerDim != dim + 1)
        {
            continue;
        }

        bool isHingeMove = centerIsHinge && (centerDeg > 3);
        bool isBistellarMove = !isHingeMove;

        if (isHingeMove)
        {
            int triangIndx = uniform(0, numNgonTriangs(centerDeg));
            auto rimVertices = manifold.orderedHingeLinkVertices(center);
            chosenMove = HingeMove!(dim, Vertex)(center, rimVertices, triangIndx);
        }

        if (isBistellarMove)
        {
            if (centerDim == dim)
            { 
                chosenMove = BistellarMove!(dim, Vertex)(center, newVertex.only);
            }
            else
            {
                auto coCenter = manifold.coCenter(center, facet);
                chosenMove = BistellarMove!(dim, Vertex)(center, coCenter);
            }
        }

        if (uniform01 > 2.0 / centerDeg)
        {
            continue;
        }

        if (!manifold.hasValidMove(chosenMove))
        {
            continue;
        }
        
        return chosenMove;
    }
    return chosenMove;
}
///
@Name("chooseRandomMove") @safe unittest
{
    // http://page.math.tu-berlin.de/~lutz/stellar/RP3
    auto rp3 = Manifold!3([[1, 2, 3, 7], [1, 2, 3, 11], [1, 2, 6, 9], [1,
            2, 6, 11], [1, 2, 7, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 7,
            10], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 8, 9], [1, 4, 8, 10], [1,
            5, 6, 8], [1, 5, 6, 11], [1, 5, 8, 10], [1, 6, 8, 9], [2, 3, 4, 8],
            [2, 3, 4, 11], [2, 3, 7, 8], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4,
            8, 10], [2, 5, 7, 8], [2, 5, 7, 9], [2, 5, 8, 10], [2, 5, 9, 10],
            [2, 6, 9, 10], [3, 4, 5, 9], [3, 4, 5, 11], [3, 4, 8, 9], [3, 5, 9,
            10], [3, 6, 7, 8], [3, 6, 7, 10], [3, 6, 8, 9], [3, 6, 9, 10], [4,
            5, 6, 7], [4, 5, 6, 11], [4, 5, 7, 9], [4, 6, 7, 10], [5, 6, 7, 8]]);

    // foreach (i; iota(10))
    // {
    //     auto mv = rp3.chooseRandomMove(716, Yes.includeHingeMoves, Yes.listAllMoves);
    //     writelnUt(mv);
    // }
}

void updateUnusedVertices(int dim, Vertex)(ref Vertex[] unusedVertices, SumType!(BistellarMove!(dim,Vertex), HingeMove!(dim,Vertex)) move)
{
    alias BM = BistellarMove!(dim, Vertex);
    move.match!(
        (BM bistellarMove)
        {
            if (bistellarMove.coCenter.length == 1)
            {
                unusedVertices.popBack;
            }
            if (bistellarMove.center.length == 1)
            {
                unusedVertices ~= bistellarMove.center;
            }                   
        },
        (_) {});
}

void undoUpdateUnusedVertices(int dim, Vertex)(ref Vertex[] unusedVertices, SumType!(BistellarMove!(dim,Vertex), HingeMove!(dim,Vertex)) move)
{
    alias BM = BistellarMove!(dim, Vertex);
    move.match!(
        (BM bistellarMove)
        {
            if (bistellarMove.coCenter.length == 1)
            {
                unusedVertices ~= bistellarMove.coCenter;
            }
            if (bistellarMove.center.length == 1)
            {
                assert(bistellarMove.center.front == unusedVertices.back);
                unusedVertices.popBack;
            }                   
        },
        (_) {});
}

void incrementMoveCounts(int dim, Vertex)(
    size_t[] bistellarTries,
    size_t[] hingeTries,
    SumType!(BistellarMove!(dim,Vertex), HingeMove!(dim,Vertex)) move)
{
    alias BM = BistellarMove!(dim, Vertex);
    alias HM = HingeMove!(dim, Vertex);

    move.match!(
        (BM bistellarMove)
        {
            ++bistellarTries[bistellarMove.coCenter.length - 1];
        },
        (HM hingeMove)
        {
            ++hingeTries[hingeMove.rim.length - 4];
        });
}

void decrementMoveCounts(int dim, Vertex)(
    size_t[] bistellarTries,
    size_t[] hingeTries,
    SumType!(BistellarMove!(dim,Vertex), HingeMove!(dim,Vertex)) move)
{
    alias BM = BistellarMove!(dim, Vertex);
    alias HM = HingeMove!(dim, Vertex);

    move.match!(
        (BM bistellarMove)
        {
            --bistellarTries[bistellarMove.coCenter.length - 1];
        },
        (HM hingeMove)
        {
            --hingeTries[hingeMove.rim.length - 4];
        });
}

void saveSample(M, S, T, P)(M manifold, ulong sampleNumber,
                S startTime, T timer, ulong[] bistellarTries, ulong[] bistellarAccepts,
                ulong[] hingeTries, ulong[] hingeAccepts, P parameters)
{
    string prefix = parameters.savedFilesPrefix ~ "_sample_"
        ~ sampleNumber.to!string;
    
    manifold.saveTo(prefix ~ ".mfd");
        
    auto saveFile = File(prefix ~ ".mfd", "a");       
    saveFile.writeln;
    saveFile.writeReports(manifold, startTime, timer,
        bistellarTries[], bistellarAccepts[], hingeTries[], hingeAccepts[], parameters);
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