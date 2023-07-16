/// Samples from manifolds (TO DO: better description)
module applications.manifold_sampler;

import algorithms : eulerCharacteristic;
import manifold;
import manifold_examples : standardSphere, octahedron;
import manifold_moves : BistellarMove, HingeMove, allHingeMoves;
import simplicial_complex : fVector;
import utility : binomial, dump, flatDegreeInDim, parseParameterFile, prettyTime;

import unit_threaded : Name, writelnUt;

import core.memory : GC;

import std.algorithm : all, count, each, filter, findSplit, joiner, map, max, maxElement, sort, sum;
import std.array : array, replace, split;
import std.conv : to;
import std.datetime.stopwatch : Duration, msecs, StopWatch;
import std.datetime.date : DateTime;
import std.datetime.systime : Clock;
import std.format : format;
import std.getopt : getopt, defaultGetoptPrinter;
import std.math : exp, sqrt, isNaN, modf;
import std.range;
import std.random : choice, uniform, uniform01;
import std.stdio : File, write, writef, writefln, writeln, stdout;
import std.string : startsWith;
import std.sumtype : SumType;
import std.typecons : Flag, Yes, No;

version (unittest) {} else {
int main(string[] args)
{
    string mfdFileName;
    string paramFileName;
    auto helpInformation = getopt(
        args,
        "initialManifoldFile", &mfdFileName,
        "parameterFile", &paramFileName);
   
    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        return 0; // exit with return code zero (success)
    }

    StopWatch timer;
    timer.start;

    static immutable parametersUsed = [
            ["bool", "disableGC"],
            ["bool", "useHingeMoves"],
            ["bool", "checkForProblems"],
            ["int", "numFacetsTarget"],
            ["int", "maxHingeMoveDeg"],
            ["int", "sweepsPerProblemCheck"],
            ["int", "triesPerStdoutReport"],
            ["int", "triesPerCollect"],
            ["int", "dim"],
            ["int", "maxCoDim2Bins"],
            ["int", "maxCoDim3Bins"],
            ["int", "maxBar2"],
            ["int", "maxBar3"],
            ["real", "dt"],
            ["real", "dtPerFileReport"],
            ["real", "dtPerSave"],
            ["real", "hingeDegreeTarget"],
            ["real", "numFacetsCoef"],
            ["real", "numHingesCoef"],
            ["real", "hingeDegreeVarCoef"],
            ["real", "cd3DegVarCoef"],
            ["real", "maxSweeps"],
            ["string", "saveFilePrefix"]
        ];
    auto params = parseParameterFile!parametersUsed(paramFileName);
    
    enum dim = 3; // TO DO: Make this selectable from parameter file
    Manifold!dim mfd;

    if (!mfdFileName.empty)
    {
        "Loading %s-manifold from file: %s".writefln(dim, mfdFileName);
        "... ".write;
        mfd = loadManifold!dim(mfdFileName);
        "done. Took %s.".format(timer.prettyTime).writeln;
        timer.reset;
    }
    else
    {
        // TO DO: User chosen options from library of triangs
        mfd = standardSphere!dim;
    }

    auto initialVertices = mfd.simplices(0).joiner.array.dup.sort.array;
    assert(initialVertices.length > 0, "initial manifold is empty");
    int[] unusedVertices = getUnusedVertices(mfd, initialVertices);

    // bistellarTries[j] counts the (j + 1) -> (dim + 1 - j) moves attempted
    ulong[dim + 1] bistellarTries;

    // bistellarAccepts[j] counts the (j + 1) -> (dim + 1 - j) moves accepted
    ulong[dim + 1] bistellarAccepts;

    // hingeTries[j] counts the hinge-moves with hinge degree j + 4 tried
    assert(params.maxHingeMoveDeg >= 4, "Hinge move degrees must be at least 4");
    ulong[] hingeTries;
    hingeTries.length = params.maxHingeMoveDeg - 3;

    // hingeAccepts[j] counts the hinge-moves with hinge degree j + 4 accepted
    ulong[] hingeAccepts;
    hingeAccepts.length = params.maxHingeMoveDeg - 3;

    ulong dtElapsed; // number of dt intervals elapsed (in sweeps)
    auto dtIncThreshold = (params.dt* params.numFacetsTarget).to!ulong;
    assert(dtIncThreshold > 0);

    typeof(timer.peek()) timePerMove;
    int sampleNumber;
    auto startTime = Clock.currTime;
    auto currentObjective = mfd.objective(params);
    timer.reset;

    if(params.disableGC)
    {
        GC.disable;
    }

    auto doneSampling = false;
    while (!doneSampling)
    {
        doneSampling = (dtElapsed * params.dt >= params.maxSweeps);
        ulong numMovesTried = bistellarTries[].sum;
        ulong numMovesAccepted = bistellarAccepts[].sum;
        if(params.useHingeMoves)
        {
            numMovesTried += hingeTries[].sum;
            numMovesAccepted += hingeAccepts[].sum;
        }

        double acceptFrac = double(numMovesAccepted) / numMovesTried;
        timePerMove = timer.peek / params.triesPerStdoutReport;

        bool dtIncreased = false;
        if (numMovesAccepted == dtIncThreshold)
        {
            ++dtElapsed;
            dtIncreased = true;
            dtIncThreshold = ((dtElapsed + 1) * params.dt* params.numFacetsTarget).to!ulong;
            assert(dtIncThreshold > numMovesAccepted);
        }

        if (unusedVertices.empty)
        {
            /* If there are no unused vertices, then the next unused vertex
            label is just the number of vertices. */
            unusedVertices ~= mfd.fVector[0].to!int;
        }
        assert(unusedVertices.all!(v => !mfd.contains(v.only)));

        auto bistellarMoves = mfd.allBistellarMoves;
        auto numBistellarMoves = bistellarMoves.walkLength;
        auto numNewVertexMoves = mfd.fVector[dim];
        auto numberOfMoves = numBistellarMoves + numNewVertexMoves;

        HingeMove!dim[] hingeMoves;
        size_t numHingeMoves = 0;
        if(params.useHingeMoves)
        {
            hingeMoves = mfd.allHingeMoves;
            numHingeMoves = hingeMoves.length;
            numberOfMoves += numHingeMoves;
        }

        auto indxOfChosenMove = uniform(0,numberOfMoves);

        real deltaObjective = 0;
        if(indxOfChosenMove < numNewVertexMoves)
        {
            // Chosen move is a 1->(dim+1) bistellar move
            auto center = mfd.randomFacetOfDim(dim);
            auto coCenter = unusedVertices.back.only;
            auto chosenMove = BistellarMove!dim(center, coCenter);

            deltaObjective = changeInObjective(mfd, chosenMove, params);
            writeln("chosen move: ", chosenMove);
        }
        else if(indxOfChosenMove < numNewVertexMoves + numBistellarMoves)
        {
            // Chosen move is a bistellar move that isn't 1->(dim+1)
            auto indx = indxOfChosenMove - numNewVertexMoves;
            auto chosenMove = bistellarMoves[indx];
            deltaObjective = changeInObjective(mfd, chosenMove, params);
            writeln("chosen move: ", chosenMove);
        }
        else
        {
            assert(params.useHingeMoves);
            auto indx = indxOfChosenMove - numNewVertexMoves - numBistellarMoves;
            auto chosenMove = hingeMoves[indx];
            deltaObjective = changeInObjective(mfd, chosenMove, params);
            writeln("chosen move: ", chosenMove);
        }
  
        doneSampling = true;

        // auto chosenMove = mfd.getRandomMove;



        // if (chosenMove.isPachner)
        // {
        //     mfd.doPachner(chosenMove.center, chosenMove.coCenter);
        //     ++bistellarTries[dim + 1 - chosenMove.center.length];
        //     // Also increment accepted #, we will fix if move rejected
        //     ++bistellarAccepts[dim + 1 - chosenMove.center.length];

        //     if (chosenMove.center.length == 1)
        //     {
        //         unusedVertices ~= chosenMove.center.front;
        //     }
        //     else if (chosenMove.center.length == dim + 1)
        //     {
        //         unusedVertices.popBack;
        //     }
        // }
        // else
        // {
        //     assert(useHingeMoves);
        //     assert(chosenMove.isHinge);
        //     mfd.doHingeMove(chosenMove.center, chosenMove.coCenter, chosenMove.triangIndx);
        //     ++hingeTries[chosenMove.coCenter.length - 4];
        //     // Also increment accepted #, we will fix if move rejected
        //     ++hingeAccepts[chosenMove.coCenter.length - 4];
        // }

        // real newObjective = mfd.objective;
        // real deltaObj = newObjective - currentObjective;
        // bool rejectMove = false;
 
        // // REJECT and UNDO move, if appropriate
        // if ((deltaObj > 0) && (uniform01 > exp(-deltaObj)))
        // {
        //     if (chosenMove.isPachner)
        //     {
        //         mfd.doPachner(chosenMove.coCenter, chosenMove.center);
        //         --bistellarAccepts[dim + 1 - chosenMove.center.length];
        //         if (chosenMove.center.length == 1)
        //         {
        //             unusedVertices.popBack;
        //         }
        //         else if (chosenMove.center.length == dim + 1)
        //         {
        //             unusedVertices ~= chosenMove.coCenter.front;
        //         }
        //     }
        //     else
        //     {
        //         --hingeAccepts[dim + 1 - chosenMove.center.length];
        //         mfd.undoHingeMove(chosenMove.center, chosenMove.coCenter, chosenMove.triangIndx);
        //     }
        // }

        //--------------------------- MAKE REPORT ----------------------------
        if ((numMovesTried % params.triesPerStdoutReport == 0) || doneSampling)
        {
            stdout.writeTimingAndTargetsReport(mfd, dtElapsed, startTime, timePerMove, acceptFrac, params);
            stdout.writeSimplexReport(mfd);
            stdout.writeObjectiveReport(mfd, params);
            stdout.writeMoveReport(bistellarTries, bistellarAccepts, hingeTries, hingeAccepts, params);
            stdout.writeHistogramReport(mfd, params);
        }

        //----------------------- WRITE TO DATA FILE --------------------------
        if ((dtIncreased && (dtElapsed % params.dtPerFileReport == 0)) || doneSampling)
        {
            auto file = File(params.saveFilePrefix ~ ".dat", "a");
            // TO DO: Implement this
            //columnReport(file);
        }

        //----------------------- SAVE CURRENT MANIFOLD ----------------------
        if ((dtIncreased && (dtElapsed % params.dtPerSave == 0)) || doneSampling)
        {
            string prefix;
            if (doneSampling)
            {
                prefix = params.saveFilePrefix ~ "_final";
            }
            else
            {
                prefix = params.saveFilePrefix ~ "_save"
                    ~ sampleNumber.to!string;
            }

            auto savedMfdFileName = prefix ~ ".mfd";
            auto graphFileName = prefix ~ ".edge_graph";

            mfd.saveTo(savedMfdFileName);
            auto saveFile = File(savedMfdFileName, "a");
            saveFile.writeln;

            saveFile.writeTimingAndTargetsReport(mfd, dtElapsed, startTime, timePerMove, acceptFrac, params);
            saveFile.writeSimplexReport(mfd);
            saveFile.writeObjectiveReport(mfd, params);
            saveFile.writeMoveReport(bistellarTries, bistellarAccepts, hingeTries, hingeAccepts, params);
            saveFile.writeHistogramReport(mfd, params);

            mfd.saveEdgeGraphTo(graphFileName);
            ++sampleNumber;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if(params.disableGC)
        {
            if (numMovesTried % params.triesPerCollect == 0)
            {
                GC.enable;
                GC.collect;
                GC.disable;
            }
        }

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        if(params.checkForProblems)
        {
            if (dtIncreased
                && ((dtElapsed * params.dt) % params.sweepsPerProblemCheck == 0))
            {
                "checking for problems ... ".write;
                auto problems = mfd.findProblems;
                if (!problems.empty)
                {
                    problems.each!writeln;
                    assert(0);
                }
                "done".writeln;
            }
        }
    }

    "Finished! Time elapsed: %s".writefln(Clock.currTime.to!DateTime - startTime.to!DateTime);
    return 0; // exit with return code zero (success)
}}

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
    immutable nCoDim3 = mfd.fVector[dim - 3];
    immutable totSqDeg = mfd.totalSquareDegree(dim - 2);
    immutable totSAsqDeg = mfd.totalSquareDegree(dim - 3);

    immutable nHingesTarget = hingesPerFacet * nFacets / params.hingeDegreeTarget;
    immutable degTarget = hingesPerFacet * nFacets / nHinges.to!real;
    immutable coDim3DegTarget = coDim3PerFacet * nFacets / nCoDim3.to!real;

    Penalty penalty;
    penalty.volumePenalty = (nFacets - params.numFacetsTarget) ^^ 2;
    penalty.globalCurvPenalty = (nHinges - nHingesTarget) ^^ 2;

    real _; // dummy for integer part
    real x = modf(degTarget, _); // fractional part
    real minPenalty = (x - x ^^ 2) * nHinges;

    // TO DO: Refer to external paper for this!
    penalty.localCurvPenalty = (
        degTarget ^^ 2 * nHinges - 2 * degTarget * hingesPerFacet * nFacets + totSqDeg) - minPenalty;

    x = modf(coDim3DegTarget, _); // fractional part
    minPenalty = (x - x ^^ 2) * nCoDim3;

    // TO DO: Refer to arxiv paper for this!
    penalty.localSolidAngleCurvPenalty = (
        coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + totSAsqDeg) - minPenalty;

    return penalty;
}

real objective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    auto pen = mfd.penalties(params);
    return params.numFacetsCoef * pen.volumePenalty
        + params.numHingesCoef * pen.globalCurvPenalty
        + params.hingeDegreeVarCoef * pen.localCurvPenalty
        + params.cd3DegVarCoef * pen.localSolidAngleCurvPenalty;
}

void writeTimingAndTargetsReport(M, S, T, W, P)(W sink, M mfd, ulong dtElapsed, S startTime, T timePerMove, double acceptFrac, P params)
{
    typeof(timePerMove) timePerSweep;
    if (acceptFrac > 0.0)
    {
        timePerSweep = timePerMove * (mfd.fVector[params.dim] / acceptFrac).to!ulong;
    }

    "-".repeat(80).joiner.writeln;
    string timeInfo = "%.1f / %s".format(dtElapsed * params.dt, params.maxSweeps);
    auto prettyStartTime = startTime.to!string.split('.').front;

    "sweeps         : %24-s|".writefln(timeInfo);
    "started at     : %24-s".writef(prettyStartTime);
    if (params.numFacetsCoef > 0.0)
    {
        writef("| # facets target     : %s", params.numFacetsTarget);
    }
    writeln;

    "end (estimate) : %23-s ".writef("*** TO DO ***");
    if (params.numHingesCoef > 0.0)
    {
        writef("| hinge degree target : %.5f", params.hingeDegreeTarget);
    }
    writeln;

    "Δt/move        : %23-s |".writefln(timePerMove.to!string);
    "Δt/sweep       : %23-s |".writefln(timePerMove.getFirstPart);
    "move accept %%  : %23.3-f |".writefln(acceptFrac);
}

void writeSimplexReport(int dim, Vertex, W)(W sink, const ref Manifold!(dim, Vertex) mfd)
{
    "-".repeat(80).joiner.writeln;
    "dimension      : ".write;
    (dim + 1).iota.each!(d => sink.writef("%12s ", d));
    sink.writeln;
    "# simplices    : ".write;
    (dim + 1).iota.each!(d => sink.writef("%12s ", mfd.fVector[d]));
    sink.writeln;
    "degree mean    : ".write;
    (dim + 1).iota.each!(d => sink.writef("%12.4f ", mfd.meanDegree(d)));
    sink.writeln;
    "degree std dev : ".write;
    (dim + 1).iota.each!(d => sink.writef("%12.4f ", mfd.degreeVariance(d).sqrt));
    sink.writeln;
}

void writeObjectiveReport(int dim, Vertex, W, P)(W sink, const ref Manifold!(dim, Vertex) mfd, P params)
{
    "-".repeat(80).joiner.writeln;
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

    if (params.hingeDegreeVarCoef > 0.0)
    {
        auto lcp = mfd.penalties(params).localCurvPenalty;
        sink.writefln("hinge deg var  : %.6e  *  %.4f  =  %.6e",
            lcp, params.hingeDegreeVarCoef, params.hingeDegreeVarCoef * lcp);
    }

    if (params.cd3DegVarCoef > 0.0)
    {
        auto lsacp = mfd.penalties(params).localSolidAngleCurvPenalty;
        sink.writefln("codim-3 deg var: %.6e  *  %.4f  =  %.6e",
            lsacp, params.cd3DegVarCoef, params.cd3DegVarCoef * lsacp);
    }
    sink.writefln("total penalty  :                             %.6e",
        mfd.objective(params));
}

void writeMoveReport(S, T, W, P)(W sink, S bistellarTries, S bistellarAccepts, T hingeTries, T hingeAccepts, P params)
{
    sink.writeln("-".repeat(80).joiner);
    sink.writeln("Bistellar Moves      # Accepted          # Tried    %    ");
    foreach (i; 0 .. params.dim + 1)
    {
        auto accepts = "%10,d".format(bistellarAccepts[i]);
        auto tries = "%10,d".format(bistellarTries[i]);
        sink.writefln("    %2s → %2-s    : %14s / %14s (%5.3f)",
            i + 1, params.dim + 1 - i, accepts, tries,
            double(bistellarAccepts[i]) / bistellarTries[i]);
    }
    if(params.useHingeMoves)
    {
        sink.writeln("Hinge Moves");
        foreach (i; 0 .. hingeAccepts.length)
        {
            auto accepts = "%10,d".format(hingeAccepts[i]);
            auto tries = "%10,d".format(hingeTries[i]);
            sink.writefln("    %2s → %2-s    : %14s / %14s (%5.3f)", i + 4,
                (i + 2) * (params.dim - 1), accepts, tries,
                double(hingeAccepts[i]) / hingeTries[i]);
        }
    }
}

void writeHistogramReport(int dim, Vertex, W, P)(W sink, const ref Manifold!(dim, Vertex) mfd, P params)
{
    "-".repeat(80).joiner.writeln;
    auto maxDeg2 = 2 + params.maxCoDim2Bins;
    auto maxDeg3 = 2 + 2 * params.maxCoDim3Bins;

    auto hist2 = mfd.degreeHistogram(mfd.dimension - 2);
    auto tail2 = (hist2.length >= maxDeg2) ? hist2[maxDeg2 .. $].sum : 0;
    auto maxDegBin2 = max(hist2.maxElement, tail2);
    auto normedHist2 = hist2.map!(freq => real(freq) / maxDegBin2);

    auto hist3 = mfd.degreeHistogram(mfd.dimension - 3);
    auto tail3 = (hist3.length >= maxDeg3) ? hist3[maxDeg3 .. $].sum : 0;
    auto maxDegBin3 = max(hist3.maxElement, tail3);
    auto normedHist3 = hist3.map!(freq => real(freq) / maxDegBin3);

    writeln(' '.repeat(27), "Codimension-2 Degree Histogram");

    static immutable bars = [
        '▏', '▎', '▍', '▌', '▋', '▊', '▉', '█'
    ];

    foreach (bin; iota(0, params.maxCoDim2Bins))
    {
        writef("%2s ", bin + 3);
        if (bin + 2 < normedHist2.length)
        {
            auto nEighths = (params.maxBar2 * normedHist2[bin + 2] * 8.0).to!int;
            auto bar = bars.back.repeat(nEighths / 8).array;
            if (nEighths % 8 > 0)
            {
                bar ~= bars[nEighths % 8];
            }
            writef("%-*s", params.maxBar2, bar);
        }
        writeln;
    }
    auto tailFreq2 = real(tail2) / maxDegBin2;
    if (tailFreq2 > 0)
    {
        auto nEighths2 = (params.maxBar2 * tailFreq2 * 8.0).to!int;
        auto bar2 = bars.back.repeat(nEighths2 / 8).array;
        if (nEighths2 % 8 > 0)
        {
            bar2 ~= bars[nEighths2 % 8];
        }
        writefln(" > %-*s", params.maxBar2, bar2);
    }

    static if (dim > 2)
    {
        "-".repeat(80).joiner.writeln;
        writeln(' '.repeat(27), "Codimension-3 Degree Histogram");
        foreach (bin; iota(0, params.maxCoDim3Bins))
        {
            writef("%2s ", 4 + bin * 2);
            if ((3 + bin * 2) < normedHist3.length)
            {
                auto nEighths = (params.maxBar3 * normedHist3[3 + bin * 2] * 8.0).to!int;
                auto bar = bars.back.repeat(nEighths / 8).array;
                if (nEighths % 8 > 0)
                {
                    bar ~= bars[nEighths % 8];
                }
                writef("%-*s", params.maxBar3, bar);
            }
            writeln;
        }
        auto tailFreq3 = real(tail3) / maxDegBin3;
        if (tailFreq3 > 0)
        {
            write(" > ");
            auto nEighths3 = (params.maxBar3 * tailFreq3 * 8.0).to!int;
            auto bar3 = bars.back.repeat(nEighths3 / 8).array;
            if (nEighths3 % 8 > 0)
            {
                bar3 ~= bars[nEighths3 % 8];
            }
            writefln("%-*s", params.maxBar3, bar3);
        }
    }
}

auto getFirstPart(T)(T time)
{
    return time.to!string.findSplit(",")[0].findSplit("and")[0]
        .replace("minutes", "mins");
}

real changeInObjective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, BistellarMove!dim move, P params)
{
    auto currentObjective = mfd.objective(params);
    writeln("current objective: ", currentObjective);
    
    return 0;
}

real changeInObjective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, HingeMove!dim move, P params)
{
    auto currentObjective = mfd.objective(params);
    writeln("current objective: ", currentObjective);
    
    return 0;
}

auto chooseRandomMove(int dim, Vertex)(Manifold!(dim, Vertex) manifold, Vertex newVertex, Flag!"includeHingeMoves" includeHingeMoves = Yes.includeHingeMoves, Flag!"listAllMoves" listAllMoves = No.listAllMoves)
{
    SumType!(BistellarMove!(dim, Vertex), HingeMove!(dim, Vertex)) chosenMove;




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

    foreach(i; iota(10))
    {
        auto mv = rp3.chooseRandomMove(716, Yes.includeHingeMoves, Yes.listAllMoves);
        writelnUt(mv);
    }
}
