module applications.manifold_sampler;

import algorithms : eulerCharacteristic;
import manifold;
import manifold_examples : standardSphere;
import moves : Move;
import simplicial_complex : fVector;
import utility : binomial, flatDegreeInDim, subsetsOfSize, subsets, toStackArray, StackArray, prettyTime;

import std.algorithm;
import std.array : split, array, replace, staticArray;
import std.conv : to;
import std.datetime.stopwatch : Duration, msecs, StopWatch;
import std.datetime.date : DateTime;
import std.datetime.systime : Clock;
import std.format : format;
import std.getopt : getopt, defaultGetoptPrinter;
import std.math : exp, sqrt, isNaN, modf;
import std.range;
import std.random : choice, dice, rndGen, uniform, uniform01;
import std.stdio : File, write, writef, writefln, writeln, stdout;
import std.traits : isFloatingPoint;
import std.variant : Algebraic, visit, tryVisit;
import std.uuid : randomUUID;

import core.memory : GC;

import unit_threaded : Name;

mixin(import("manifold_sampler.config").parseConfig);

version (unittest) {} else {
void main(string[] args)
{
    
    string mfdFile;
    auto helpInformation = getopt(args, "manifoldFile", &mfdFile);
    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
    }

    auto filePrefix = mfdFile.findSplit(".mfd")[0];
    if (filePrefix.empty)
    {
        auto runID = randomUUID();
        filePrefix = runID.to!string[0 .. 8];
    }

    StopWatch timer;
    timer.start;

    Manifold!dim mfd;

    if (!mfdFile.empty)
    {
        "Loading %s-manifold from file: %s".writefln(dim, mfdFile);
        "... ".write;
        mfd = loadManifold!dim(mfdFile);
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
    int[] unusedVertices = mfd.getUnusedVertices(initialVertices);
    
    // bistellarTries[j] counts the (j + 1) -> (dim + 1 - j) moves attempted
    ulong[dim + 1] bistellarTries;

    // bistellarAccepts[j] counts the (j + 1) -> (dim + 1 - j) moves accepted
    ulong[dim + 1] bistellarAccepts;

    // hingeTries[j] counts the hinge-moves with hinge degree j + 4 tried
    static assert(maxHingeMoveDeg >= 4, "Hinge move degrees must be at least 4");
    ulong[maxHingeMoveDeg - 3] hingeTries;

    // hingeAccepts[j] counts the hinge-moves with hinge degree j + 4 accepted
    ulong[maxHingeMoveDeg - 3] hingeAccepts;

    ulong dtElapsed; // number of dt intervals elapsed (in sweeps)

    auto numMovesTried = bistellarTries[].sum + hingeTries[].sum;
    auto numMovesAccepted = bistellarAccepts[].sum + hingeAccepts[].sum;
    double acceptFrac = double(numMovesAccepted) / numMovesTried;

    auto timePerMove = timer.peek / triesPerStdoutReport;
    typeof(timePerMove) timePerSweep;
    if (acceptFrac > 0.0)
    {
        timePerSweep = timePerMove * (mfd.fVector.back / acceptFrac).to!ulong;
    }

    auto startTime = Clock.currTime;
    timer.reset;

    // **************************** PRINT REPORTS ****************************
    timingAndTargetsReport(dtElapsed, startTime, timePerMove, acceptFrac);
    mfd.simplexReport;
    mfd.objectiveReport;
    moveReport(bistellarTries, bistellarAccepts, hingeTries, hingeAccepts);
    mfd.histogramReport;
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

string parseConfig()(string str)
{
    string ret;
    foreach (line; str.split("\n"))
    {
        if (line.walkLength > 0 && line[0] != '#')
        {
            auto parts = line.split(" = ");
            ret ~= `enum ` ~ parts[0] ~ ` = ` ~ parts[1] ~ `;`;
        }
    }
    return ret;
}

private struct Penalty
{
    real volumePenalty;
    real globalCurvPenalty;
    real localCurvPenalty;
    real localSolidAngleCurvPenalty;
}

Penalty penalties(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    immutable hingesPerFacet = dim * (dim + 1) / 2;
    immutable coDim3PerFacet = binomial(dim + 1, dim - 2);

    immutable nFacets = mfd.fVector[dim];
    immutable nHinges = mfd.fVector[dim - 2];
    immutable nCoDim3 = mfd.fVector[dim - 3];
    immutable totSqDeg = mfd.totalSquareDegree(dim - 2);
    immutable totSAsqDeg = mfd.totalSquareDegree(dim - 3);

    immutable nHingesTarget = hingesPerFacet * nFacets / hingeDegreeTarget;
    immutable degTarget = hingesPerFacet * nFacets / nHinges.to!real;
    immutable coDim3DegTarget = coDim3PerFacet * nFacets / nCoDim3.to!real;

    Penalty penalty;
    penalty.volumePenalty = (nFacets - numFacetsTarget) ^^ 2;
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

real objective(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    auto pen = mfd.penalties;
    return pen.volumePenalty;
}

void timingAndTargetsReport(S,T,U)(ulong dtElapsed, S startTime, T timePerMove, U acceptFrac)
{
    "-".repeat(80).joiner.writeln;
    string timeInfo = "%.1f / %s".format(dtElapsed * dt, maxSweeps);
    auto prettyStartTime = startTime.to!string.split('.').front;

    "sweeps         : %24-s".writef(timeInfo);
    "| tracking # moves    : %s".writefln(trackNumMoves);
  
    "started at     : %24-s".writef(prettyStartTime);
    if (numFacetsCoef > 0.0)
    {
        writef("| # facets target     : %s", numFacetsTarget);
    }
    writeln;

    "end (estimate) : %23-s ".writef("*** TO DO ***");
    if (numHingesCoef > 0.0)
    {
        writef("| hinge degree target : %.5f", hingeDegreeTarget);
    }
    writeln;

    "Δt/move        : %23-s |".writefln(timePerMove.to!string);
    "Δt/sweep       : %23-s |".writefln(timePerMove.getFirstPart);
    "move accept %%  : %23.3-f |".writefln(acceptFrac);
}

void simplexReport(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    "-".repeat(80).joiner.writeln;
    "dimension      : ". write;
    (dim + 1).iota.each!(d => "%12s ".writef(d));
    writeln;
    "# simplices    : ". write;
    (dim + 1).iota.each!(d => "%12s ".writef(mfd.fVector[d]));
    writeln;
    "degree mean    : ". write;
    (dim + 1).iota.each!(d => "%12.4f ".writef(mfd.meanDegree(d)));
    writeln;
    "degree std dev : ". write;
    (dim + 1).iota.each!(d => "%12.4f ".writef(mfd.degreeVariance(d).sqrt));
    writeln;
}

void objectiveReport(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    "-".repeat(80).joiner.writeln;
    writeln("    Penalty", ' '.repeat(10), "Raw",
        ' '.repeat(11), "Coef", ' '.repeat(9), "Value");

    if (numFacetsCoef > 0.0)
    {
        auto vp = mfd.penalties.volumePenalty;
        writefln("# facets       : %.6e  *  %.4f  =  %.6e",
           vp, numFacetsCoef, numFacetsCoef * vp);
    }

    if (numHingesCoef > 0.0)
    {
        auto gcp = mfd.penalties.globalCurvPenalty;
        writefln("# hinges       : %.6e  *  %.4f  =  %.6e",
            gcp, numHingesCoef, numHingesCoef * gcp);
    }

    if (hingeDegreeVarCoef > 0.0)
    {
        auto lcp = mfd.penalties.localCurvPenalty;
        writefln("hinge deg var  : %.6e  *  %.4f  =  %.6e",
            lcp, hingeDegreeVarCoef, hingeDegreeVarCoef * lcp);
    }

    if (cd3DegVarCoef > 0.0)
    {
        auto lsacp = mfd.penalties.localSolidAngleCurvPenalty;
        writefln("codim-3 deg var: %.6e  *  %.4f  =  %.6e",
            lsacp, cd3DegVarCoef, cd3DegVarCoef * lsacp);
    }
    writefln("total penalty  :                             %.6e",
        mfd.objective);
}

void moveReport(S, T)(S bistellarTries, S bistellarAccepts, T hingeTries, T hingeAccepts)
{
   "-".repeat(80).joiner.writeln;
    writeln("Bistellar Moves      # Accepted          # Tried    %    ");
    foreach (i; 0 .. dim + 1)
    {
        auto accepts = "%10,d".format(bistellarAccepts[i]);
        auto tries = "%10,d".format(bistellarTries[i]);
        writefln("    %2s → %2-s    : %14s / %14s (%5.3f)",
            i + 1, dim + 1 - i, accepts, tries,
            double(bistellarAccepts[i]) / bistellarTries[i]);
    }
    if(useHingeMoves)
    {
        writeln("Hinge Moves");
        foreach(i; 0 .. hingeAccepts.length)
        {
            auto accepts = "%10,d".format(hingeAccepts[i]);
            auto tries = "%10,d".format(hingeTries[i]);
            writefln("    %2s → %2-s    : %14s / %14s (%5.3f)", i + 4,
                (i + 2) * (dim - 1), accepts, tries,
                double(hingeAccepts[i]) / hingeTries[i]);
        }
    }
}

void histogramReport(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    "-".repeat(80).joiner.writeln;
    auto maxDeg2 = 2 + maxCoDim2Bins;
    auto maxDeg3 = 2 + 2 * maxCoDim3Bins;

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


    foreach (bin; iota(0, maxCoDim2Bins))
    {
        writef("%2s ", bin + 3);
        if (bin + 2 < normedHist2.length)
        {
            auto nEighths = (maxBar2 * normedHist2[bin + 2] * 8.0).to!int;
            auto bar = bars.back.repeat(nEighths / 8).array;
            if (nEighths % 8 > 0)
            {
                bar ~= bars[nEighths % 8];
            }
            writef("%-*s", maxBar2, bar);
        }
        writeln;
    }
    auto tailFreq2 = real(tail2) / maxDegBin2;
    if (tailFreq2 > 0)
    {
        auto nEighths2 = (maxBar2 * tailFreq2 * 8.0).to!int;
        auto bar2 = bars.back.repeat(nEighths2 / 8).array;
        if (nEighths2 % 8 > 0)
        {
            bar2 ~= bars[nEighths2 % 8];
        }
        writefln(" > %-*s", maxBar2, bar2);
    }

    static if(dim > 2)
    {
        "-".repeat(80).joiner.writeln;
        writeln(' '.repeat(27), "Codimension-3 Degree Histogram");
        foreach (bin; iota(0, maxCoDim3Bins))
        {
            writef("%2s ", 4 + bin * 2);
            if ((3 + bin * 2) < normedHist3.length)
            {
                auto nEighths = (maxBar3 * normedHist3[3 + bin * 2] * 8.0).to!int;
                auto bar = bars.back.repeat(nEighths / 8).array;
                if (nEighths % 8 > 0)
                {
                    bar ~= bars[nEighths % 8];
                }
                writef("%-*s", maxBar3, bar);
            }
            writeln;
        }
        auto tailFreq3 = real(tail3) / maxDegBin3;
        if (tailFreq3 > 0)
        {
            write(" > ");
            auto nEighths3 = (maxBar3 * tailFreq3 * 8.0).to!int;
            auto bar3 = bars.back.repeat(nEighths3 / 8).array;
            if (nEighths3 % 8 > 0)
            {
                bar3 ~= bars[nEighths3 % 8];
            }
            writefln("%-*s", maxBar3, bar3);
        }
    }
}
// struct Sampler(Vertex, int dim)
// {
// private:



//     StopWatch timer;
//     DateTime startTime;
// public:

//     void report(Writer)(auto ref Writer w)
//     {






//     }

//     void columnReport(Writer)(auto ref Writer w, string separator = ",")
//     {
//         static wroteColumnLabels = false;

//         auto degHists = (dim - 1).iota.map!(d => manifold.degreeHistogram(d));

//         if (!wroteColumnLabels)
//         {
//             w.write("num_moves_accepted", separator);
//             w.write((dim + 1).iota.map!(d => "num_simps_dim%s".format(d)).joiner(separator));
//             w.write(separator);
//             w.write((dim - 1).iota.map!(d => "tot_sqr_deg_dim%s".format(d)).joiner(separator));
//             w.write(separator);
//             w.write((dim - 1).iota.map!(d => "mean_deg_dim%s".format(d)).joiner(separator));
//             w.write(separator);
//             w.write((dim - 1).iota.map!(d => "stddev_deg_dim%s".format(d)).joiner(separator));
//             w.write(separator);
//             w.writeln((dim - 1).iota.map!(d => "deg_histogram_dim%s".format(d))
//                     .joiner(
//                         separator ~ "XXX" ~ separator));
//             wroteColumnLabels = true;
//         }

//         w.write(bistellarAccepts[].sum + hingeAccepts[].sum, separator);
//         w.write(manifold.fVector.map!(to!string).joiner(separator));
//         w.write(separator);
//         w.write((dim - 1).iota.map!(d =>
//                 manifold.totalSquareDegree(d)
//                 .to!string).joiner(separator));
//         w.write(separator);
//         w.write((dim - 1).iota.map!(d =>
//                 manifold.meanDegree(d).to!string).joiner(separator));
//         w.write(separator);
//         w.write((dim - 1).iota.map!(d =>
//                 manifold.degreeVariance(d)
//                 .sqrt.to!string).joiner(separator));
//         foreach (d; 0 .. dim - 1)
//         {
//             w.write(separator);
//             w.write(degHists[d].map!(to!string).joiner(separator));
//             if (d < dim - 2)
//             {
//                 w.write(separator, "XXX");
//             }
//         }
//         w.writeln;
//     }
// }

Objective objectiveParts(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    auto volPen = numFacetsCoef * volumePenalty;
    auto gcPen = numHingesCoef * globalCurvaturePenalty;
    auto lcPen = hingeDegreeVarCoef * localCurvaturePenalty;
    auto lsacPen = cd3DegVarCoef * localSACurvaturePenalty;

    return Objective(volPen, gcPen, lcPen, lsacPen);
}



void sample(Vertex, int dim)(ref Sampler!(Vertex, dim) s)
{
    if (disableGC)
    {
        GC.disable;
    }
    startTime = Clock.currTime.to!DateTime;
    timer.start;

    int sampleNumber;
    auto dtIncThreshold = (dt * numFacetsTarget).to!ulong;
    assert(dtIncThreshold > 0);
    auto doneSampling = false;
    while (!doneSampling)
    {
        if (dtElapsed * dt >= maxSweeps)
        {
            doneSampling = true;
        }

        size_t numMovesTried = bistellarTries[].sum + hingeTries[].sum;
        size_t numMovesAccepted = bistellarAccepts[].sum + hingeAccepts[].sum;

        bool dtIncreased = false;
        if (numMovesAccepted == dtIncThreshold)
        {
            ++s.dtElapsed;
            dtIncreased = true;
            dtIncThreshold = ((dtElapsed + 1)
                    * dt * numFacetsTarget).to!ulong;
            assert(dtIncThreshold > numMovesAccepted);
        }

        assert(unusedVertices.all!(v => !s.manifold.contains(v.only)));
        if (unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            unusedVertices ~= manifold_.fVector[0].to!int;
        }

        real oldObj = objective;
        long oldNumMoves = manifold.numValidMoves;

        manifold_.Move chosenMove;

        real totMoves = manifold.numValidMoves + manifold.fVector[$ - 1];
        if (uniform(0, totMoves) < manifold.numValidMoves)
        {
            chosenMove = manifold.moves.choice;
            while (manifold.contains(chosenMove.coCenter))
            {
                chosenMove = manifold.moves.choice;
            }
        }
        else
        {
            auto facet = manifold.randomFacetOfDim(dim);
            chosenMove = manifold_.Move(facet, [s.unusedVertices[$ - 1]]);
        }

        manifold_.doPachner(chosenMove.center, chosenMove.coCenter);

        ++s.bistellarTries[dim + 1 - chosenMove.center.length];
        ++s.bistellarAccepts[dim + 1 - chosenMove.center.length];

        if (chosenMove.center.length == 1)
        {
            unusedVertices ~= chosenMove.center.front;
        }

        if (chosenMove.center.length == dim + 1)
        {
            unusedVertices.popBack;
        }

        long newNumMoves = manifold.numValidMoves;
        real nMovesFactor = min(1.0, real(oldNumMoves) / newNumMoves);
        real deltaObj = objective - oldObj;

        // REJECT MOVE if appropriate
        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj) * nMovesFactor))
        {
            --s.bistellarAccepts[dim + 1 - chosenMove.center.length];

            chosenMove = manifold.Move(chosenMove.coCenter, chosenMove.center);

            manifold_.doPachner(chosenMove.center, chosenMove.coCenter);
            if (chosenMove.center.length == 1)
            {
                unusedVertices ~= chosenMove.center.front;
            }

            if (chosenMove.center.length == dim + 1)
            {
                unusedVertices.popBack;
            }
        }

        //--------------------------- MAKE REPORT ----------------------------
        if ((numMovesTried % triesPerStdoutReport == 0)
            || doneSampling)
        {
            report(stdout);
        }

        //----------------------- WRITE TO DATA FILE --------------------------
        if ((dtIncreased && (dtElapsed % dtPerFileReport == 0))
            || doneSampling)
        {
            auto file = File(saveFilePrefix ~ ".dat", "a");
            columnReport(file);
        }

        //----------------------- SAVE CURRENT MANIFOLD ----------------------
        if ((dtIncreased && (dtElapsed % dtPerSave == 0))
            || doneSampling)
        {
            string prefix;
            if (doneSampling)
            {
                prefix = saveFilePrefix ~ "_final";
            }
            else
            {
                prefix = saveFilePrefix ~ "_save"
                    ~ sampleNumber.to!string;
            }

            auto mfdFileName = prefix ~ ".mfd";
            auto graphFileName = prefix ~ ".edge_graph";

            manifold.saveTo(mfdFileName);
            auto saveFile = File(mfdFileName, "a");
            saveFile.writeln;
            report(saveFile);

            manifold.saveEdgeGraphTo(graphFileName);
            ++sampleNumber;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if (disableGC && (numMovesTried % triesPerCollect == 0))
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        static if (checkForProblems)
        {
            if (dtIncreased
                && ((dtElapsed * dt) % sweepsPerProblemCheck == 0))
            {
                "checking for problems ... ".write;
                auto problems = manifold.findProblems;
                if (!problems.empty)
                {
                    problems.each!writeln;
                    assert(0);
                }
                "done".writeln;
            }
        }
    }

    "Finished! Time elapsed: %s".writefln(
        Clock.currTime.to!DateTime - startTime);
}

auto getFirstPart(T)(T time)
{
    return time.to!string.findSplit(",")[0].findSplit("and")[0]
        .replace("minutes", "mins");
}
