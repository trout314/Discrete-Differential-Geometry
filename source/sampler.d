module sampler;

import algorithms : eulerCharacteristic;
import core.memory : GC;
import manifold : coCenter, degreeHistogram, degreeVariance, doHingeMove, findCoCenter, findProblems, getRandomHingeMove, hasValidHingeMove, 
    linkVerticesAtHinge, movesAtFacet, Manifold, saveEdgeGraphTo, saveTo, totalSquareDegree, doPachner, meanDegree, undoHingeMove;
import simplicial_complex : fVector;
import std.algorithm : all, each, filter, findSplit, joiner, map, max, mean, min, maxElement, sort, sum;
import std.array : replace;
import std.conv : to;
import std.datetime : Duration, msecs;
import std.datetime.date : DateTime;
import std.datetime.stopwatch : StopWatch;
import std.datetime.systime : Clock;
import std.format : format;
import std.math : exp, sqrt, isNaN, modf;
import std.random : choice, dice, rndGen, uniform, uniform01;
import std.range : array, back, chain, empty, enumerate, front, iota, only, popBack, popFront, repeat, retro,
    save, walkLength;
import std.stdio : File, write, writef, writefln, writeln, stdout;
import std.traits : isFloatingPoint;
import std.variant : Algebraic, visit, tryVisit;
import unit_threaded : Name;
import utility : binomial, flatDegreeInDim, subsetsOfSize, subsets, toStackArray, StackArray;

import manifold_examples : standardSphere;

import moves : Move;

import std.array : staticArray;

/// struct containing the sampling parameters
struct Parameters
{
    // TO DO: getter and setters to check things...

    int numFacetsTarget;
    real hingeDegreeTarget = flatDegreeInDim[3];

    real numFacetsCoef = 0.01;
    real numHingesCoef = 0.0;
    real hingeDegreeVarCoef = 0.0;
    real cd3DegVarCoef = 0.0;
    real maxSweeps = 0.0;   // 1 sweep = 1 accepted move / facet

    // Time increment (in units of sweeps) used for finer-graned intervals
    real dt = 0.1;

    int dtPerFileReport = int.max;
    int dtPerSave = int.max;

    string saveFilePrefix;

    bool disableGC = true;
    int triesPerCollect = 500;
    int triesPerStdoutReport = 10;

    enum checkForProblems = false; 
    int sweepsPerProblemCheck = int.max;

    bool useHingeMoves = true;

    int dtPerHistory = 20;
    enum historyLength = 20;
}

private struct Objective
{
    real volPen;
    real gcPen;
    real lcPen;
    real lsacPen;
}

struct Sampler(Vertex, int dim)
{
private:
    Parameters params;

    Manifold!(dim, Vertex) manifold_;
    const(Vertex)[] unusedVertices;

    // bistellarTries[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] bistellarTries;
    ulong[dim + 1] bistellarAccepts;

    // hingeTries[j] counts the hinge moves with hinge degree j + 4;
    // TO DO: remove magic constant here 3 + 4 = 7 = max hinge degree
    ulong[4] hingeTries;
    ulong[4] hingeAccepts;

    ulong dtElapsed;    // number of dt intervals elapsed (in sweeps)
    StopWatch timer;
    DateTime startTime;
public:
    void setParameters(Parameters params_)
    {
        params = params_;
        // TO DO: put in some checking here (positive, reasonable?
        // no change in dt?)
    }

    this(Manifold!(dim, Vertex) initialManifold)
    {
        manifold_ = initialManifold;
        auto vertices = manifold_.simplices(0).joiner.array.dup.sort;

        // all gaps in list of vertices should be unusued vertices
        if(vertices.front != 0)
        {
            unusedVertices ~= vertices.front.iota.array;
        }
        foreach(i; 0 .. vertices.length - 1)
        {
            if(vertices[i] + 1 != vertices[i+1])
            {
                unusedVertices ~= iota(vertices[i]+1, vertices[i+1]).array;
            }
        }
        assert(unusedVertices.all!(v => !manifold.contains(v.only)));
    }

    ref const(Manifold!(dim, Vertex)) manifold()() const
    {
        return manifold_;
    }

    void report(Writer)(auto ref Writer w)
    {
        auto numMovesTried = bistellarTries[].sum + hingeTries[].sum;
        auto numMovesAccepted = bistellarAccepts[].sum + hingeAccepts[].sum;
        double acceptFrac = double(numMovesAccepted) / numMovesTried;       

        enum width = 103;
        auto maxBins = 24;
        auto maxBar2 = 45;
        auto maxBar3 = 45;

        auto timePerMove = timer.peek / params.triesPerStdoutReport;
        typeof(timePerMove) timePerSweep;
        if(acceptFrac > 0.0)
        {
            timePerSweep = timePerMove 
                * (manifold.fVector.back / acceptFrac).to!ulong;
        }
        timer.reset;

        w.writeln('╔', '═'.repeat(24), '╦', '═'.repeat(30), '╦', '═'.repeat(22),
            '╦', '═'.repeat(22), '╗');

        string timeInfo = "%.1f / %s".format(dtElapsed * params.dt, params.maxSweeps);
        w.writefln("║ sweeps : %13-s ║ move accept fraction : %.3-f ║"
            ~ " Δt/move : %9-s  ║ Δt/sweep : %9-s ║",
            timeInfo, acceptFrac, timePerMove.getFirstPart,
            timePerSweep.getFirstPart);

        w.writeln('╠', '═'.repeat(24), '╩', '═'.repeat(9), '╦', '═'.repeat(20), '╩',
            '═'.repeat(20), "╦═╩", '═'.repeat(22), '╣');

        w.writefln("║ start : %24-s ║ estimated end : %23-s ║",
            this.startTime.to!string, "TO DO!");
        w.writeln('╠', '═'.repeat(34), '╩', '═'.repeat(41), '╩', '═'.repeat(24), '╣');

        alias formatStrings = (string rest) => 
            r"%(%" ~ ((width - dim - 21)/(dim + 1)).to!string ~ rest ~ r"┊ %) ║";
        w.writefln("║ dimension      : " ~ formatStrings("s"), (dim+1).iota);
        w.writefln("║ num simplices  : " ~ formatStrings(",d"), manifold.fVector);
        w.writefln("║ mean degree    : " ~ formatStrings(".3f"),
            (dim+1).iota.map!(dim => manifold.meanDegree(dim)));

        w.write("║                : ");
        w.writeln(" ");

        w.writefln("║ std dev degree : " ~ formatStrings(".3f"),
            (dim+1).iota.map!(dim => manifold.degreeVariance(dim).sqrt));

        w.write("║                : ");
        w.writeln(" ");

        auto vp = this.volumePenalty;
        auto gcp = this.globalCurvaturePenalty;
        auto lcp = this.localCurvaturePenalty;
        auto lsacp = this.localSACurvaturePenalty;

        w.writeln('╠', '═'.repeat(width-2), '╣');
        w.write("║", "  Component", ' '.repeat(7), "Raw",
            ' '.repeat(8), "Coef", ' '.repeat(7), "Value",
            ' '.repeat(6), "Δ/Sweep", ' '.repeat(3), "P-val ║");
        w.writeln("         Targets         ║");
        
        if (params.numFacetsCoef > 0.0)
        {
            w.writef("║ facets      : %.4e * %.4f = %.4e ┊",
                vp, params.numFacetsCoef, params.numFacetsCoef * vp);
            w.writefln(" num facets   : %8,d ║", params.numFacetsTarget);
        }

        if (params.numHingesCoef > 0.0)
        {
            w.writef("║ hinges      : %.4e * %.4f = %.4e ┊",
                gcp, params.numHingesCoef, params.numHingesCoef * gcp);
            w.writefln(" hinge degree : %8.5f ║", params.hingeDegreeTarget);
        }

        if (params.hingeDegreeVarCoef > 0.0)
        {
            w.writef("║ cd2 var deg : %.4e * %.4f = %.4e ┊",
                lcp, params.hingeDegreeVarCoef, params.hingeDegreeVarCoef * lcp);
            w.writeln;
        }

        if (params.cd3DegVarCoef > 0.0)
        {
            w.writef("║ cd3 var deg : %.4e * %.4f = %.4e ┊",
                lsacp, params.cd3DegVarCoef,  params.cd3DegVarCoef * lsacp);
            w.writeln;
        }
        w.writef("║                             TOTAL : %.4e ┊",
            this.objective);

        w.writeln;

        w.writeln('╠', '═'.repeat((width-2)/2), '╦', '═'.repeat((width-2)/2),'╣');
        w.write("║ Bistellar Move   Accept              Try    %    ║");
        w.writeln(" Hinge Move       Accept              Try    %    ║");
        foreach(i; 0 .. dim + 1)
        {
            auto accepts = "%14,d".format(bistellarAccepts[i]);
            auto tries = "%14,d".format(bistellarTries[i]);
            w.writef("║ %2s→ %2-s : %14s / %14s (%5.3f) ║",
                i + 1, dim + 1 - i, accepts, tries,
                double(bistellarAccepts[i])/bistellarTries[i]);
            if (i < hingeAccepts.length)
            {
                accepts = "%14,d".format(hingeAccepts[i]);
                tries = "%14,d".format(hingeTries[i]);
                w.writef(" %2s→ %2-s : %14s / %14s (%5.3f) ║", i + 4,
                        (i + 2) * (dim - 1), accepts, tries,
                        double(hingeAccepts[i]) / hingeTries[i]);
            }
            w.writeln;
        }
        w.writef("║ TOTAL  : %14,d / %14,d (%5.3f) ║",
            bistellarAccepts[].sum, bistellarTries[].sum,
            double(bistellarAccepts[].sum)/bistellarTries[].sum);
        w.writefln(" TOTAL  : %14,d / %14,d (%5.3f) ║",
            hingeAccepts[].sum, hingeTries[].sum,
            double(hingeAccepts[].sum)/hingeTries[].sum);
        w.writeln('╠', '═'.repeat((width-2)/2), '╬', '═'.repeat((width-2)/2),'╣');

        auto maxDeg2 = 2 + maxBins;
        auto maxDeg3 = 2 + 2*maxBins;

        auto hist2 = this.manifold.degreeHistogram(this.manifold.dimension - 2);
        auto tail2 = (hist2.length >= 2 + maxDeg2) ? hist2[maxDeg2 .. $].sum : 0;
        auto maxDegBin2 = max(hist2.maxElement, tail2);
        auto normedHist2 = hist2.map!(freq => real(freq) / maxDegBin2);

        auto hist3 = this.manifold.degreeHistogram(this.manifold.dimension - 3);
        auto tail3 = (hist3.length >= maxDeg3) ? hist3[maxDeg3 .. $].sum : 0;
        auto maxDegBin3 = max(hist3.maxElement, tail3);
        auto normedHist3 = hist3.map!(freq => real(freq) / maxDegBin3);

        auto space = ' '.repeat((width - 67)/4);
        w.writefln("║ %s%s%s ║ %s%s%s ║",
            space, "Codimension-2 Degree Histogram", space,
            space, "Codimension-3 Degree Histogram", space);
        static immutable bars = ['▏','▎','▍','▌','▋','▊','▉','█'];
        foreach(bin; iota(0, maxBins))
        {
            w.writef("║ %2s ", bin + 3);
            if(bin + 2 < normedHist2.length)
            {
                auto nEighths = (maxBar2 * normedHist2[bin + 2] * 8.0).to!int;
                auto bar = bars.back.repeat(nEighths / 8).array;
                if (nEighths % 8 > 0)
                {
                    bar ~= bars[nEighths % 8];
                }
                w.writef("%-*s", maxBar2, bar);
            }
            else
            {
                w.write(' '.repeat(maxBar2));
            }
            w.write(" ║ ");
            w.writef("%2s ", 4 + bin * 2);
            if((3 + bin * 2) < normedHist3.length)
            {
                auto nEighths = (maxBar3 * normedHist3[3 + bin * 2] * 8.0).to!int;
                auto bar = bars.back.repeat(nEighths / 8).array;
                if (nEighths % 8 > 0)
                {
                    bar ~= bars[nEighths % 8];
                }
                w.writef("%-*s", maxBar3, bar);
            }
            else
            {
                w.write(' '.repeat(maxBar3));
            }
            w.writeln(" ║");
        }

        w.write("║  > ");
        auto tailFreq2 = real(tail2) / maxDegBin2;
        auto nEighths2 = (maxBar2 * tailFreq2 * 8.0).to!int;
        auto bar2 = bars.back.repeat(nEighths2 / 8).array;
        if (nEighths2 % 8 > 0)
        {
            bar2 ~= bars[nEighths2 % 8];
        }

        auto tailFreq3 = real(tail3) / maxDegBin3;
        auto nEighths3 = (maxBar3 * tailFreq3 * 8.0).to!int;
        auto bar3 = bars.back.repeat(nEighths3 / 8).array;
        if (nEighths3 % 8 > 0)
        {
            bar3 ~= bars[nEighths3 % 8];
        }

        if (tailFreq2 > 0)
        {
            w.writef("%-*s", maxBar2, bar2);
        }
        else
        {
            w.write(' '.repeat(maxBar2));
        }
        w.write(" ║  > ");

        if (tailFreq3 > 0)
        {
            w.writef("%-*s", maxBar3, bar3);
        }
        else
        {
            w.write(' '.repeat(maxBar3));
        }
        w.writeln(" ║");
        w.writeln('╚', '═'.repeat((width-2)/2), '╩', '═'.repeat((width-2)/2), '╝');
    }

    void columnReport(Writer)(auto ref Writer w, string separator = ",")
    {
        static wroteColumnLabels = false;

        auto degHists = (dim - 1).iota.map!(d => manifold.degreeHistogram(d));

        if (!wroteColumnLabels)
        {
            w.write("num_moves_accepted", separator);
            w.write((dim + 1).iota.map!(d => "num_simps_dim%s".format(d)).joiner(separator));
            w.write(separator);
            w.write((dim - 1).iota.map!(d => "tot_sqr_deg_dim%s".format(d)).joiner(separator));
            w.write(separator);
            w.write((dim - 1).iota.map!(d => "mean_deg_dim%s".format(d)).joiner(separator));
            w.write(separator);
            w.write((dim - 1).iota.map!(d => "stddev_deg_dim%s".format(d)).joiner(separator));
            w.write(separator);
            w.writeln((dim - 1).iota.map!(d => "deg_histogram_dim%s".format(d)).joiner(
                separator ~ "XXX" ~ separator));
            wroteColumnLabels = true;
        }        

        w.write(bistellarAccepts[].sum + hingeAccepts[].sum, separator);
        w.write(manifold.fVector.map!(to!string).joiner(separator));
        w.write(separator);
        w.write((dim - 1).iota.map!(d =>
            manifold.totalSquareDegree(d).to!string).joiner(separator));
        w.write(separator);
        w.write((dim - 1).iota.map!(d =>
            manifold.meanDegree(d).to!string).joiner(separator));
        w.write(separator);
        w.write((dim - 1).iota.map!(d =>
            manifold.degreeVariance(d).sqrt.to!string).joiner(separator));
        foreach(d; 0 .. dim - 1)
        {
            w.write(separator);
            w.write(degHists[d].map!(to!string).joiner(separator));
            if (d < dim - 2)
            {
                w.write(separator, "XXX");
            }
        }
        w.writeln;
    }
}

real volumePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    return (numF - s.params.numFacetsTarget) ^^ 2;
}

real globalCurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;
    immutable numHingesTarget = hPerF * numF / s.params.hingeDegreeTarget;
    return (numH - numHingesTarget) ^^ 2;
}

real localCurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;

    immutable sqDeg = s.manifold.totalSquareDegree(dim - 2);
    immutable real degTarget = hPerF * numF / numH.to!real;

    real _; // dummy for integer part
    immutable real x = degTarget.modf(_);   // fractional part
    immutable minPenalty = (x - x^^2) * numH;

    // TO DO: Refer to arxiv paper for this!
    return (degTarget^^2 * numH - 2*degTarget*hPerF*numF + sqDeg) - minPenalty;
}
// TO DO: unittests!


real localSACurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    immutable numCD3 = s.manifold.fVector[dim - 3];
    immutable cd3PerF = binomial(dim + 1, dim - 2);

    immutable sqDeg = s.manifold.totalSquareDegree(dim - 3);
    immutable real cd3DegTarget = cd3PerF * numF / numCD3.to!real;

    real _; // dummy for integer part
    immutable real x = cd3DegTarget.modf(_);   // fractional part
    immutable minPenalty = (x - x^^2) * numCD3;

    // TO DO: Refer to arxiv paper for this!
    return (cd3DegTarget^^2 * numCD3 - 2*cd3DegTarget*cd3PerF*numF + sqDeg) - minPenalty;
}

Objective objectiveParts(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    auto volPen = s.params.numFacetsCoef * s.volumePenalty;
    auto gcPen = s.params.numHingesCoef * s.globalCurvaturePenalty;
    auto lcPen = s.params.hingeDegreeVarCoef * s.localCurvaturePenalty;
    auto lsacPen = s.params.cd3DegVarCoef * s.localSACurvaturePenalty;

    return Objective(volPen, gcPen, lcPen, lsacPen);
}

real objective(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    auto obj = s.objectiveParts;
    return obj.volPen + obj.gcPen + obj.lcPen + obj.lsacPen;
}

void sample(Vertex, int dim)(ref Sampler!(Vertex, dim) s)
{
    if (s.params.disableGC)
    {
        GC.disable;
    }
    s.startTime = Clock.currTime.to!DateTime;
    s.timer.start;
    
    int sampleNumber;
    auto dtIncThreshold = (s.params.dt * s.params.numFacetsTarget).to!ulong;
    assert(dtIncThreshold > 0);
    auto doneSampling = false;
    while (!doneSampling)
    {
        if(s.dtElapsed * s.params.dt >= s.params.maxSweeps)
        {
            doneSampling = true;
        }

        size_t numMovesTried = s.bistellarTries[].sum + s.hingeTries[].sum;
        size_t numMovesAccepted = s.bistellarAccepts[].sum + s.hingeAccepts[].sum;
        
        bool dtIncreased = false;
        if (numMovesAccepted == dtIncThreshold)
        {
            ++s.dtElapsed;
            dtIncreased = true;
            dtIncThreshold = ((s.dtElapsed + 1) 
                * s.params.dt * s.params.numFacetsTarget).to!ulong;
            assert(dtIncThreshold > numMovesAccepted);
        }

        assert(s.unusedVertices.all!(v => !s.manifold.contains(v.only)));
        if(s.unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            s.unusedVertices ~= s.manifold_.fVector[0].to!int;
        }
        
        real oldObj = s.objective;

        s.manifold_.Move chosenMove;

        real totMoves = s.manifold.numValidMoves + s.manifold.fVector[$-1];
        if (uniform(0, totMoves) < s.manifold.numValidMoves)
        {
            chosenMove = s.manifold.moves.choice;
            while (s.manifold.contains(chosenMove.coCenter))
            {
                chosenMove = s.manifold.moves.choice;
            }
        }
        else
        {
            auto facet = s.manifold.randomFacetOfDim(dim);
            chosenMove = s.manifold_.Move(facet, [s.unusedVertices[$-1]]);
        }
        
        
        s.manifold_.doPachner(chosenMove.center, chosenMove.coCenter);

        ++s.bistellarTries[dim + 1 - chosenMove.center.length];
        ++s.bistellarAccepts[dim + 1 - chosenMove.center.length];

        if (chosenMove.center.length == 1)
        {
            s.unusedVertices ~= chosenMove.center.front;
        }

        if (chosenMove.center.length == dim + 1)
        {
            s.unusedVertices.popBack;
        }
        
        real deltaObj = s.objective - oldObj;        
        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj))) // REJECT MOVE
        {
            --s.bistellarAccepts[dim + 1 - chosenMove.center.length];

            chosenMove = s.manifold.Move(chosenMove.coCenter, chosenMove.center);

            s.manifold_.doPachner(chosenMove.center, chosenMove.coCenter);
            if (chosenMove.center.length == 1)
            {
                s.unusedVertices ~= chosenMove.center.front;
            }

            if (chosenMove.center.length == dim + 1)
            {
                s.unusedVertices.popBack;
            }
        }

        //--------------------------- MAKE REPORT ----------------------------
        if ((numMovesTried % s.params.triesPerStdoutReport == 0)
            || doneSampling)
        {
            s.report(stdout);
        }

        //----------------------- WRITE TO DATA FILE --------------------------
    //     if ((dtIncreased && (s.dtElapsed % s.params.dtPerFileReport == 0))
    //         || doneSampling)
    //     {
    //         auto file = File(s.params.saveFilePrefix ~ ".dat", "a");
    //         s.columnReport(file);
    //     }

    //     //----------------------- SAVE CURRENT MANIFOLD ----------------------
    //     if ((dtIncreased && (s.dtElapsed % s.params.dtPerSave == 0))
    //         || doneSampling)
    //     {
    //         string prefix;
    //         if (doneSampling)
    //         {
    //             prefix = s.params.saveFilePrefix ~ "_final";
    //         }
    //         else
    //         {
    //             prefix = s.params.saveFilePrefix ~ "_save" 
    //                 ~ sampleNumber.to!string;
    //         }

    //         auto mfdFileName = prefix ~ ".mfd";
    //         auto graphFileName = prefix ~ ".edge_graph";

    //         s.manifold.saveTo(mfdFileName);
    //         auto saveFile = File(mfdFileName, "a");
    //         saveFile.writeln;
    //         s.report(saveFile);

    //         s.manifold.saveEdgeGraphTo(graphFileName);
    //         ++sampleNumber;
    //     }

    //     //------------------------- COLLECT GARBAGE --------------------------
    //     if (s.params.disableGC && (numMovesTried % s.params.triesPerCollect == 0))
    //     {
    //         GC.enable;
    //         GC.collect;
    //         GC.disable;
    //     }

    //     //----------------------- CHECK FOR PROBLEMS ----------------------- 
    //     static if (s.params.checkForProblems)
    //     {
    //         if (dtIncreased
    //             && ((s.dtElapsed * s.params.dt) % s.params.sweepsPerProblemCheck == 0))
    //         {
    //             "checking for problems ... ".write;
    //             auto problems = s.manifold.findProblems;
    //             if (!problems.empty)
    //             {
    //                 problems.each!writeln;
    //                 assert(0);
    //             }
    //             "done".writeln;
    //         }
    //     }
    }

    "Finished! Time elapsed: %s".writefln(
        Clock.currTime.to!DateTime -  s.startTime);
}

@Name("sample") unittest
{
    auto mfd = Manifold!3([[11,12,13,14],[11,12,13,15],[11,12,14,15],[11,13,14,15],[12,13,14,15]]);
    auto s = Sampler!(int, 3)(mfd);
    import std.algorithm : canFind;
    assert(11.iota.all!(i => s.unusedVertices.canFind(i)));
    
    // TO DO: Tests!
}

auto getFirstPart(T)(T time)
{
    return time.to!string.findSplit(",")[0].findSplit("and")[0]
        .replace("minutes", "mins");
}
