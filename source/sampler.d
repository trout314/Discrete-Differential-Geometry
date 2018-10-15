import algorithms : eulerCharacteristic;
import core.memory : GC;
import dstats.regress : linearRegress;
import manifold : coCenter, degreeHistogram, degreeVariance, doHingeMove, findProblems, getRandomHingeMove, hasValidHingeMove, 
    linkVerticesAtHinge, movesAtFacet, Manifold, saveEdgeGraphTo, saveTo, standardSphere, totalSquareDegree, doPachner, meanDegree, undoHingeMove;
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
import utility : binomial, flatDegreeInDim, subsetsOfSize, subsets, toStackArray, toStaticArray, StackArray;


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
    enum historyLength = 30;
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

    // -------- history of quantities of interest ----------
    size_t[dim + 1][params.historyLength] fVecHistory;
    real[dim + 1][params.historyLength] meanDegHistory;
    real[dim + 1][params.historyLength] stdDevDegHistory;

    Objective[params.historyLength] objHistory;

    // update history
    void advanceHistory()
    {
        foreach(i; 1 .. params.historyLength)
        {
            fVecHistory[$ - i] = fVecHistory[$ - i - 1];
            meanDegHistory[$ - i] = meanDegHistory[$ - i - 1];
            stdDevDegHistory[$ - i] = stdDevDegHistory[$ - i - 1];
            objHistory[$ - i] = objHistory[$ - i - 1];
        }

        fVecHistory.front = manifold.fVector;
        meanDegHistory.front = (dim + 1).iota.map!(
            d => manifold.meanDegree(d)).toStaticArray!(dim + 1);
        stdDevDegHistory.front = (dim + 1).iota.map!(
            d => manifold.degreeVariance(d).sqrt).toStaticArray!(dim + 1);
        objHistory.front = this.objectiveParts;
    }

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
        foreach(i; 0 .. vertices.length - 1)
        {
            // all gaps in list of vertices should be unusued vertices
            if(vertices[i] + 1 != vertices[i+1])
            {
                foreach(v; iota(vertices[i]+1, vertices[i+1]))
                {
                    unusedVertices ~= v;
                }
            }
        }
        assert(unusedVertices.all!(v => !manifold.contains(v.only)));

        this.advanceHistory;
    }

    ref const(Manifold!(dim, Vertex)) manifold()() const
    {
        return manifold_;
    }

    // return historic average of mean degree
    auto historicMeanDeg(int dim)
    {
        return meanDegHistory[].map!(entry => entry[dim]).mean;
    }

    auto historicStdDevDeg(int dim)
    {
        return stdDevDegHistory[].map!(entry => entry[dim]).mean;
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

        foreach(d; 0 .. dim - 1)
        {

            w.writef("mean degree in dimension %2s : %s",
                d, manifold.meanDegree(d));
            auto degHist = meanDegHistory[d];
            w.writeRegression(degHist, params.dt);
            w.writeln;

            w.writef("std dev degree in dimension %2s : %s",
                d, manifold.degreeVariance(d).sqrt);
            auto degStdDev = meanDegHistory[].map!(x => x[d].sqrt).array;
            w.writeRegression(degStdDev, params.dt);
            w.writeln;
        }


        w.writeln('╔', '═'.repeat(24), '╦', '═'.repeat(30), '╦', '═'.repeat(22),
            '╦', '═'.repeat(22), '╗');

        string timeInfo = "%.1f / %s".format(dtElapsed * params.dt, params.maxSweeps);
        w.writefln("║ sweeps : %13-s ║ move accept fraction : %.3-f ║"
            ~ " Δt/move : %9-s  ║ Δt/sweep : %9-s ║",
            timeInfo, acceptFrac, timePerMove.getFirstPart,
            timePerSweep.getFirstPart);

        w.writeln('╠', '═'.repeat(24), '╩', '═'.repeat(9), '╦', '═'.repeat(20), '╩',
            '═'.repeat(20), "╦═╩", '═'.repeat(22), '╣');

        auto nHistStored = objHistory[].filter!(x => !x.volPen.isNaN).walkLength;
        w.writefln("║ start : %24-s ║ estimated end : %23-s ║ # regression pts : %3s ║",
            this.startTime.to!string, "TO DO!", nHistStored);
        w.writeln('╠', '═'.repeat(34), '╩', '═'.repeat(41), '╩', '═'.repeat(24), '╣');

        alias formatStrings = (string rest) => 
            r"%(%" ~ ((width - dim - 21)/(dim + 1)).to!string ~ rest ~ r"┊ %) ║";
        w.writefln("║ dimension      : " ~ formatStrings("s"), (dim+1).iota);
        w.writefln("║ num simplices  : " ~ formatStrings(",d"), manifold.fVector);
        w.writefln("║ mean degree    : " ~ formatStrings(".3f"),
            (dim+1).iota.map!(dim => manifold.meanDegree(dim)));

        w.write("║                : ");
        foreach(d; 0 .. dim - 1)
        {
            auto degHist = meanDegHistory[d];
            w.write(" ");
            w.writeRegression(degHist, params.dt);
        }
        w.writeln(" ");

        w.writefln("║ std dev degree : " ~ formatStrings(".3f"),
            (dim+1).iota.map!(dim => manifold.degreeVariance(dim).sqrt));

        w.write("║                : ");
        foreach(d; 0 .. dim - 1)
        {
            auto degStdDev = meanDegHistory[].map!(x => x[d].sqrt).array;
            w.write(" ");
            w.writeRegression(degStdDev, params.dt);
        }
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
            auto hist = objHistory[0 .. nHistStored].map!(
                x => params.numFacetsCoef * x.volPen).array;
            w.writeRegression(hist, params.dt);
            w.writefln(" num facets   : %8,d ║", params.numFacetsTarget);
        }

        if (params.numHingesCoef > 0.0)
        {
            w.writef("║ hinges      : %.4e * %.4f = %.4e ┊",
                gcp, params.numHingesCoef, params.numHingesCoef * gcp);
            auto hist = objHistory[0 .. nHistStored].map!(
                x => params.numHingesCoef * x.gcPen).array;
            w.writeRegression(hist, params.dt);
            w.writefln(" hinge degree : %8.5f ║", params.hingeDegreeTarget);
        }

        if (params.hingeDegreeVarCoef > 0.0)
        {
            w.writef("║ cd2 var deg : %.4e * %.4f = %.4e ┊",
                lcp, params.hingeDegreeVarCoef, params.hingeDegreeVarCoef * lcp);
            auto hist = objHistory[0 .. nHistStored].map!(
                x => params.hingeDegreeVarCoef * x.lcPen).array;
            w.writeRegression(hist, params.dt);
            w.writeln;
        }

        if (params.cd3DegVarCoef > 0.0)
        {
            w.writef("║ cd3 var deg : %.4e * %.4f = %.4e ┊",
                lsacp, params.cd3DegVarCoef,  params.cd3DegVarCoef * lsacp);
            auto hist = objHistory[0 .. nHistStored].map!(
                x => params.cd3DegVarCoef * x.lsacPen).array;
            w.writeRegression(hist, params.dt);
            w.writeln;
        }
        w.writef("║                             TOTAL : %.4e ┊",
            this.objective);

        auto totObjHist = objHistory[0 .. nHistStored].map!(
            x => params.numFacetsCoef * x.volPen 
                + params.numHingesCoef * x.gcPen
                + params.hingeDegreeVarCoef * x.lcPen 
                + params.cd3DegVarCoef * x.lsacPen).array;
        w.writeRegression(totObjHist, params.dt);
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

    void columnReport(Writer)(auto ref Writer w)
    {

    }
}

void writeRegression(H, Writer)(auto ref Writer w, H history, real dt)
{
    static if (isFloatingPoint!(typeof(history[].front)))
    {
        auto filteredHist = history[].filter!(
            val => !val.isNaN).array;
    }
    else
    {
        auto filteredHist = history[].filter!(
            val => val != 0).array;
    }
    
    auto nSweepsHistory = filteredHist.length.iota.retro.map!(i => dt * i);
    auto result = linearRegress(filteredHist, nSweepsHistory, 1.repeat);
    w.writef(" %8.1e ┊ %5.3f ║", result.betas.front, result.p.front);
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
    while (s.dtElapsed * s.params.dt < s.params.maxSweeps)
    {
        size_t numMovesTried = s.bistellarTries[].sum + s.hingeTries[].sum;
        size_t numMovesAccepted = s.bistellarAccepts[].sum + s.hingeAccepts[].sum;
        
        bool dtIncreased = false;
        if (dtIncThreshold == numMovesAccepted)
        {
            ++s.dtElapsed;
            dtIncreased = true;
            dtIncThreshold = numMovesAccepted +
                (s.params.dt * s.manifold.fVector.back).to!ulong;
            assert(dtIncThreshold >= numMovesAccepted);
        }

        assert(s.unusedVertices.all!(v => !s.manifold.contains(v.only)));
        if(s.unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            s.unusedVertices ~= s.manifold_.fVector[0].to!int;
        }

        auto facet_ = s.manifold_.randomFacetOfDim(dim).toStaticArray!(dim + 1);
        auto facet = facet_[];

        struct BistellarMove
        {
            double weight;
            typeof(facet.subsets.front()) center;
            typeof(s.manifold.coCenter(center, facet)) coCenter;            
        }

        struct HingeMove
        {
            double weight;
            typeof(facet.subsets.front()) hinge;
            typeof(s.manifold.linkVerticesAtHinge(hinge, facet)) hingeLink;
            // used for storing hinge move used (if any)
            // we start it at -2 because this is invalid, and not
            // the invalid index -1 used by getRandomHingeMove
            int diskIndex_ = -2;
        }

        alias Move_ = Algebraic!(BistellarMove, HingeMove);
        StackArray!(Move_, 2^^(dim + 1) - 1) moves;

        foreach(f; facet.subsets)
        {
            auto face_ = f.toStackArray!(Vertex, dim + 1);
            auto face = face_[];
            auto fDim = face.length - 1;
            auto fDeg = s.manifold.degree(face);
            if ((fDim == dim) || (fDeg == dim + 1 - fDim))
            {
                moves ~= Move_(BistellarMove());
                auto m_ = moves.back.peek!BistellarMove;
                m_.weight = 1.0 / (dim - fDim + 1);
                m_.center = f;
            }
            // TO DO: Get rid of magic constant here. It's max deg hinge move.
            else if ((fDim == dim - 2) && fDeg <= 7 && s.params.useHingeMoves)
            {
                moves ~= Move_(HingeMove());
                auto m_ = moves.back.peek!HingeMove;
                m_.weight = 1.0L / fDeg;
                m_.hinge = f;
            }
        }
        
        real oldObj = s.objective;
        size_t mIndx_;
        bool done = false;
        while(!done)     // Choose a move involving the facet
        {
            mIndx_ = moves[].map!(_ => _.visit!(m => m.weight)).dice;
            bool moveInvalid = false;

            // TO DO: Why can't I get ref access to work here?!
            moves[mIndx_].visit!(
                (BistellarMove m)
                {
                    auto m_ = moves[mIndx_].peek!BistellarMove;
                    if(m_.center.walkLength == dim + 1)
                    {
                        // (dim + 1) -> 1 moves are always valid
                        done = true;
                        m_.coCenter = [s.unusedVertices.back];
                    }
                    else
                    {
                        auto coMove = s.manifold_.coCenter(m_.center, facet);
                        assert(!coMove.empty);
                        if(!s.manifold.contains(coMove))
                        {
                            m_.coCenter = coMove;
                            done = true;
                        }
                        else
                        {
                            moveInvalid = true;
                        }
                    }

                    if (done)
                    {
                        s.manifold_.doPachner(m_.center.array, m_.coCenter.array);

                        ++s.bistellarTries[dim + 1 - m_.center.walkLength];
                        ++s.bistellarAccepts[dim + 1 - m_.center.walkLength];
                        if (m_.center.walkLength == 1)
                        {
                            s.unusedVertices ~= m.center.front;
                        }
                        else if (m_.center.walkLength == dim + 1)
                        {
                            s.unusedVertices.popBack;
                        }
                    }
                },
                (HingeMove m)
                {
                    auto m_ = moves[mIndx_].peek!HingeMove;

                    auto lnkVerts = s.manifold.linkVerticesAtHinge(m_.hinge, facet);
                    auto indx = s.manifold.getRandomHingeMove(m_.hinge.array.dup, lnkVerts.array.dup);
                    if (indx >= 0)  // Some hinge move was valid
                    {
                        m_.hingeLink = lnkVerts;
                        m_.diskIndex_ = indx;
                        s.manifold_.doHingeMove(m_.hinge.array.dup, m_.hingeLink.array.dup, m_.diskIndex_);
                        ++s.hingeTries[m_.hingeLink.walkLength - 4];
                        ++s.hingeAccepts[m_.hingeLink.walkLength - 4];
                        done = true;
                    }
                    else            // No hinge move was valid
                    {
                        assert(indx == -1);
                        moveInvalid = true;
                    }
                },
                () {assert(0, "move type not implemented");}
            );

            if(moveInvalid)
            {
                assert(!done);
                moves[mIndx_] = moves.back;
                assert(moves.length > 0);
                moves.length = moves.length - 1;
            }
        }


        real deltaObj = s.objective - oldObj;        

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj))) // REJECT MOVE
        {
            moves[mIndx_].visit!(
                (BistellarMove m)
                {
                    auto m_ = moves[mIndx_].peek!BistellarMove;            
                    s.manifold_.doPachner(m.coCenter.array, m.center.array);
                    --s.bistellarAccepts[dim + 1 - m_.center.walkLength];
                    if (m.center.walkLength == 1)
                    {
                        s.unusedVertices.popBack;
                    }
                    else if (m.center.walkLength == dim + 1)
                    {
                        s.unusedVertices ~= m.coCenter.front;
                    }
                },
                (HingeMove m)
                {
                    assert(!m.hinge.empty);
                    assert(!m.hingeLink.empty);
                    assert(!m.diskIndex_ >= 0);
                    
                    auto m_ = moves[mIndx_].peek!HingeMove;

                    s.manifold_.undoHingeMove(m_.hinge, m_.hingeLink, m_.diskIndex_);
                    --s.hingeAccepts[m_.hingeLink.walkLength - 4];
                },
                () {assert(0, "move type not implemented");}
            );
        }

        //--------------------------- MAKE REPORT ----------------------------
        if (numMovesTried % s.params.triesPerStdoutReport == 0)
        {
            s.report(stdout);
        }

        //----------------------- WRITE TO DATA FILE --------------------------
        if (dtIncreased && (s.dtElapsed % s.params.dtPerFileReport == 0))
        {
            // TO DO: Implement this!
        }

        //----------------------- SAVE CURRENT MANIFOLD ----------------------
        if (dtIncreased && (s.dtElapsed % s.params.dtPerSave == 0))
        {
            auto prefix = s.params.saveFilePrefix 
                ~ sampleNumber.to!string;
            auto mfdFileName = prefix ~ ".mfd";
            auto graphFileName = prefix ~ ".edge_graph";

            s.manifold.saveTo(mfdFileName);
            auto saveFile = File(mfdFileName, "a");
            saveFile.writeln;
            s.report(saveFile);

            s.manifold.saveEdgeGraphTo(graphFileName);
            ++sampleNumber;
        }

        //----------------------- RECORD HISTORY --------------------------
        if (dtIncreased && (s.dtElapsed % s.params.dtPerHistory == 0))
        {
            // TO DO: Implement this!
            s.advanceHistory;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if (s.params.disableGC && (numMovesTried % s.params.triesPerCollect == 0))
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        static if (s.params.checkForProblems)
        {
            if (dtIncreased
                && ((s.dtElapsed * s.params.dt) % s.params.sweepsPerProblemCheck == 0))
            {
                "checking for problems ... ".write;
                auto problems = s.manifold.findProblems;
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
        Clock.currTime.to!DateTime -  s.startTime);
}

@Name("sample") unittest
{
    // TO DO: Tests!
}

auto getFirstPart(T)(T time)
{
    return time.to!string.findSplit(",")[0].findSplit("and")[0]
        .replace("minutes", "mins");
}
