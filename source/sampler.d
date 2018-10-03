import algorithms : eulerCharacteristic;
import manifold : coCenter, degreeHistogram, degreeVariance, doHingeMove, findProblems, getRandomHingeMove, hasValidHingeMove, 
    linkVerticesAtHinge, movesAtFacet, Manifold, standardSphere, totalSquareDegree, doPachner, meanDegree, undoHingeMove;
import simplicial_complex : fVector;
import std.algorithm : all, each, filter, joiner, map, max, min, maxElement, sum;
import std.conv : to;
import std.datetime : Duration, msecs;
import std.datetime.date : DateTime;
import std.datetime.stopwatch : StopWatch;
// import std.format : format, formattedWrite;
import std.math : exp, sqrt;
import std.random : choice, dice, rndGen, uniform, uniform01;
import std.range : array, back, chain, empty, enumerate, front, iota, only, popBack, popFront, repeat,
    save, walkLength;
import std.stdio : write, writef, writefln, writeln, stdout;
import unit_threaded : Name;
import utility : subsetsOfSize, subsets, toStackArray, toStaticArray, StackArray;
import std.format : format;

import core.memory : GC;
import std.datetime.systime : Clock;

import std.variant : Algebraic, visit, tryVisit;


// I wish we could compute these, but acos dosn't work at compile time. These
// come from wolfram alpha input: Table[N[2 Pi/ArcCos[1/k],30],{k, 2, 16}]
static immutable real[17] flatDegreeInDim = [
    real.nan,   // Not applicable in dimension zero!
    2.00000000000000000000000000000,
    6.00000000000000000000000000000,
    5.10429931211954041318017937918,
    4.76679212271567770016331231084,
    4.58814743301323607389345550261,
    4.47728161419381561316532718870,
    4.40168886795573189523776294354,
    4.34681580829256787810763853238,
    4.30515772121519709317292314208,
    4.27244785078781511448809296727,
    4.24607958792781091933915226732,
    4.22436998865935854222871451330,
    4.20618365430421015310353357902,
    4.19072666439811610044839625288,
    4.17742710626470998673691504053,
    4.16586250565979517934736897387];


//-------------------------------- SETTINGS ------------------------------------
// TO DO: make setting setting a coef to 0.0 disable un-needed code

enum int numFacetsTarget = 10_000;
enum real hingeDegreeTarget = flatDegreeInDim[3];

enum real numFacetsCoef = 0.01;
enum real numHingesCoef = 0.01;
enum real hingeDegreeVarCoef = 0.4;

enum int triesPerReport = 100_000;
enum int maxTries = 10_000_000;

enum int triesPerCollect = 500;

enum bool checkForProblems = false; 
enum int triesPerProblemCheck = 200_000_000;

enum useHingeMoves = true;

struct Sampler(Vertex, int dim)
{
private:
    Manifold!(dim, Vertex) manifold_;
    const(Vertex)[] unusedVertices;

    // bistellarTries[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] bistellarTries;
    ulong[dim + 1] bistellarAccepts;

    // hingeTries[j] counts the hinge moves with hinge degree j + 4;
    // TO DO: remove magic constant here 3 + 4 = 7 = max hinge degree
    ulong[4] hingeTries;
    ulong[4] hingeAccepts;

    StopWatch timer;

    enum historyLength = 40;
    static struct History
    {
        real[dim + 1][historyLength] meanDegree;
        real[dim + 1][historyLength] stdDevDegree;
        real[historyLength] objective;

        void advance() {
            foreach(list; History.tupleof)
            {
                iota(1, list.length).each!(
                    i => list[$ - i] = list[$ - i - 1]);
            }
        }
    }

    DateTime startTime;
    History history;

public:
    this(Manifold!(dim, Vertex) initialManifold)
    {
        manifold_ = initialManifold;
    }

    ref const(Manifold!(dim, Vertex)) manifold()() const
    {
        return manifold_;
    }

    void report(Writer)(auto ref Writer w)
    {
        enum width = 103;
        auto maxBins = 29;
        auto maxBar2 = 45;
        auto maxBar3 = 45;


        auto tt = timer.peek.total!"usecs" / real(triesPerReport);
        timer.reset;

        w.writeln('╔', '═'.repeat(32), '╦', '═'.repeat(23), '╦', '═'.repeat(44), '╗');
        w.writefln("║  move accept fraction : %5.3f  ║  usec/move : %7.0-f  ║  started : %30-s  ║",
            double(bistellarAccepts[].sum + hingeAccepts[].sum)/(bistellarTries[].sum + hingeTries[].sum),
            tt, this.startTime.to!string);
        w.writeln('╠', '═'.repeat(32), '╩', '═'.repeat(23), '╩', '═'.repeat(44), '╣');

        alias formatStrings = (string rest) => 
            r"%(%" ~ ((width - dim - 21)/(dim + 1)).to!string ~ rest ~ r"┊ %) ║";
        w.writefln("║ dimension      : " ~ formatStrings("s"), (dim+1).iota);
        w.writefln("║ num simplices  : " ~ formatStrings(",d"), manifold.fVector);
        w.writefln("║ mean degree    : " ~ formatStrings(".2f"),
            (dim+1).iota.map!(dim => manifold.meanDegree(dim)));
        w.writefln("║ std dev degree : " ~ formatStrings(".2f"),
            (dim+1).iota.map!(dim => manifold.degreeVariance(dim).sqrt));
        w.writefln("║ mean history   : " ~ formatStrings("s"), 5.repeat(dim+1));
        w.writefln("║ std dev hist.  : " ~ formatStrings("s"), 5.repeat(dim+1));

        auto vp = this.volumePenalty;
        auto gcp = this.globalCurvaturePenalty;
        auto lcp = this.localCurvaturePenalty;

        w.writeln('╠', '═'.repeat(width-2), '╣');
        w.writefln("║ facets  : %.4e = %.4f * %.4e",
            numFacetsCoef * vp, numFacetsCoef, vp);
        w.writefln("║ hinges  : %.4e = %.4f * %.4e",
            numHingesCoef * gcp, numHingesCoef, gcp);
        w.writefln("║ var deg : %.4e = %.4f * %.4e",
            hingeDegreeVarCoef * lcp, hingeDegreeVarCoef, lcp);
        w.writef("║ TOTAL   : %.4e history: ", this.objective);
        w.writeHistory(this.objective, this.history.objective);
        w.writeln;

        w.writeln('╠', '═'.repeat(width-2), '╣');
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
            w.writef(" %2s→ %2-s : %14s / %14s (%5.3f) ║",
                i + 4, (i + 2)*(dim - 1), accepts, tries,
                double(hingeAccepts[i])/hingeTries[i]);
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

        history.advance;
        history.meanDegree.front = (dim + 1).iota.map!(
            d => manifold.meanDegree(d)).toStaticArray!(dim + 1);

        history.stdDevDegree.front = (dim + 1).iota.map!(
            d => manifold.degreeVariance(d).sqrt).toStaticArray!(dim + 1);
        history.objective.front = this.objective;
    }

    void columnReport(Writer)(auto ref Writer w)
    {

    }
}

void writeHistory(Writer, Value, size_t historyLength)(auto ref Writer w, Value currentVal,
    Value[historyLength] records)
{
    foreach(i; 0 .. historyLength)
    {
        if (records[i] < currentVal)
        {
            w.write("+");
        }
        else if (records[i] > currentVal)
        {
            w.write("-");
        }
        else if (records[i] == currentVal)
        {
            w.write("=");
        }
        else
        {
            w.write("X");
        }
    }
}


real volumePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    return (numF - numFacetsTarget) ^^ 2;
}

real globalCurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;
    immutable numHingesTarget = hPerF * numF / hingeDegreeTarget;
    return (numH - numHingesTarget) ^^ 2;
}

real localCurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;

    immutable sqDeg = s.manifold.totalSquareDegree(dim - 2);
    immutable real degTarget = hPerF * numF / numH.to!real;

    import std.math : modf;

    real _; // dummy for integer part
    immutable real x = degTarget.modf(_);   // fractional part
    immutable real minPenalty = (x - x^^2) * numH;

    // TO DO: Refer to arxiv paper for this!
    return (degTarget^^2 * numH - 2*degTarget*hPerF*numF + sqDeg) - minPenalty;
}
// TO DO: unittests!


real objective(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    return numFacetsCoef * s.volumePenalty
        + numHingesCoef * s.globalCurvaturePenalty
        + hingeDegreeVarCoef * s.localCurvaturePenalty;
}

void sample(Vertex, int dim)(ref Sampler!(Vertex, dim) s)
{
    GC.disable;
    s.startTime = Clock.currTime.to!DateTime;
    s.timer.start;
    
    size_t numMovesTried;
    while (numMovesTried < maxTries)
    {
        numMovesTried = s.bistellarTries[].sum + s.hingeTries[].sum;

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
            else if ((fDim == dim - 2) && fDeg <= 7 && useHingeMoves)
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

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        static if (checkForProblems)
        {
            if (numMovesTried % triesPerProblemCheck == 0)
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

        //------------------------- COLLECT GARBAGE --------------------------
        if (numMovesTried % triesPerCollect == 0)
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //--------------------------- MAKE REPORT ----------------------------
        if (numMovesTried % triesPerReport == 0)
        {
            s.report(stdout);
            s.timer.reset;
        }

    }
}

@Name("sample") unittest
{
    // TO DO: Tests!
}