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

enum int numFacetsTarget = 16000;
enum real hingeDegreeTarget = flatDegreeInDim[3];

enum real numFacetsCoef = 0.01;
enum real numHingesCoef = 0.0;
enum real hingeDegreeVarCoef = 0.0;

enum int triesPerReport = 100_000;
enum int maxTries = 100_000_000;

enum int triesPerCollect = 1_000;
enum int triesPerProblemCheck = 1_000_000;

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
        auto tt = timer.peek.total!"usecs" / real(triesPerReport);
        timer.reset;

        w.writeln('#'.repeat(80));
        w.writefln("# %s  |  usec/move : %s", Clock.currTime.to!DateTime, tt);

        w.writeln("#", '-'.repeat(79));
        w.writefln("# fVector           : %s", manifold.fVector);
        w.writefln("# average degree    : %s",
            (dim+1).iota.map!(dim => manifold.meanDegree(dim)));
        w.writefln("# std dev in degree : %s",
            (dim+1).iota.map!(dim => manifold.degreeVariance(dim).sqrt));

        auto vp = this.volumePenalty;
        auto gcp = this.globalCurvaturePenalty;
        auto lcp = this.localCurvaturePenalty;

        w.writeln("#", '-'.repeat(79));
        w.writefln("# number of facets penalty (A): %e = %f * %e",
            numFacetsCoef * vp, numFacetsCoef, vp);
        w.writefln("# mean hinge-degee penalty (B): %e = %f * %e",
            numHingesCoef * gcp, numHingesCoef, gcp);
        w.writefln("# var hinge-degree penalty (C): %e = %f * %e",
            hingeDegreeVarCoef * lcp, hingeDegreeVarCoef, lcp);
        w.writefln("#             TOTAL OBJECTIVE : %e %s %e", this.objective,
            "Raw A + Raw B :", vp + gcp);

        w.write("#", '-'.repeat(11), " Bistellar Moves ", '-'.repeat(11));
        w.writeln("-|", '-'.repeat(12), " Hinge Moves ", '-'.repeat(13));

        // w.writeln("# MOVE   :       DONE /      TRIED");
        foreach(i; 0 .. dim + 1)
        {
            w.writef("# %2s→ %2-s : %13,s / %13,s |",
                i + 1, dim + 1 - i, bistellarAccepts[i], bistellarTries[i]);
            if (i < hingeAccepts.length)
            {
                w.writefln(" %2s→ %2-s : %12,s / %13,s ",
                i + 4, (i + 2)*(dim - 1), hingeAccepts[i], hingeTries[i]);
            }
            else
            {
                w.writeln;
            }
        }

        // iota(dim + 1).each!(i => w.writefln("# %2s→ %2-s : %10,s / %10,s",
        //     i + 1, dim + 1 - i, bistellarAccepts[i], bistellarTries[i]));
        // w.writefln("# TOTALS : %10,s / %10,s", bistellarAccepts[].sum, bistellarTries[].sum);

        // w.writeln("#", '-'.repeat(79));
        // w.writeln("# MOVE   : DONE  /  TRIED");
        // iota(4).each!(i => w.writefln("# deg %s : %10s / %10s",
        //     i + 4, hingeAccepts[i], hingeTries[i]));
        // w.writefln("# TOTALS : %10,s / %10,s", hingeAccepts[].sum, hingeTries[].sum);


        auto hist = this.manifold.degreeHistogram(this.manifold.dimension - 2);

        auto maxBar = 68;
        auto maxDeg = 16;
        auto maxDegBin = hist.maxElement;
        auto normedHist = hist.map!(freq => real(freq) / maxDegBin);
        writeln("# degree | frequency");
        foreach(bin; 2 .. maxDeg)
        {
            if(bin < normedHist.length)
            {
                auto nAst = (maxBar * normedHist[bin]).to!int;
                "# %6s | %s".writef(bin + 1, '*'.repeat(nAst));

                if ((normedHist[bin] > 0) && (nAst == 0))
                {
                    ".".writeln;
                }
                else
                {
                    writeln;

                }
            }
            else
            {
                "# %6s |".writefln(bin + 1);
            }
        }
        if (normedHist.length >= maxDeg)
        {
            auto tailFreq = real(hist[maxDeg .. $].sum) / maxDegBin;
            auto nAst = (maxBar * tailFreq).to!int;
            "#  >= %s | %s".writef(maxDeg + 1, '*'.repeat(nAst));
            if ((tailFreq > 0) && (nAst == 0))
            {
                ".".writeln;
            }
            else
            {
                writeln;
            }
        }
        else
        {
            "#  >= %s | ".writefln(maxDeg + 1);
        }
    }

    void columnReport(Writer)(auto ref Writer w)
    {

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
    s.timer.start;

    while ((s.bistellarTries[].sum + s.hingeTries[].sum) < maxTries)
    {
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
        StackArray!(Move_, 2^^(dim + 1) - 1) movesNew_;

        foreach(f; facet.subsets)
        {
            auto face_ = f.toStackArray!(Vertex, dim + 1);
            auto face = face_[];
            auto fDim = face.length - 1;
            auto fDeg = s.manifold.degree(face);
            if ((fDim == dim) || (fDeg == dim + 1 - fDim))
            {
                movesNew_ ~= Move_(BistellarMove());
                movesNew_.back.peek!BistellarMove.weight = 1.0 / (dim - fDim + 1);
                movesNew_.back.peek!BistellarMove.center = f;
            }
            // TO DO: Get rid of magic constant here. It's max deg hinge move.
            else if ((fDim == dim - 2) && fDeg <= 7 && useHingeMoves)
            {
                movesNew_ ~= Move_(HingeMove());
                movesNew_.back.peek!HingeMove.weight = 1.0L / fDeg;
                movesNew_.back.peek!HingeMove.hinge = f;
            }
        }
        
        real oldObj = s.objective;
        size_t mIndx_;
        bool done = false;
        while(!done)     // Choose a move involving the facet
        {
            mIndx_ = movesNew_[].map!(_ => _.visit!(m => m.weight)).dice;
            bool moveInvalid = false;

            // TO DO: Why can't I get ref access to work here?!
            movesNew_[mIndx_].visit!(
                (BistellarMove m)
                {
                    auto m_ = movesNew_[mIndx_].peek!BistellarMove;
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
                    auto m_ = movesNew_[mIndx_].peek!HingeMove;

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
                movesNew_[mIndx_] = movesNew_.back;
                movesNew_.length = movesNew_.length - 1;
            }
        }


        real deltaObj = s.objective - oldObj;        

        if ((deltaObj > 0) && (uniform01 > exp(-deltaObj))) // REJECT MOVE
        {
            movesNew_[mIndx_].visit!(
                (BistellarMove m)
                {
                    auto m_ = movesNew_[mIndx_].peek!BistellarMove;            
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
                    
                    auto m_ = movesNew_[mIndx_].peek!HingeMove;

                    s.manifold_.undoHingeMove(m.hinge.array.dup, m.hingeLink.array.dup, m.diskIndex_);
                    --s.hingeAccepts[m_.hingeLink.walkLength - 4];
                },
                () {assert(0, "move type not implemented");}
            );
        }

        // (s.bistellarTries[].sum + s.hingeTries[].sum).writeln;
        //--------------------------- MAKE REPORT ----------------------------
        if ((s.bistellarTries[].sum + s.hingeTries[].sum) % triesPerReport == 0)
        {
            s.report(stdout);
            s.timer.reset;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if ((s.bistellarTries[].sum + s.hingeTries[].sum) % triesPerCollect == 0)
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //----------------------- CHECK FOR PROBLEMS ----------------------- 
        if ((s.bistellarTries[].sum + s.hingeTries[].sum) % triesPerProblemCheck == 0)
        {
            assert(s.manifold.findProblems.empty, 
                s.manifold.findProblems.joiner(", ").array.to!string);
        }
    }
}

@Name("sample") unittest
{
    // TO DO: Tests!
}