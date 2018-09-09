import algorithms : eulerCharacteristic;
import manifold : coCenter, degreeVariance, findProblems, movesAtFacet, Manifold, standardSphere, totalSquareDegree, doPachner, meanDegree;
import simplicial_complex : fVector;
import std.algorithm : all, each, filter, joiner, map, max, maxElement, sum;
import std.conv : to;
import std.datetime : Duration, msecs;
import std.datetime.date : DateTime;
import std.datetime.stopwatch : StopWatch;
// import std.format : format, formattedWrite;
import std.math : exp, sqrt;
import std.random : choice, dice, rndGen, uniform, uniform01;
import std.range : array, back, chain, empty, front, iota, popBack, popFront, repeat,
    save, walkLength;
import std.stdio : write, writefln, writeln, stdout;
import unit_threaded : Name;
import utility : subsetsOfSize, subsets, toStaticArray;

import core.memory : GC;
import std.datetime.systime : Clock;



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

enum int numFacetsTarget = 100;
enum real hingeDegreeTarget = flatDegreeInDim[4];

enum real numFacetsCoef = 0.01;
enum real numHingesCoef = 0.05;
enum real hingeDegreeVarCoef = 0.05;

enum int triesPerReport = 100;
enum int maxTries = 1000000;

enum int triesPerCollect = 1000;
enum int triesPerProblemCheck = 200;

enum useHingeMoves = true;

struct Sampler(Vertex, int dim)
{
private:
    Manifold!(dim, Vertex) manifold_;
    const(Vertex)[] unusedVertices;

    // tryCount[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] tryCount;
    ulong[dim + 1] acceptCount;

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

        w.writeln("#", '-'.repeat(79));
        w.writeln("# MOVE   : DONE  /  TRIED");
        iota(dim + 1).each!(i => w.writefln("# %s -> %s : %s / %s",
            i + 1, dim + 1 - i, acceptCount[i], tryCount[i]));
        w.writefln("# TOTALS : %s / %s", acceptCount[].sum, tryCount[].sum);        
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
    while (s.tryCount[].sum < maxTries)
    {
        assert(s.unusedVertices.all!(v => !s.manifold_.contains([v])));
        if(s.unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            s.unusedVertices ~= s.manifold_.fVector[0].to!int;
        }

        auto facet_ = s.manifold_.randomFacetOfDim(dim).toStaticArray!(dim + 1);
        auto facet = facet_[];
        auto moves = s.manifold_.movesAtFacet(facet).map!array.array;

        auto choiceIndx = moves.map!(m => 1.0L / (dim - m.walkLength + 2)).array.dice;
        auto move = moves[choiceIndx];

        auto coMove = (move.walkLength == dim + 1) 
            ? [s.unusedVertices.back] 
            : s.manifold_.coCenter(move, facet);

        real oldObj = s.objective;

        s.manifold_.doPachner(move, coMove);
        ++s.tryCount[dim + 1 - move.walkLength];

        real deltaObj = s.objective - oldObj;
        if ((deltaObj < 0) || (uniform01 < exp(-deltaObj))) // ACCEPT MOVE
        {
            ++s.acceptCount[dim + 1 - move.walkLength];
            if (move.walkLength == 1)
            {
                s.unusedVertices ~= move.front;
            }
            else if (move.walkLength == dim + 1)
            {
                s.unusedVertices.popBack;
            }
        }
        else                                                // REJECT MOVE
        {
            s.manifold_.doPachner(coMove, move);
        }

        //--------------------------- MAKE REPORT ----------------------------
        if (s.tryCount[].sum % triesPerReport == 0)
        {
            s.report(stdout);
            s.timer.reset;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if (s.tryCount[].sum % triesPerCollect == 0)
        {
            GC.enable;
            GC.collect;
            GC.disable;
        }

        //------------------------- COLLECT GARBAGE --------------------------
        if (s.tryCount[].sum % triesPerProblemCheck == 0)
        {
            write("checking for problems! ...");
            assert(s.manifold.findProblems.empty, 
                s.manifold.findProblems.joiner(", ").array.to!string);
            writeln("ok!");
        }
    }
}

@Name("sample") unittest
{
    // TO DO: Tests!
}