import algorithms : eulerCharacteristic;
import manifold : coCenter, movesAtFacet, Manifold, standardSphereFacets, doPachnerImpl;
import simplicial_complex : fVector;
import std.algorithm : all, each, filter, joiner, map, max, maxElement, sum;
import std.conv : to;
import std.datetime : Duration, msecs;
import std.datetime.stopwatch : StopWatch;
// import std.format : format, formattedWrite;
import std.math : exp, sqrt;
import std.random : choice, rndGen, uniform01;
import std.range : array, back, chain, empty, front, iota, popBack, popFront, repeat,
    save, walkLength;
import std.stdio : write, writefln, writeln, stdout;
import unit_threaded : Name;
import utility : subsetsOfSize, subsets;

import core.memory : GC;
import std.datetime.systime : Clock;



// I wish we could compute these, but acos dosn't work at compile time. These
// come from wolfram alpha input: Table[N[2 Pi/ArcCos[1/k],30],{k, 2, 16}]
enum real[17] flatDegreeInDim = [
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
enum real numHingesCoef = 0.1;
enum real hingeDegreeVarCoef = 0.0;

enum int triesPerReport = 10000;
enum int maxTries = 1000000;

real meanHingeDegree(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    immutable real numHinges = manifold.fVector[dim - 2];
    immutable real numFacets = manifold.fVector[dim];

    // See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf
    return numFacets * (dim + 1) * dim / (2 * numHinges);
}

struct Sampler(Vertex, int dim)
{
private:
    Manifold!(dim, Vertex) manifold;
    const(Vertex)[] unusedVertices;

    // tryCount[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] tryCount;
    ulong[dim + 1] acceptCount;

    StopWatch moveTimer;
    StopWatch outputTimer;
    StopWatch timer;
public:
    this(Manifold!(dim, Vertex) initialManifold)
    {
        manifold = initialManifold;
    }
    
    void report(Writer)(auto ref Writer w)
    {
        auto tt = timer.peek.total!"usecs" / real(triesPerReport);
        timer.reset;

        auto mt = moveTimer.peek.total!"usecs" / real(triesPerReport);
        moveTimer.reset;

        // auto gt = gcTimer.peek.total!"usecs" / real(triesPerReport);
        // gcTimer.reset;

        w.writeln('#'.repeat(80));
        w.writefln("# %s", Clock.currTime);
        w.writefln("# usec/move (moves) : %s", mt);
        // w.writefln("# usec/move (GC)    : %s", gt);
        w.writefln("# usec/move total   : %s", tt);

        w.writeln("#", '-'.repeat(79));
        w.writefln("# fVector              : %s", manifold.fVector);
        w.writefln("# number of facets     : %s (target: %s)",
            manifold.fVector.back, numFacetsTarget);
        w.writefln("# average hinge-degree : %s (target: %s)",
            manifold.meanHingeDegree, hingeDegreeTarget);

        auto vp = this.volumePenalty;
        auto gcp = this.globalCurvaturePenalty;
        auto lcp = this.localCurvaturePenalty;

        w.writeln("#", '-'.repeat(79));
        w.writefln("# number of facets penalty : %e = %f * %e",
            numFacetsCoef * vp, numFacetsCoef, vp);
        w.writefln("# mean hinge-degee penalty : %e = %f * %e",
            numHingesCoef * gcp, numHingesCoef, gcp);
        w.writefln("# var hinge-degree penalty : %e = %f * %e",
            hingeDegreeVarCoef * lcp, hingeDegreeVarCoef, lcp);
        w.writefln("#          TOTAL OBJECTIVE : %e", this.objective);

        w.writeln("#", '-'.repeat(79));
        w.writeln("# MOVE   : DONE  /  TRIED");
        iota(dim + 1).each!(i => w.writefln("# %s -> %s : %s / %s",
            i + 1, dim + 1 - i, acceptCount[i], tryCount[i]));
        w.writefln("# TOTALS : %s / %s", acceptCount[].sum, tryCount[].sum);        
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
    // TO DO: IMPORTANT! Make this take difference between 
    //
    // ACTUAL = sum_over_s of (deg(s) - targetDeg)^2
    //
    // IDEAL = sum_over_s of (deg(s) - targetDeg)^2 
    //         (if all s had degrees closest to targetDeg)

    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;

    immutable sqDeg = s.manifold.totSqrDegree;
    immutable real degTarget = hPerF * numF / numH.to!real;

    import std.math : modf;

    real _; // dummy for integer part
    immutable real x = degTarget.modf(_);   // fractional part
    immutable real minPenalty = (x - x^^2) * numH;

    // TO DO: Refer to arxiv paper for this!
    return (degTarget^^2 * numH - 2*degTarget*hPerF*numF + sqDeg) - minPenalty;
}

real objective(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    return numFacetsCoef * s.volumePenalty
        + numHingesCoef * s.globalCurvaturePenalty
        + hingeDegreeVarCoef * s.localCurvaturePenalty;
}

void sample(Vertex, int dim)(ref Sampler!(Vertex, dim) s) @safe
{
    // GC.disable;
    s.timer.start;
    // s.moveTimer.start;
    while (s.tryCount[].sum < maxTries)
    {
        assert(s.unusedVertices.all!(v => !s.manifold.contains([v])));
        if(s.unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            s.unusedVertices ~= s.manifold.fVector[0].to!int;
        }

        auto facet = s.manifold.randomFacetOfDim(dim).array;
        auto move = s.manifold.movesAtFacet(facet).map!array.array.choice;

        auto coMove = (move.walkLength == dim + 1) 
            ? [s.unusedVertices.back] 
            : s.manifold.coCenter(move, facet);

        real oldObj = s.objective;

        s.manifold.doPachnerImpl(move, coMove);
        ++s.tryCount[dim + 1 - move.walkLength];

        real newObj = s.objective;

        real deltaObj = newObj - oldObj;
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
            s.manifold.doPachnerImpl(coMove, move);
        }

        //--------------------------- MAKE REPORT ----------------------------
        if (s.tryCount[].sum % triesPerReport == 0)
        {
            s.moveTimer.stop;
            scope(exit) s.moveTimer.start;

            () @trusted {s.report(stdout);} ();
            s.timer.reset;
        }
    }
}

@Name("sample") unittest
{
    // TO DO: Tests!
}