import algorithms : eulerCharacteristic;
import manifold : degreeHistogram, coCenter, movesAtFacet, Manifold, standardSphereFacets, doPachnerImpl;
import simplicial_complex : fVector;
import std.algorithm : all, each, filter, joiner, map, max, maxElement, sum;
import std.conv : to;
import std.datetime : Duration, msecs;
import std.datetime.stopwatch : StopWatch;
import std.format : format;
import std.math : exp, sqrt;
import std.random : choice, rndGen, uniform01;
import std.range : array, back, chain, empty, front, iota, popBack, popFront, repeat,
    save, walkLength;
import std.stdio : writefln, writeln;
import unit_threaded : Name;
import utility : subsetsOfSize, subsets;

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

enum dim = 3;

enum int numFacetsTarget = 1000;
enum real hingeDegreeTarget = flatDegreeInDim[dim];


enum real numFacetsCoef = 0.1;
enum real numHingesCoef = 0.05;
enum real hingeDegreeVarCoef = 0.0;

//
enum triesPerReport = 1000;
enum maxTries = 50000;

real meanHingeDegree(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    immutable real numHinges = manifold.fVector[dim - 2];
    immutable real numFacets = manifold.fVector[dim];

    // See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf
    return numFacets * (dim + 1) * dim / (2 * numHinges);
}

struct Sampler(Vertex, int dim)
{
    Manifold!(dim, Vertex) manifold;
    const(Vertex)[] unusedVertices;

    // tryCount[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] tryCount;
    ulong[dim + 1] acceptCount;

    ulong totalSquaredHingeDegree;

    this(Manifold!(dim, Vertex) initialManifold)
    {
        manifold = initialManifold;
        totalSquaredHingeDegree = initialManifold.simplices(dim - 2)
            .map!(s => manifold.degree(s)^^2).sum;
    }

    private void updateTotSqDeg(
        const(Vertex)[] center,
        const(Vertex)[] coCenter)
    {
        // TO DO: Implement this!
    }
}

real volumePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    return numFacetsCoef * (numF - numFacetsTarget) ^^ 2;
}

real globalCurvaturePenalty(Vertex, int dim)(const ref Sampler!(Vertex, dim) s)
{
    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;

    immutable numHingesTarget = hPerF * numF / hingeDegreeTarget;
    return numHingesCoef * (numH - numHingesTarget) ^^ 2;
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

    immutable sqDeg = s.totalSquaredHingeDegree;
    immutable real degTarget = hPerF * numF / numH.to!real;

    // TO DO: Refer to arxiv paper for this!
    return hingeDegreeVarCoef * (
        (degTarget^^2 * numH - 2*degTarget*hPerF*numF + sqDeg)
        -numH)^^2; // TO DO: Fix this fudge! See TO DO above...
}

void sample(Vertex, int dim)(ref Sampler!(Vertex, dim) s)
{
    StopWatch timer;
    timer.start;

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

        real oldObjective = s.volumePenalty
            + s.globalCurvaturePenalty 
            + s.localCurvaturePenalty;

        s.manifold.doPachnerImpl(move, coMove);
        ++s.tryCount[dim + 1 - move.walkLength];

        real newObjective = s.volumePenalty
            + s.globalCurvaturePenalty 
            + s.localCurvaturePenalty;

        real deltaObj = newObjective - oldObjective;
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

        // --------------------------- MAKE REPORT ----------------------------
        if (s.tryCount[].sum % triesPerReport == 0)
        {
            writeln;
            '-'.repeat(80).writeln;
            writeln(" MOVE  :  DONE  /  TRIED  |  msec/move : ",
                    timer.peek.total!"msecs" / real(triesPerReport));
            '-'.repeat(80).writeln;
            iota(dim + 1).each!(indx => writeln(indx + 1, " -> ", dim + 1 - indx,
                    " : ", s.acceptCount[indx], " / ", s.tryCount[indx]));

            writeln("TOTALS : ", s.acceptCount[].sum, " / ", s.tryCount[].sum);
            '-'.repeat(80).writeln;
            writeln("fVector              : ", s.manifold.fVector);
            writeln("number of facets     : ", s.manifold.fVector.back,
                    " (target: ", numFacetsTarget, ")");
            writeln("average hinge-degree : ", s.manifold.meanHingeDegree,
                    " (target: ", hingeDegreeTarget, ")");
            '-'.repeat(80).writeln;
            writeln("number of facets penalty : ", s.volumePenalty);
            writeln("mean hinge-degee penalty : ", s.globalCurvaturePenalty);
            writeln("var hinge-degree penalty : ", s.localCurvaturePenalty);
            writeln("         TOTAL OBJECTIVE : ", s.volumePenalty
            + s.globalCurvaturePenalty 
            + s.localCurvaturePenalty);
            '-'.repeat(80).writeln;

            timer.reset;
        }
    }
}
 

@Name("sample") unittest
{

}