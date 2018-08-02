import algorithms : eulerCharacteristic;
import manifold : degreeHistogram, doPachner, getCoCenter, Manifold, pachnerMoves,
    standardSphereFacets;
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


enum real numFacetsCoef = 0.01;
enum real numHingesCoef = 0.0;
enum real hingeDegreeVarCoef = 0.0;

//
enum triesPerReport = 500;
enum maxTries = 200000;

//------------------------------------------------------------------------------

real[3] objectiveParts(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    // TO DO: Why does this allocate a closure? Fix it?

    immutable numHinges = manifold.fVector[dim - 2];
    immutable numFacets = manifold.fVector[dim];

    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numHingesTarget = dim * (dim + 1) * numFacets / (2 * hingeDegreeTarget);

    immutable real degStdDev = manifold.degreeHistogram(dim - 2)
        .byKeyValue.map!(p => p.key * (p.value - manifold.meanHingeDegree) ^^ 2).sum.sqrt;

    // TO DO: Put in a more realistic value here...
    immutable real minDegStdDev = 0.5 * numHinges;

    return [numFacetsCoef * (numFacets - numFacetsTarget) ^^ 2,
            numHingesCoef * (numHinges - numHingesTarget) ^^ 2,
            hingeDegreeVarCoef * (degStdDev - minDegStdDev) ^^ 2];
}

real meanHingeDegree(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    immutable real numHinges = manifold.fVector[dim - 2];
    immutable real numFacets = manifold.fVector[dim];

    // See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf
    return numFacets * (dim + 1) * dim / (2 * numHinges);
}

void sampleOld()
{

    // tryCount[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] tryCount;
    ulong[dim + 1] acceptCount;

    auto manifold = Manifold!dim(standardSphereFacets(dim));
    auto oldManifold = manifold;
    auto oldObjective = manifold.objectiveParts[].sum;
    manifold.Vertex[] unusedVertices;
    assert(unusedVertices.all!(v => !manifold.contains([v])));

    StopWatch timer;
    timer.start;
    while (tryCount[].sum < maxTries)
    {
        if(unusedVertices.empty)
        {
            /* If no unused vertices then next unused vertex label
            is just the number of vertices. */
            unusedVertices ~= manifold.fVector[0].to!int;
        }
        assert(unusedVertices.all!(v => !manifold.contains([v])));

        oldObjective = manifold.objectiveParts[].sum;
        auto moves = manifold.pachnerMoves;
        auto chosenMove = moves.choice;

        manifold.Vertex vertexToRemember;
        if (chosenMove.length == 1)
        {
            unusedVertices ~= chosenMove;
            manifold.doPachner(chosenMove);
        }
        else if (chosenMove.length == manifold.dimension + 1)
        {
            manifold.doPachner(chosenMove, unusedVertices.back);
            vertexToRemember = unusedVertices.back;
            unusedVertices.popBack;
        }
        else
        {
            manifold.doPachner(chosenMove);
        }
        ++tryCount[dim + 1 - chosenMove.length];

        if (manifold.objectiveParts[].sum < oldObjective)
        {
            oldManifold = manifold;
            ++acceptCount[dim + 1 - chosenMove.length];
        }
        else
        {
            immutable acceptProb = exp(oldObjective - manifold.objectiveParts[].sum);
            assert(acceptProb >= 0.0);
            assert(acceptProb <= 1.0);

            // TO DO: Put #pachner-moves effect into accept/reject (important!)

            if (uniform01 <= acceptProb)
            {
                oldManifold = manifold;
                ++acceptCount[dim + 1 - chosenMove.length];
            }
            else
            {
                manifold = oldManifold;

                // Make sure to undo any changes to list of unused vertices
                if (chosenMove.length == 1)
                {
                    unusedVertices.popBack;
                }
                else if (chosenMove.length == manifold.dimension + 1)
                {
                    unusedVertices ~= vertexToRemember;
                }
            }
        }

        // --------------------------- MAKE REPORT ----------------------------
        if (tryCount[].sum % triesPerReport == 0)
        {
            writeln;
            '-'.repeat(80).writeln;
            writeln(" MOVE  :  DONE  /  TRIED  |  msec/move : ",
                    timer.peek.total!"msecs" / real(triesPerReport));
            '-'.repeat(80).writeln;
            iota(dim + 1).each!(indx => writeln(indx + 1, " -> ", dim + 1 - indx,
                    " : ", acceptCount[indx], " / ", tryCount[indx]));

            writeln("TOTALS : ", acceptCount[].sum, " / ", tryCount[].sum);
            '-'.repeat(80).writeln;
            writeln("fVector              : ", manifold.fVector);
            writeln("number of facets     : ", manifold.fVector.back,
                    " (target: ", numFacetsTarget, ")");
            writeln("average hinge-degree : ", manifold.meanHingeDegree,
                    " (target: ", hingeDegreeTarget, ")");
            '-'.repeat(80).writeln;
            writeln("number of facets penalty : ", manifold.objectiveParts[0]);
            writeln("mean hinge-degee penalty : ", manifold.objectiveParts[1]);
            writeln("var hinge-degree penalty : ", manifold.objectiveParts[2]);
            writeln("         TOTAL OBJECTIVE : ", manifold.objectiveParts[].sum);
            '-'.repeat(80).writeln;

            timer.reset;
        }
    }
}

struct Sampler(Vertex, int dim)
{
    Manifold!(dim, Vertex) manifold;

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

    }
}

real volumePenalty(Vertex, int dim)(Sampler!(Vertex, dim) s)
{
    immutable numF = s.manifold.fVector[dim];
    return numFacetsCoef * (numF - numFacetsTarget) ^^ 2;
}

real globalCurvaturePenalty(Vertex, int dim)(Sampler!(Vertex, dim) s)
{
    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numF = s.manifold.fVector[dim];
    immutable numH = s.manifold.fVector[dim - 2];
    immutable hPerF = dim * (dim + 1) / 2;

    immutable numHingesTarget = hPerF * numF / hingeDegreeTarget;
    return numHingesCoef * (numH - numHingesTarget) ^^ 2;
}

real localCurvaturePenalty(Vertex, int dim)(Sampler!(Vertex, dim) s)
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


void sample(Vertex, int dim)(Sampler!(Vertex, dim) s)
{
    auto simp = s.manifold.randomFacetOfDim(dim).array;

    auto possibleMoves = simp.subsets.map!array.filter!(
        face => s.manifold.degree(face) == dim + 2 - face.walkLength);

    import std.stdio : writeln;
    simp.writeln;
    possibleMoves.writeln;
}
 

@Name("sample") unittest
{
    // octahedron    
    auto m = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    auto s = Sampler!(int, 2)(m);

        real tot = s.volumePenalty 
        + s.globalCurvaturePenalty + s.localCurvaturePenalty;

    writeln("volume           : ", s.volumePenalty);
    writeln("global curvature : ", s.globalCurvaturePenalty);
    writeln("local curvature  : ", s.localCurvaturePenalty);
    writeln("-----------------:-------------");
    writeln("total objective  : ", tot, "\n");

    s.manifold.writeln;
    s.manifold.doPachner([1,2]);
    s.manifold.writeln;

    s.sample;


    writeln("volume           : ", s.volumePenalty);
    writeln("global curvature : ", s.globalCurvaturePenalty);
    writeln("local curvature  : ", s.localCurvaturePenalty);
    writeln("-----------------:-------------");
    writeln("total objective  : ", tot, "\n");

    s.manifold.getCoCenter([3,4]).writeln;
    s.manifold.getCoCenter([3,4], [3,4,5]).writeln;

}