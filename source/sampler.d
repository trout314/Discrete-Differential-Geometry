import algorithms : eulerCharacteristic;
import manifold : degreeHistogram, doPachner, Manifold, pachnerMoves,
    standardSphereFacets;
import simplicial_complex : fVector;
import std.algorithm : all, each, joiner, map, max, maxElement, sum;
import std.conv : to;
import std.datetime : AutoStart, Duration, msecs, StopWatch;
import std.format : format;
import std.math : exp, sqrt;
import std.random : back, choice, rndGen, uniform01;
import std.range : array, chain, empty, front, iota, popBack, popFront, repeat,
    save;
import std.stdio : writefln, writeln;
import utility : subsetsOfSize;

//-------------------------------- SETTINGS ------------------------------------
    
enum int numFacetsTarget = 1000;
enum real numFacetsCoef = 0.1;

enum real meanHingeDegreeTarget = 5.1;
enum real numHingesCoef = 0.5;

// TO DO: make setting this to 0.0 disable tracking of hinge degrees
enum real degreeStdDevCoef = 0.0;

enum dim = 3;
enum triesPerReport = 200;
enum maxVertices = 1000;
enum maxTries = 20000;

//------------------------------------------------------------------------------

real[3] objectiveParts(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    // TO DO: Why does this allocate a closure? Fix it?

    immutable numHinges = manifold.fVector[dim - 2];
    immutable numFacets = manifold.fVector[dim];

    /* The target mean hinge-degree and current number of facets imply a target
    number of hinges. See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf */
    immutable numHingesTarget = dim*(dim + 1) * numFacets 
        / (2 * meanHingeDegreeTarget);

    immutable real degStdDev = manifold.degreeHistogram(dim - 2).byKeyValue.map!(
        p => p.key * (p.value - manifold.meanHingeDegree)^^2).sum.sqrt;

    // TO DO: Put in a more realistic value here...
    immutable real minDegStdDev = 0.5 * numHinges;

    return [numFacetsCoef * (numFacets - numFacetsTarget)^^2,
        numHingesCoef * (numHinges - numHingesTarget)^^2,
        degreeStdDevCoef * (degStdDev - minDegStdDev)^^2];
}

real meanHingeDegree(Vertex, int dim)(const ref Manifold!(dim, Vertex) manifold)
{
    immutable real numHinges = manifold.fVector[dim - 2];
    immutable real numFacets = manifold.fVector[dim];

    // See Eq. 7, p. 5 of https://arxiv.org/pdf/1208.1514.pdf
    return numFacets * (dim + 1) * dim / (2 *  numHinges);
}

void sample()
{

     // tryCount[j] counts j + 1 -> dim + 1 - j moves tried
    ulong[dim + 1] tryCount;
    ulong[dim + 1] acceptCount;

    auto manifold = Manifold!dim(standardSphereFacets(dim));
    auto oldManifold = manifold;
    auto oldObjective = manifold.objectiveParts[].sum;

    auto unusedVertices = iota(dim + 2, maxVertices).array;
    assert(unusedVertices.all!(v => !manifold.contains([v])));

    auto timer = StopWatch(AutoStart.yes);  
    while(!unusedVertices.empty && tryCount[].sum < maxTries)
    {
        assert(unusedVertices.all!(v => !manifold.contains([v])));

        oldObjective = manifold.objectiveParts[].sum;
        auto moves = manifold.pachnerMoves;
        auto chosenMove = moves.choice;

        manifold.Vertex vertexToRemember;
        if(chosenMove.length == 1)
        {
            unusedVertices ~= chosenMove;
            manifold.doPachner(chosenMove);
        }
        else if(chosenMove.length == manifold.dimension + 1)
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

        if(manifold.objectiveParts[].sum < oldObjective)
        {
            oldManifold = manifold;
            ++acceptCount[dim + 1 - chosenMove.length];
        }
        else
        {
            immutable acceptProb 
                = exp(oldObjective - manifold.objectiveParts[].sum);
            assert(acceptProb >= 0.0);
            assert(acceptProb <= 1.0);

            // TO DO: Put #pachner-moves effect into accept/reject (important!)

            if(uniform01 <= acceptProb)
            {
                oldManifold = manifold;
                ++acceptCount[dim + 1 - chosenMove.length];
            }
            else
            {
                manifold = oldManifold;
                
                // Make sure to undo any changes to list of unused vertices
                if(chosenMove.length == 1)
                {
                    unusedVertices.popBack;
                }
                else if(chosenMove.length == manifold.dimension + 1)
                {
                    unusedVertices ~= vertexToRemember;
                }
            }
        }

        // --------------------------- MAKE REPORT ----------------------------
        if(tryCount[].sum % triesPerReport == 0)
        {
            writeln;
            '-'.repeat(80).writeln;   
            writeln(" MOVE  :  DONE  /  TRIED  |  msec/move : ",
                timer.peek.msecs / real(triesPerReport));
            '-'.repeat(80).writeln;   
            iota(dim + 1).each!(indx => writeln(
                indx + 1, " -> ", dim + 1 - indx, " : ", 
                acceptCount[indx], " / ", tryCount[indx]));
            
            writeln("TOTALS : ", acceptCount[].sum, " / ", tryCount[].sum);
            '-'.repeat(80).writeln;   
            writeln("fVector              : ", manifold.fVector);
            writeln("number of facets     : ", manifold.fVector.back,
                " (target: ", numFacetsTarget, ")");
            writeln("average hinge-degree : ", manifold.meanHingeDegree,
                " (target: ", meanHingeDegreeTarget, ")");
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