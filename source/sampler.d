import std.algorithm : all, each, map, sum, max, maxElement, joiner;
import std.range : array, iota;
import std.conv : to;
import std.random : choice, back, uniform01, rndGen;
import std.range : array, chain, empty, front, popBack, repeat;
import manifold_small : SmallManifold, pachnerMoves, doPachner;
import utility : subsetsOfSize;

import std.math : exp;
import std.format : format;
import std.stdio : writeln, writefln;

import simplicial_complex : fVector;
import simplicial_complex_algorithms : eulerCharacteristic;

import std.datetime : StopWatch, AutoStart, Duration, msecs;

// static this()
// {
//     rndGen.seed(8379);
// }

void sample()
{
    auto timer = StopWatch(AutoStart.yes);

    enum dim = 3;
    enum triesPerReport = 10;

    immutable numFacetsTarget = 100;
    real numFacetsCoef = 0.1;

    real meanHingeDegreeTarget = 4.8;
    real numHingesCoef = 0.1; 

    immutable maxVertices = 100;

     // attemptCount[j] counts j + 1 -> dim + 1 - j moves
     ulong[dim + 1] tryCount;
     ulong[dim + 1] acceptCount;

    auto manifold = SmallManifold!dim((dim + 2).iota.subsetsOfSize(dim + 1));
    auto oldManifold = SmallManifold!dim(manifold.facets);

    auto unusedVertices = iota(dim + 2,maxVertices).array;

    auto fVec = manifold.fVector;
    auto oldObjective =
            numFacetsCoef * (fVec[dim] - real(numFacetsTarget))^^2
            + numHingesCoef * (fVec[dim-2] - dim*(dim + 1)*fVec[dim] / (2 * meanHingeDegreeTarget))^^2;

    assert(unusedVertices.all!(v => !manifold.contains([v])));
  
    while(!unusedVertices.empty)
    {
        assert(unusedVertices.all!(v => !manifold.contains([v])));

        fVec = manifold.fVector;
        oldObjective =
            numFacetsCoef * (fVec[dim] - real(numFacetsTarget))^^2
            + numHingesCoef * (fVec[dim-2] - dim*(dim + 1)*fVec[dim] / (2 * meanHingeDegreeTarget))^^2;

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

        fVec = manifold.fVector;
        auto objective =
            numFacetsCoef * (fVec[dim] - real(numFacetsTarget))^^2
            + numHingesCoef * (fVec[dim-2] - dim*(dim + 1)*fVec[dim] / (2 * meanHingeDegreeTarget))^^2;

        if(objective < oldObjective)
        {
            oldManifold = SmallManifold!dim(manifold.facets);
            ++acceptCount[dim + 1 - chosenMove.length];
        }
        else
        {
            auto acceptProb = exp(oldObjective - objective);
            assert(acceptProb >= 0.0);
            assert(acceptProb <= 1.0);

            if(uniform01 <= acceptProb)
            {
                oldManifold = SmallManifold!dim(manifold.facets);
                ++acceptCount[dim + 1 - chosenMove.length];
            }
            else
            {
                manifold = SmallManifold!dim(oldManifold.facets);
                objective = oldObjective;

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
            fVec = manifold.fVector;

            writeln;
            "----------------------------------".writeln;   
            " MOVE  :  DONE  /  TRIED".writeln;
            "----------------------------------".writeln;   
            iota(dim + 1).each!(indx => writeln(
                indx + 1, " -> ", dim + 1 - indx, " : ", 
                acceptCount[indx], " / ", tryCount[indx]));
            
            writeln("totals : ", acceptCount[].sum, " / ", tryCount[].sum);
            "----------------------------------".writeln;   
            writeln("objective : ", objective);
            writeln("fVector   : ", fVec);
            writeln("avg h-deg : ", real(fVec[$-1])/fVec[$-3]*((dim+1)*dim)/2);
            "----------------------------------".writeln;   
            writeln("msec/move : ", timer.peek.msecs / real(tryCount[].sum));

        }
    }
}