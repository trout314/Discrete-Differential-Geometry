import std.algorithm : all, each, map, sum, max, maxElement, joiner;
import std.range : array, iota;
import std.conv : to;
import std.random : choice, back, uniform01, rndGen;
import std.range : array, empty, front, popBack, repeat;
import manifold_small : SmallManifold, pachnerMoves, doPachner;
import utility : subsetsOfSize;

import std.math : exp;
import std.format : format;
import std.stdio : writeln, writefln;

import simplicial_complex : fVector;

// static this()
// {
//     rndGen.seed(8379);
// }

void sample()
{
    enum dim = 2;
    enum triesPerReport = 1;
    immutable numFacetsTarget = 100;
    real nFacetCoef = 0.1;
    immutable maxVertices = 100;

     // attemptCount[j] counts j + 1 -> dim + 1 - j moves
     ulong[dim + 1] tryCount;
     ulong[dim + 1] acceptCount;

    auto manifold = SmallManifold!dim((dim + 2).iota.subsetsOfSize(dim + 1));
    auto oldManifold = SmallManifold!dim(manifold.facets);

    auto unusedVertices = iota(dim + 2,maxVertices).array;
    auto objectiveBefore = nFacetCoef * (manifold.numFacets - real(numFacetsTarget))^^2;

    assert(unusedVertices.all!(v => !manifold.contains([v])));
  
    while(!unusedVertices.empty)
    {
        assert(unusedVertices.all!(v => !manifold.contains([v])));

        objectiveBefore = nFacetCoef * (manifold.numFacets - real(numFacetsTarget))^^2;
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

        auto objectiveAfter = nFacetCoef * (manifold.numFacets - real(numFacetsTarget))^^2;
        
        if(objectiveAfter < objectiveBefore)
        {
            oldManifold = SmallManifold!dim(manifold.facets);
            ++acceptCount[dim + 1 - chosenMove.length];
        }
        else
        {
            auto acceptProb = exp(objectiveBefore - objectiveAfter);
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
            auto indxWidth = "%,s".format(dim + 1).length;

            auto acceptWidth = "%,s".format(acceptCount[].sum).length;
            auto tryWidth = "%,s".format(tryCount[].sum).length;
            auto totalMoveWidth = max(
                2*indxWidth + acceptWidth + tryWidth + 10, 21);

            writeln;            
            "-".repeat(totalMoveWidth).joiner.writeln;
            "%*s : %*s / %*s".writefln(
                2*indxWidth + 4, "MOVES",
                acceptWidth, "DONE",
                tryWidth, "TRIED");
            "-".repeat(totalMoveWidth).joiner.writeln;

            foreach(indx; 0 .. dim + 1)
            {
                  "%*s -> %*s : %*s / %*s".writefln(
                    indxWidth, indx + 1,
                    indxWidth, dim + 1 - indx,
                    acceptWidth, "%,s".format(acceptCount[indx]),
                    tryWidth, "%,s".format(tryCount[indx]));
            }
          
            "%*s : %*,s / %*,s".writefln(
                2*indxWidth + 4, "total",
                acceptWidth, acceptCount[].sum,
                tryWidth, tryCount[].sum);       
        }
    }
}