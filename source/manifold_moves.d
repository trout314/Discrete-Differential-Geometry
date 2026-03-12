/// Types and functions for working with bistellar moves on a combinatorial n-manifold
module manifold_moves;

import std.algorithm, std.array, std.conv, std.random, std.range, std.traits;
import unit_threaded;
import manifold, simplicial_complex, utility;

alias isIRof = isInputRangeOf;
alias isIRofIRof = isInputRangeOfInputRangeOf;

/**
Represents a bistellar move, also called a Pachner move. See:
Udo Pachner, "P.L. Homeomorphic Manifolds are Equivalent by Elementary Shellings",
European Journal of Combinatorics, Volume 12, Issue 2, 1991, Pages 129-145,
ISSN 0195-6698, https://doi.org/10.1016/S0195-6698(13)80080-7
**/
struct BistellarMove(int dim, Vertex=int)
{
public:
    /// Constructs a bistellar move with specified center and co-center
    this(R1, R2)(R1 centerVertices, R2 coCenterVertices)
    {
        auto center_ = replaceEmptyLiteral!Vertex(centerVertices);
        auto coCenter_ = replaceEmptyLiteral!Vertex(coCenterVertices);

        static assert(isIRof!(typeof(center_), Vertex));
        static assert(isIRof!(typeof(coCenter_), Vertex));

        assert(center_.walkLength > 0, "center must have at least one vertex");
        assert(coCenter_.walkLength > 0, "co-center must have at least one vertex");
        assert(center_.walkLength + coCenter_.walkLength == dim + 2,
            "total number of vertices in center and co-center must be dim + 2");
        assert(setIntersection(center_, coCenter_).empty,
            "center and co-center must not have any common vertices");

        this.lenCenter = center_.walkLength;
        copy(center_, this.center);
        copy(coCenter_, this.coCenter);
    }

    /// Returns the center of the move
    inout(Vertex)[] center()() inout
    {
        return vertices[0..lenCenter];
    }

    /// Returns the co-center of the move
    inout(Vertex)[] coCenter()() inout
    {
        return vertices[lenCenter..$];
    }

    /// Returns a string representation of the move
    string toString()() const
    {
        return "bistellar move, center=" ~ center.to!string 
            ~ " coCenter=" ~ coCenter.to!string;
    }

    alias dimension = dim;
private:
    size_t lenCenter;
    Vertex[dim + 2] vertices;
}
///
@Name("BistellarMove") pure /* nothrow */ @safe unittest
{
    auto mv = BistellarMove!4([1,2,3],[4,5,6]);
    mv.center.shouldBeSameSetAs([1,2,3]);
    mv.coCenter.shouldBeSameSetAs([4,5,6]);
}
///
@Name("BistellarMove (errors)") pure @system unittest
{
    BistellarMove!4([1,2,3],[4,5]).throwsWithMsg(
        "total number of vertices in center and co-center must be dim + 2");

    BistellarMove!4([1,2,3],[4,5,6,7]).throwsWithMsg(
        "total number of vertices in center and co-center must be dim + 2");

    BistellarMove!4([1,2,3,4,5,6],[]).throwsWithMsg(
        "co-center must have at least one vertex");

    BistellarMove!4([],[1,2,3,4,5,6]).throwsWithMsg(
        "center must have at least one vertex");

    BistellarMove!4([1,2,3],[3,4,5]).throwsWithMsg(
        "center and co-center must not have any common vertices");
}

/******************************************************************************
Takes as input a bistellar move and a slice containing the starting fVector
and modifies the fVector as if the given move were performed
*/
auto modifyFVector(Move)(size_t[] fVector_, Move move)
{
    static assert(isInstanceOf!(BistellarMove, Move),
        "must be a bistellar move");

    // We modify the fVector for the removal of the original star
    auto centerDim = move.center.length - 1;
    auto dim = fVector_.length.to!int - 1;
    foreach(d; centerDim .. dim + 1)
    {
        fVector_[d] -= binomial(dim + 1 - centerDim, d - centerDim);
    }

    // Now modify it for the addition of the new star
    auto coDim = dim + 1 - move.center.length;
    foreach(d; coDim .. dim + 1)
    {
        fVector_[d] += binomial(dim + 1 - coDim, d - coDim);
    }
}
///
@Name("modifyFVector (dim 3)") pure @safe unittest
{
    alias BMove = BistellarMove!3;

    size_t[] fVector = [0,0,0,0];

    /* A 1 -> 4 move should give net: +1 vertices, +4 edges, +6 triangles,
    +3 tetrahdra */
    fVector.modifyFVector(BMove([1,2,3,4],[5]));
    fVector.shouldEqual([1, 4, 6, 3]);

    /* A 2 -> 3 move should give net: +0 vertices, +1 edges, +2 triangles,
    +1 tetrahedra */
    fVector.modifyFVector(BMove([1,2,3],[4,5]));
    fVector.shouldEqual([1,5,8,4]);

    fVector.modifyFVector(BMove([4,5],[1,2,3]));    // 3 -> 2 move
    fVector.modifyFVector(BMove([5],[1,2,3,4]));    // 4 -> 1 move

    // Should be back where we started
    fVector.shouldEqual([0,0,0,0]);
}

///
@Name("modifyFVector (dim 2)") pure @safe unittest
{
    alias BMove = BistellarMove!2;
    
    size_t[] fVector = [0,0,0];

    // A 2 -> 2 move should leave the f-vector unchanged 
    fVector.modifyFVector(BMove([1,2],[3,4]));
    fVector.shouldEqual([0, 0, 0]);

    // A 1 -> 3 move should give net: +1 vertices, +3 edges, +2 triangles
    fVector.modifyFVector(BMove([1,2,3],[4]));
    fVector.shouldEqual([1, 3, 2]);

    // Doing a 3 -> 1 move should return the f-vector to its original value
    fVector.modifyFVector(BMove([4],[1,2,3]));
    fVector.shouldEqual([0, 0, 0]);
}

bool hasValidMove(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const ref BistellarMove!(dim, Vertex) move)
{
    return !manifold.contains(move.coCenter);
}
