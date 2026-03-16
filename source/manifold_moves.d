/// Types and functions for working with bistellar moves on a combinatorial n-manifold
module manifold_moves;

import std.algorithm, std.array, std.conv, std.random, std.range, std.traits;
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
pure /* nothrow */ @safe unittest
{
    auto mv = BistellarMove!4([1,2,3],[4,5,6]);
    mv.center.shouldBeSameSetAs([1,2,3]);
    mv.coCenter.shouldBeSameSetAs([4,5,6]);
}
///
pure @system unittest
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
auto modifyFVector(T, Move)(T[] fVector_, Move move)
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
pure @safe unittest
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
pure @safe unittest
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

/******************************************************************************
Represents a 4-4 hinge move on a 3-manifold. Replaces the star of a degree-4
edge (4 tetrahedra) with an alternative triangulation of the same region
(also 4 tetrahedra) using a different diagonal of the quadrilateral link.

This is not a Pachner move (involves 6 vertices, not 5) but preserves the
f-vector.
**/
struct HingeMove(Vertex = int)
{
    Vertex[2] removedEdge;  // sorted
    Vertex[2] addedEdge;    // sorted
    Vertex[4] linkCycle;    // vertices in cyclic order around removedEdge

    /// Returns the 4 old facets (tetrahedra containing the removed edge).
    /// Each tet is the join of removedEdge with a consecutive pair from linkCycle.
    Vertex[4][4] oldFacets()() const
    {
        Vertex[4][4] result;
        foreach (i; 0 .. 4)
        {
            auto next = (i + 1) % 4;
            Vertex[4] tet = [removedEdge[0], removedEdge[1],
                             linkCycle[i], linkCycle[next]];
            tet[].sort();
            result[i] = tet;
        }
        return result;
    }

    /// Returns the 4 new facets (tetrahedra containing the added edge).
    /// The added edge diagonal splits the quadrilateral link into two triangles;
    /// each triangle is joined with each endpoint of removedEdge.
    Vertex[4][4] newFacets()() const
    {
        // Determine which diagonal: 0 connects cycle[0,2], 1 connects cycle[1,3]
        Vertex[3][2] triangles;
        Vertex[2] diag0 = [min(linkCycle[0], linkCycle[2]),
                           max(linkCycle[0], linkCycle[2])];
        if (addedEdge == diag0)
        {
            triangles[0] = [linkCycle[0], linkCycle[1], linkCycle[2]];
            triangles[1] = [linkCycle[0], linkCycle[2], linkCycle[3]];
        }
        else
        {
            triangles[0] = [linkCycle[0], linkCycle[1], linkCycle[3]];
            triangles[1] = [linkCycle[1], linkCycle[2], linkCycle[3]];
        }

        Vertex[4][4] result;
        foreach (t; 0 .. 2)
            foreach (e; 0 .. 2)
            {
                Vertex[4] tet = [removedEdge[e],
                                 triangles[t][0], triangles[t][1], triangles[t][2]];
                tet[].sort();
                result[t * 2 + e] = tet;
            }
        return result;
    }

    /// Returns the HingeMove that undoes this one.
    HingeMove inverse()() const
    {
        HingeMove inv;
        inv.removedEdge = addedEdge;
        inv.addedEdge = removedEdge;

        // Compute the link cycle of addedEdge in the post-move configuration.
        // The inverse cycle interleaves removedEdge endpoints with the
        // "other" pair of link cycle vertices (those not in addedEdge).
        Vertex[2] diag0 = [min(linkCycle[0], linkCycle[2]),
                           max(linkCycle[0], linkCycle[2])];
        if (addedEdge == diag0)
            inv.linkCycle = [removedEdge[0], linkCycle[1],
                             removedEdge[1], linkCycle[3]];
        else
            inv.linkCycle = [removedEdge[0], linkCycle[0],
                             removedEdge[1], linkCycle[2]];
        return inv;
    }

    string toString()() const
    {
        return "hinge move: remove " ~ removedEdge[].to!string
            ~ " add " ~ addedEdge[].to!string
            ~ " cycle " ~ linkCycle[].to!string;
    }
}
///
pure @safe unittest
{
    alias HM = HingeMove!int;

    // Construct a hinge move on edge [0,1] with link cycle [2,4,3,5], diagonal 1
    HM hm;
    hm.removedEdge = [0, 1];
    hm.addedEdge = [4, 5];
    hm.linkCycle = [2, 4, 3, 5];

    // Old facets: join [0,1] with consecutive pairs of cycle
    auto old = hm.oldFacets;
    old[].sort();
    int[4][4] expectedOld = [[0,1,2,4], [0,1,3,4], [0,1,3,5], [0,1,2,5]];
    expectedOld[].sort();
    assert(old == expectedOld);

    // New facets: diagonal [4,5] splits quad [2,4,3,5] into triangles
    auto nf = hm.newFacets;
    nf[].sort();
    int[4][4] expectedNew = [[0,2,4,5], [1,2,4,5], [0,3,4,5], [1,3,4,5]];
    expectedNew[].sort();
    assert(nf == expectedNew);
}
///
pure @safe unittest
{
    alias HM = HingeMove!int;

    // Test diagonal 0: addedEdge connects cycle[0] and cycle[2]
    HM hm;
    hm.removedEdge = [0, 1];
    hm.addedEdge = [2, 3];
    hm.linkCycle = [2, 4, 3, 5];

    auto nf = hm.newFacets;
    nf[].sort();
    // Triangles: {2,4,3} and {2,3,5}
    int[4][4] expectedNew = [[0,2,3,4], [1,2,3,4], [0,2,3,5], [1,2,3,5]];
    expectedNew[].sort();
    assert(nf == expectedNew);
}
///
pure @safe unittest
{
    alias HM = HingeMove!int;

    // Test inverse().inverse() roundtrip
    HM hm;
    hm.removedEdge = [0, 1];
    hm.addedEdge = [4, 5];
    hm.linkCycle = [2, 4, 3, 5];

    auto inv = hm.inverse;
    assert(inv.removedEdge == [4, 5]);
    assert(inv.addedEdge == [0, 1]);

    // inv.inverse should recover the original old/new facets (possibly reordered)
    auto roundtrip = inv.inverse;
    auto origOld = hm.oldFacets;  origOld[].sort();
    auto rtOld = roundtrip.oldFacets;  rtOld[].sort();
    assert(origOld == rtOld);
    auto origNew = hm.newFacets;  origNew[].sort();
    auto rtNew = roundtrip.newFacets;  rtNew[].sort();
    assert(origNew == rtNew);
}

/// Check that the added edge doesn't already exist in the manifold.
bool hasValidHingeMove(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const ref HingeMove!Vertex move)
{
    return !manifold.contains(move.addedEdge[]);
}
