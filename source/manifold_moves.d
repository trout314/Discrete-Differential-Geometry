/// Types and functions for working with bistellar and hinge moves on a combinatorial n-manifold
module manifold_moves;

import std.conv : to;
import utility : isInputRangeOf, isInputRangeOfInputRangeOf, throwsWithMsg, replaceEmptyLiteral;
import std.range : array, empty, walkLength;
import std.algorithm : copy, map, setIntersection;
import unit_threaded : Name, shouldBeSameSetAs, shouldEqual, shouldBeTrue, shouldBeFalse, writelnUt;
import polygons : nGonTriangs;

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

/// Represents a hinge-move
struct HingeMove(int dim, Vertex = int, int maxHingeDeg=7)
{
public:   
    /// Returns the move's hinge
    inout(Vertex)[] hinge()() inout
    {
        return vertices[0..dim-1];
    }

    /// Returns the move's rim
    inout(Vertex)[] rim()() inout
    {
        return vertices[dim-1..dim-1+rimLength];
    }

    /// Returns the new triangulated surface normal to the hinge
    auto normalSurface()()
    {
        return rimLength.nGonTriangs[triangIndx].map!(
            facet => facet.map!(indx => rim[indx]));
        // return simplicialComplex!2([[1,2,3],[2,3,4]]);
    }

    /// Returns a string representation of the move
    string toString()() const
    {
        return "hinge move, hinge=" ~ hinge.to!string 
            ~ " rim=" ~ rim.to!string
            ~ " triangIndx=" ~ triangIndx.to!string;
    }

    // Returns the index of the normal surface triangulation 
    ulong triangIndx()() const
    {
        return triangIndx_;
    }

    this(R1, R2)(R1 hingeVertices, R2 rimVertices, int triangIndx = -1)
    {
        auto hinge_ = replaceEmptyLiteral!Vertex(hingeVertices);
        auto rim_ = replaceEmptyLiteral!Vertex(rimVertices);

        static assert(isIRof!(typeof(hinge_), Vertex));
        static assert(isIRof!(typeof(rim_), Vertex));

        assert(hinge_.walkLength == dim-1, "hinge must have dim - 1 vertices");
        assert(rim_.walkLength >= 4, "rim must have at least 4 vertices");
        assert(rim_.walkLength <= maxHingeDeg, "rim must have at most "
            ~ maxHingeDeg.to!string ~ " vertices");
        assert(setIntersection(hinge_, rim_).empty,
            "hinge and rim must not have any common vertices");


        this.rimLength = rim_.walkLength;
        this.triangIndx_ = triangIndx;
        copy(hinge_, this.hinge);
        copy(rim_, this.rim);
    }
private:
    Vertex[(dim-1) + maxHingeDeg] vertices;
    size_t rimLength;

    // indicates which triangulation of a disk to use
    // (where the boundary of this disk is the coCenter)
    int triangIndx_;
}
///
@Name("HingeMove") pure @safe unittest
{
    auto mv = HingeMove!3([1,2],[3,4,5,6],0);

    mv.hinge.shouldBeSameSetAs([1,2]);
    mv.rim.shouldBeSameSetAs([3,4,5,6]);
    mv.normalSurface.map!array.array.shouldBeSameSetAs([[3,4,5],[3,5,6]]);

}
///
@Name("HingeMove (errors)") pure @system unittest
{
    HingeMove!4([1,2],[3,4,5,6]).throwsWithMsg(
        "hinge must have dim - 1 vertices");

    HingeMove!4([1,2,3],[4,5,6]).throwsWithMsg(
        "rim must have at least 4 vertices");

    HingeMove!4([1,2,3],[3,4,5,6]).throwsWithMsg(
        "hinge and rim must not have any common vertices");

    enum maxRim = 5;
    HingeMove!(4, int, maxRim)([1,2,3],[4,5,6,7,8,9,10]).throwsWithMsg(
        "rim must have at most 5 vertices");
}
