/// Types and functions for working with bistellar and hinge moves on a combinatorial n-manifold
module manifold_moves;

import std.algorithm, std.array, std.conv, std.random, std.range, std.sumtype, std.traits;
import unit_threaded;
import polygons, manifold, simplicial_complex, utility;

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

/** Represents a hinge-move. These moves remove the star of a hinge
H and replace it with the join of (boundary of H) and a triangulated
disk T whose boundary is the link of H. We think of T as normal to H.
**/
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

    /// Returns the new triangulated disk normal to the hinge
    auto normalDisk()()
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
    size_t triangIndx_;
}
///
@Name("HingeMove") pure @safe unittest
{
    auto mv = HingeMove!3([1,2],[3,4,5,6],0);
    mv.hinge.shouldBeSameSetAs([1,2]);
    mv.rim.shouldBeSameSetAs([3,4,5,6]);

    mv.normalDisk.map!array.array.shouldBeSameSetAs([[3,4,5],[3,5,6]]);

    mv.triangIndx_ = 1; 
    mv.normalDisk.map!array.array.shouldBeSameSetAs([[3,4,6],[4,5,6]]);
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

/******************************************************************************
Takes as input a move and a slice containing the starting fVector and modifies
the fVector as if the given move were performed
*/
auto modifyFVector(Move)(size_t[] fVector_, Move move)
{
    enum moveIsBistellar = isInstanceOf!(BistellarMove, Move);
    enum moveisHinge = isInstanceOf!(HingeMove, Move);
    static assert(moveIsBistellar || moveisHinge,
        "must be a bistellar or hinge move");
    static if(moveIsBistellar)
    {
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
    else static if(moveisHinge)
    {
        
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

bool hasValidHingeMove(Vertex, int dim)(
    ref const(Manifold!(dim, Vertex)) manifold,
    ref const(HingeMove!(dim, Vertex)) hingeMove)
{
    auto linkVertices_ = hingeMove.rim;

    // TO DO: Replace the use of triangIndex in this function with
    // the built-in normalDisk member function of HingeMove
    size_t diskIndx = hingeMove.triangIndx;

    static if (dim < 3)
    {
        return false;
    }
    else
    {
        // TO DO: Decide what to do about this magic constantS 7, 14
        // (It comes from nGonTriangs only supporting up to 7-gon.)
        auto deg = linkVertices_.walkLength.to!int;
        Unqual!Vertex[7] linkVerticesBuff = linkVertices_.staticArray!7;
        auto linkVertices = linkVerticesBuff[0 .. deg];

        auto linkEdgeBuffer = chain(linkVertices[], linkVertices.front.only)
            .slide(2).take(deg).joiner.toStackArray!(Unqual!Vertex, 14);
        auto linkEdges = linkEdgeBuffer[];
        iota(deg).each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

        assert(diskIndx < deg.nGonTriangs.walkLength);
        auto diskFacetsBuffer = deg.nGonTriangs[diskIndx]
            .joiner.map!(i => linkVertices[i])
            .toStackArray!(Unqual!Vertex, (7 - 2) * 3);
        foreach(indx; 0 .. deg - 2)
        {
            diskFacetsBuffer[][3*indx .. 3*indx + 3].sort;
        }
        auto diskFacets = diskFacetsBuffer[].chunks(3);

        // TO DO: rework this to avoid creating simplicial complex
        // (and allocating of course!)
        auto disk = SimplicialComplex!(Vertex, 2)(diskFacets);

        // None of the "internal" edges can already be in manifold
        return disk.simplices(1)
            .filter!(edge => disk.star(edge).walkLength == 2)
            .all!(edge => !manifold.contains(edge));
    }
}
// TO DO: unittests!

auto orderedHingeLinkVertices(Vertex, int dim, S)(
    const ref Manifold!(dim, Vertex) mfd,
    S hinge)
if (isIRof!(S, const(Vertex)))
{
    auto lnk = mfd.link(hinge).map!array.array;
    auto rim = lnk.front.array;
    auto hingeDeg = mfd.degree(hinge);

    while (rim.length < hingeDeg)
    {
        auto lastVert = rim[$-2];
        auto thisVert = rim[$-1];
        auto lnkLnk = SimplicialComplex!Vertex(lnk).link([thisVert])
            .map!array.array;
        if (lnkLnk[0][0] == lastVert)
        {
            rim ~= lnkLnk[1][0];
        }
        else
        {
            rim ~= lnkLnk[0][0];
        }
    }
    return rim;
}


/*******************************************************************************
Returns a list of all the hinge moves in this manifold.
*/
HingeMove!(dim, Vertex)[] allHingeMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    static if(dim < 2)
    {
        return [];
    }
    else
    {
        alias MV = HingeMove!(dim, Vertex);
        alias SimpComp = SimplicialComplex!Vertex; 
        MV[] result;
        foreach(hinge; mfd.simplices(dim-2))
        {
            auto hingeDeg = mfd.degree(hinge); 
            if (hingeDeg > 3 && hingeDeg <= 7)
            {
                auto lnk = mfd.link(hinge).map!array.array;
                auto rim = lnk.front.array;
                while (rim.length < hingeDeg)
                {
                    auto lastVert = rim[$-2];
                    auto thisVert = rim[$-1];
                    auto lnkLnk = SimpComp(lnk).link([thisVert])
                        .map!array.array;
                    if (lnkLnk[0][0] == lastVert)
                    {
                        rim ~= lnkLnk[1][0];
                    }
                    else
                    {
                        rim ~= lnkLnk[0][0];
                    }
                }

                foreach(diskIndx; numNgonTriangs(hingeDeg).iota)
                {
                    auto move = MV(hinge, rim, diskIndx);
                    if (mfd.hasValidHingeMove(move))
                    {
                        result ~= move;
                    }
                }
            }
        }
        return result;        
    }
}

///
unittest
{
    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]];  
    auto twoPts = [[6], [7]];
    auto mfd = Manifold!3(productUnion(octahedron, twoPts));

    // 48 = (num hinges with moves) * (moves per hinge)
    //    = (12 deg4 edges in octahedron + 2*6 deg4 edges of
    //       the form [pt in twoPts, pt in octahedron])
    //      * 2 triangulations of each orthogonal disk
    //    = 24 * 2
    mfd.allHingeMoves.length.shouldEqual(48);

    // All of these moves are valid (i.e. not blocked by an
    // existing simplex)
    mfd.allHingeMoves.length.shouldEqual(48);

    // Spot check done by hand.
    alias MV = HingeMove!3;
    mfd.allHingeMoves.shouldBeSameSetAs([
        // Moves where center is edge in octaheron
        MV([0, 1],[2, 6, 4, 7],0), MV([0, 1],[2, 6, 4, 7],1),
        MV([0, 2],[1, 6, 3, 7],0), MV([0, 2],[1, 6, 3, 7],1),
        MV([0, 3],[2, 6, 4, 7],0), MV([0, 3],[2, 6, 4, 7],1),
        MV([0, 4],[3, 6, 1, 7],0), MV([0, 4],[3, 6, 1, 7],1),
        MV([1, 2],[0, 6, 5, 7],0), MV([1, 2],[0, 6, 5, 7],1),
        MV([1, 4],[0, 6, 5, 7],0), MV([1, 4],[0, 6, 5, 7],1),
        MV([1, 5],[2, 6, 4, 7],0), MV([1, 5],[2, 6, 4, 7],1),
        MV([2, 3],[0, 6, 5, 7],0), MV([2, 3],[0, 6, 5, 7],1),
        MV([2, 5],[1, 6, 3, 7],0), MV([2, 5],[1, 6, 3, 7],1),
        MV([3, 4],[0, 6, 5, 7],0), MV([3, 4],[0, 6, 5, 7],1),
        MV([3, 5],[2, 6, 4, 7],0), MV([3, 5],[2, 6, 4, 7],1),
        MV([4, 5],[3, 6, 1, 7],0), MV([4, 5],[3, 6, 1, 7],1),
        // Moves where center is edge with a vertex in twoPts
        MV([0, 6],[1, 2, 3, 4],0), MV([0, 6],[1, 2, 3, 4],1),
        MV([1, 6],[0, 2, 5, 4],0), MV([1, 6],[0, 2, 5, 4],1),
        MV([2, 6],[0, 1, 5, 3],0), MV([2, 6],[0, 1, 5, 3],1),
        MV([3, 6],[0, 2, 5, 4],0), MV([3, 6],[0, 2, 5, 4],1),
        MV([4, 6],[0, 3, 5, 1],0), MV([4, 6],[0, 3, 5, 1],1),
        MV([5, 6],[1, 2, 3, 4],0), MV([5, 6],[1, 2, 3, 4],1),
        MV([0, 7],[1, 2, 3, 4],0), MV([0, 7],[1, 2, 3, 4],1),
        MV([1, 7],[0, 2, 5, 4],0), MV([1, 7],[0, 2, 5, 4],1),
        MV([2, 7],[0, 1, 5, 3],0), MV([2, 7],[0, 1, 5, 3],1),
        MV([3, 7],[0, 2, 5, 4],0), MV([3, 7],[0, 2, 5, 4],1),
        MV([4, 7],[0, 3, 5, 1],0), MV([4, 7],[0, 3, 5, 1],1),
        MV([5, 7],[1, 2, 3, 4],0), MV([5, 7],[1, 2, 3, 4],1),
     ]);
}

bool hasValidMove(int dim, Vertex, Move)(const ref Manifold!(dim, Vertex) manifold, Move move)
{
    alias BM = BistellarMove!(dim, Vertex);
    alias HM = HingeMove!(dim, Vertex);
    return move.match!(
        (BM bistellarMove) => !manifold.contains(bistellarMove.coCenter),
        (HM hingeMove) => manifold.hasValidHingeMove(hingeMove));
}
