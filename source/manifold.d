module manifold;

import algorithms : eulerCharacteristic, is2Sphere, isCircle,
    isConnected, isPureOfDim, join;
import simplicial_complex : fVector, simplicialComplex, SimplicialComplex,
    assertValidSimplex, loadSimplicialComplex;
import std.algorithm : all, any, copy, canFind, each, equal, filter, find,
    joiner, map, maxElement, setDifference, merge, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown;
import std.range : array, back, chain, chunks, cycle, ElementType, empty,
    enumerate, front, iota, isInputRange, isOutputRange, only, popBack,
    popFront, put, replace, retro, save, slide, take, walkLength, zip;
import unit_threaded : Name, shouldEqual, shouldBeSameSetAs, shouldBeEmpty, should;
import utility : binomial, isInputRangeOf, isInputRangeOfInputRangeOf, isSubsetOf, numNgonTriangs, nGonTriangs, productUnion, SmallMap,
    StackArray, staticIota, subsets, subsetsOfSize, throwsWithMsg, toStackArray, swapPop;
import std.stdio : File, writeln;
import std.typecons : Flag, No, Yes;
import std.math : isClose;

import std.stdio : writeln, writefln, write;

import std.traits : Unqual;
import std.random : uniform, choice;

import moves : Move_ = Move;

import std.array : staticArray;

import std.algorithm : countUntil, copy, map, min;

import manifold_examples : trigonalBipyramid, octahedron, standardSphere;

alias isIRof = isInputRangeOf;
alias isIRofIRof = isInputRangeOfInputRangeOf;

//dfmt off
/*******************************************************************************
Manifold type... TO DO: More info here.

Manifold models a compact combinatorial n-manifold. Supports doing pachner moves
quickly
*/
struct Manifold(int dimension_, Vertex_ = int)
{
private:
    alias SimpComp = SimplicialComplex!Vertex_;
    SimpComp simpComp;

    alias NSimplex = StackArray!(Vertex, dimension + 1);
    size_t[NSimplex] degreeMap;
    
    alias Ridge = Vertex[dimension_];
    alias RidgeLink = StackArray!(Vertex, 2);
    RidgeLink[Ridge] ridgeLinks;

    size_t[dimension_ + 1] numSimplices;
    ulong[dimension - 1] totSqrDegrees;

    //--------------- Helper Functions ---------------
    static NSimplex toNSimp(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.toStackArray!(Vertex_, dimension + 1, R);
    }

    static Ridge toRidge(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.staticArray!dimension;
    }

    alias Facet = Vertex[dimension_ + 1];
    static Facet toFacet(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.staticArray!(dimension_ + 1);
    }
    size_t numValidMoves_;
    Move[] moves_;
    size_t[Vertex[]] indxOfCenter;
    size_t[][Vertex[]] indicesOfCoCenter;
public:
    /// Dimension of the manifold
    static immutable dimension = dimension_;

    /// Vertex type used in the manifold
    alias Vertex = Vertex_;

    /// Move type used for the manifold
    alias Move = Move_!(dimension_, Vertex_);

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    /// We can initialize the manifold from an input range of input ranges
    /// of vertices
    this(F)(F initialFacets) if(isIRofIRof!(F, const(Vertex)))
    {
        initialFacets.each!(f => this.insertFacet(f));
        numSimplices[] = simpComp.fVector[];
        foreach(d; 0 .. dimension - 1)
        {
            totSqrDegrees[d] = simplices(d).map!(s => this.degree(s)^^2).sum;
        }

        foreach(simp, deg; degreeMap)
        {
            if(simp[].walkLength < dimension + 1)
            {
                if(deg == dimension + 2 - simp[].walkLength)
                {
                    auto coCenter = this.findCoCenter(simp[]);
                    if(!this.contains(coCenter))
                    {
                        ++numValidMoves_;
                    }
                    moves_ ~= Move(simp[], coCenter.sort);
                }
            }
        }

        foreach(indx, mv; moves_.enumerate)
        {
            this.indxOfCenter[mv.center.idup] = indx;
            this.indicesOfCoCenter[mv.coCenter.idup] ~= indx;
        }
    
        // this.findProblems.shouldBeEmpty;
    }

    this(this) pure @safe
    {
        ridgeLinks = ridgeLinks.dup;
        degreeMap = degreeMap.dup;
        moves_ = moves_.dup;
    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        bool notInDegreeMap = (toNSimp(simplex) in this.degreeMap) is null;
        assert(this.simpComp.contains(simplex) == !notInDegreeMap,
            "simplices in degreeMap and internal simplicial complex disagree");
        return !notInDegreeMap;
    }

    /***************************************************************************
    Returns the degree of a simplex in the manifold.
    */
    size_t degree(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        assert(this.contains(simplex),
            "called degree on a simplex not in the manifold");
        assert(degreeMap[toNSimp(simplex)] == simpComp.star(simplex).walkLength,
            "degree in degreeMap and internal simplicial complex disagree");
        return degreeMap[toNSimp(simplex)];
    }

    /***************************************************************************
    Returns an array containing the number of simplices in each dimension.
    This is called the "fVector" of the manifold.
    */
    const(size_t)[] fVector()() const 
    {
        assert(numSimplices[] == this.simpComp.fVector);
        return numSimplices[];
    }

    // Special version of insertFacet to update tracked info
    // (EXCEPT numSimplices, which is updated after each pachner move)
    private void insertFacet(F)(F facet_) if (isIRof!(F, Vertex))
    {       
        assert(facet_.walkLength == dimension + 1,
            "facet has wrong dimension");

        auto facetBuffer = toFacet(facet_);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);
        assert(toNSimp(facet) !in degreeMap,
            "tried to insert a facet already in the manifold");

        this.simpComp.insertFacet!(No.checkForFacetFaces)(facet);

        foreach(simplex_; facet.subsets)
        {
            auto simplex = toNSimp(simplex_);
            ++degreeMap[simplex];

        
            if(simplex.length <= dimension - 1)
            {
                int i = simplex.length.to!int - 1;
                assert(i >= 0);
                assert(i < totSqrDegrees.length);
                totSqrDegrees[i] += 2*degreeMap[simplex] - 1;
            }
        
            if(simplex.length == dimension)
            {
                auto oppVerts = facet.filter!(v => !simplex_.canFind(v));
                assert(oppVerts.walkLength == 1);
                auto oppVert = oppVerts.front;
               
                auto ridge = toRidge(simplex_);
                auto ptrToLink = ridge in ridgeLinks;
                if(!ptrToLink)
                {
                    assert(simplex in degreeMap);
                    assert(degreeMap[simplex] == 1);
                    auto rL = RidgeLink();
                    rL ~= oppVert;
                    ridgeLinks[ridge] = rL;
                }
                else
                {
                    assert(simplex in degreeMap);
                    assert(degreeMap[simplex] == 2);
                    (*ptrToLink) ~= oppVert;
                    (*ptrToLink)[].sort;
                }
                assert(ridge in ridgeLinks);
            }
        }
        assert(toNSimp(facet) in degreeMap);
    }

    // Special version of removeFacet to update tracked info.
    private void removeFacet(F)(F facet_) if (isIRof!(F, Vertex))
    {
        assert(facet_.walkLength == dimension + 1);

        auto facetBuffer = toFacet(facet_);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);      
        assert(toNSimp(facet) in degreeMap);

        this.simpComp.removeFacet(facet);
        
        foreach(simplex_; facet.subsets)
        {
            auto simplex = toNSimp(simplex_);
            assert(simplex in degreeMap);

            --degreeMap[simplex];

            if(simplex.length <= dimension - 1)
            {
                totSqrDegrees[simplex.length - 1] -= 2*degreeMap[simplex] + 1;
            }

            if((simplex.length == dimension) && (degreeMap[simplex] == 1))
            {
                auto ridge = toRidge(simplex_);
                assert(ridge in ridgeLinks);

                auto oppVerts = facet.filter!(v => !simplex_.canFind(v));
                assert(oppVerts.walkLength == 1);
                auto oppVert = oppVerts.front;

                assert(ridgeLinks[ridge][].canFind(oppVert));

                if(ridgeLinks[ridge].length == 2)
                {
                    if(ridgeLinks[ridge][0] == oppVert)
                    {
                        ridgeLinks[ridge][0] = ridgeLinks[ridge][1];
                    }
                    else
                    {
                        assert(ridgeLinks[ridge][1] == oppVert);
                    }
                    ridgeLinks[ridge].length = 1;
                }
                else
                {
                    assert(ridgeLinks[ridge].length == 1);
                }

            }

            if(degreeMap[simplex] == 0)
            {
                degreeMap.remove(simplex);
                if(simplex.length == dimension)
                {
                    auto ridge = toRidge(simplex_);
                    assert(ridgeLinks[ridge].length == 1);
                    ridgeLinks.remove(ridge);
                }
            }

        }
        assert(toNSimp(facet) !in degreeMap);
        assert(!this.simpComp.containsFacet(facet));
    }

    /// We provide access to the manifold as a simplicial complex
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const pure nothrow @nogc @safe 
    {
        return simpComp; 
    }

    alias asSimplicialComplex this;

    const(Move)[] moves()() const
    {
        return moves_[];
    }

    size_t numValidMoves()() const
    {
        return numValidMoves_;
    }

}

/*******************************************************************************
Returns a list of all the possible pachner moves except for the 1->(dim+1) moves
that are valid in this manifold. (Note that 1->(dim+1) moves are always valid.)
*/
Move_!(dim, Vertex)[] computePachnerMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    auto result = mfd.computeMBPMoves.filter!(
        mv => !mfd.contains(mv.coCenter)).array;
    return result;
}


/*******************************************************************************
Returns a list of all the (maybe blocked) hinge moves in this manifold.
*/
Move_!(dim, Vertex)[] computeMBHMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    alias MV = Move_!(dim, Vertex); 
    MV[] result;
    foreach(hinge; mfd.simplices(dim-2))
    {
        auto hingeDeg = mfd.degree(hinge); 
        if (hingeDeg >= 3 && hingeDeg <= 7)
        {
            auto lnk = mfd.link(hinge).map!array.array;
            auto coCen = lnk.front.array;
            while(coCen.length < hingeDeg)
            {
                auto lastVert = coCen[$-2];
                auto thisVert = coCen[$-1];
                auto lnkLnk = mfd.SimpComp(lnk).link([thisVert])
                    .map!array.array;
                if (lnkLnk[0][0] == lastVert)
                {
                    coCen ~= lnkLnk[1][0];
                }
                else
                {
                    coCen ~= lnkLnk[0][0];
                }
            }

            foreach(i; numNgonTriangs(hingeDeg).iota)
            {
                result ~= MV(hinge, coCen, i);
            }
        }
    }
    return result;
}

/*******************************************************************************
Returns a list of all the valid hinge moves in this manifold.
*/
Move_!(dim, Vertex)[] computeHingeMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    return mfd.computeMBHMoves.filter!(mv => 
        mfd.hasValidHingeMove(mv.coCenter, mv.triangIndx)).array;

    // alias MV = Move_!(dim, Vertex); 
    // alias SC = SimplicialComplex!(Vertex, dim);
    // MV[] result;

    // foreach(mv; mfd.computeMBHMoves)
    // {
    //     auto deg = mv.coCenter.length.to!int;
    //     auto diskIndx = mv.triangIndx;
    //     auto diskFacets = deg.nGonTriangs[diskIndx];

    //     if (mfd.hasValidHingeMove(mv.coCenter, diskIndx))
    //     {
    //         result ~= mv;
    //     }
    // }
    // return result;
}

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
    mfd.computeMBHMoves.length.shouldEqual(48);

    // All of these moves are valid (i.e. not blocked by an
    // existing simplex)
    mfd.computeHingeMoves.length.shouldEqual(48);

    // Spot check done by hand.
    alias MV = mfd.Move;
    mfd.computeMBHMoves.shouldBeSameSetAs([
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

///
@Name("computePachnerMoves") pure @safe unittest
{
    static foreach(d; 2 .. 8)
    {
        {
            auto m = standardSphere!d;
            m.computePachnerMoves.shouldBeEmpty;
        }
    }

    // trigonal bipyramid
    auto tb = Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,4]]);
    tb.computePachnerMoves.shouldBeSameSetAs([
        Move_!2([0],[1,2,3]),
        Move_!2([4],[1,2,3]),
        Move_!2([1,2],[0,4]),
        Move_!2([1,3],[0,4]),
        Move_!2([2,3],[0,4])
    ]);

    // octahedron
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    oct.computePachnerMoves.shouldBeSameSetAs([
        Move_!2([0,1],[2,4]),
        Move_!2([0,2],[1,3]),
        Move_!2([0,3],[2,4]),
        Move_!2([0,4],[1,3]),
        Move_!2([1,5],[2,4]),
        Move_!2([2,5],[1,3]),
        Move_!2([3,5],[2,4]),
        Move_!2([4,5],[1,3]),
        Move_!2([1,2],[0,5]),
        Move_!2([2,3],[0,5]),
        Move_!2([3,4],[0,5]),
        Move_!2([1,4],[0,5])
    ]);

    auto m = Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,5],
        [2,4,5],[3,4,5]]);    
    m.computePachnerMoves.shouldBeSameSetAs([
        Move_!2([0],[1,2,3]),
        Move_!2([5],[2,3,4]),
        Move_!2([1,2],[0,4]),
        Move_!2([1,3],[0,4]),
        Move_!2([2,3],[0,5]),
        Move_!2([2,4],[1,5]),
        Move_!2([3,4],[1,5])
    ]);

    // two-point suspension over boundary of 3-simplex
    auto tps = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],
        [1,2,3,5],[1,2,4,5],[1,3,4,5],[2,3,4,5]]);
}


/*******************************************************************************
Returns a list of all the (perhaps blocked) pachner moves except for the
1->(dim+1) moves in this manifold. Note that 1->(dim+1) moves are always valid.
*/
Move_!(dim, Vertex)[] computeMBPMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    Move_!(dim, Vertex)[] result;
    foreach(simp_, deg; mfd.degreeMap)
    {
        if(simp_.length < dim + 1)
        {
            auto simp = simp_[];
            if(deg == mfd.dimension + 2 - simp.walkLength)
            {
                auto coCenter = mfd.findCoCenter(simp);
                result ~= Move_!(dim, Vertex)(simp, coCenter);
            }
        }
    }
    return result;
}

///
@Name("Manifold doc tests") pure @safe unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));

    octahedron.fVector.shouldEqual([6UL,12,8]);
    octahedron.eulerCharacteristic.shouldEqual(2);

    auto tetrahedron = Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    octahedron.star([1,2]).shouldBeSameSetAs([[0,1,2], [1,2,5]]);    

    octahedron.computePachnerMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[0,1], [0,2], [0,3], [0,4], [1,2], [1,4],  // 1-simplices
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5]]);      

    tetrahedron.computePachnerMoves.shouldBeEmpty;
    
    // TO DO: FINISH CHECKS!
    // octahedron.doPachner([1,2]);
    // octahedron.doPachner([0,5]);
}
// NOTE: The following unittest cannot be @safe since throwsWithMsg 
// catches an Error
///
@Name("Manifold (errors)") pure @system unittest
{
    Manifold!2([[1,2,3,4]]).throwsWithMsg("facet has wrong dimension");

    auto sphere3 = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];    
    Manifold!3(chain(sphere3, sphere3)).throwsWithMsg(
        "tried to insert a facet already in the manifold");
    
}

///
@Name("computePachnerMoves") pure @safe unittest
{
    auto m = Manifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    m.computePachnerMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[1], [5], [2, 3], [2, 4], [3, 4]]);

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.computePachnerMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[2, 3], [0, 1], [1, 5], [4, 5], [0, 3], [1, 4],
        [1, 2], [0, 4], [0, 2], [2, 5], [3, 5], [3, 4]]);
}

// NOTE: The following unittest cannot be @safe since throwsWithMsg 
// catches an Error
///
@Name("doPachner (errors)") pure @system unittest
{
    // Can't do 2->2 move on the boundary of a 3-simplex
    auto m = Manifold!2([[1,2,3],[1,2,4], [1,3,4], [2,3,4]]);   
    m.doPachner([1,2], [3,4]).throwsWithMsg("coCenter of move in manifold");
}
///
@Name("doPachner") pure unittest
{
    auto m = standardSphere!2;
    m.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    m.doPachner([1,2,3], [4]);
    m.facets.shouldBeSameSetAs(
        [[0,1,2],[0,1,3], [0,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    m.doPachner([0,2,3], [7]);
    m.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,7],
        [0,3,7], [1,2,4], [1,3,4], [2,3,4], [2,3,7]]);

    m.doPachner([7],[0,2,3]);
    m.doPachner([4],[1,2,3]);
    m.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);


    auto n = standardSphere!2;
    n.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    n.doPachner([1,2,3], [4]);
    n.doPachner([1,2], [0,4]);
    n.facets.shouldBeSameSetAs(
        [[0,1,3], [0,1,4], [0,2,3], [0,2,4], [1,3,4], [2,3,4]]);
    n.doPachner([0,4], [1,2]);
    n.doPachner([4], [1,2,3]);
    n.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);      

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.doPachner([1,2],[0,5]);
    octahedron.degree([0]).shouldEqual(5);
    octahedron.degree([5]).shouldEqual(5);
    octahedron.degree([1]).shouldEqual(3);
    octahedron.degree([2]).shouldEqual(3);

    // We can undo the 2->2 move
    octahedron.doPachner([0,5], [1,2]);
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.doPachner([0,1,2], [99]);
    octahedron.doPachner([99], [0,1,2]);
}

/******************************************************************************
Do a bistellar move replacing the star(center) with star(coCenter)
TO DO: Better docs here!
*/
void doPachner(Vertex, int dim, C, D)(
    ref Manifold!(dim, Vertex) manifold,
    C center,
    D coCenter
)
if (isIRof!(C, const(Vertex)) && isIRof!(D, const(Vertex)))
{
    assert(manifold.contains(center), "center of move not in manifold");
    assert(!manifold.contains(coCenter), "coCenter of move in manifold");
    if(coCenter.walkLength > 1)
    {
        auto pm = Move_!(dim, Vertex)(center, coCenter);
        assert(manifold.computePachnerMoves.canFind(pm), "not a valid pachner move");
    }
 
    // Buffer for holding vertices in center (followed by coCenter)
    const(Vertex)[dim + 2] cBuffer = chain(center, coCenter)
        .staticArray!(dim + 2);
 
    immutable cenLen = center.walkLength;
    immutable coCenLen = (dim + 2) - cenLen;
    assert(coCenLen == coCenter.walkLength);

    auto center_ = cBuffer[0 .. cenLen];
    auto coCenter_ = cBuffer[cenLen .. $];

    immutable centerDim = cast(int) cenLen - 1;
    immutable coCenterDim = cast(int) coCenLen - 1;

    auto oldPiece = productUnion(coCenter_.subsetsOfSize(coCenterDim), center_.only);
    auto newPiece = productUnion(center_.subsetsOfSize(centerDim), coCenter_.only);

    alias SC = SimplicialComplex!(Vertex, dim);
    alias MFD = Manifold!(dim, Vertex);

    assert(SC(oldPiece).isPureOfDim(dim));
    assert(SC(newPiece).isPureOfDim(dim));
    assert(MFD(chain(oldPiece, newPiece)).numFacets == dim + 2);
    assert(manifold.star(center).map!array.array.sort
        .equal!equal(oldPiece.map!array.array.sort));

    oldPiece.each!(f => manifold.removeFacet(f));
    newPiece.each!(f => manifold.insertFacet(f));
    manifold.numSimplices.modifyFVector(cenLen);

    auto pm = Move_!(dim, Vertex)(center_, coCenter_);
    manifold.modifyMoveDataOnMove(pm, oldPiece, newPiece);
}

auto findCoCenter(Vertex, int dim, C)(
    const ref Manifold!(dim, Vertex) manifold,
    C center
)
if (isIRof!(C, const(Vertex)))
{
    assert(manifold.contains(center));
    return manifold.link(center).joiner.array.dup.sort.uniq.array;
}

/*******************************************************************************
Returns the coCenter for the Pachner move with given center. For efficiency we
also need to know a facet with face center.
*/
auto coCenter(Vertex, int dim, C, F)(
    const ref Manifold!(dim, Vertex) mfd,
    C center,
    F facet
)
if (isIRof!(C, const(Vertex)) && isIRof!(F, const(Vertex)))
{
    assert(mfd.contains(facet));
    assert(mfd.contains(center));
    assert(center.isSubsetOf(facet));
    
    // The coCenter of a facet is a new vertex not in the manifold.
    assert(center.walkLength < dim + 1);

    // TO DO: Clean this up!
    auto ridges = facet.subsetsOfSize(dim)
        .filter!(r => center.isSubsetOf(r)).map!(r => mfd.toRidge(r));
    auto coCenterVerts = ridges.map!(r => mfd.ridgeLinks[r][])
        .joiner.array.dup.sort.uniq.array;

    assert(coCenterVerts.equal(mfd.findCoCenter(center.array)));
    return coCenterVerts;
}

// TO DO: Separate unittesting for findCoCenter


///
@Name("facets(dim) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge02 = [0,2];
    int[2] edge12 = [1,2];

    () pure nothrow @nogc @safe {
        auto facetsRange = sc.facets(1);
        auto savedRange = facetsRange.save;

        int[2][3] facs;
        facs[0] = facetsRange.front;
        facetsRange.popFront;
        facs[1] = facetsRange.front;
        facetsRange.popFront;
        facs[2] = facetsRange.front;
        facetsRange.popFront;
        
        facs[].sort();

        assert(facs[0] == edge01);
        assert(facs[1] == edge02);
        assert(facs[2] == edge12);

        assert(facetsRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("facets() (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge02 = [0,2];
    int[2] edge12 = [1,2];

    () pure nothrow @nogc @safe {
        auto facetsRange = sc.facets;
        auto savedRange = facetsRange.save;

        int[2][3] facs;
        facs[0] = facetsRange.front;
        facetsRange.popFront;
        facs[1] = facetsRange.front;
        facetsRange.popFront;
        facs[2] = facetsRange.front;
        facetsRange.popFront;
        
        facs[].sort();

        assert(facs[0] == edge01);
        assert(facs[1] == edge02);
        assert(facs[2] == edge12);

        assert(facetsRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("star(range) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge12 = [1,2];
    int[1] vertex1 = [1];

    () pure nothrow @nogc @safe {
        auto starRange = sc.star(vertex1[]);
        auto savedRange = starRange.save;

        int[2][2] edges;
        edges[0][] = starRange.front[];
        starRange.popFront;
        edges[1][] = starRange.front[];
        starRange.popFront;

        edges[].sort();        
        assert(edges[0] == edge01);
        assert(edges[1] == edge12);

        assert(starRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("link(range) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    immutable(int[1]) v = [1];

    () pure nothrow @nogc @safe {
        auto linkRange = sc.link(v[]);
        auto savedRange = linkRange.save;

        int[2] vertices;
        vertices[0] = linkRange.front.front;
        linkRange.popFront;
        vertices[1] = linkRange.front.front;
        linkRange.popFront;

        vertices[].sort();        
        assert(vertices[0] == 0);
        assert(vertices[1] == 2);

        assert(linkRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("contains (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge07 = [0,7];

    () pure nothrow @nogc @safe {
        assert(sc.contains(edge01[]));
        assert(!sc.contains(edge07[]));
    }();  
}

/******************************************************************************
Modifies an fVector for a pachner move with center simplex that contains `n`
vertices.
*/
auto modifyFVector(size_t[] fVector_, size_t centerLength)
{
    assert(centerLength >= 1, "center length must be at least one");
    assert(fVector_.length >= centerLength, "fVector must have length at least the center length");

    // We modify the fVector for the removal of the original star
    auto centerDim = centerLength - 1;
    auto dim = fVector_.length.to!int - 1;
    foreach(d; centerDim .. dim + 1)
    {
        fVector_[d] -= binomial(dim + 1 - centerDim, d - centerDim);
    }

    // Now modify it for the addition of the new star
    auto coDim = dim + 1 - centerLength;
    foreach(d; coDim .. dim + 1)
    {
        fVector_[d] += binomial(dim + 1 - coDim, d - coDim);
    }
}
///
@Name("modifyFVector") pure @safe unittest
{
    size_t[4] fVec = [0,0,0,0];

    /* A 1 -> 4 move should give net: +1 vertices, +4 edges, +6 triangles,
    +3 tetrahdra */
    fVec.modifyFVector(4);
    fVec.shouldEqual([1, 4, 6, 3]);

    /* A 2 -> 3 move should give net: +0 vertices, +1 edges, +2 triangles,
    +1 tetrahedra */
    fVec.modifyFVector(3);
    fVec.shouldEqual([1,5,8,4]);

    fVec.modifyFVector(2);    // 3 -> 2 move
    fVec.modifyFVector(1);    // 4 -> 1 move

    // Should be back where we started
    fVec.shouldEqual([0,0,0,0]);

    size_t[3] fVec2 = [0,0,0];
    fVec2.modifyFVector(2); // 2 -> 2 move
    fVec2.shouldEqual([0,0,0]);

    fVec2.modifyFVector(3); // 1 -> 3 move

    // a 1 -> 3 move should give net: +1 vertices, +3 edges, +2 triangles
    fVec2.shouldEqual([1, 3, 2]);

    fVec2.modifyFVector(1); // 3 -> 1 move
    fVec2.shouldEqual([0,0,0]);
}

/******************************************************************************
Returns a range containing the valid pachner moves whose center is a face of
the given facet
*/
auto movesAtFacet(Vertex, int dim, F)(
    const ref Manifold!(dim, Vertex) mfd,
    F facet)
if (isIRof!(F, const(Vertex)))
{
    return facet.subsets.map!array.filter!(
        center => (center.walkLength == dim + 1) 
            || (mfd.degree(center) == dim + 2 - center.walkLength 
                && !mfd.contains(mfd.coCenter(center, facet))));
}
///
@Name("movesAtFacet") pure @safe unittest
{
    /* If the manifold is the boundary of a simplex (i.e. a sphere with the
    minimum number of facets) then only the type 1 -> (dim + 1) Pachner moves
    are valid. */    
    foreach(dim; staticIota!(1,4))
    {
        immutable sphere = standardSphere!dim;
        foreach(f; sphere.facets)
        {
            sphere.movesAtFacet(f).shouldBeSameSetAs([f.array]);
        }
    }

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.movesAtFacet([0,1,2]).shouldBeSameSetAs([
        [0,1], [0,2], [1,2],    // 2 -> 2 moves
        [0,1,2]                 // 1 -> 3 move
    ]);
    
    foreach(f; octahedron.facets)
    {
        octahedron.movesAtFacet(f).shouldBeSameSetAs(
            chain([f.array], f.subsetsOfSize(2).map!array).array);
    }
    
    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);
    pyramid.movesAtFacet([0,1,2]).shouldBeSameSetAs([
        [0],                    // 3 -> 1 move
        [1,2],                  // 2 -> 2 move
        [0,1,2]                 // 1 -> 3 move
    ]);
    pyramid.movesAtFacet([0,2,3]).shouldBeSameSetAs([
        [0],                    // 3 -> 1 move
        [2,3],                  // 2 -> 2 move
        [0,2,3]                 // 1 -> 3 move
    ]);
    pyramid.movesAtFacet([1,3,4]).shouldBeSameSetAs([
        [4],                    // 3 -> 1 move
        [1,3],                  // 2 -> 2 move
        [1,3,4]                 // 1 -> 3 move
    ]);
}

/******************************************************************************
* Returns a manifold loaded from the file specified by fileName. If fileName
* is the empty string, the returned manifold is the standard sphere.
*/
Manifold!(dim, Vertex) loadManifold(int dim, Vertex = int)(string fileName)
{
    SimplicialComplex!Vertex sc = loadSimplicialComplex!Vertex(fileName);
    return Manifold!(dim, Vertex)(sc.facets);
}

///
@Name("loadManifold") @system unittest
{
    auto m = loadManifold!2(
        "data/manifold_sampler_unittest_load.dat");
    auto expected = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5],
        [0, 5, 6], [0, 1, 6], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
    assert(m == expected);

    assertThrown(loadManifold!2(
        "data/manifold_sampler_unittest_bad_load.dat"));
}

/******************************************************************************
Saves a manifold to a file specified by fileName.
*/
void saveTo(int dimension, Vertex = int)(
    const ref Manifold!(dimension, Vertex) mfd,
    string fileName)
{
    auto saveFile = File(fileName, "w"); // Open in write-only mode
    saveFile.write(mfd.asSimplicialComplex);
}

///
@Name("saveTo") @system unittest
{
    auto fileName = "data/manifold_sampler_unittest_save.dat";
    auto sphere = standardSphere!4;
    sphere.saveTo(fileName);
    Manifold!4 loadedSphere = loadManifold!4(fileName);
    assert(loadedSphere == sphere);
}

/******************************************************************************
Returns the mean degree of simplices of the given dimension 'dim'
*/
real meanDegree(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfdDim);
    immutable nSimps = mfd.fVector[dim];
    immutable nFacets = mfd.fVector.back;
    immutable simpsPerFacet = binomial(mfdDim + 1, dim + 1);
    return simpsPerFacet * nFacets / real(nSimps);
}
///
@Name("meanDegree") pure @safe unittest
{
    foreach(dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        assert((dim + 1).iota.all!(d => m.meanDegree(d) == (dim - d + 1)));
    }

    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    // pyramid has two vertices with deg(v)=3 and three vertices with deg(v)=4
    assert(pyramid.meanDegree(0).isClose((2 * 3 + 3 * 4)/5.0));
    assert(pyramid.meanDegree(1).isClose(2.0));
    assert(pyramid.meanDegree(2).isClose(1.0));
}

/******************************************************************************
Returns total (degree(s))^^2 for all simplices s of given dimension 'dim'.
*/
ulong totalSquareDegree(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfd.dimension);
    if(dim == mfdDim - 1)
    {
        // All codimension-1 simplices (ridges) have degree 2
        return 4 * mfd.numSimplices[dim];
    }
    else if(dim == mfdDim)
    {
        // All facets have degree 1
        return mfd.numSimplices[dim];
    }
    return mfd.totSqrDegrees[dim];
}
///
@Name("totalSquareDegree") pure @safe unittest
{
    foreach(dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        foreach(d; 0 .. dim + 1)
        {
            // (dim + 2) choose (d + 1) faces of dimension d
            // each with degree (dim - d + 1)
            m.totalSquareDegree(d).shouldEqual(
                binomial(dim + 2, d + 1) * (dim - d + 1)^^2);
        }
    }

    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    pyramid.totalSquareDegree(0).shouldEqual(2 * 3^^2 + 3 * 4^^2);
    pyramid.totalSquareDegree(1).shouldEqual(9 * 2^^2);
    pyramid.totalSquareDegree(2).shouldEqual(6 * 1^^2);
}

/******************************************************************************
Returns the variance in the degree of simplices of the given dimension 'dim'
*/
real degreeVariance(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfdDim);
    if(dim >= mfdDim - 1)
    {
        return 0;   // No variance in ridge and facet dimension
    }

    immutable meanDeg = mfd.meanDegree(dim);
    immutable meanSqrDeg = mfd.totalSquareDegree(dim) / real(mfd.fVector[dim]);
    return meanSqrDeg - meanDeg^^2;
}
///
@Name("degreeVariance") pure @safe unittest
{
    foreach(dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        assert((dim+1).iota.all!(d => m.degreeVariance(d) == 0));
    }

    // pyramid has two vertices with deg(v)=3 and three vertices with deg(v)=4
    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    immutable real md0 = pyramid.meanDegree(0);
    assert(pyramid.degreeVariance(0)
        .isClose((2 * (3 - md0)^^2 + 3 * (4 - md0)^^2)/5.0));

    assert(pyramid.degreeVariance(1).isClose(0.0));
    assert(pyramid.degreeVariance(2).isClose(0.0));

}

/******************************************************************************
Does the 'diskIndx'-th hinge move associated at 'hinge'. Must give the
link of this hinge as `hingeLink`
*/
void doHingeMove(Vertex, int dim, H, K)(
    ref Manifold!(dim, Vertex) manifold,
    H hinge_,
    K linkVertices_,
    int diskIndx)
if (isIRof!(H, const(Vertex)) && isIRof!(K, const(Vertex)))
{
    // TO DO: Clean this up!
    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);

    static assert(dim >= 3,
        "no hinge moves in dimension less than 3");
    auto hingeBuffer = hinge_.staticArray!(dim - 1);
    auto hinge = hingeBuffer[];

    // TO DO: Decide what to do about this magic constant 7
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    Unqual!Vertex[7] linkVerticesBuff = linkVertices_.staticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    auto linkEdgeBuffer = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.toStackArray!(Unqual!Vertex, 14);
    auto linkEdges = linkEdgeBuffer[];
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

    assert(diskIndx < deg.nGonTriangs.walkLength);
    auto diskFacets = deg.nGonTriangs[diskIndx];

    auto oldPiece = productUnion(hinge.only, linkEdges.chunks(2));
    auto newPiece = productUnion(
        hinge.subsetsOfSize(dim - 2),
        diskFacets.map!(f => f.map!(i => linkVertices[i]).array.dup.sort));

    alias SC = SimplicialComplex!(Vertex, dim);
    assert(SC(newPiece).isPureOfDim(dim));   
    assert(SC(oldPiece).isPureOfDim(dim));
    assert(manifold.star(hinge).map!array.array.sort
        .equal!equal(oldPiece.map!array.array.sort));

    foreach(f; oldPiece)
    {
        manifold.removeFacet(f);    
    }
    // oldPiece.each!(f => manifold.removeFacet(f));

    // Modify fVector, start with simplices removed
    manifold.numSimplices[dim] -= deg;      // (# 1-simplices in Lk(hinge)) * (1 hinge) 
    manifold.numSimplices[dim - 1] -= deg;  // (# 0-simplices in Lk(hinge)) * (1 hinge)
    manifold.numSimplices[dim - 2] -= 1;    // 1 hinge

    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);

    newPiece.each!(f => manifold.insertFacet(f));

    // Now, add in the new simplices to the stored fVector
    manifold.numSimplices[1] += (deg - 3);  // # edges in disk
    manifold.numSimplices[2] += (deg - 2);  // # triangles in disk

    // (# 0-simplices in bdry(hinge)) * (# 1-simplices in disk)
    manifold.numSimplices[2] += (dim - 1) * (deg - 3);

    foreach(d; 3 .. dim)
    {
        // (# (d - 3)-simplices in bdry(hinge)) * (# 2-simplices in disk)
        manifold.numSimplices[d] += binomial(dim - 1, d - 2) * (deg - 2);            

        // (# (d - 2)-simplices in bdry(hinge)) * (# 1-simplices in disk)
        manifold.numSimplices[d] += binomial(dim - 1, d - 1) * (deg - 3);
    }

    // (# (dim - 2)-simplices in bdry(hinge)) * (# 2-simplices in disk)
    manifold.numSimplices[dim] += (dim - 1) * (deg - 2);

    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);
}


/******************************************************************************
Does the 'diskIndx'-th hinge move associated at 'hinge'. Must give the
link of this hinge as `hingeLink`
*/
void undoHingeMove(Vertex, int dim, H, K)(
    ref Manifold!(dim, Vertex) manifold,
    H hinge_,
    K linkVertices_,
    int diskIndx)
if (isIRof!(H, const(Vertex)) && isIRof!(K, const(Vertex)))
{
    // TO DO: Clean this up!
    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);

    static assert(dim >= 3,
        "no hinge moves in dimension less than 3");
    auto hingeBuffer = hinge_.staticArray!(dim - 1);
    auto hinge = hingeBuffer[];

    // TO DO: Decide what to do about this magic constant 7
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    Unqual!Vertex[7] linkVerticesBuff = linkVertices_.staticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    auto linkEdgeBuffer = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.toStackArray!(Unqual!Vertex, 14);
    auto linkEdges = linkEdgeBuffer[];
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

    assert(diskIndx < deg.nGonTriangs.walkLength);

    auto diskFacetsBuffer = deg.nGonTriangs[diskIndx]
        .joiner.map!(i => linkVertices[i])
        .toStackArray!(Unqual!Vertex, (7 - 2) * 3);
    foreach(indx; 0 .. deg - 2)
    {
        diskFacetsBuffer[][3*indx .. 3*indx + 3].sort;
    }

    auto diskFacets = diskFacetsBuffer[].chunks(3);

    auto newPiece = productUnion(hinge.only, linkEdges.chunks(2));
    auto oldPiece = productUnion(hinge.subsetsOfSize(dim - 2), diskFacets);

    // All of the disk should be in manifold
    assert(diskFacets.all!(f => manifold.contains(f)));
    alias SC = SimplicialComplex!(Vertex, dim);
    assert(SC(newPiece).isPureOfDim(dim));   
    assert(SC(oldPiece).isPureOfDim(dim));

    oldPiece.each!(f => manifold.removeFacet(f));

    // STUFF REMOVED
    manifold.numSimplices[1] -= (deg - 3);  // # edges in disk
    manifold.numSimplices[2] -= (deg - 2);  // # triangles in disk
    
    // (# 0-simplices in bdry(hinge)) * (# 1-simplices in disk)
    manifold.numSimplices[2] -= (dim - 1) * (deg - 3);

    foreach(d; 3 .. dim)
    {
        // (# (d - 3)-simplices in bdry(hinge)) * (# 2-simplices in disk)
        manifold.numSimplices[d] -= binomial(dim - 1, d - 2) * (deg - 2);            

        // (# (d - 2)-simplices in bdry(hinge)) * (# 1-simplices in disk)
        manifold.numSimplices[d] -= binomial(dim - 1, d - 1) * (deg - 3);
    }

    // (# (dim - 2)-simplices in bdry(hinge)) * (# 2-simplices in disk)
    manifold.numSimplices[dim] -= (dim - 1) * (deg - 2);

    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);

    newPiece.each!(f => manifold.insertFacet(f));

    // STUFF ADDED
    manifold.numSimplices[dim] += deg;      // (# 1-simplices in Lk(hinge)) * (1 hinge) 
    manifold.numSimplices[dim - 1] += deg;  // (# 0-simplices in Lk(hinge)) * (1 hinge)
    manifold.numSimplices[dim - 2] += 1;    // 1 hinge

    assert(manifold.fVector == manifold.asSimplicialComplex.fVector);
}

///
@Name("hinge moves") pure /* @safe */ unittest
{
    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]];
    
    auto twoPts = [[6], [7]];
    auto mfd = Manifold!3(productUnion(octahedron, twoPts));

    assert(mfd.hasValidHingeMove([1,2,3,4], 1));
    mfd.doHingeMove([0,6], [1,2,3,4], 1);
    mfd.facets.shouldBeSameSetAs([[1, 4, 5, 7], [0, 1, 2, 7], [1, 4, 5, 6],
        [0, 2, 3, 7], [3, 4, 5, 7], [0, 3, 4, 7], [3, 4, 5, 6], [0, 1, 4, 7],
        [1, 2, 5, 6], [1, 2, 5, 7], [2, 3, 5, 6], [2, 3, 5, 7],[0, 1, 2, 4],
        [0, 2, 3, 4], [1, 2, 4, 6], [2, 3, 4, 6]]);

    mfd.undoHingeMove([0,6], [1,2,3,4], 1);
    mfd.facets.shouldBeSameSetAs(productUnion(octahedron, twoPts).map!array);

    auto twoOtherPts = [[8], [9]];
    auto mfd4 = Manifold!4(
        productUnion(productUnion(octahedron, twoPts), twoOtherPts));

    mfd4.doHingeMove([0,6,8], [1,2,3,4], 1);
    mfd4.undoHingeMove([0,6,8], [1,2,3,4], 1);
    mfd4.facets.shouldBeSameSetAs(productUnion(
        productUnion(octahedron, twoPts), twoOtherPts).map!array);    
}

bool hasValidHingeMove(Vertex, int dim, K)(
    ref const(Manifold!(dim, Vertex)) manifold,
    K linkVertices_,
    int diskIndx)
if (isIRof!(K, const(Vertex)))
{
    static assert(dim >= 3,
        "no hinge moves in dimension less than 3");

    // TO DO: Decide what to do about this magic constantS 7, 14
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    Unqual!Vertex[7] linkVerticesBuff = linkVertices_.staticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    auto linkEdgeBuffer = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.toStackArray!(Unqual!Vertex, 14);
    auto linkEdges = linkEdgeBuffer[];
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

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
// TO DO: unittests!


@Name("fVector") @system unittest
{
    auto hyperbolicDodecahedral = loadManifold!3(
        "data/manifold_sampler_unittest_dodecahedral.dat");

    () pure @safe {
        // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
        assert(standardSphere!7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
                9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
                (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

        // TO DO: Get reference for this...
        assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);
    } ();
}

string[] findProblems(Vertex, int dim)(const ref Manifold!(dim, Vertex) mfd)
{
    alias SC = SimplicialComplex!Vertex;
    string[] problems;

    // TO DO: Create a findProblems function for simplicial complexes
    // and check that here

    if(!mfd.simpComp.isPureOfDim(dim))
    {
        problems ~= "not all facets have the correct dimension";
    }
    
    if(!mfd.simpComp.isConnected)
    {
        problems ~= "facets do not define a connected simplicial complex";
    }

    static if(dim >= 1)
    {
        if(!mfd.simplices(dim - 1).all!(
            s => mfd.simpComp.star(s).walkLength == 2))
        {
            problems ~= "found a ridge with degree not equal to 2";
        }
    }

    static if(dim >= 2)
    {
        if(!mfd.simplices(dim - 2).all!(s => SC(mfd.simpComp.link(s)).isCircle))
        {
            problems ~= "found a hinge whose link is not a circle";
        }
    }

    static if(dim >= 3)
    {
        if(!mfd.simplices(dim - 3).all!(s => SC(mfd.simpComp.link(s)).is2Sphere))
        {
            problems ~= "found a codimension-3 simplex whose link is not a 2-sphere";
        }
    }

    if(mfd.numSimplices[] != mfd.simpComp.fVector[])
    {
        problems ~= "number of simplices incorrect";        
    }

    foreach(d; 0 .. dim - 1)
    {
        auto correct = mfd.simplices(d).map!(s => mfd.degree(s)^^2).sum;
        if(mfd.totSqrDegrees[d] != correct)
        {
            problems ~= "found incorrect total squared degree";
            break;
        }
    }

    foreach(d; 0 .. dim + 1)
    {
        foreach(s; mfd.simpComp.simplices(d))
        {
            if(mfd.toNSimp(s) !in mfd.degreeMap)
            {
                problems ~= "found a simplex in simpComp that is not in degreeMap";
                goto done;
            }
        }
    }
    done:

    foreach(simplex; mfd.degreeMap.byKey)
    {
        if(!mfd.simpComp.contains(simplex[]))
        {
            problems ~= "found a simplex in degreeMap that is not in simpComp";
            break;
        }
    }

    foreach(pair; mfd.degreeMap.byKeyValue)
    {
        auto simplex = pair.key[];
        auto deg = pair.value;

        if(deg != mfd.simpComp.star(simplex).walkLength)
        {
            problems ~= "found a simplex in degreeMap with incorrect degree";
            break;
        }
    }

    foreach(pair; mfd.ridgeLinks.byKeyValue)
    {
        auto ridge = pair.key[];
        auto link = pair.value[];

        if(link.walkLength != mfd.simpComp.star(ridge).walkLength)
        {
            problems ~= "found a ridge in ridgeLinks whose link has incorrect number of vertices";
            break;
        }
    }

    foreach(pair; mfd.ridgeLinks.byKeyValue)
    {
        auto ridge = pair.key[];
        auto link = pair.value[];

        // TO DO: Clean this up!
        if(link.array != mfd.simpComp.link(ridge).joiner.array.dup.sort.array)
        {
            problems ~= "found a ridge in ridgeLinks whose link has the wrong vertices";
            break;
        }
    }

    foreach(ridge; mfd.ridgeLinks.byKey)
    {
        if(!mfd.simpComp.contains(ridge[]))
        {
            problems ~= "found a ridge in ridgeLinks that is not in simpComp";
            break;
        }
    }

    foreach(ridge; mfd.simpComp.simplices(dim - 1))
    {
        if(mfd.toRidge(ridge) !in mfd.ridgeLinks)
        {
            problems ~= "found a ridge in simpComp that is not in ridgeLinks";
            break;
        }
    }

    return problems;
}
///
@Name("findProblems") pure @safe unittest
{
    auto m3 = standardSphere!3;
    assert(m3.findProblems.empty);

    m3 = Manifold!3([[1,2,3,4]]);
    m3.findProblems.shouldBeSameSetAs([
        "found a hinge whose link is not a circle",
        "found a ridge with degree not equal to 2",
        "found a codimension-3 simplex whose link is not a 2-sphere"
    ]);

    // These facets separately define 2-spheres, but share the vertex 1
    auto sphere2a = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]];
    auto sphere2b = [[1,5,6], [1,5,7], [1,6,7], [5,6,7]];
    auto m2 = Manifold!2(chain(sphere2a, sphere2b));
    m2.findProblems.shouldBeSameSetAs([
        "found a hinge whose link is not a circle"
    ]);

    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5],
        [3,4,5], [1,4,5]];
    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];
    m2 = Manifold!2(chain(octahedron, sphere));
    m2.findProblems.shouldBeSameSetAs([
        "facets do not define a connected simplicial complex"
    ]);

    // These facets separately define 2-spheres, but share the vertex 1
    auto sphere3a = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];
    auto sphere3b = [[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]];
    m3 = Manifold!3(chain(sphere3a, sphere3b));
    m3.findProblems.shouldBeSameSetAs([
        "found a codimension-3 simplex whose link is not a 2-sphere"
    ]);

    // TO DO: This test trips an internal assert. I don't want to remove
    // the internal assert. Decide what to do...

    // auto m1 = Manifold!1([[0,1],[1,2],[0,2],[0,3]]);
    // m3.findProblems.shouldBeSameSetAs([
    //     "found ...."
    // ]);
}

/******************************************************************************
Rerturns the vertices in the link. [v_1, v_2, v_3 ... v_deg] where the edges
in the link are [v_1, v_2], [v_2, v_3], ... [v_deg, v_1]
*/
auto linkVerticesAtHinge(int dim, Vertex, S, F)(
    const ref Manifold!(dim, Vertex) mfd,
    S hinge,
    F facet)
if (isIRof!(S, const(Vertex)) && isIRof!(F, const(Vertex)))
{
    assert(mfd.contains(hinge));
    assert(mfd.contains(facet));
    assert(hinge.isSubsetOf(facet));

    auto verticesFound = setDifference(facet, hinge).array;
    auto finalVertex = verticesFound.front;

    auto facet_ = facet.array;
    
    auto done = false;
    do
    {
        Vertex nextVertex;

        auto ridgeLink = mfd.ridgeLinks[
            mfd.toRidge(merge(hinge, verticesFound.back.only))][];
        assert(ridgeLink.length == 2);

        if(ridgeLink.front == verticesFound[$ - 2])
        {
            nextVertex = ridgeLink.back;
        }
        else
        {
            nextVertex = ridgeLink.front;
        }

        if (nextVertex != finalVertex)
        {
            verticesFound ~= nextVertex;
        }
        else
        {
            done = true;
        }
    }
    while (!done);

    return verticesFound;
}
///
pure @safe unittest
{
    auto m = standardSphere!3;
    assert(m.linkVerticesAtHinge([0,1],[0,1,2,3]).equal([2,3,4]));

    // TO DO: More tests...
}

size_t[] degreeHistogram(Vertex, int mfdDim)(
    auto ref const(Manifold!(mfdDim, Vertex)) mfd,
    size_t dim)
{
    auto degData = mfd.degreeMap.byKeyValue
        .filter!(p => p.key.length == dim + 1);

    size_t[] histogram;
    foreach(p; degData)
    {
        if (p.value > histogram.length)
        {
            histogram.length = p.value;
        }
        ++histogram[p.value - 1];
    }

    return histogram;
}
///
unittest
{
    // TO DO: More tests...

    auto m = standardSphere!3;
    assert(4.iota.map!(k => m.degreeHistogram(k)).equal!equal(
        [[0,0,0,5], [0,0,10], [0,10], [5]]));
}

/******************************************************************************
Choose a random hinge move. If none is possible, returns -1. Otherwise, returns
the index of the triangulation in nGonTriangs for the chosen move. The vertices
of the disk are given by linkVertices_ and these should be the link of hinge_
in the manifold `mfd`
*/
int getRandomHingeMove(Vertex, int dim, K)(
    const ref Manifold!(dim, Vertex) mfd,
    K linkVertices_)
if (isIRof!(K, const(Vertex)))
{
    auto deg = linkVertices_.walkLength.to!int;
    auto possibleDisks = deg.nGonTriangs.array;

    auto indx = uniform(0, possibleDisks.length.to!int);
    while(!mfd.hasValidHingeMove(linkVertices_, indx))
    {
        possibleDisks[indx] = possibleDisks.back;
        possibleDisks.popBack;
        if(possibleDisks.empty)
        {
            return -1;
        }
        indx = uniform(0, possibleDisks.length.to!int);
    }

    return indx;
}
///
unittest
{
    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]];
    
    auto twoPts = [[6], [7]];
    auto mfd = Manifold!3(productUnion(octahedron, twoPts));

    auto x = mfd.getRandomHingeMove([1,2,3,4]);
    assert(x != -1); // both possible hinge moves should be valid
    assert(x == 0 || x == 1); // only two moves possible
    mfd.doHingeMove([0,6], [1,2,3,4], x);
    mfd.undoHingeMove([0,6], [1,2,3,4], x);
    mfd.facets.map!array.shouldBeSameSetAs(
        productUnion(octahedron, twoPts).map!array);

    // TO DO: More tests, test for -1 return
}

void saveEdgeGraphTo(int dimension, Vertex = int)(
    const ref Manifold!(dimension, Vertex) mfd,
    string fileName)
{
    auto saveFile = File(fileName, "w"); // Open in write-only mode
    foreach(edge; mfd.asSimplicialComplex.simplices(1))
    {
        saveFile.writeln(edge.front, " ", edge.back);
    }
}

private void report(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    writeln("   current manifold: ", mfd);
    "   moves:".writeln;
    foreach(indx, move; mfd.moves_.enumerate)
    {
        write("      ", indx, ": ", move);
        if (mfd.contains(move.coCenter))
        {
            writeln(", blocked");
        }
        else
        {
            "".writeln;
        }
    }
    "   indxOfCenter:".writeln;
    foreach(cen, indx; mfd.indxOfCenter)
    {
        writeln("      ", cen, ": ", indx);            
    }
    "   indicesOfCoCenter:".writeln;
    foreach(coCen, indices; mfd.indicesOfCoCenter)
    {
        writeln("      ", coCen, ": ", indices);            
    }
    writeln("   num valid moves: ", mfd.numValidMoves_);
    writeln("   computed num valid moves: ", mfd.computePachnerMoves.length);

}

private void checkMoveData(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    // writeln("   checking move records...");
    foreach(i, mv; mfd.moves_.enumerate)
    {
        // writeln("      checking move #", i);
        if (mfd.indxOfCenter[mv.center] != i)
        {
            // writeln("bad index of center for move. actual index i=", i,
            //     " indxOfCenter[", mv.center, "]=", indxOfCenter[mv.center]);
            debug mfd.report();
            assert(0);
        }

        if (mv.coCenter !in mfd.indicesOfCoCenter || 
            !mfd.indicesOfCoCenter[mv.coCenter].canFind(i))
        {
            // writeln("bad or missing indicesOfCoCenter[", mv.coCenter, "]");
            debug mfd.report();
            assert(0);
        }
    }

    foreach(cc, ccIndxList; mfd.indicesOfCoCenter)
    {
        // writeln("      checking coCenter ", cc);
        foreach(i; ccIndxList)
        {
            if (mfd.moves_[i].coCenter != cc)
            {
                // writeln("bad coCenter indices. indicesOfCocenter[", cc, "]=", ccIndxList,
                //     " but moves[", i, "].coCenter=", moves[i].coCenter);
                debug mfd.report();
                assert(0);
            }
        }
    }

    if (mfd.numValidMoves_ != mfd.computePachnerMoves.length)
    {
        // writeln("numValidMoves_ = ", numValidMoves_, " but...");
        // writeln("this.computePachnerMoves.length = ", this.computePachnerMoves.length);
        debug mfd.report();
        assert(0);
    }

    import unit_threaded : shouldBeSameSetAs;
    import std.algorithm : each, sort;
    import std.range : array;
    mfd.moves_.array.shouldBeSameSetAs(mfd.computeMBPMoves.dup.sort);

    // report();
}

// TO DO: figure out how to make oldPiece and newPiece const
void modifyMoveDataOnMove(Vertex, int dim, S)(
    ref Manifold!(dim, Vertex) mfd,
    Move_!(dim, Vertex) mv,
    ref S oldPiece,
    ref S newPiece)
{
    // Note that this function is designed to run *after* any changes to the manifold's
    // underlying simplicial complex (and stored fVector) have already happened.

    auto cen = mv.center;
    auto coCen = mv.coCenter;

    auto cenLen = cen.walkLength.to!int;
    assert(cenLen>0);

    auto coCenLen = coCen.walkLength.to!int;
    assert(coCenLen>0);

    auto oldP = mfd.SimpComp(oldPiece);
    auto newP = mfd.SimpComp(newPiece);
    auto bdry = mfd.SimpComp(productUnion(
        cen.subsetsOfSize(cenLen-1),
        coCen.subsetsOfSize(coCenLen-1)));
    
    const(Vertex)[][] oldInterior;
    const(Vertex)[][] newInterior;
    const(Vertex)[][] boundary;
    static foreach(d; 0 .. dim + 1)
    {
        oldInterior ~= oldP.simplices(d).filter!(s => cen.isSubsetOf(s)).array;
        newInterior ~= newP.simplices(d).filter!(s => coCen.isSubsetOf(s)).array;
        boundary ~= bdry.simplices(d);
    }

    // writeln(mv);
    // writeln("      newP: ", newP);
    // writeln("      oldP: ", oldP);
    // writeln("      bdry: ", boundary);
    // writeln("      oldInterior: ", oldInterior);
    // writeln("      newInterior: ", newInterior);

    size_t[] toRemove = [];
    size_t[] toUpdate = [];
    mfd.Move[] newMoves = [];

    if (coCenLen > 1)
    {
        assert(cen in mfd.indxOfCenter, "pachner move center indx not found");
        assert(coCen in mfd.indicesOfCoCenter, "pachner move cocenter indx not found");
    }

    // writeln("   oldInterior simplices ");
    foreach(simp; oldInterior)
    {
        // writeln("      simp = ", simp);

        // any moves blocked by an old simplex are now valid
        if(simp in mfd.indicesOfCoCenter)
        {
            // writeln("         adding ",indicesOfCoCenter[simp].length,
            //      " to numValidMoves_ because of unblocked coCenter ", simp); 
            mfd.numValidMoves_ += mfd.indicesOfCoCenter[simp].length;
        }

        if(simp.walkLength < dim + 1)
        {
            auto simpMoveIndx = mfd.indxOfCenter[simp];
            // writeln("         remove index ", simpMoveIndx);
            toRemove ~= simpMoveIndx;
        }
    }

    // writeln("   newInterior simplices ");
    if (coCenLen < dim + 1)
    {
        // new center (old coCenter) is valid (unless it is trivial)
        assert(!coCen.empty);
        // writeln("      adding 1 to numValidMoves_ since new center ", coCen, " is not a facet");
        mfd.numValidMoves_++;
    }
    
    foreach(simp; newInterior)
    {
        // writeln("      simp = ", simp);

        // any moves blocked by a new simplex are now not valid
        if(simp in mfd.indicesOfCoCenter)
        {
            auto simpCoCen = mfd.findCoCenter(simp);
            // writeln("         subtracting ",indicesOfCoCenter[simp].length,
            //      " from numValidMoves_ because of newly blocked coCenter ", simp); 
            mfd.numValidMoves_ -= mfd.indicesOfCoCenter[simp].length;
        }

        if (simp.walkLength < dim + 1)
        {
            auto simpCoCen = mfd.findCoCenter(simp);
            // TO DO: use known facet to do faster coCenter lookup
            newMoves ~= mfd.Move(simp, simpCoCen);
            // writeln("         adding ", newMoves[$-1]);
        }
    }

    // writeln("   boundary simplices ");
    foreach(simp; boundary)
    {
        auto d = simp.walkLength.to!int - 1;
        assert(d >= 0);
        assert(d < dim);

        // writeln("      simp = ", simp);

        bool wasMove = (simp in mfd.indxOfCenter) !is null;
        bool wasValidMove = false;
        if(wasMove)
        {
            auto simpMoveIndx = mfd.indxOfCenter[simp];
            auto simpMove = mfd.moves_[simpMoveIndx];
            auto simpCen = simpMove.center;
            auto simpCoCen = simpMove.coCenter;
            wasValidMove = !mfd.contains(simpCoCen);
        }

        bool isMove = (mfd.degree(simp) == dim - d + 1);

        // TO DO: Use faster version of findCoCenter that takes a facet!
        bool isValidMove = isMove && !mfd.contains(mfd.findCoCenter(simp).sort);
        
        if(wasMove && !isMove)
        {
            // writeln("         remove index ", indxOfCenter[simp.array]);
            toRemove ~= mfd.indxOfCenter[simp.array];
        }
        else if (!wasMove && isMove)
        {
            newMoves ~= typeof(mv)(simp.array, mfd.findCoCenter(simp).sort);
            // writeln("         add move: ", newMoves[$-1]);
        }
        else if (wasMove && isMove)
        {
            // writeln("         update index ", indxOfCenter[simp.array]);
            toUpdate ~= mfd.indxOfCenter[simp.array];
        }

        if (wasValidMove && !isValidMove)
        {
            // writeln("         subtracting 1 from numValidMoves_ because of update at ", simp);
            --mfd.numValidMoves_;
        }
        
        if (!wasValidMove && isValidMove)
        {
            // writeln("         adding 1 from numValidMoves_ because of update at ", simp);
            ++mfd.numValidMoves_;
        }
    }
    

    toRemove.sort;
    // writeln("   indices toRemove = ", toRemove.retro);
    // writeln("   indices toUpdate = ", toUpdate);
    // writeln("   newMoves = ");
    foreach(m; newMoves)
    {
        // writeln("      ", m);
    }

    // "   updating moves...".writeln;
    foreach(i; toUpdate)
    {
        auto thisCoCen = mfd.moves_[i].coCenter;
        auto thisOldCoCen = thisCoCen.idup;
        auto thisCen = mfd.moves_[i].center;
        // writeln("      updating #", i, " ", moves_[i]);
        
        // Update vertices in cocenter of move
        Vertex[] newVerts;
        Vertex[] oldVerts;
        if (newP.contains(thisCen))
        {
            newVerts = newP.link(thisCen).joiner.array.dup.sort.uniq.array;
        }
        if (oldP.contains(thisCen))
        {
            oldVerts = oldP.link(thisCen).joiner.array.dup.sort.uniq.array;
        }

        if (oldVerts.empty)
        {
            oldVerts = newVerts;
        }

        if (newVerts.empty)
        {
            newVerts = oldVerts;
        }

        assert(newVerts.walkLength == oldVerts.walkLength);
        foreach(j; 0 .. newVerts.walkLength)
        {               
            thisCoCen = thisCoCen.replace(oldVerts[j].only, newVerts[j].only);
        }
        
        mfd.moves_[i] = typeof(mv)(thisCen, thisCoCen.sort);
        // writeln("      new move #", i, " ", moves_[i]);


        mfd.indicesOfCoCenter[mfd.moves_[i].coCenter.idup] ~= i;

        // Update indices pointing to updated cocenter
        if (thisOldCoCen in mfd.indicesOfCoCenter)
        {
            auto newIndices = mfd.indicesOfCoCenter[thisOldCoCen].filter!(k => k != i).array;
            if (!newIndices.empty)
            {
                mfd.indicesOfCoCenter[thisOldCoCen] = newIndices;
            }
            else
            {
                // writeln("      indicesOfCoCenter.remove(", thisOldCoCen, ")");
                mfd.indicesOfCoCenter.remove(thisOldCoCen);
            }
        }
    }


    // "   removing moves...".writeln;
    foreach(i; toRemove.retro)
    {
        assert(mfd.moves_.length > 0);
        auto lastIndx = mfd.moves_.length - 1;
        auto lastMove = mfd.moves_[lastIndx];
        auto goneMove = mfd.moves_[i];
        mfd.indxOfCenter[lastMove.center.idup] = i;

        // writeln("      removing index ", i, " by swapping ", i, " <-> ", lastIndx);

        if (lastMove.coCenter in mfd.indicesOfCoCenter)
        {
            // writeln("         for lastMove: ", lastMove);
            // writeln("         old indicesOfCoCenter[lastMove.coCenter]: ",
                // indicesOfCoCenter[lastMove.coCenter]);
            mfd.indicesOfCoCenter[lastMove.coCenter.idup]
                = mfd.indicesOfCoCenter[lastMove.coCenter.idup].replace(lastIndx.only, i.only).
                    sort.uniq.filter!(indx => indx != lastIndx).array;
            // writeln("         new indicesOfCoCenter[lastMove.coCenter]: ",
                // indicesOfCoCenter[lastMove.coCenter]);
        }
        else
        {
            // writeln("         old coCenterIndices: --removed--");
            // writeln("         new coCenterIndices: --removed--");
        }

        mfd.moves_.swapPop(i);
        if (goneMove.coCenter in mfd.indicesOfCoCenter)
        {
            // writeln("         for goneMove: ", goneMove);
            // writeln("         old indicesOfCoCenter[goneMove.coCenter]: ",
                // indicesOfCoCenter[goneMove.coCenter]);

            if (goneMove.coCenter != lastMove.coCenter)
            {
                mfd.indicesOfCoCenter[goneMove.coCenter.idup] = 
                    mfd.indicesOfCoCenter[goneMove.coCenter].filter!(k => k != i).array;
            }

            // writeln("         new indicesOfCoCenter[goneMove.coCenter]: ",
                // indicesOfCoCenter[goneMove.coCenter]);

            if (mfd.indicesOfCoCenter[goneMove.coCenter].all!(k => k == i))
            {
                if (lastMove.coCenter != goneMove.coCenter)
                {
                    // writeln("         indicesOfCoCenter.remove(goneMove.coCenter)");
                    mfd.indicesOfCoCenter.remove(goneMove.coCenter);
                }
            }
        }
        mfd.indxOfCenter.remove(goneMove.center);
    }

    // "   adding new moves_...".writeln;
    foreach(newMove; newMoves)
    {
        // writeln("      adding ", newMove);
        mfd.moves_ ~= newMove;
        mfd.indxOfCenter[newMove.center.idup] = mfd.moves_.length - 1;
        mfd.indicesOfCoCenter[newMove.coCenter.idup] ~= mfd.moves_.length - 1;
    }
}

@Name("overall move") unittest
{
    auto m = standardSphere!2;
    auto testMoves =  [
        [[0, 2, 3], [4]],
        [[0, 2], [1, 4]],
        [[0, 1, 4], [6]],
        [[3, 4], [0, 2]],
        [[3], [0, 1, 2]],
        [[0, 4], [2, 6]],
        [[2, 4, 6], [10]],
        [[4, 6], [1, 10]],
        [[0], [1, 2, 6]],
        [[4], [1, 2, 10]],
        [[1, 2, 6], [14]]
        ];
    foreach(mv; testMoves)
    {
        m.doPachner(mv[0], mv[1]);
        m.checkMoveData();
    }

    m = standardSphere!2;
    auto testMoves2 = [
        [[1,2,3],[4]],
        [[0],[1,2,3]],
        [[1,2,3],[6]]
    ];

    foreach(mv; testMoves2)
    {
        m.doPachner(mv[0], mv[1]);
        m.checkMoveData();
    }

    static foreach (d; iota(2,6))
    {{
        // writeln("dimension = ", d);

        int numMoves = 20;

        auto mfd = standardSphere!d;
        mfd.computePachnerMoves.shouldBeEmpty;
        mfd.numValidMoves_.shouldEqual(0);
        mfd.moves_.length.shouldEqual(2^^(d+2) - d - 4);
        m.checkMoveData();

        mfd.Move[] movesDone;
        foreach(i; numMoves.iota)
        {
            real totMoves = mfd.numValidMoves_ + mfd.fVector[$-1];
            if (uniform(0, totMoves) < mfd.numValidMoves_)
            {
                mfd.Move chosenMove = mfd.moves_.choice;
                while (mfd.contains(chosenMove.coCenter))
                {
                    chosenMove = mfd.moves_.choice;
                }
                movesDone ~= chosenMove;
            }
            else
            {
                auto facet = mfd.randomFacetOfDim(d).array;
                movesDone ~= mfd.Move(facet, [i+d+2]);
            }
            // mfd.writeln;
            mfd.doPachner(movesDone[$-1].center, movesDone[$-1].coCenter);
            m.checkMoveData();
        }
        // Do inverse moves, in reverse order
        foreach(move_; movesDone.retro)
        {
            mfd.doPachner(move_.coCenter, move_.center);
            m.checkMoveData();
        }
    }}
}

@Name("value semantics") pure /* @safe */ unittest
{
    auto m1 = octahedron;
    auto m2 = m1;

    m2.doPachner([0,1], [2,4]);
    m2.facets.should.not ~ octahedron.facets;
    m2.facets.should.not ~ m1.facets;
    m1.facets.should ~ octahedron.facets;
}

@Name("value semantics for contained simplices") pure @safe unittest
{
    auto vertices = [0,1,2];
    auto edge1 = vertices[0..2];
    auto edge2 = vertices[1..3];
    // edge1 and edge2 slices now overlap at index 1

    auto edge3 = [0,2];

    auto m = Manifold!1([edge1, edge2, edge3]);
    m.facets.shouldBeSameSetAs([[0,1],[1,2],[0,2]]);
    vertices[1] = 42;
    [42].should.not in m.simplices(0);
    m.facets.shouldBeSameSetAs([[0,1],[1,2],[0,2]]);
}