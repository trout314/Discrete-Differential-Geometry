import algorithms : canFind, eulerCharacteristic, is2Sphere, isCircle,
    isConnected, isPureOfDim, join;
import simplicial_complex : fVector, simplicialComplex, SimplicialComplex,
    assertValidSimplex, loadSimplicialComplex;
import std.algorithm : all, any, copy, canFind, each, equal, filter, find,
    joiner, map, maxElement, setDifference, merge, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown;
import std.range : array, back, chain, chunks, cycle, ElementType, empty,
    enumerate, front, iota, isInputRange, isOutputRange, only, popBack, popFront,
    put, replace, retro, save, slide, take, walkLength, zip;
import unit_threaded : Name, shouldEqual, shouldBeSameSetAs;
import utility : binomial, isInputRangeOf, isInputRangeOfInputRangeOf,
    isSubsetOf, nGonTriangs, productUnion, SmallMap, StackArray,
    staticIota, subsets, subsetsOfSize, throwsWithMsg, toStackArray, swapPop;
import std.stdio : File, writeln;
import std.typecons : Flag, No, Yes;
import std.math : approxEqual;

import std.stdio : writeln, writefln, write;

import std.traits : Unqual;
import std.random : uniform;

import moves : Move;

import std.array : staticArray;

import std.algorithm : countUntil, copy, map, min;

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

    alias Move_ = Move!(dimension_, Vertex_);

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
public:
    // TO DO: Make private
    Move_[] moves;
    size_t[Vertex[]] indxOfCenter;
    size_t[][Vertex[]] indicesOfCoCenter;

    // TO DO: figure out how to make oldPiece and newPiece const
    void modifyMoveDataOnMove(S)(Move!(dimension, Vertex) mv, ref S oldPiece, ref S newPiece)
    {
        auto oldP = SimpComp(oldPiece);
        auto newP = SimpComp(newPiece);
        writeln("   newP: ", newP);
        writeln("   oldP: ", oldP);


        auto cen = mv.center;
        auto coCen = mv.coCenter;

        auto cenLen = cen.walkLength.to!int;
        assert(cenLen>0);

        auto coCenLen = coCen.walkLength.to!int;
        assert(coCenLen>0);

        size_t[] toRemove = [];
        size_t[] toUpdate = [];
        Move_[] newMoves = [];

        if (coCenLen > 1)
        {
            assert(cen in indxOfCenter, "pachner move center indx not found");
            assert(coCen in indicesOfCoCenter, "pachner move cocenter indx not found");
        }

        static foreach (d; 0 .. this.dimension)
        {
            writeln("   checking oldPiece simplices of dimension ", d);
            foreach(simp; newP.simplices(d))
            {
                write("     checking simp = ", simp);
                bool wasMove = (simp in indxOfCenter) !is null;
                bool isNowMove = (degree(simp) == dimension - d + 1);

                if(wasMove && !isNowMove)
                {
                    writeln(", remove index ", indxOfCenter[simp.array]);
                    toRemove ~= indxOfCenter[simp.array];
                }
                else if (!wasMove && isNowMove)
                {
                    newMoves ~= typeof(mv)(simp.array, this.findCoCenter(simp).sort);
                    writeln(", add move: ", newMoves[$-1]);
                }
                else if (wasMove && isNowMove)
                {
                    writeln(", update index ", indxOfCenter[simp.array]);
                    toUpdate ~= indxOfCenter[simp.array];
                }
                else
                {
                    "".writeln;
                }
            }
        }

        toRemove.sort;
        writeln("   indices toRemove = ", toRemove.retro);
        writeln("   indices toUpdate = ", toUpdate);
        writeln("   newMoves = ");
        foreach(m; newMoves)
        {
            writeln("      ", m);
        }

        "   updating moves...".writeln;
        foreach(i; toUpdate)
        {
            auto thisCoCen = moves[i].coCenter;
            auto thisOldCoCen = thisCoCen.idup;
            auto thisCen = moves[i].center;
            writeln("      updating index ", i);
            
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
            writeln("         newVerts: ", newVerts);
            writeln("         oldVerts: ", oldVerts);

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
            moves[i] = typeof(mv)(thisCen, thisCoCen.sort);
            indicesOfCoCenter[thisCoCen.idup] ~= i;

            // Update indices pointing to updated cocenter
            auto newIndices = indicesOfCoCenter[thisOldCoCen].filter!(k => k != i).array;
            if (newIndices.walkLength > 0)
            {
                indicesOfCoCenter[thisOldCoCen] = newIndices;
            }
            else
            {
                indicesOfCoCenter.remove(thisOldCoCen);
            }
        }


        "   removing moves...".writeln;
        foreach(i; toRemove.retro)
        {
            writeln("      removing index ", i);
            auto lastIndx = moves.length - 1;
            auto lastMove = moves[lastIndx];
            auto goneMove = moves[i];
            indxOfCenter[lastMove.center.idup] = i;
            indicesOfCoCenter[lastMove.coCenter.idup]
                = indicesOfCoCenter[lastMove.coCenter].replace(lastIndx.only, i.only);
            
            this.moves.swapPop(i);
            indicesOfCoCenter.remove(goneMove.coCenter);
            indxOfCenter.remove(goneMove.center);
        }


        "   adding new moves...".writeln;
        foreach(newMove; newMoves)
        {
            writeln("      adding move: ", newMove);
            moves ~= newMove;
            indxOfCenter[newMove.center.idup] = moves.length - 1;
            indicesOfCoCenter[newMove.coCenter.idup] ~= moves.length - 1;
        }


        // TO DO: remove or add (simplex, coCenter) to list
        // of pachner moves if needed, due to degree change

    }

    /// Dimension of the manifold
    static immutable dimension = dimension_;

    /// Vertex type used in the manifold
    alias Vertex = Vertex_;

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

        int numValidMoves;
        foreach(simp, deg; degreeMap)
        {
            if(simp[].walkLength < dimension + 1)
            {
                if(deg == dimension + 2 - simp[].walkLength)
                {
                    auto coCenter = this.findCoCenter(simp[]);
                    if(!this.contains(coCenter))
                    {
                        ++numValidMoves;
                    }
                    moves ~= Move_(simp[], coCenter.sort);
                }
            }
        }

        foreach(indx, mv; moves.enumerate)
        {
            this.indxOfCenter[mv.center.idup] = indx;
            this.indicesOfCoCenter[mv.coCenter.idup] ~= indx;
        }
    }

    this(this) pure @safe
    {
        ridgeLinks = ridgeLinks.dup;
        degreeMap = degreeMap.dup;
        moves = moves.dup;
    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        bool notInDegreeMap = (toNSimp(simplex) in degreeMap) is null;
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
}

/*******************************************************************************
Returns a list of all the possible pachner moves except for the 1->(dim+1) moves
that are valid in this manifold. (Note that 1->(dim+1) moves are always valid.)
*/
Move!(dim, Vertex)[] computePachnerMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    Move!(dim, Vertex)[] result;
    foreach(simp_, deg; mfd.degreeMap)
    {
        if(simp_.length < dim + 1)
        {
            auto simp = simp_[];
            if(deg == mfd.dimension + 2 - simp.walkLength)
            {
                auto coCenter = mfd.findCoCenter(simp);
                if(coCenter.empty || (!mfd.contains(coCenter)))
                {
                    result ~= Move!(dim, Vertex)(simp, coCenter);
                }
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

    tetrahedron.computePachnerMoves.map!(mv => mv.center.array).empty.shouldEqual(true);

    // octahedron.doPachner([0,1], [2,4]);


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
@Name("doPachner") pure @safe unittest
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
    if(coCenter.walkLength > 1)
    {
        // auto pm = Move!(dim, Vertex)(center, coCenter);
        // assert(manifold.hasValidPachnerMove(pm), "not a valid pachner move");
    }
    assert(manifold.contains(center), "center of move not in manifold");
    assert(!manifold.contains(coCenter), "coCenter of move in manifold");
 
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
    // auto bdry = productUnion(center_.subsetsOfSize(centerDim), coCenter_.subsetsOfSize(coCenterDim));

    alias SC = SimplicialComplex!(Vertex, dim);
    alias MFD = Manifold!(dim, Vertex);
    assert(SC(oldPiece).isPureOfDim(dim));
    assert(SC(newPiece).isPureOfDim(dim));
    assert(MFD(chain(oldPiece, newPiece)).numFacets == dim + 2);
    assert(manifold.star(center).map!array.array.sort
        .equal!equal(oldPiece.map!array.array.sort));
    // assert(SC(bdry).isPureOfDim(dim-1));
    // static foreach(d; dim.iota)
    // {
    //     // assert(SC(bdry).simplices(d).length == dim + 2);
    // }

    oldPiece.each!(f => manifold.removeFacet(f));
    newPiece.each!(f => manifold.insertFacet(f));
    manifold.numSimplices.modifyFVector(cenLen);

    auto pm = Move!(dim, Vertex)(center_, coCenter_);
    manifold.modifyMoveDataOnMove(pm, oldPiece, newPiece);

    // TO DO: Do sanity checking for manifold
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

/******************************************************************************
Returns the "standard" triangulation of a sphere of the given dimension. This
is just the boundary of the simplex of one higher dimension.
*/
Manifold!(dim, int) standardSphere(int dim)()
{
    static assert(dim >= 1);
    return Manifold!(dim, int)(iota(dim + 2).subsetsOfSize(dim + 1));
}
///
@Name("standardSphere") unittest
{
    standardSphere!1.facets.shouldBeSameSetAs(
        [[0,1], [1,2], [0,2]]); 
    standardSphere!2.facets.shouldBeSameSetAs(
        [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    standardSphere!3.facets.shouldBeSameSetAs(
        [[0,1,2,3], [0,1,2,4], [0,1,3,4], [0,2,3,4],[1,2,3,4]]);
}

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
    real epsilon = 1e-20;
    assert(pyramid.meanDegree(0).approxEqual((2 * 3 + 3 * 4)/5.0, epsilon));
    assert(pyramid.meanDegree(1).approxEqual(2.0, epsilon));
    assert(pyramid.meanDegree(2).approxEqual(1.0, epsilon));
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
    immutable real epsilon = 1e-20;
    assert(pyramid.degreeVariance(0)
        .approxEqual((2 * (3 - md0)^^2 + 3 * (4 - md0)^^2)/5.0, epsilon));

    assert(pyramid.degreeVariance(1).approxEqual(0.0, epsilon));
    assert(pyramid.degreeVariance(2).approxEqual(0.0, epsilon));

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
@Name("hinge moves") pure @safe unittest
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