import algorithms : canFind, eulerCharacteristic, is2Sphere, isCircle,
    isConnected, isPureOfDim, join;
import fluent.asserts;
import simplicial_complex : fVector, simplicialComplex, SimplicialComplex, assertValidSimplex, loadSimplicialComplex;
import std.algorithm : all, any, copy, canFind, each, equal, filter, find, joiner,
    map, maxElement, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown;
import std.range : array, back, chain, chunks, cycle, ElementType, empty, enumerate, front, iota,
    isInputRange, only, popFront, save, slide, take, walkLength;
import unit_threaded : Name;
import utility : binomial, isInputRangeOf, isInputRangeOfInputRangeOf, isSubsetOf, nGonTriangs, productUnion, SmallMap,
    StackArray, staticIota, subsets, subsetsOfSize, throwsWithMsg, toImmutStaticArray, toStackArray, toStaticArray;
import std.stdio : File, writeln;
import std.typecons : Flag, No, Yes;
import std.math : approxEqual;

import std.stdio : writeln, writefln;


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
        return range.toStaticArray!dimension;
    }

    alias Facet = Vertex[dimension_ + 1];
    static Facet toFacet(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.toStaticArray!(dimension_ + 1);
    }
public:
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

        assert(this.findProblems.empty, 
            this.findProblems.joiner(", ").array.to!string);
    }

    this(this) pure @safe
    {
        ridgeLinks = ridgeLinks.dup;
        degreeMap = degreeMap.dup;
    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        assert(this.simpComp.facets.any!(f => simplex.isSubsetOf(f))
            == !((toNSimp(simplex) in degreeMap) is null));

        return !((toNSimp(simplex) in degreeMap) is null);
    }

    /***************************************************************************
    Returns the degree of a simplex in the simplicial complex.
    */
    size_t degree(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        assert(this.contains(simplex));
        assert(degreeMap[toNSimp(simplex)] == star(simplex).walkLength);
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
        assert(toNSimp(facet) !in degreeMap);

        this.simpComp.insertFacet!(No.checkForFacetFaces)(facet);

        foreach(simplex_; facet.subsets)
        {
            auto simplex = toNSimp(simplex_);
            ++degreeMap[simplex];
        
            if(simplex.length <= dimension - 1)
            {
                totSqrDegrees[simplex.length - 1] += 2*degreeMap[simplex] - 1;
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

    private void printDiagnosticReport() const
    {
        import std.stdio : writeln;
        debug writeln("dimension: ", dimension);
        debug writeln("\tSComp: ", this.asSimplicialComplex.facets);
        debug writeln("numSimplices: ", numSimplices);
        debug writeln("\tdegreeMap: ", degreeMap);
        debug writeln("\tridgeLinks: ", ridgeLinks);        
    }

    /// We provide access to the manifold as a simplicial complex
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const pure nothrow @nogc @safe 
    {
        return simpComp; 
    }

    alias asSimplicialComplex this;
}

/*******************************************************************************
Returns a list of all the possible pachner moves for the given manifold
*/
const(Vertex)[][] pachnerMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    const(Vertex)[][] result;
    foreach(simp_, deg; mfd.degreeMap)
    {
        auto simp = simp_[];
        if(deg == mfd.dimension + 2 - simp.walkLength)
        {
            auto coCenter = mfd.findCoCenter(simp);
            if(coCenter.empty || (!mfd.contains(coCenter)))
            {
                result ~= simp.dup;
            }
        }           
    }
    return result;
}

///
@Name("Manifold doc tests") @safe unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));

    octahedron.fVector.should.equal([6UL,12,8]);
    octahedron.eulerCharacteristic.should.equal(2);

    auto tetrahedron = Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    octahedron.star([1,2]).should.containOnly([[0,1,2], [1,2,5]]);    

    octahedron.pachnerMoves.should.containOnly(
        [[0,1], [0,2], [0,3], [0,4], [1,2], [1,4],  // 1-simplices
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5],
         [0,1,2], [0,1,4], [0,2,3], [0,3,4],        // 2-simplices
         [1,2,5], [1,4,5], [2,3,5], [3,4,5]]);      

    tetrahedron.pachnerMoves.should.containOnly(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]);

    tetrahedron.doPachner([1,2,3], 5);

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

    Manifold!2([[1,2,3]]).throwsWithMsg(
        "found a ridge with degree not equal to 2, " 
        ~ "found a hinge whose link is not a circle");

    Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,5,6], [1,5,7], [1,6,7],
        [5,6,7]]).throwsWithMsg("found a hinge whose link is not a circle");

    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5],
        [3,4,5], [1,4,5]];

    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];

    Manifold!2(chain(octahedron, sphere)).throwsWithMsg(
        "facets do not define a connected simplicial complex");

    // These facets separately define 3-spheres, but share the vertex 1
    auto sphere3 = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];
    auto sphere3A =[[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]];

    Manifold!3(chain(sphere3, sphere3A)).throwsWithMsg(
        "found a codimension-3 simplex whose link is not a 2-sphere");
}

// More tests
@Name("Manifold (additional)") pure @safe unittest
{
    static assert(!__traits(compiles, Manifold!2([["a", "bubba", "gump"]])));

    auto m1 = Manifold!(1, string)([["a", "b"], ["b", "c"], ["a", "c"]]);
    static assert(m1.dimension == 1);
    static assert(is(m1.Vertex == string));
}

///
@Name("pachnerMoves") @safe unittest
{
    auto m = Manifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    m.pachnerMoves.should.containOnly([[1], [2, 3], [2, 3, 5],
        [1, 3, 4], [1, 2, 3], [5], [2, 4, 5], [3, 4, 5], [1, 2, 4], [2, 4],
        [3, 4]]);

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.pachnerMoves.should.containOnly(
        [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5], [3,4,5], [1,4,5],
        [2, 3], [0, 1], [1, 5], [4, 5], [0, 3], [1, 4],
        [1, 2], [0, 4], [0, 2], [2, 5], [3, 5], [3, 4]]);
}

/*******************************************************************************
Does a pachner move of type 1 -> (dim + 1) with the new vertex given by
the user.
*/
void doPachner(Vertex, int dim, C)(
    ref Manifold!(dim, Vertex) manifold,
    C centerFacet,
    const(Vertex) newVertex)
if (isIRof!(C, const(Vertex)))
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}
///
@Name("doPachner 1 -> (dim + 1)") @safe unittest
{
    auto m = standardSphere!2;
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    m.doPachner([1,2,3], 4);
    m.facets.should.containOnly(
        [[0,1,2],[0,1,3], [0,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    m.doPachner([0,2,3], 7);
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,7],
        [0,3,7], [1,2,4], [1,3,4], [2,3,4], [2,3,7]]);

    m.doPachner([7]);
    m.doPachner([4]);
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
}

/*******************************************************************************
Does a pachner move which replaces the star of the given simplex. This is a
move of type (dim + 2 - center.length) -> center.length
*/
void doPachner(Vertex, int dim, C)(
    ref Manifold!(dim, Vertex) manifold,
    C center)
if (isIRof!(C, const(Vertex)))
{
    // TO DO: Better error
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");

    auto coCenter = manifold.findCoCenter(center);
    assert(!coCenter.empty);

    manifold.doPachnerImpl(center, coCenter);
}
///
@Name("doPachner n -> (dim + 2 - n), 1 < n < dim + 2") @safe unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.doPachner([1,2]);
    octahedron.degree([0]).should.equal(5);
    octahedron.degree([5]).should.equal(5);
    octahedron.degree([1]).should.equal(3);
    octahedron.degree([2]).should.equal(3);

    // We can undo the 2->2 move
    octahedron.doPachner([0,5]);
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.doPachner([0,1,2], 99);
    octahedron.doPachner([99]);
}
// NOTE: The following unittest cannot be @safe since throwsWithMsg 
// catches an Error
///
@Name("doPachner (errors)") @system unittest
{
    // Can't do 2->2 move on the boundary of a 3-simplex
    auto m = Manifold!2([[1,2,3],[1,2,4], [1,3,4], [2,3,4]]);   
    m.doPachner([1,2]).throwsWithMsg("bad pachner move");
}
///
@Name("doPachner (general)") unittest
{
    auto m = standardSphere!2;
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    m.doPachner([1,2,3], 4);
    m.doPachner([1,2]);
    m.facets.should.containOnly(
        [[0,1,3], [0,1,4], [0,2,3], [0,2,4], [1,3,4], [2,3,4]]);
    m.doPachner([0,4]);
    m.doPachner([4]);
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);      
}

// Factor out the common code for the two types of doPachner
void doPachnerImpl(Vertex, int dim, C, D)(
    ref Manifold!(dim, Vertex) manifold,
    C center,
    D coCenter
)
if (isIRof!(C, const(Vertex)) && isIRof!(D, const(Vertex)))
{
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");
 
    // Buffer for holding vertices in center (followed by coCenter)
    const(Vertex)[(dim + 2)] cBuffer = chain(center, coCenter)
        .toStaticArray!(dim + 2);
 
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

    assert(coCenterVerts.equal(mfd.findCoCenter(center)));
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
    standardSphere!1.facets.should.containOnly(
        [[0,1], [1,2], [0,2]]); 
    standardSphere!2.facets.should.containOnly(
        [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    standardSphere!3.facets.should.containOnly(
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
@Name("modifyFVector") unittest
{
    size_t[4] fVec = [0,0,0,0];

    /* A 1 -> 4 move should give net: +1 vertices, +4 edges, +6 triangles,
    +3 tetrahdra */
    fVec.modifyFVector(4);
    fVec.should.equal([1, 4, 6, 3]);

    /* A 2 -> 3 move should give net: +0 vertices, +1 edges, +2 triangles,
    +1 tetrahedra */
    fVec.modifyFVector(3);
    fVec.should.equal([1,5,8,4]);

    fVec.modifyFVector(2);    // 3 -> 2 move
    fVec.modifyFVector(1);    // 4 -> 1 move

    // Should be back where we started
    fVec.should.equal([0,0,0,0]);

    size_t[3] fVec2 = [0,0,0];
    fVec2.modifyFVector(2); // 2 -> 2 move
    fVec2.should.equal([0,0,0]);

    fVec2.modifyFVector(3); // 1 -> 3 move

    // a 1 -> 3 move should give net: +1 vertices, +3 edges, +2 triangles
    fVec2.should.equal([1, 3, 2]);

    fVec2.modifyFVector(1); // 3 -> 1 move
    fVec2.should.equal([0,0,0]);
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
@Name("movesAtFacet") unittest
{
    /* If the manifold is the boundary of a simplex (i.e. a sphere with the
    minimum number of facets) then only the type 1 -> (dim + 1) Pachner moves
    are valid. */    
    foreach(dim; staticIota!(1,4))
    {
        immutable sphere = standardSphere!dim;
        foreach(f; sphere.facets)
        {
            sphere.movesAtFacet(f).should.containOnly([f.array]);
        }
    }

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.movesAtFacet([0,1,2]).should.containOnly([
        [0,1], [0,2], [1,2],    // 2 -> 2 moves
        [0,1,2]                 // 1 -> 3 move
    ]);
    
    foreach(f; octahedron.facets)
    {
        octahedron.movesAtFacet(f).should.containOnly(
            chain([f.array], f.subsetsOfSize(2).map!array).array);
    }
    
    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);
    pyramid.movesAtFacet([0,1,2]).should.containOnly([
        [0],                    // 3 -> 1 move
        [1,2],                  // 2 -> 2 move
        [0,1,2]                 // 1 -> 3 move
    ]);
    pyramid.movesAtFacet([0,2,3]).should.containOnly([
        [0],                    // 3 -> 1 move
        [2,3],                  // 2 -> 2 move
        [0,2,3]                 // 1 -> 3 move
    ]);
    pyramid.movesAtFacet([1,3,4]).should.containOnly([
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
@Name("loadManifold") unittest
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
@Name("saveTo") unittest
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
@Name("meanDegree") unittest
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
@Name("totalSquareDegree") unittest
{
    foreach(dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        foreach(d; 0 .. dim + 1)
        {
            // (dim + 2) choose (d + 1) faces of dimension d
            // each with degree (dim - d + 1)
            m.totalSquareDegree(d).should.equal(
                binomial(dim + 2, d + 1) * (dim - d + 1)^^2);
        }
    }

    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    pyramid.totalSquareDegree(0).should.equal(2 * 3^^2 + 3 * 4^^2);
    pyramid.totalSquareDegree(1).should.equal(9 * 2^^2);
    pyramid.totalSquareDegree(2).should.equal(6 * 1^^2);
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
@Name("degreeVariance") unittest
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
    auto hingeBuffer = hinge_.toStaticArray!(dim - 1);
    auto hinge = hingeBuffer[];

    // TO DO: Decide what to do about this magic constant 7
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    auto linkVerticesBuff = linkVertices_.toStaticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    // TO DO: Get rid of allocations here!
    auto linkEdges = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.array;
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

    assert(diskIndx < deg.nGonTriangs.walkLength);
    auto diskFacets = deg.nGonTriangs[diskIndx];

    auto oldPiece = productUnion(hinge.only, linkEdges.chunks(2));
    auto newPiece = productUnion(
        hinge.subsetsOfSize(dim - 2),
        diskFacets.map!(f => f.map!(i => linkVertices[i])));

    version(unittest)
    {
        alias SC = SimplicialComplex!(Vertex, dim);
        auto disk = SC(diskFacets);

        // None of the "internal" edges can already be in manifold
        assert(disk.simplices(1)
            .filter!(edge => disk.star(edge).walkLength == 2)
            .all!(edge => !manifold.contains(edge)));

        assert(SC(newPiece).isPureOfDim(dim));   
        assert(SC(oldPiece).isPureOfDim(dim));
        assert(manifold.star(hinge).map!array.array.sort
            .equal!equal(oldPiece.map!array.array.sort));
    }

    oldPiece.each!(f => manifold.removeFacet(f));

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
    auto hingeBuffer = hinge_.toStaticArray!(dim - 1);
    auto hinge = hingeBuffer[];

    // TO DO: Decide what to do about this magic constant 7
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    auto linkVerticesBuff = linkVertices_.toStaticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    // TO DO: Get rid of allocations here!
    auto linkEdges = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.array;
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

    assert(diskIndx < deg.nGonTriangs.walkLength);
    auto diskFacets = deg.nGonTriangs[diskIndx];

    auto newPiece = productUnion(hinge.only, linkEdges.chunks(2));
    auto oldPiece = productUnion(
        hinge.subsetsOfSize(dim - 2),
        diskFacets.map!(f => f.map!(i => linkVertices[i])));

    version(unittest)
    {
        alias SC = SimplicialComplex!(Vertex, dim);
        auto disk = SC(diskFacets.map!(f => f.map!(i => linkVertices[i])));

        // All of the disk should be in manifold
        assert(disk.facets.all!(f => manifold.contains(f)));

        assert(SC(newPiece).isPureOfDim(dim));   
        assert(SC(oldPiece).isPureOfDim(dim));
    }

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

unittest
{
    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]];
    
    auto twoPts = [[6], [7]];
    auto mfd = Manifold!3(productUnion(octahedron, twoPts));

    assert(mfd.hasValidHingeMove([1,2,3,4], 1));
    mfd.doHingeMove([0,6], [1,2,3,4], 1);
    mfd.facets.should.containOnly([[1, 4, 5, 7], [0, 1, 2, 7], [1, 4, 5, 6],
        [0, 2, 3, 7], [3, 4, 5, 7], [0, 3, 4, 7], [3, 4, 5, 6], [0, 1, 4, 7],
        [1, 2, 5, 6], [1, 2, 5, 7], [2, 3, 5, 6], [2, 3, 5, 7],[0, 1, 2, 4],
        [0, 2, 3, 4], [1, 2, 4, 6], [2, 3, 4, 6]]);

    mfd.undoHingeMove([0,6], [1,2,3,4], 1);
    mfd.facets.should.containOnly(productUnion(octahedron, twoPts).map!array);

    auto twoOtherPts = [[8], [9]];
    auto mfd4 = Manifold!4(
        productUnion(productUnion(octahedron, twoPts), twoOtherPts));

    mfd4.doHingeMove([0,6,8], [1,2,3,4], 1);
    mfd4.undoHingeMove([0,6,8], [1,2,3,4], 1);
    mfd4.facets.should.containOnly(productUnion(
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

    // TO DO: Decide what to do about this magic constant 7
    // (It comes from nGonTriangs only supporting up to 7-gon.)
    auto deg = linkVertices_.walkLength.to!int;
    auto linkVerticesBuff = linkVertices_.toStaticArray!7;
    auto linkVertices = linkVerticesBuff[0 .. deg];

    // TO DO: Get rid of allocations here!
    auto linkEdges = chain(linkVertices[], linkVertices.front.only)
        .slide(2).take(deg).joiner.array;
    deg.iota.each!(indx => linkEdges[2*indx .. 2*(indx + 1)].sort);

    assert(diskIndx < deg.nGonTriangs.walkLength);
    auto diskFacets = deg.nGonTriangs[diskIndx];

    auto disk = SimplicialComplex!(Vertex, 2)(
            diskFacets.map!(f => f.map!(i => linkVertices[i])));

    // None of the "internal" edges can already be in manifold
    return disk.simplices(1)
        .filter!(edge => disk.star(edge).walkLength == 2)
        .all!(edge => !manifold.contains(edge));
}

@Name("fVector") unittest
{
    // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
    assert(standardSphere!7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
            9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
            (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

    auto hyperbolicDodecahedral = loadManifold!3(
        "data/manifold_sampler_unittest_dodecahedral.dat");

    // TO DO: Get reference for this...
    assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);
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
unittest
{
    auto m = standardSphere!3;
    assert(m.findProblems.empty);
}
