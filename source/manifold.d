import algorithms : canFind, eulerCharacteristic, is2Sphere, isCircle,
    isConnected, isPureOfDim, join;
import fluent.asserts;
import simplicial_complex : fVector, SimplicialComplex, assertValidSimplex;
import std.algorithm : all, any, copy, canFind, each, equal, filter, find, joiner,
    map, maxElement, sort, sum, uniq;
import std.conv : to;
import std.range : array, back, chain, ElementType, empty, enumerate, front, iota,
    isInputRange, popFront, save, walkLength;
import unit_threaded : Name;
import utility : binomial, isSubsetOf, SmallMap, staticIota, subsets, subsetsOfSize,
    throwsWithMsg, toImmutStaticArray;

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

    static struct Simplex(int maxDim)
    {
        Vertex_[maxDim + 1] vertices;
        size_t nValid;
        auto opSlice()
        {
            return vertices[0 .. nValid];
        }
    }

    static Simplex!dimension_ toSimplex(R)(R range)
    if (isInputRange!R && is(ElementType!R : Vertex_))
    {
        assert(range.walkLength <= dimension_ + 1);
        Simplex!dimension_ simplexToReturn;
        copy(range, simplexToReturn.vertices[]);
        simplexToReturn.nValid = range.walkLength;
        return simplexToReturn;
    }

    alias Ridge = Vertex[dimension_];

    SimpComp simpComp;
    size_t[Simplex!dimension_] degreeMap;

    alias RidgeLink = Simplex!1;
    RidgeLink[Ridge] ridgeLinks;

    size_t[dimension_ + 1] numSimplices;   
    ulong totalSquaredHingeDegree;
public:
    /// Dimension of the manifold
    static immutable dimension = dimension_;

    /// Vertex type used in the manifold
    alias Vertex = Vertex_;

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    /// We can initialize the manifold from a range of ranges of vertices
    this(F)(F initialFacets)
    if (isInputRange!F && isInputRange!(ElementType!F)
        && is(ElementType!(ElementType!F) : Vertex))
    {
        initialFacets.each!(f => this.insertFacet(f));

        assert(this.isPureOfDim(dimension),
            "not all facets have the correct dimension");

        assert(this.isConnected,
            "facets do not define a connected simplicial complex");

        static if(dimension >= 1)
        {
            assert(simplices(dimension - 1).all!(s => degree(s) == 2),
                "found a ridge with degree not equal to 2");
        }

        static if(dimension >= 2)
        {
            // TO DO: Figure out why we need "this.link" below (low priority)
            assert(simplices(dimension - 2).all!(
                s => SimpComp(this.link(s)).isCircle),
                "found a hinge whose link is not a circle");
        }

        static if(dimension >= 3)
        {
            // TO DO: Figure out why we need "this.link" below (low priority)
            assert(simplices(dimension - 3).all!(
                s => SimpComp(this.link(s)).is2Sphere),
                "found a codimension-3 simplex whose link is not a 2-sphere");
        }

        numSimplices[] = simpComp.fVector[];

        static if(dimension >= 2)
        {
            totalSquaredHingeDegree = simplices(dimension - 2)
                .map!(s => this.degree(s)^^2).sum;
        }


        // Sanity checking. TO DO: More checking!
        // foreach(ridge; simplices(dimension - 1))
        // {
        //     auto ridgeVertices = ridge.to!(Vertex[dimension]);
        //     Vertex[2] linkInManifoldAA = ridgeLinks[ridgeVertices];
        //     Vertex[] linkFromSimpComp = simpComp.link(ridge).map!(s => s.front).array.dup.sort;
        // }

    }

    this(this) pure @safe
    {
        ridgeLinks = ridgeLinks.dup;
        degreeMap = degreeMap.dup;

    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(S)(S simplex) const
    if (isInputRange!S && is(ElementType!S : Vertex))
    {
        assert(this.simpComp.facets.any!(f => simplex.isSubsetOf(f))
            == !((toSimplex(simplex) in degreeMap) is null));

        return !((toSimplex(simplex) in degreeMap) is null);
    }

    /***************************************************************************
    Returns the degree of a simplex in the simplicial complex.
    */
    size_t degree(S)(S simplex) const
    if (isInputRange!S && is(ElementType!S : Vertex))
    {
        assert(this.contains(simplex));
        assert(degreeMap[toSimplex(simplex)] == star(simplex).walkLength);
        return degreeMap[toSimplex(simplex)];
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
    private void insertFacet(F)(F facet_)
    if (isInputRange!F && is(ElementType!F : Vertex))
    {       
        assert(facet_.walkLength == dimension + 1,
            "facet has wrong dimension");

        Vertex[dimension + 1] facetBuffer;
        copy(facet_, facetBuffer[]);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);

        assert(toSimplex(facet) !in degreeMap);

        this.simpComp.insertFacet(facet);

        foreach(simplex_; facet.subsets)
        {
            auto simplex = toSimplex(simplex_);

            ++degreeMap[simplex];
        
            if(simplex_.walkLength == dimension - 1)
            {
                totalSquaredHingeDegree += 2*degreeMap[simplex] - 1;
            }
        
            if(simplex_.walkLength == dimension)
            {
                auto oppVerts = facet.filter!(v => !simplex_.canFind(v));
                assert(oppVerts.walkLength == 1);
                auto oppVert = oppVerts.front;
               
                auto ridge = simplex_.toImmutStaticArray!dimension_();
                auto ptrToLink = ridge in ridgeLinks;
                if(!ptrToLink)
                {
                    assert(simplex in degreeMap);
                    assert(degreeMap[simplex] == 1);

                    // If the degree of the ridge is one we make sure it is
                    // the *first* vertex in the Vertex[2] pair that valid
                    Vertex[2] linkVertices;
                    linkVertices.front = oppVert;
                    size_t nValid = 1;
                    ridgeLinks[ridge] = RidgeLink(linkVertices, nValid);
                }
                else
                {                   
                    assert(simplex in degreeMap);
                    assert(degreeMap[simplex] == 2);

                    (*ptrToLink).vertices.back = oppVert;
                    (*ptrToLink).nValid = 2;
                    (*ptrToLink).vertices[].sort;
                }
                assert(ridge in ridgeLinks);

            }
        }
        assert(toSimplex(facet) in degreeMap);
    }

    // Special version of removeFacet to update tracked info.
    private void removeFacet(F)(F facet_)
    if (isInputRange!F && is(ElementType!F : Vertex))
    {
        assert(facet_.walkLength == dimension + 1);

        Vertex[dimension + 1] facetBuffer;
        copy(facet_, facetBuffer[]);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);      
        assert(toSimplex(facet) in degreeMap);

        this.simpComp.removeFacet(facet);
        
        foreach(simplex_; facet.subsets)
        {
            auto simplex = toSimplex(simplex_);

            assert(simplex in degreeMap);
            --degreeMap[simplex];

            if(simplex_.walkLength == dimension - 1)
            {
                totalSquaredHingeDegree -= 2*degreeMap[simplex] - 1;
            }

            // If the degree of a ridge is one, make sure it is the first
            // vertex in the Vertex[2] pair that is the valid one
            if((simplex_.walkLength == dimension) && (degreeMap[simplex] == 1))
            {
                auto ridge = simplex_.toImmutStaticArray!dimension();
                assert(ridge in ridgeLinks);

                auto oppVerts = facet.filter!(v => !simplex_.canFind(v));
                assert(oppVerts.walkLength == 1);
                auto oppVert = oppVerts.front;

                assert(ridgeLinks[ridge].vertices[].canFind(oppVert));
                if(ridgeLinks[ridge].nValid == 2)
                {
                    if(ridgeLinks[ridge].vertices[].front == oppVert)
                    {
                        ridgeLinks[ridge].vertices[].front = ridgeLinks[ridge].vertices.back;
                    }
                    else
                    {
                        assert(ridgeLinks[ridge].vertices[].back == oppVert);
                    }
                    ridgeLinks[ridge].vertices[].back = Vertex.init;
                    ridgeLinks[ridge].nValid = 1;
                }
                else
                {
                    assert(ridgeLinks[ridge].nValid == 1);
                }

            }

            if(degreeMap[simplex] == 0)
            {
                degreeMap.remove(simplex);
                
                if(simplex_.walkLength == dimension)
                {
                    auto ridge = simplex_.toImmutStaticArray!dimension();
                    assert(ridgeLinks[ridge].nValid == 1);
                    ridgeLinks.remove(ridge);
                }
            }

        }
        assert(toSimplex(facet) !in degreeMap);
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

    auto totSqrDegree() const
    {
        return totalSquaredHingeDegree;
    }

}

/*******************************************************************************
Returns a list of all the possible pachner moves for the given manifold
*/
const(Vertex)[][] pachnerMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold)
{
    const(Vertex)[][] result;
    foreach(simp_, deg; manifold.degreeMap)
    {
        auto simp = simp_[];
        if(deg == manifold.dimension + 2 - simp.walkLength)
        {
            auto coCenter = manifold.findCoCenter(simp);
            if(coCenter.empty || (!manifold.contains(coCenter)))
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

    // TO DO: Improve this test! The simplicial complex defined by given facets
    // has additional things wrong with it than simply a non deg-2 ridge 
    Manifold!2([[1,2,3]]).throwsWithMsg(
        "found a ridge with degree not equal to 2");

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
void doPachner(Vertex, int dim)(
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] centerFacet,
    const(Vertex) newVertex)
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}
///
@Name("doPachner 1 -> (dim + 1)") @safe unittest
{
    auto m = Manifold!2(standardSphereFacets(2));
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
void doPachner(Vertex, int dim)(
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center)
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
    auto m = Manifold!2(standardSphereFacets(2));
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
void doPachnerImpl(Vertex, int dim)(
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center,
    const(Vertex)[] coCenter_
)
{
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");
    immutable centerDim = center.walkLength.to!int - 1;
    immutable coCenterDim = coCenter_.walkLength.to!int - 1;
    alias SC = SimplicialComplex!Vertex;

    auto oldPiece = (coCenterDim == 0)
        ? SC([center])
        : join(SC(coCenter_.subsetsOfSize(coCenterDim)), SC([center]));
    assert(oldPiece.isPureOfDim(dim));
    assert(manifold.star(center).map!array.array.sort
        .equal!equal(oldPiece.facets.map!array.array.sort));

    auto newPiece = (centerDim == 0)
        ? SC([coCenter_])
        : join(SC(center.subsetsOfSize(centerDim)), SC([coCenter_]));
    assert(newPiece.isPureOfDim(dim));
    assert(newPiece.facets.walkLength + oldPiece.facets.walkLength == dim + 2);

    oldPiece.facets.each!(f => manifold.removeFacet(f));
    newPiece.facets.each!(f => manifold.insertFacet(f));
    manifold.numSimplices.modifyFVector(center.length);

    // TO DO: Do sanity checking for manifold
}

// TO DO: Implement 4->4 moves (Important!)

auto findCoCenter(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center
)
{
    assert(manifold.contains(center));
    return manifold.link(center).joiner.array.dup.sort.uniq.array;
}

/*******************************************************************************
Returns the coCenter for the Pachner move with given center. For efficiency we
also need to know a facet with face center.
*/
auto coCenter(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center,
    const(Vertex)[] facet
)
{
    assert(manifold.contains(facet));
    assert(manifold.contains(center));
    assert(center.isSubsetOf(facet));
    
    // The coCenter of a facet is a new vertex not in the manifold.
    assert(center.walkLength < dim + 1);


    // TO DO: Clean this up!
    auto ridges = facet.subsetsOfSize(dim)
        .filter!(r => center.isSubsetOf(r)).map!(r => r.toImmutStaticArray!dim());
    auto coCenterVerts = ridges.map!(r => manifold.ridgeLinks[r].vertices[])
        .array.joiner.array.dup.sort.uniq.array;

    assert(coCenterVerts.equal(manifold.findCoCenter(center)));
    return coCenterVerts;
}

// TO DO: Separate unittesting for findCoCenter

// TO DO: Use standardSphereFacets below in appropriate places

/*******************************************************************************
Returns a lazy range that goes through the facets of the "standard" sphere
of dimension `dim`. This manifold is just the boundary if the "standard"
(dim+1)-simplex [0,1,2, ... (dim+1)].
*/
auto standardSphereFacets(int dim)
{
    assert(dim >= 1);
    return iota(dim + 2).subsetsOfSize(dim + 1).map!array;
}
///
@Name("standardSphereFacets") @safe unittest
{
    standardSphereFacets(2).should.containOnly(
        [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);

    auto s = Manifold!2(standardSphereFacets(2));
    s.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);

    foreach(dim; staticIota!(1,9))
    {
        immutable m = Manifold!dim(standardSphereFacets(dim));
        assert(m.numFacets == dim + 2);
    }
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
@Name("contains (pure nothrow /* @nogc */ @safe)") pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge07 = [0,7];

    () pure nothrow /* @nogc */ @safe {
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
auto movesAtFacet(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] facet
)
{
    alias m = manifold;  

    return facet.subsets.map!array.filter!(
        center => (center.walkLength == dim + 1) 
            || (m.degree(center) == dim + 2 - center.walkLength 
                && !m.contains(m.coCenter(center, facet))));
}
///
@Name("movesAtFacet") unittest
{
    auto emptyMfd = Manifold!3();
    assert(emptyMfd.movesAtFacet([]).empty);

    /* If the manifold is the boundary of a simplex (i.e. a sphere with the
    minimum number of facets) then only the type 1 -> (dim + 1) Pachner moves
    are valid. */    
    foreach(dim; staticIota!(1,4))
    {
        immutable sphere = Manifold!dim(standardSphereFacets(dim));
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


// TO DO: Adapt old code below for new manifold type!

// /******************************************************************************
// * Returns a manifold loaded from the file specified by fileName. If fileName
// * is the empty string, the returned manifold is the standard sphere.
// */
// Manifold!dim loadManifold(size_t dim = dimManifold)(string fileName)
// {
//     auto manifoldFile = File(fileName, "r"); // Open file in read-only mode

//     string facets;
//     foreach (string line; manifoldFile.lines)
//     {
//         // We allow comments starting with '#'
//         if (line.front != '#')
//         {
//             facets ~= line.strip;
//         }
//     }

//     try
//     {
//         auto loaded = Manifold!dim(facets.parse!(Simplex[]));
//         return loaded;
//     }
//     catch (Exception ex)
//     {
//         ex.msg ~= "\n\n\t ERROR: encountered malformed facet list in initial manifold triangulation file: "
//             ~ fileName ~ "\n";
//         throw ex;
//     }
// }
// ///
// unittest
// {

//     auto manifold = loadManifold!2("data/manifold_sampler_unittest_load.dat");
//     auto expectedManifold = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4,
//             5], [0, 5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
//     assert(manifold == expectedManifold);

//     assertThrown(loadManifold!2("data/manifold_sampler_unittest_bad_load.dat"));
// }

// /******************************************************************************
// * Saves a manifold to file specified by fileName.
// */
// void saveTo(size_t dimension)(Manifold!dimension manifold, string fileName)
// {
//     auto saveFile = File(fileName, "w"); // Open in write-only mode

//     saveFile.writeln("# created ", Clock.currTime.to!DateTime);
//     saveFile.write("[");
//     foreach (record; manifold.facetRecords[0 .. $ - 1])
//     {
//         saveFile.write(record.facet, ",");
//     }
//     saveFile.writeln(manifold.facetRecords.back.facet, "]");
// }
// ///
// unittest
// {

//     auto fileName = "data/manifold_sampler_unittest_save.dat";
//     auto sphere = Manifold!4(standardSphere(4));
//     sphere.saveTo(fileName);
//     Manifold!4 loadedSphere = loadManifold!4(fileName);
//     assert(loadedSphere == sphere);
// }

// unittest
// {

//     auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
//             [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);

//     assert(octahedron.fVector == [6, 12, 8]);

//     auto sphere7 = Manifold!7(standardSphere(7));

//     // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
//     assert(sphere7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
//             9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
//             (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

//     auto hyperbolicDodecahedral = loadManifold!3("data/manifold_sampler_unittest_dodecahedral.dat");

//     // TO DO: Get reference for this...
//     assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);

// }
