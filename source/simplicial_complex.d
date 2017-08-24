import simplex : facesOfDim, hasFace, simplex, Simplex, assertValidSimplex;
import std.algorithm : all, any, canFind, chunkBy, copy, countUntil, each,
    equal, filter, find, joiner, map, maxElement, remove, setDifference,
    setIntersection, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown;
import std.range : array, chunks, empty, enumerate, ElementType, front, iota,
    isInputRange, popFront, refRange, walkLength, zip;
import std.traits : isArray, Unqual;
import std.typecons : Tuple, tuple;
import utility : isSubsetOf, SmallMap, staticIota, throwsWithMsg, subsetsOfSize, subsets;
import std.stdio : writeln;
import unit_threaded : Name;
import fluent.asserts : should, Assert;

import simplicial_complex_algorithms : connectedComponents, eulerCharacteristic,
    isCircle, isConnected, isPureOfDim, isOrientableSurfaceOfGenus, join;

import utility : capture;

/// Basic Functionality
@Name("doc tests") /* @safe */ pure unittest
{
    // create an empty simplicial complex with vertices of default type `int`
    SimplicialComplex!() sc;
    assert(sc.numFacets == 0);
    assert(sc.facets.empty);

    // get the type of vertex used
    static assert(is(sc.VertexType == int));

    /* while you can specify the vertex type used as a template parameter,
    we will stick to `int` for now */
    static assert(is(SimplicialComplex!int == typeof(sc)));
    
    // for clarity
    alias s = simplex;

    // insert some facets
    sc.insertFacet(s(1,2,3));
    sc.insertFacet(s(2,3,4));
    sc.insertFacet(s(4,5));
    sc.insertFacet(s(5,6));

    /* can also create simplicial complexes directly from an array of arrays of
    vertices which specify the facets */
    auto sc2 = SimplicialComplex!()([[4,5], [5,6], [1,2,3], [2,3,4]]);
    assert(sc == sc2);

    // a helper function template that returns a newly constructed complex
    auto sc3 = simplicialComplex([[4,5], [5,6], [1,2,3], [2,3,4]]);
    assert(sc3 == sc);

    // simplicial complexes are value types
    auto sc4 = sc3;
    sc4.insertFacet(s(6,7));
    assert(sc3 == sc);

    /* get the vertices in the facets, returned in order of increasing dimension
    and dictionary order within a dimension */
    assert(sc.facets.equal([[4,5], [5,6], [1,2,3], [2,3,4]]));
    assert(sc.numFacets == 4);

    // Can get a sorted array of all the facet simplices of a given dimension
    assert(sc.facets!0.empty);
    assert(sc.facets!1.equal([s(4,5), s(5,6)]));
    assert(sc.facets!2.equal([s(1,2,3), s(2,3,4)]));
    assert(sc.facets!3.empty);

    // Can get a sorted array of all the simplices of a given dimension
    assert(sc.simplices!0.equal([s(1), s(2), s(3), s(4), s(5), s(6)]));
    assert(sc.simplices!1.equal([s(1,2), s(1,3), s(2,3), s(2,4), s(3,4), 
        s(4,5), s(5,6)]));
    assert(sc.simplices!2.equal(sc.facets!2));
    assert(sc.simplices!3.empty);

    // get the f-vector of the simplicial complex
    assert(sc.fVector == [6,7,2]);

    // check for the presence of simplices
    assert(sc.contains(s(4,5)));
    assert(sc.contains(s(2,3,4)));
    assert(sc.contains(s(6)));
    assert(!sc.contains(s(1,2,5)));
    assert(!sc.contains(s(4,6)));

    // get the star of a simplex an array of facets
    assert(sc.star(s(5)).equal([[4,5], [5,6]]));
    assert(sc.star(s(4)).equal([[4,5], [2,3,4]]));
    assert(sc.star(s(2,3)).equal([[1,2,3], [2,3,4]]));

    // get link of a simplex as list of facets
    assert(sc.link(s(5)).equal!equal([[4], [6]]));
    assert(sc.link(s(4)).equal!equal([[5], [2,3]]));
    assert(sc.link(s(2,3)).equal!equal([[1], [4]]));

    // can print out a simplicial complex
    assert(sc.toString == "[[4, 5], [5, 6], [1, 2, 3], [2, 3, 4]]");

    // can construct from a range of ranges of vertices
    auto sc5 = SimplicialComplex!()(iota(3).map!(i => iota(3*i, 3*i + 3)));
    assert(sc5.facets.equal([[0, 1, 2], [3, 4, 5], [6, 7, 8]]));
}

/// Some restrictions
@Name("doc tests (errors)") pure @system unittest
{
    // -------------------------------------------------------------------------
    // Restrictions
    // -------------------------------------------------------------------------

    auto sc = simplicialComplex([[4,5], [5,6], [1,2,3], [2,3,4]]);

    // cannot insert a simplex if it is already the face of an existing facet
    sc.insertFacet([4,5]).throwsWithMsg!Error(
        "expected a simplex not already in the simplicial complex");

    sc.insertFacet([2]).throwsWithMsg!Error(
        "expected a simplex not already in the simplicial complex");

    // cannot remove a facet that isn't part of the simplicial complex
    simplicialComplex([[3,4]]).removeFacet([2,3]).throwsWithMsg!Error(
        "tried to remove a facet not in the simplicial complex");

    // not allowed to ask for link of a simplex not in the simplicial complex
    simplicialComplex([[1,3,4], [2,3,4]]).link([1,2]).throwsWithMsg!Error(
        "expected a simplex in the simplicial complex");

    sc.facets(-1).throwsWithMsg!Error("expected a non-negative dimension");
}

@Name("additional tests") unittest
{
    alias s = simplex;
    alias sComp = simplicialComplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimplicialComplex!()();
    static assert(is(sc.VertexType == int)); // Default vertex type is int

    sc.facets.empty.should.equal(true);
    sc.numFacets.should.equal(0);
    foreach(d; staticIota!(0, 16))
    {
        sc.facets!d.empty.should.equal(true);
        sc.simplices!d.empty.should.equal(true);
    }
   
    // ---------------------------- ADD SOME FACETS ----------------------------
    sc.insertFacet(s(1,2));
    sc.insertFacet([2,3]);
    sc.insertFacet([3, 4, 5]);
    
    sc.should.equal(sComp([[1,2], [2,3], [3,4,5]]));

    // -------------------------- TEST FUNCTIONALITY ---------------------------
    sc.facets.walkLength.should.equal(3);
    sc.numFacets.should.equal(3);

    sc.facets!0.empty.should.equal(true);
    sc.facets!1.should.equal([s(1,2), s(2,3)]);
    sc.facets!2.should.equal([s(3,4,5)]);

    sc.simplices!0.should.containOnly([s(1), s(2), s(3), s(4), s(5)]);
    sc.simplices!1.should.containOnly([s(1,2), s(2,3), s(3,4), s(3,5), s(4,5)]);
    sc.simplices!2.should.containOnly([s(3,4,5)]);

    foreach(d; staticIota!(3, 16))
    {
        sc.facets!d.empty.should.equal(true);
        sc.facets(d).empty.should.equal(true);
        sc.simplices!d.empty.should.equal(true);
        sc.simplices(d).empty.should.equal(true);
    }
}

@Name("insertFacet") unittest
{
    auto sc = SimplicialComplex!()();
    sc.insertFacet([1,2]);
    assert(sc.facets.equal([[1,2]]));
    sc.insertFacet([2,3,4]);
    assert(sc.facets.equal([[1,2], [2,3,4]]));
    sc.insertFacet([1,5,6]);
    assert(sc.facets.equal([[1,2], [1,5,6], [2,3,4]]));
    sc.insertFacet([1,5,6,7]);
    assert(sc.facets.equal([[1,2], [2,3,4], [1,5,6,7]]));
}

@Name("insertFacet/removeFacet") unittest
{
    auto sc = simplicialComplex([[1,2], [1,3], [1,4], [4,5,6]]);

    assert(sc.facets.equal([[1,2], [1,3], [1,4], [4,5,6]]));
    sc.removeFacet([1,3]);
    assert(sc.facets.equal([[1,2], [1,4], [4,5,6]]));
    sc.removeFacet([1,2]);
    assert(sc.facets.equal([[1,4], [4,5,6]]));
    sc.removeFacet([4,5,6]);
    assert(sc.facets.equal([[1,4]]));
    sc.removeFacet([1,4]);
    assert(sc.facets.empty);

    auto sc2 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    sc2.removeFacet([1,2,3]);
    assert(sc2.facets.equal([[1,2,4], [1,3,4], [2,3,4]]));
    sc2.insertFacet([1,2,5]);
    sc2.insertFacet([1,3,5]);
    sc2.insertFacet([2,3,5]);
    assert(sc2.facets.equal([[1,2,4], [1,2,5], [1,3,4], [1,3,5], [2,3,4],
        [2,3,5]]));

    SimplicialComplex!int sc3;
    sc3.insertFacet([1,2,3]);
    assert(sc3.facets.equal([[1,2,3]]));
    sc3.removeFacet(simplex(1,2,3));
    assert(sc3.facets.empty);
}

/*******************************************************************************
A simplicial complex type whose vertices are of type `Vertex`.
*/
struct SimplicialComplex(Vertex = int)
{
private:
    /* We store arrays containing all the vertices in the facets of a given 
    dimension. These arrays are flattened, so the array storing the 2-dim
    facets has the form [v1, v2, v3, w1, w2, w3, ... ] where [v1, v2, v3] give
    the first simplex, [w1, w2, w3] the second, etc.
    */
    SmallMap!(int, Vertex[]) facetVertices;
public:
    /***************************************************************************
    The type of the vertices in this simplicial complex.
    */
    alias VertexType = Vertex;

    /***************************************************************************
    Construct a simplicial complex from an array of vertex arrays, specifying
    the facets.
    */
    this(const(Vertex)[][] facets)
    {
        // Check that none of the arrays reference the same memory
        assert(facets.map!(f => f.ptr).array.sort().uniq.walkLength == facets.length);

        // TO DO: Maybe check that no facets are faces of any others?
        facets.each!(f => insertFacet(f));
    }

    this(R)(R facets) if (isInputRange!R 
        && isInputRange!(ElementType!R) 
        && is(ElementType!(ElementType!R) : VertexType))
    {
        facets.each!(f => insertFacet(f.array));
    }

    // Postblit makes sure copy doesn't share data
    this(this)
    {
        facetVertices = facetVertices.dup;
        facetVertices.keys.each!(k => facetVertices[k] = facetVertices[k].dup);
    }

    /***************************************************************************
    Inserts the simplex s as a facet in the simplicial complex.
    */
    void insertFacet(int dim)(const Simplex!(dim, Vertex) s)
    {
        s.vertices.assertValidSimplex(dim);
        insertFacet(s.vertices.dup);
    }

    /***************************************************************************
    Inserts a facet (given as an input range of vertices) into the simplicial
    complex.
    */
    void insertFacet(V)(V vertices) if (isInputRange!V)
    {
        // TO DO: Improve this! Why does this allocate a closure?

        int dim = vertices.walkLength.to!int - 1;
        assert(dim >= 0);
        vertices.assertValidSimplex(dim);

        static assert(is(Unqual!(ElementType!V) == Vertex));
        assert(!this.contains(vertices),
            "expected a simplex not already in the simplicial complex");

        // First we remove any existing facets which are faces of inserted facet
        // TO DO: Improve this?
        vertices.subsets
            .filter!(vSet => this.facets.canFind(vSet.array))
            .each!(vSet => this.removeFacet(vSet.array));


        if (dim in facetVertices)
        {
            facetVertices[dim] ~= vertices.array.dup;            
        }
        else
        {
            facetVertices.insert(dim, vertices.array.dup);
        }
        
        // TO DO: Improve this sorting function? Seems yucky!
        facetVertices[dim][] = facetVertices[dim].chunks(dim + 1)
            .array.sort().joiner.array[];
    }

    /***************************************************************************
    Removes the simplex `s` as a facet in the simplicial complex.
    */
    void removeFacet(int dim)(const Simplex!(dim, Vertex) s)
    {
        removeFacet(s.vertices.dup);
    }

    /***************************************************************************
    Removes an existing facet of the simplicial complex a facet (given as an 
    input range of vertices)
    */
    void removeFacet(V)(V vertices) if (isInputRange!V)
    {
        assert(this.contains(vertices), "tried to remove a facet not in the simplicial complex");
        
        int dim = vertices.walkLength.to!int - 1;

        auto indx = facetVertices[dim].chunks(dim + 1).countUntil(vertices);
        copy(facetVertices[dim][(dim + 1) * (indx + 1)  .. $],
             facetVertices[dim][(dim + 1) * indx .. $ - (dim + 1)]);
        
        facetVertices[dim] = facetVertices[dim][0 .. $ - (dim + 1)];
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension.
    */
    auto facets(int dim)() const
    {
        static assert(dim >= 0, "expected a non-negative dimension but got "
            ~ "dimension " ~ dim.to!string);
        return facets(dim).map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension as an 
    array of arrays of vertices.
    */
    auto facets(int dim) const
    {
        assert(dim >= 0, "expected a non-negative dimension");

        if (dim !in facetVertices)
        {
            typeof(facetVertices[dim]) empty;
            return empty.chunks(dim + 1);
        }
        else
        {
            return facetVertices[dim].chunks(dim + 1);
        }
    }

    /***************************************************************************
    Returns the facets of the simplicial complex. These are simplicies that 
    are not the face of another simplex. They are returned in increasing order 
    of dimension and in lexicographic order within dimensions.
    */
    auto facets() const
    {
        // TO DO: Needed capture here to make this @nogc. Clean up capture!
        return facetVertices.keys.capture(&this)
            .map!(d => d.d0.facets(d)).joiner;
    }

    /***************************************************************************
    Returns the number of facets
    */
    size_t numFacets() const pure nothrow @nogc @safe
    {
        return this.facets.walkLength;
    }

    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()`
    */
    auto link(int dim, V)(const Simplex!(dim, V) s) const
    {
        assert(contains(s), "expected a simplex in the simplicial complex");

        // TO DO: Make this @nogc
        // Why can't I use f.d0 in place of the simplex captured by star?

        // These work...
        // return this.star(s).map!(f => setDifference(f, f.d0.vertices.array));
        return this.star(s).map!(f => setDifference(f, s.vertices.array));

        // But this one causes lots of unittest failtures (aliasing?)
        // return this.star(s).map!(f => setDifference(f, f.d0.vertices));

        // return this.facets.capture(s)
        //     .filter!(f => f.d0.vertices.isSubsetOf(f))
        //     .map!(f => setDifference(f, f.d0.vertices));
    }
    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()`
    */
    auto link(V)(V vertices) const if (isInputRange!V)
    {
        assert(contains(vertices), "expected a simplex in the simplicial complex");

        // Here f.d0 is vertices since star has captured it. TO DO: Clean this up!
        return this.star(vertices).map!(f => setDifference(f, f.d0));
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets()
    */ 
    auto star(int dim, V)(const Simplex!(dim, V) s) const
    {
        return this.facets.capture(s).filter!(f => f.d0.vertices.isSubsetOf(f));
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets()
    */ 
    auto star(V)(V v) const if (isInputRange!V)
    {
        return this.facets.capture(v).filter!(f => f.d0.isSubsetOf(f));
    }

    /***************************************************************************
    Returns true if the simplex `s` is in this simplicial complex and false 
    otherwise
    */
    bool contains(int dim, V)(const Simplex!(dim, V) s) const
    {
        return this.contains(s.vertices);
    }

    /***************************************************************************
    Returns true if the simplex with vertices given by the input range 
    `vertices` is in this simplicial complex and false otherwise
    */
    bool contains(V)(V vertices) const if (isInputRange!V)
    {
        static assert(is(Unqual!(ElementType!V) == Vertex));       
        return this.facets.any!(f => vertices.isSubsetOf(f));
    }

    /***************************************************************************
    Get simplices of dimension `dim`
    */
    auto simplices(int dim)() const
    {
        // TO DO: Reduce gc presure here. (Can't make @nogc I think.)

        static assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);
        return simplices(dim).map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /***************************************************************************
    Get the simplices of dimension `dim` as lists of vertices.
    */
    auto simplices(int dim) const
    {
        // TO DO: Reduce gc presure here. (Can't make @nogc I think.)

        assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);

        Vertex[][] simplicesSeen;
        auto dims = facetVertices.keys;
        foreach(d; dims.filter!(d => d >= dim))
        {
            foreach(f; this.facets(d))
            {
                simplicesSeen ~= f.subsetsOfSize(dim + 1)
                    .map!(s => s.array.dup).array;
            }
        }
        return simplicesSeen.sort().uniq;
    }
}

/*******************************************************************************
Returns a nice looking representation of the simplicial complex as a string.
*/
string toString(Vertex)(const ref SimplicialComplex!Vertex sc)
{
    return sc.facets.to!string;
}

/*******************************************************************************
Get the f-vector of the simplicial complex. The returned array lists the
number of simplices in each dimension.
*/
int[] fVector(Vertex)(const ref SimplicialComplex!Vertex sc)
{   
    immutable maxDim = sc.facetVertices.keys.maxElement.to!int;
    return iota(maxDim + 1).map!(dim => sc.simplices(dim).walkLength.to!int)
        .array; 
}
///
@Name("fVector")
unittest
{
    SimplicialComplex!() sc;
    sc.insertFacet([1,2]);
    assert(sc.fVector == [2,1]);        
}

/*******************************************************************************
Helper function template returning a newly constructed simplicial complex from
an array of facets (given as arrays of vertices.)
*/
SimplicialComplex!Vertex simplicialComplex(Vertex)(const(Vertex)[][] initialFacets)
{
    return SimplicialComplex!Vertex(initialFacets);
}
///
@Name("simplicialComplex")
unittest
{
    auto sc = simplicialComplex([[1,2], [2,3], [3,4,5], [6,7,8]]);
    assert(sc.facets.equal([[1,2], [2,3], [3,4,5], [6,7,8]]));

    Assert.equal(sc, simplicialComplex(sc.facets.array));

    const(int)[][] noFacets;
    Assert.equal(simplicialComplex(noFacets).facets.empty, true);
}

///
@Name("facets(dim) (pure nothrow @nogc @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];
    int[3] triangle = [2,4,5];

    () pure nothrow @nogc @safe {
        auto facs1 = sc.facets(1);
        auto saved1 = facs1.save;

        auto facs2 = sc.facets(2);
        auto saved2 = facs2.save;

        assert(facs1.front == edge[]);
        facs1.popFront;
        assert(facs1.empty);    
        assert(!saved1.empty);    
        
        assert(facs2.front == triangle[]);
        facs2.popFront;
        assert(facs2.empty);    
        assert(!saved2.empty);    
    }();
}

///
@Name("facets (pure nothrow @nogc @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];
    int[3] triangle = [2,4,5];

    () pure nothrow @nogc @safe {
        auto facs = sc.facets;
        auto saved = facs.save;
        assert(facs.front == edge[]);
        facs.popFront;
        assert(facs.front == triangle[]);
        facs.popFront;
        assert(facs.empty);    
        assert(!saved.empty);    
    }();
}

///
@Name("star(range) (pure nothrow @nogc @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];

    () pure nothrow @nogc @safe {
        auto starRange = sc.star(edge[]);
        auto saved = starRange.save;

        assert(starRange.front == edge[]);

        // Can still access the captured simplex
        assert(starRange.front.d0 == edge[]);

        starRange.popFront;
        assert(starRange.empty);
        assert(!saved.empty);
    }();
}

///
@Name("star(simplex) (pure nothrow @nogc @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] facetAns = [1,2];

    () pure nothrow @nogc @safe {
        auto starRange = sc.star(simplex(1,2));
        auto saved = starRange.save;

        assert(starRange.front == facetAns[]);

        // Can still access the captured simplex
        assert(starRange.front.d0 == simplex(1,2));

        starRange.popFront;
        assert(starRange.empty);
        assert(!saved.empty);
    }();
}

///
@Name("link(range) (pure nothrow @nogc @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[1] vertex = [1];

    () pure nothrow @nogc @safe {
        auto linkRange = sc.link(vertex[]);
        auto saved = linkRange.save;

        assert(linkRange.front.front == 2);

        // Can still access the captured simplex
        // assert(linkRange.front.d0 == vertex[]);

        linkRange.popFront;
        assert(linkRange.empty);
        assert(!saved.empty);
    }();  
}

///
@Name("link(simplex) (pure nothrow @safe)") unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    () pure nothrow /* @nogc */ @safe {
        assert(sc.link(simplex(1)).front.front == 2);
    }();
}
