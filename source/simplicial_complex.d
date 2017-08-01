import simplex : facesOfDim, hasFace, simplex, Simplex;
import std.algorithm : all, any, canFind, chunkBy, copy, countUntil, each,
    equal, filter, find, joiner, map, maxElement, remove, setDifference,
    setIntersection, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown, enforce;
import std.range : array, chunks, empty, enumerate, ElementType, front, iota,
    isInputRange, refRange, walkLength, zip;
import std.traits : isArray, Unqual;
import std.typecons : Tuple, tuple;
import utility : isSubsetOf, SmallMap, subsets, staticIota, subsetsOfSize, throwsWithMsg;
import std.stdio : writeln;
import unit_threaded : Name;
import fluent.asserts : should, Assert;

import simplicial_complex_algorithms : connectedComponents, eulerCharacteristic,
    isCircle, isConnected, isPureOfDim, isOrientableSurfaceOfGenus, join;

/// Basic Functionality
@Name("doc tests")
@safe pure unittest
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
    assert(sc.facets == [[4,5], [5,6], [1,2,3], [2,3,4]]);
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
    assert(sc.star(s(5)) == [[4,5], [5,6]]);
    assert(sc.star(s(4)) == [[4,5], [2,3,4]]);
    assert(sc.star(s(2,3)) == [[1,2,3], [2,3,4]]);

    // get link of a simplex as list of facets
    assert(sc.link(s(5)) == [[4], [6]]);
    assert(sc.link(s(4)) == [[5], [2,3]]);
    assert(sc.link(s(2,3)) == [[1], [4]]);

    // can print out a simplicial complex
    assert(sc.toString == "[[4, 5], [5, 6], [1, 2, 3], [2, 3, 4]]");

    // -------------------------------------------------------------------------
    // Restrictions
    // -------------------------------------------------------------------------

    // cannot insert a simplex if it is already the face of an existing facet
    sc.insertFacet(s(4,5)).assertThrown;
    sc.insertFacet(s(2)).throwsWithMsg("expected a simplex not already in the "
        ~ "simplicial complex, but got [2] and already have facet [1, 2, 3]");

    throwsWithMsg(simplicialComplex([[3,4]]).removeFacet([2,3]), "tried to "
        ~ "remove a facet [2, 3] not in the simplicial complex");

    throwsWithMsg(simplicialComplex([[1,3,4], [2,3,4]]).link(s(1,2)),
        "expected a simplex in the simplicial complex");

    throwsWithMsg(sc.facets(-1),"expected a non-negative dimension but got "
        ~ "dimension -1"); 


}

@Name("additional tests")
unittest
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

@Name("insertFacet")
unittest
{
    auto sc = SimplicialComplex!()();
    sc.insertFacet([1,2]);
    assert(sc.facets == [[1,2]]);
    sc.insertFacet([2,3,4]);
    assert(sc.facets == [[1,2], [2,3,4]]);
    sc.insertFacet([1,5,6]);
    assert(sc.facets == [[1,2], [1,5,6], [2,3,4]]);
    sc.insertFacet([1,5,6,7]);
    assert(sc.facets == [[1,2], [2,3,4], [1,5,6,7]]);
}

@Name("insertFacet/removeFacet")
unittest
{
    auto sc = simplicialComplex([[1,2], [1,3], [1,4], [4,5,6]]);

    assert(sc.facets == [[1,2], [1,3], [1,4], [4,5,6]]);
    sc.removeFacet([1,3]);
    assert(sc.facets == [[1,2], [1,4], [4,5,6]]);
    sc.removeFacet([1,2]);
    assert(sc.facets == [[1,4], [4,5,6]]);
    sc.removeFacet([4,5,6]);
    assert(sc.facets == [[1,4]]);
    sc.removeFacet([1,4]);
    assert(sc.facets == []);

    auto sc2 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    sc2.removeFacet([1,2,3]);
    assert(sc2.facets == [[1,2,4], [1,3,4], [2,3,4]]);
    sc2.insertFacet([1,2,5]);
    sc2.insertFacet([1,3,5]);
    sc2.insertFacet([2,3,5]);
    assert(sc2.facets == [[1,2,4], [1,2,5], [1,3,4], [1,3,5], [2,3,4],
        [2,3,5]]);

    SimplicialComplex!int sc3;
    sc3.insertFacet([1,2,3]);
    assert(sc3.facets == [[1,2,3]]);
    sc3.removeFacet(simplex(1,2,3));
    assert(sc3.facets == []);
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
    this(Vertex[][] facets)
    {
        // TO DO: Maybe check that no facets are faces of any others?
        facets.each!(f => insertFacet(f));
    }

    // Postblit makes sure copy doesn't share data
    this(this)
    {
        facetVertices = facetVertices.dup;
    }

    /***************************************************************************
    Inserts the simplex s as a facet in the simplicial complex.
    */
    void insertFacet(int dim)(const Simplex!(dim, Vertex) s)
    {
        insertFacet(s.vertices);
    }

    /***************************************************************************
    Inserts a facet (given as an input range of vertices) into the simplicial
    complex.
    */
    void insertFacet(V)(V vertices) if (isInputRange!V)
    {
        /* TO DO: Validate vertex lists. How can we do this without duplicating
        functionality from simplex.d */

        static assert(is(Unqual!(ElementType!V) == Vertex));
        enforce(!contains(vertices), "expected a simplex not already in the "
            ~ "simplicial complex, but got " ~ vertices.to!string 
            ~ " and already have facet " 
            ~ facets.find!(f => vertices.isSubsetOf(f)).front.to!string);

        // First we remove any existing facets which are faces of inserted facet
        // TO DO: Improve this?
        vertices.subsets
            .filter!(vSet => this.facets.canFind(vSet))
            .each!(vSet => this.removeFacet(vSet));

        int dim = vertices.walkLength.to!int - 1;

        if (dim in facetVertices)
        {
            facetVertices[dim] ~= vertices.dup;            
        }
        else
        {
            facetVertices.insert(dim, vertices.dup);
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
        removeFacet(s.vertices);
    }

    /***************************************************************************
    Removes an existing facet of the simplicial complex a facet (given as an 
    input range of vertices)
    */
    void removeFacet(V)(V vertices) if (isInputRange!V)
    {
        enforce(this.contains(vertices), "tried to remove a facet "
            ~ vertices.to!string ~ " not in the simplicial complex");

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
        enforce(dim >= 0, "expected a non-negative dimension but got dimension "
            ~ dim.to!string);
        if (dim !in facetVertices)
        {
            Vertex[][] empty;
            return empty;
        }
        else
        {
            return facetVertices[dim].chunks(dim + 1).array.to!(Vertex[][]);
        }
    }

    /***************************************************************************
    Returns the facets of the simplicial complex. These are simplicies that 
    are not the face of another simplex. They are returned in increasing order 
    of dimension and in lexicographic order within dimensions.
    */
    VertexType[][] facets() const pure @safe
    {
        auto dims = facetVertices.keys.array;
        return dims.map!(d => facets(d)).joiner.array;
    }

    /***************************************************************************
    Returns the number of facets
    */
    size_t numFacets() const pure @safe
    {
        return this.facets.walkLength;
    }

    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()`
    */
    Vertex[][] link(int dim)(const Simplex!(dim, Vertex) s) const
    {
        enforce(contains(s), "expected a simplex in the simplicial complex");
        return this.star(s).map!(f => setDifference(f, s.vertices).array).array;
    }
    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()`
    */
    Vertex[][] link(V)(V vertices) const if (isInputRange!V)
    {
        static assert(is(Unqual!(ElementType!V) == Vertex));
        return this.star(vertices).map!(
            f => setDifference(f, vertices).array).array;
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets()
    */ 
    VertexType[][] star(int dim)(const Simplex!(dim, Vertex) simplex) const
    {
        return star(simplex.vertices);
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets()
    */ 
    VertexType[][] star(V)(V vertices) const if (isInputRange!V || isArray!V)
    {
        return this.facets.filter!(f => vertices.isSubsetOf(f)).array;
    }

    /***************************************************************************
    Returns true if the simplex `s` is in this simplicial complex and false 
    otherwise
    */
    bool contains(int dim)(const Simplex!(dim, Vertex) s) const
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
        static assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);
        return simplices(dim).map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /***************************************************************************
    Get the simplices of dimension `dim` as lists of vertices.
    */
    auto simplices(int dim) const
    {
        assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);
        Vertex[][] simplicesSeen;
        auto dims = facetVertices.keys.array;
        foreach(key; dims.filter!(d => d >= dim))
        {
            foreach(facet; facetVertices[key].chunks(key + 1))
            {
                simplicesSeen ~= facet.subsetsOfSize(dim + 1)
                    .map!(f => f.dup).array;
            }
        }
        return simplicesSeen.sort().uniq;
    }
}

/*******************************************************************************
Returns a nice looking representation of the simplicial complex as a string.
*/
string toString(Vertex)(const SimplicialComplex!Vertex sc)
{
    return sc.facets.to!string;
}

/*******************************************************************************
Get the f-vector of the simplicial complex. The returned array lists the
number of simplices in each dimension.
*/
int[] fVector(Vertex)(const SimplicialComplex!Vertex sc)
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
SimplicialComplex!Vertex simplicialComplex(Vertex)(Vertex[][] initialFacets)
{
    return SimplicialComplex!Vertex(initialFacets);
}
///
@Name("simplicialComplex")
unittest
{
    auto sc = simplicialComplex([[1,2], [2,3], [3,4,5], [6,7,8]]);
    assert(sc.facets == [[1,2], [2,3], [3,4,5], [6,7,8]]);

    Assert.equal(sc, simplicialComplex(sc.facets));

    int[][] noFacets;
    Assert.equal(simplicialComplex(noFacets).facets.empty, true);
}

// TO DO: Orientability tester! Fix genus stuff once that's done...