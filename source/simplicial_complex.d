import simplex : facesOfDim, hasFace, simplex, Simplex;

import std.algorithm : all, any, canFind, chunkBy, countUntil, each, filter, find, joiner, map, 
    maxElement, remove, setDifference, sort, sum, uniq;
import std.conv : to;
import std.exception : assertThrown, enforce;
import std.range : array, chunks, empty, enumerate, ElementType, front, iota, isInputRange, refRange, walkLength, zip;
import std.traits : isArray, Unqual;
import std.typecons : Tuple, tuple;

import utility : isSubsetOf, SmallMap, subsets, subsetsOfSize, throwsWithMsg;

import std.stdio : writeln;

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
        facets.each!(f => insertFacet(f));
    }
    ///
    unittest
    {
        auto sc = SimplicialComplex!int([[1,2], [2,3], [3,4,5]]);
        assert(sc.facets == [[1,2], [2,3], [3,4,5]]);
    }

    /***************************************************************************
    Inserts the simplex s as a facet in the simplicial complex.
    */
    void insertFacet(int dim)(const Simplex!(dim, Vertex) s)
    {
        insertFacet(s.vertices);
    }
    ///
    unittest
    {
        alias s = simplex;

        SimplicialComplex!int sc;
        sc.insertFacet(s(1,2,3));
        assert(sc.facets!2.array == [s(1,2,3)]);
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
        enforce(!contains(vertices), "insertFacet expects an inserted simplex "
            ~ "not already in the simplicial complex, but got " 
            ~ vertices.to!string ~ " and already have facet " 
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
    ///
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

    /***************************************************************************
    Removes the simplex `s` as a facet in the simplicial complex.
    */
    void removeFacet(int dim)(const Simplex!(dim, Vertex) s)
    {
        removeFacet(s.vertices);
    }
    ///
    unittest
    {
        SimplicialComplex!int sc;
        sc.insertFacet([1,2,3]);
        assert(sc.facets == [[1,2,3]]);
        sc.removeFacet(simplex(1,2,3));
        assert(sc.facets == []);
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
        facetVertices[dim][indx * (dim + 1) .. $ - (dim + 1)]
            = facetVertices[dim][(indx + 1) * (dim + 1) .. $];
        facetVertices[dim] = facetVertices[dim][0 .. $ - (dim + 1)];
    }
    ///
    unittest
    {
        auto sc = SimplicialComplex!()([[1,2], [1,3], [1,4], [4,5,6]]);

        assert(sc.facets == [[1,2], [1,3], [1,4], [4,5,6]]);
        sc.removeFacet([1,3]);
        assert(sc.facets == [[1,2], [1,4], [4,5,6]]);
        sc.removeFacet([1,2]);
        assert(sc.facets == [[1,4], [4,5,6]]);
        sc.removeFacet([4,5,6]);
        assert(sc.facets == [[1,4]]);
        sc.removeFacet([1,4]);
        assert(sc.facets == []);

        throwsWithMsg(sc.removeFacet([2,3]), "tried to remove a facet [2, 3] "
            ~ "not in the simplicial complex");
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension.
    */
    auto facets(int dim)() const
    {
        static assert(dim >= 0, "facets expected a non-negative dimension but "
            ~ "got dimension " ~ dim.to!string);
        return facets(dim).map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension as an 
    array of arrays of vertices.
    */
    auto facets(int dim) const
    {
        enforce(dim >= 0, "facets expected a non-negative dimension but "
            ~ "got dimension " ~ dim.to!string);
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
        static assert(is(ElementType!V == Vertex));
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

    /***************************************************************************
    Get the f-vector of the simplicial complex. The returned array lists the
    number of simplices in each dimension.
    */
    auto fVector() const
    {   
        immutable maxDim = facetVertices.keys.maxElement.to!int;
        return iota(maxDim + 1).map!(dim => simplices(dim).walkLength).array; 
    }
    ///
    unittest
    {
        SimplicialComplex!() sc;
        sc.insertFacet([1,2]);
        assert(sc.fVector == [2,1]);        
    }

    /***************************************************************************
    Returns the Euler characteristic of the simplicial complex
    */
    auto eulerCharacteristic() const
    {
        return fVector.enumerate.map!(f => (-1)^^f.index * f.value).sum;
    }
    ///
    unittest
    {
        SimplicialComplex!() sc;
        sc.insertFacet([1,2]);
        assert(sc.eulerCharacteristic == 1);
        sc.insertFacet([2,3]);
        sc.insertFacet([2,4]);
        assert(sc.eulerCharacteristic == 1);
        sc.insertFacet([3,4]);
        assert(sc.eulerCharacteristic == 0);
        sc.insertFacet([1,3]);
        assert(sc.eulerCharacteristic == -1);        
    }


    // TO DO: Finish this....
    auto connectedComponents()
    {
        static struct FacetRecord
        {
            Vertex[] facet;
            int label;      // labels which connected component facet is in
            bool seenNear;  // looked at facets in the star of the vertices?
        }

        /* We will label each facet with an integer 1, 2, 3 ... etc, marking its
        connected component */
        // auto labels = new int[numFacets];
        // auto allFacets = facets.array;


        /* We keep track of the facets we've yet to examine, the facets we've 
        labeled (but are not yet done with) and those that are completed (which
        means both labeled and all the facets in the star of its vertices are
        also labeled. */
        FacetRecord[] records = facets.map!(f => FacetRecord(f)).array;

        int currentLabel = 1;

        // Deal with the facets that are (necessarily isolated) vertices
        auto nVertFac = facets!0.walkLength;
        auto vertexFacets = records[0 .. nVertFac];

        foreach(indx, ref f; vertexFacets)
        {
            f.label = currentLabel;
            f.seenNear = true;
            ++currentLabel;
        }

        while(!records.find!(r => r.label == 0).empty)
        {
            records.find!(r => r.label == 0).front.label = currentLabel;

            // Now we propagate the current label as far as we can
            while(records.canFind!(r => r.label != 0 && !r.seenNear))
            {
                auto toDo = records.find!(r => r.label != 0 && !r.seenNear);

                foreach(f; toDo.front.facet.map!(v => star(simplex(v))).joiner)
                {
                    auto findIt = records.find!(r => r.facet == f);
                    assert(!findIt.empty);
                    findIt.front.label = currentLabel;
                }

                toDo.front.seenNear = true;           
            }
            currentLabel += 1;
        }

        return records.chunkBy!((r1, r2) => r1.label == r2.label)
            .map!(rList => SimplicialComplex(rList.map!(r => r.facet).array));           
    }
    ///
    unittest
    {
        auto sc = SimplicialComplex!()([[1],[9], [2,3],[3,4],[5,6],[6,7,8]]);

        auto c1 = SimplicialComplex!()([[1]]);
        auto c2 = SimplicialComplex!()([[9]]);
        auto c3 = SimplicialComplex!()([[2,3], [3,4]]);
        auto c4 = SimplicialComplex!()([[5,6], [6,7,8]]);

        assert(sc.connectedComponents.array == [c1, c2, c3, c4]);       
    }


    /***************************************************************************
    Returns a nice looking representation of the simplicial complex as a string.
    */
    auto toString() const
    {
        return this.facets.to!string;
    }
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a 1-sphere (circle)
*/
bool isCircle(Vertex)(SimplicialComplex!Vertex sc)
{
    // TO DO: What about connected components?

    // All facets must be edges and the star of each vertex must contain 2 edges
    return sc.facets.all!(f => f.walkLength == 2)
        && sc.simplices!0.all!(v => sc.star(v).walkLength == 2);
}
///
unittest
{
    // Start with a circle with 4 edges
    auto s = SimplicialComplex!()([[1,2], [2,3], [3,4], [1,4]]);
    assert(s.isCircle);

    s.removeFacet([2,3]);
    assert(!s.isCircle);

    // Simplicial complex must be homeomorphic to a circle, not just homotopic
    s.insertFacet([2,3,5]);
    assert(!s.isCircle);

    // WARNING: Currently returns true for disjoint copies of the circle.
    // TO DO: Fix this
    auto disjointCircles = SimplicialComplex!()([
        [1,2], [2,3], [1,3],
        [4,5], [5,6], [4,6]]);

    assert(disjointCircles.isCircle);
}



///
unittest
{
    import simplicial_complex_test : test;
    assert(test!SimplicialComplex);
}

///
@safe pure unittest
{
    // create an empty simplicial complex
    SimplicialComplex!() sc;
    assert(sc.numFacets == 0);
    assert(sc.facets == []);

    // for clarity
    alias s = simplex;

    // insert some facets
    sc.insertFacet(s(1,2,3));
    sc.insertFacet(s(2,3,4));
    sc.insertFacet(s(4,5));
    sc.insertFacet(s(5,6));

    /* get the vertices in the facets, returned in order of increasing dimension
    and dictionary order within a dimension */ 
    assert(sc.facets == [[4,5], [5,6], [1,2,3], [2,3,4]]);
    assert(sc.numFacets == 4);

    // Can get a sorted array of all the facet simplices of a given dimension
    assert(sc.facets!1.array == [s(4,5), s(5,6)]);

    // Can get a sorted array of all the simplices of a given dimension
    assert(sc.simplices!1.array == [s(1,2), s(1,3), s(2,3), s(2,4), s(3,4), 
        s(4,5), s(5,6)]);

    // get the f-vector of the simplicial complex
    assert(sc.fVector == [6,7,2]);

    // check for the presence of simplices
    assert(sc.contains(s(4,5)));
    assert(sc.contains(s(2,3,4)));
    assert(sc.contains(s(6)));
    assert(!sc.contains(s(1,2,5)));
    assert(!sc.contains(s(4,6)));

    // get the star of a simplex as list of facets
    assert(sc.star(s(5)) == [[4,5], [5,6]]);
    assert(sc.star(s(4)) == [[4,5], [2,3,4]]);
    assert(sc.star(s(2,3)) == [[1,2,3], [2,3,4]]);

    // get link of a simplex as list of facets
    assert(sc.link(s(5)) == [[4], [6]]);
    assert(sc.link(s(4)) == [[5], [2,3]]);
    assert(sc.link(s(2,3)) == [[1], [4]]);

    assert(sc.toString == "[[4, 5], [5, 6], [1, 2, 3], [2, 3, 4]]");

    // -------------------------------------------------------------------------
    // Restrictions
    // -------------------------------------------------------------------------

    // cannot insert a simplex if it is already the face of an existing facet
    sc.insertFacet(s(4,5)).assertThrown;
    sc.insertFacet(s(2)).throwsWithMsg(
        "insertFacet expects an inserted simplex not already in the simplicial "
        ~ "complex, but got [2] and already have facet [1, 2, 3]");
}