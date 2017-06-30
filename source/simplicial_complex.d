import simplex : facesOfDim, hasFace, simplex, Simplex;

import std.algorithm : all, any, canFind, filter, find, joiner, map, maxElement,
    setDifference, sort, uniq;
import std.conv : to;
import std.exception : assertThrown, enforce;
import std.range : array, ElementType, front, iota, isInputRange, walkLength;
import std.traits : isArray, Unqual;
import std.typecons : Tuple, tuple;

import utility : isSubsetOf, SmallMap, subsetsOfSize, throwsWithMsg;

/// A simplicial complex type whose vertices are of type `Vertex`.
struct SimplicialComplex(Vertex = int)
{
    /***************************************************************************
    The type of the vertices in this simplicial complex. */
    alias VertexType = Vertex;

    /***************************************************************************
    Inserts the simplex s as a facet in the simplicial complex. */
    void insertFacet(int dim)(Simplex!(dim, Vertex) s)
    {
        insertFacet(s.vertices);
    }

    /***************************************************************************
    Inserts a facet (given as an input range of vertices) into the simplicial
    complex. */
    void insertFacet(V)(V vertices) if (isInputRange!V)
    {
        static assert(is(Unqual!(ElementType!V) == Vertex));
        enforce(!this.contains(vertices), "insertFacet expects an inserted "
            ~ "simplex not already in the simplicial complex, but got " 
            ~ vertices.to!string ~ " and already have facet " 
            ~ facets.find!(f => vertices.isSubsetOf(f)).front.to!string);

        auto numVertices = vertices.walkLength;
        if (numVertices in facetLists)
        {
            facetLists[numVertices] ~= vertices.dup;            
        }
        else
        {
            facetLists.insert(numVertices, [vertices.dup]);
        }
        
        facetLists[numVertices].sort();
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension. */
    auto facets(int dim)() const
    {
        static assert(dim > 0);
        return facetLists[dim + 1].map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension as an 
    array of arrays of vertices. */
    auto facets(int dim) const
    {
        assert(dim > 0);
        return facetLists[dim + 1];
    }

    /***************************************************************************
    Returns the facets of the simplicial complex. These are simplicies that 
    are not the face of another simplex. They are returned in increasing order 
    of dimension and in lexicographic order within dimensions. */
    VertexType[][] facets() pure @safe /* const */
    {
        auto sizes = facetLists.keys.array.dup.sort();        
        return sizes.map!(s => facetLists[s]).joiner.array;
    }

    /***************************************************************************
    Returns the number of facets */
    size_t numFacets() pure @safe
    {
        return this.facets.walkLength;
    }

    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()` */
    Vertex[][] link(int dim)(const Simplex!(dim, Vertex) s)
    {
        return this.star(s).map!(f => setDifference(f, s.vertices).array).array;
    }
    /***************************************************************************
    Returns the facets in the link of the simplex `s` as an array of arrays of 
    vertices, given in same order as they appear in `facets()` */
    Vertex[][] link(V)(V vertices) if (isInputRange!V)
    {
        static assert(is(ElementType!V == Vertex));
        return this.star(vertices).map!(
            f => setDifference(f, vertices).array).array;
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets() */ 
    VertexType[][] star(int dim)(const Simplex!(dim, Vertex) simplex)
    {
        return star(simplex.vertices);
    }

    /***************************************************************************
    Returns the star of the given simplex as an array of arrays of vertices of 
    the facets. These are given in the same order as specified facets() */ 
    VertexType[][] star(V)(V vertices) if (isInputRange!V || isArray!V)
    {
        return this.facets.filter!(f => vertices.isSubsetOf(f)).array;
    }

    /***************************************************************************
    Returns true if the simplex `s` is in this simplicial complex and false 
    otherwise */
    bool contains(int dim)(const Simplex!(dim, Vertex) s)
    {
        return this.contains(s.vertices);
    }

    /***************************************************************************
    Returns true if the simplex with vertices given by the input range 
    `vertices` is in this simplicial complex and false otherwise */
    bool contains(V)(V vertices) if (isInputRange!V)
    {
        static assert(is(Unqual!(ElementType!V) == Vertex));
        return this.facets.any!(f => vertices.isSubsetOf(f));
    }

    /***************************************************************************
    Get simplices of dimension `dim` */
    auto simplices(int dim)()
    {
        static assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);
        return simplices(dim).map!(verts => Simplex!(dim, Vertex)(verts));
    }

    /// Get the simplices of dimension `dim` as lists of vertices.
    auto simplices(int dim)
    {
        assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);
        Vertex[][] simplicesSeen;
        foreach(key; facetLists.keys.filter!(key => key >= dim+1))
        {
            foreach(facet; facetLists[key])
            {
                simplicesSeen ~= facet.subsetsOfSize(dim + 1).array;
            }
        }
        return simplicesSeen.sort().uniq;
    }

    /** Get the f-vector of the simplicial complex. The returned array lists the
    number of simplices in each dimension. */
    auto fVector()
    {   
        immutable maxDim = facetLists.keys.maxElement.to!int;
        return iota(maxDim).map!(dim => simplices(dim).walkLength).array; 
    } 

    string toString() @safe pure
    {
        return this.facets.to!string;
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    SmallMap!(size_t, Vertex[][]) facetLists;
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