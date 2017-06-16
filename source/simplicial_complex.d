import std.algorithm : all, any, canFind, filter, find, joiner, map, maxElement,
    setDifference, sort, uniq;
import std.conv : to;
import std.exception : enforce, assertThrown;
import std.range : array, ElementType, front, iota, isInputRange, walkLength;
import std.typecons : Tuple, tuple;

import simplex : Simplex, simplex, hasFace, facesOfDim;
import utility : isSubsetOf, subsetsOfSize, SmallMap, throwsWithMsg;

/// A simplicial complex type ... yada yada
struct SimplicialComplex(Vertex = int)
{
    alias VertexType = Vertex;

    // Inserts a facet into the simplicial complex.
    void insertFacet(int dim)(Simplex!(dim, Vertex) simplex)
    {
        enforce(!this.contains(simplex), "insertFacet expects an inserted "
            ~ "simplex not already in the simplicial complex, but got " 
            ~ simplex.toString ~ " and already have facet " 
            ~ facets.find!(f => simplex.vertices.isSubsetOf(f))
                    .front.to!string);

        facetLists[dim + 1] ~= simplex.vertices.dup;

        facetLists[dim + 1].sort();
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension.
    */
    auto facets(int dim)()
    {
        static assert(dim > 0);
        return facetLists[dim + 1].map!(vList => Simplex!(dim, Vertex)(vList));
    }

    /***************************************************************************
    Returns the facets of the simplicial complex. These are simplicies that 
    are not the face of another simplex. They are returned in increasing order 
    of dimension and in lexicographic order within dimensions. */
    int[][] facets() pure @safe
    {

        auto sizes = facetLists.keys.array.sort();        
        return sizes.map!(s => facetLists[s]).joiner.array;
    }

    /***************************************************************************
    Returns the number of facets */
    size_t numFacets() pure @safe
    {
        return facets.walkLength;
    }

    /***************************************************************************
    Returns the link of the given simplex as a list of facets, given in same 
    order as they appear in facets() */
    int[][] link(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.star(simplex)
            .map!(facet => setDifference(facet, simplex.vertices).array).array;
    }

    /***************************************************************************
    Returns ...*/
    int[][] link(V)(V vertices) if (isInputRange!V)
    {
        static assert(is(ElementType!V == Vertex));
        return this.star(vertice).map!(
            facet => setDifference(facet, vertices).array).array;

        return [];
    }


    /* Returns the star of the given simplex as a list of facets in same order 
    as they appear in facets() */ 
    int[][] star(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.filter!(f => simplex.vertices.isSubsetOf(f)).array;
    }
  
    // Returns true if simplex is in this simplicial complex and false otherwise
    bool contains(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.any!(f => simplex.vertices.isSubsetOf(f));
    }

    /// Get simplices of dimension `dim`
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

    ulong[] fVector()
    {   
        int maxDim = facetLists.keys.maxElement.to!int;
        return iota(maxDim).map!(dim => simplices(dim).walkLength).array; 
    } 

    string toString() @safe
    {
        return this.facets.to!string;
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    SmallMap!(size_t, int[][]) facetLists;
}

unittest
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

    // --------------------------------------------------------------
    // Restrictions
    // --------------------------------------------------------------

    // cannot insert a simplex if it is already the face of an existing facet
    sc.insertFacet(s(4,5)).assertThrown;
    sc.insertFacet(s(2)).throwsWithMsg(
        "insertFacet expects an inserted simplex not already in the simplicial "
        ~ "complex, but got [2] and already have facet [1, 2, 3]");
}