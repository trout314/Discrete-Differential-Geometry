import std.algorithm : all, any, canFind, filter, find, joiner, map, sort;
import std.conv : to;
import std.exception : enforce, assertThrown;
import std.range : array, front, walkLength;
import std.typecons : Tuple, tuple;

import simplex : Simplex, simplex, hasFace;
import utility : SmallMap, throwsWithMsg;

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
            ~ facets.find!(f => f.hasFace_Verts(simplex.vertices))
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
            .map!(facet => facet.oppositeFace_Verts(simplex.vertices)).array;
    }

    /* Returns the star of the given simplex as a list of facets in same order 
    as they appear in facets() */ 
    int[][] star(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.filter!(f => f.hasFace_Verts(simplex.vertices)).array;
    }

    // Returns true if simplex is in this simplicial complex and false otherwise
    bool contains(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.any!(f => f.hasFace_Verts(simplex.vertices));
    }

    auto simplices(int dim)()
    {
        static assert(dim >= 0);
        Simplex!(dim, Vertex)[] simplicesSeen;
        foreach(key; facetLists.keys)
        {
            if (key >= dim + 1)
            {
                foreach(facet; facetLists[key])
                {
                    //facet.to!Simplex(dim, Vertex).facesOfDim!dim;
                }
            }
        }
        return 4;
    }


    int[] fVector()
    {   
        return [];
    } 

    string toString() @safe
    {
        return this.facets.to!string;
    }



    private:
    // Lists of facets, indexed by number of vertices in the facet
    SmallMap!(size_t, int[][]) facetLists;
}

bool hasFace_Verts(const(int)[] s1, const(int)[] s2)
{
    return s2.all!(vertex => s1.canFind(vertex));
}

auto oppositeFace_Verts(const(int)[] simplex, const(int)[] face)
{
    return simplex.filter!(v => !face.canFind(v)).array.to!(int[]);
}

auto faces_Verts(const(int)[] simplex)
{
    int[][] faces;
    foreach (bitChoice; 1 .. 2 ^^ simplex.length)
    {
        int[] face;
        foreach (index, vertex; simplex)
        {
            if ((1 << index) & bitChoice)
            {
                face ~= vertex;
            }
        }
        faces ~= face.dup;
    }
    faces.sort!((a, b) => a.length < b.length);
    return faces;
}

unittest
{
    SimplicialComplex!() sc;
    assert(sc.numFacets == 0);
    assert(sc.facets == []);

    // insert some facets
    alias s = simplex;

    sc.insertFacet(s(1,2,3));
    sc.insertFacet(s(2,3,4));
    sc.insertFacet(s(4,5));
    sc.insertFacet(s(5,6));

    /* get list of facets back, in order of increasing dimension and dictionary 
    order within a dimension */ 
    assert(sc.facets == [[4,5], [5,6], [1,2,3], [2,3,4]]);
    assert(sc.numFacets == 4);

    import std.stdio : writeln;
    writeln(sc.simplices!1);

    // get the f-vector of the simplicial complex
    import std.stdio : writeln;
    sc.fVector.writeln;
    [1,2,3].faces_Verts.writeln;
//    assert(sc.fVector == [1,2,3]);

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