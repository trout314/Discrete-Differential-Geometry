import std.algorithm : all, any, canFind, filter, map, sort;
import std.exception : enforce;
import std.conv : to;
import std.range : array, walkLength;

import simplex : Simplex, simplex;

import std.range : array;
import std.algorithm : joiner, sort, map;


version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
}



/// A simplicial complex type ... yada yada
struct SimplicialComplex(Vertex = int)
{
    alias VertexType = Vertex;

    // Inserts a facet into the simplicial complex.
    void insertFacet(int dim)(Simplex!(dim, Vertex) simplex)
    {
        enforce(!facets.any!(f => f.hasFace_Verts(simplex.vertices)),
            "inserted simplex must not be a face of an existing facet");

        facetLists[dim + 1] ~= simplex.vertices.dup;

        sort(facetLists[dim + 1]);
    }

    // Returns the facets of the simplicial complex.
    // These are simplicies that are not the face of another simplex.
    // They are returned in increasing order of dimension and in
    // lexicographic order within dimensions.
    // e.g. [8],[2,3],[2,4],[1,6,7]
    int[][] facets()
    {

        auto sizes = sort(facetLists.keys);        
        return sizes.map!(s => facetLists[s]).joiner.array;
    }

    // Returns the number of facets
    size_t numFacets()
    {
        return facets.walkLength;
    }

    // Returns the link of the given simplex as a list of facets
    // in same order as they appear in facets()
    int[][] link(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.star(simplex)
            .map!(facet => facet.oppositeFace_Verts(simplex.vertices)).array;
    }

    // Returns the star of the given simplex as a list of facets
    // in same order as they appear in facets() 
    int[][] star(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.filter!(f => f.hasFace_Verts(simplex.vertices)).array;
    }

    // Returns true if simplex is in this simplicial complex and false otherwise
    bool contains(int dim)(auto ref const Simplex!(dim, Vertex) simplex)
    {
        return this.facets.any!(f => f.hasFace_Verts(simplex.vertices));
    }

    string toString()
    {
        return "SComp" ~ this.facets.to!string;
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    int[][][size_t] facetLists;

    // TO DO: Make the above into an array of (size_t, Simplex[])
    // pairs. Probably won't ever have enough distince facet
    // dimenstion to justify using an associative array
}

bool hasFace_Verts(const(int)[] s1, const(int)[] s2)
{
    return s2.all!(vertex => s1.canFind(vertex));
}

auto oppositeFace_Verts(const(int)[] simplex, const(int)[] face)
{
    return simplex.filter!(v => !face.canFind(v)).array.to!(int[]);
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

    // get list of facets back, in order of increasing dimension
    // and dictionary order within a dimension 
    assert(sc.facets == [[4,5], [5,6], [1,2,3], [2,3,4]]);
    assert(sc.numFacets == 4);

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

    // it is an error to insert a simplex that is
    // already the face of an existing facet
    assertThrown(sc.insertFacet(s(4,5)));
    assertThrown(sc.insertFacet(s(2)));
}