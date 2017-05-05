import simplex : Simplex, simplex;

version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
}

/// A simplicial complex type ... yada yada
struct SimplicialComplex
{
    // Inserts a facet into the simplicial complex.
    void insertFacet(int[] vertices) pure nothrow @safe
    in
    {
        import std.algorithm : isSorted, findAdjacent;
        assert(vertices.isSorted, "vertices must be sorted");
        assert(vertices.findAdjacent.length == 0, "repeated vertices not allowed");
    }
    body
    {
        import std.range.primitives : walkLength;
        immutable numVerts = vertices.walkLength;
        auto ptr = numVerts in facets;
        if (ptr is null)
        {
            facets[numVerts] ~= vertices;
        }
        else
        {
            import std.algorithm : canFind;
            assert(!facets[numVerts].canFind(vertices), "facet must not already exist");
            *ptr ~= vertices;
        }
    }

    // Returns the number of facets
    size_t numFacets()
    {
        import std.algorithm : map, sum;
        return facets.byValue.map!(facetList => facetList.length).sum;
    }


    private:
    // Facets indexed by number of vertices
    int[][][size_t] facets;
    
}

unittest
{
    SimplicialComplex sc;

    // insert a 2-simplex [1,5,7]
    sc.insertFacet([1,5,7]);
    sc.insertFacet([1,5,8]);
    sc.insertFacet([7,9]);

    // can get the number of facets
    assert(sc.numFacets == 3);

    sc.numFacets.writeln;

    // vertices in an inserted facet must be sorted
    assertThrown!Error(sc.insertFacet([1,5,2,3]));

    // vertices may not be repeated
    assertThrown!Error(sc.insertFacet([1,3,3]));

    // cannot insert and already existing facet again
    assertThrown!Error(sc.insertFacet([7,9]));

    writeln(sc.facets);
}