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
        auto ptr = numVerts in facets_;
        if (ptr is null)
        {
            facets_[numVerts] ~= vertices;
        }
        else
        {
            import std.algorithm : canFind;
            assert(!facets_[numVerts].canFind(vertices), "facet must not already exist");
            *ptr ~= vertices;
        }

        import std.algorithm : sort;
        sort(facets_[numVerts]);

        // TO DO: Why doesn't sort(*ptr) work?
    }

    // Returns the facets of the simplicial complex. These are simplicies that are not the face of another simplex.
    auto facets()
    {
        import std.range : array;
        import std.algorithm : joiner, sort, map;

        auto sizes = sort(facets_.keys);        
        return sizes.map!(s => facets_[s]).joiner.array;
    }

    // Returns the number of facets
    size_t numFacets()
    {
        import std.range : walkLength;
        return facets.walkLength;
    }

    // Returns the link of the given simplex as a list of facets
    // in same order as they appear in facets()
    auto link(int[] simplex)
    {
    }

    // Returns the closure of the star of the given simplex. 
    auto star(int[] simplex)
    {
    }

    private:
    // Facets indexed by number of vertices
    int[][][size_t] facets_;
    
}

unittest
{
    SimplicialComplex sc;

    // insert a 2-simplex [1,5,7]
    sc.insertFacet([1,5,8]);
    sc.insertFacet([1,5,7]);
    sc.insertFacet([7,9]);

    assert(sc.facets == [[7, 9], [1,5,7], [1,5,8]]);

    // can get the number of facets
    assert(sc.numFacets == 3);
    writeln(sc.facets);
    // vertices in an inserted facet must be sorted
    assertThrown!Error(sc.insertFacet([1,5,2,3]));

    // vertices may not be repeated
    assertThrown!Error(sc.insertFacet([1,3,3]));

    // cannot insert and already existing facet again
    assertThrown!Error(sc.insertFacet([7,9]));

    writeln(sc.facets_);
}