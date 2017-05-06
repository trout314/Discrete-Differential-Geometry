version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
}

bool hasFace(Simplex s1, Simplex s2)
{
    import std.algorithm : all, canFind;
    return s2.all!(vertex => s1.canFind(vertex));
}

private alias Simplex = int[];

/// A simplicial complex type ... yada yada
struct SimplicialComplex
{
    // Inserts a facet into the simplicial complex.
    void insertFacet(Simplex facet) pure nothrow @safe
    in
    {
        import std.algorithm : isSorted, findAdjacent;
        assert(facet.isSorted, "vertices must be sorted");
        assert(facet.findAdjacent.length == 0,
            "repeated vertices not allowed");
    }
    body
    {
        import std.range.primitives : walkLength;
        immutable numVerts = facet.walkLength;
        auto ptr = numVerts in facets_;
        if (ptr is null)
        {
            facets_[numVerts] ~= facet;
        }
        else
        {
            import std.algorithm : canFind;
            assert(!facets_[numVerts].canFind(facet),
                "facet must not already exist");
            *ptr ~= facet;
        }

        import std.algorithm : sort;
        sort(facets_[numVerts]);

        // TO DO: Why doesn't sort(*ptr) work?
    }

    // Returns the facets of the simplicial complex.
    // These are simplicies that are not the face of another simplex.
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
    auto link(Simplex simplex)
    {
    }

    // Returns the star of the given simplex as a list of facets
    // in same order as they appear in facets() 
    auto star(Simplex simplex)
    {
        import std.algorithm : filter;
        return facets.filter!(facet => facet.hasFace(simplex));
    }

    bool contains(Simplex simplex)
    {
        import std.algorithm : any;
        return facets.any!(facet => facet.hasFace(simplex));
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    Simplex[][size_t] facets_;
    
}

unittest
{
    SimplicialComplex sc;

    // insert a 2-simplex [1,5,7]
    sc.insertFacet([1,5,8]);
    sc.insertFacet([1,5,7]);
    sc.insertFacet([7,9]);

    assert(sc.facets == [[7, 9], [1,5,7], [1,5,8]]);

    // can check for the presence of simplices
    assert(sc.contains([7,9]));
    assert(!sc.contains([1,5,6]));
    assert(sc.contains([1,5]));
    assert(sc.contains([7]));
    assert(!sc.contains([8, 9]));

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