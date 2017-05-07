version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
    import std.range : array;
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
    void insertFacet(Simplex simplex)
    in
    {
        import std.algorithm : isSorted, findAdjacent, any;
        assert(simplex.isSorted, "vertices in a simplex must be sorted");
        assert(simplex.findAdjacent.length == 0,
            "repeated vertices in a simplex not allowed");
        assert(!facets.any!(f => f.hasFace(simplex)),
            "inserted facet must not be a face of an existing facet");
    }
    body
    {
        import std.range : walkLength;
        immutable numVerts = simplex.walkLength; 
        facets_[numVerts] ~= simplex;

        import std.algorithm : sort;
        sort(facets_[numVerts]);
    }

    // Returns the facets of the simplicial complex.
    // These are simplicies that are not the face of another simplex.
    // They are returned in increasing order of dimension and in
    // lexicographic order within dimensions.
    // e.g. [1],[2,3],[2,4],[5,6,7]
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
        return facets.filter!(f => f.hasFace(simplex));
    }

    bool contains(Simplex simplex)
    {
        import std.algorithm : any;
        return facets.any!(f => f.hasFace(simplex));
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    Simplex[][size_t] facets_;
    
}

unittest
{
    SimplicialComplex sc;

    // insert some facets
    sc.insertFacet([1,5,8]);
    sc.insertFacet([1,5,7]);
    sc.insertFacet([7,9]);

    // get list of facets back, in order of increasing dimension
    // and dictionary order within a dimension 
    assert(sc.facets == [[7, 9], [1,5,7], [1,5,8]]);

    // get the number of facets
    assert(sc.numFacets == 3);

    // check for the presence of simplices
    assert(sc.contains([7,9]));
    assert(!sc.contains([1,5,6]));
    assert(sc.contains([1,5]));
    assert(sc.contains([7]));
    assert(!sc.contains([8, 9]));

    // get the star of a simplex as list of facets
    assert(sc.star([9]).array == [[7,9]]);
    assert(sc.star([7]).array == [[7, 9], [1,5,7]]);
    assert(sc.star([1,5]).array == [[1,5,7], [1,5,8]]);
    
    // vertices in an inserted facet must be sorted
    assertThrown!Error(sc.insertFacet([1,5,2,3]));

    // vertices may not be repeated
    assertThrown!Error(sc.insertFacet([1,3,3]));

    // cannot insert a simplex that is already a face of an existing facet
    assertThrown!Error(sc.insertFacet([7,9]));
    assertThrown!Error(sc.insertFacet([7]));
}