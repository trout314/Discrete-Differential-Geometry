version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
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
        import std.range : walkLength;
        assert(simplex.isSorted, "vertices in a simplex must be sorted");
        assert(simplex.findAdjacent.length == 0,
            "repeated vertices in a simplex not allowed");
        assert(!facets.any!(f => f.hasFace(simplex)),
            "inserted simplex must not be a face of an existing facet");
        assert(simplex.walkLength != 0,
            "inserted simplex must have at least one vertex");
    }
    body
    {
        import std.range : walkLength;
        immutable numVerts = simplex.walkLength; 
        facetLists[numVerts] ~= simplex;

        import std.algorithm : sort;
        sort(facetLists[numVerts]);
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

        auto sizes = sort(facetLists.keys);        
        return sizes.map!(s => facetLists[s]).joiner.array;
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
        import std.range : array;
        import std.algorithm : map;
        return star(simplex).map!(facet => facet.oppositeFace(simplex)).array;
    }

    // Returns the star of the given simplex as a list of facets
    // in same order as they appear in facets() 
    auto star(Simplex simplex)
    {
        import std.range : array;
        import std.algorithm : filter;
        return facets.filter!(f => f.hasFace(simplex)).array;
    }

    // Returns true if simplex is in this simplicial complex and false otherwise
    bool contains(Simplex simplex)
    {
        import std.algorithm : any;
        return facets.any!(f => f.hasFace(simplex));
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    Simplex[][size_t] facetLists;
}

bool hasFace(Simplex s1, Simplex s2)
{
    import std.algorithm : all, canFind;
    return s2.all!(vertex => s1.canFind(vertex));
}

auto oppositeFace(Simplex simplex, Simplex face)
{
    import std.algorithm : filter, canFind;
    import std.range : array;
    return simplex.filter!(vertex => !face.canFind(vertex)).array;
}

unittest
{
    SimplicialComplex sc;

    // insert some facets
    sc.insertFacet([1,2,3]);
    sc.insertFacet([2,3,4]);
    sc.insertFacet([4,5]);
    sc.insertFacet([5,6]);

    // get list of facets back, in order of increasing dimension
    // and dictionary order within a dimension 
    assert(sc.facets == [[4,5], [5, 6], [1,2,3], [2,3,4]]);

    // get the number of facets
    assert(sc.numFacets == 4);

    // check for the presence of simplices
    assert(sc.contains([4,5]));
    assert(sc.contains([2,3,4]));
    assert(sc.contains([6]));
    assert(!sc.contains([1,2,5]));
    assert(!sc.contains([4,6]));

    // get the star of a simplex as list of facets
    assert(sc.star([5]) == [[4,5], [5,6]]);
    assert(sc.star([4]) == [[4,5], [2,3,4]]);
    assert(sc.star([2,3]) == [[1,2,3], [2,3,4]]);

    // get link of a simplex as list of facets
    assert(sc.link([4]) == [[5], [2,3]]);

    // star and link of an empty simplex give the entire complex
    assert(sc.star([]) == sc.facets);
    assert(sc.link([]) == sc.facets);

    // the empty simplex is considered to be included in any 
    // non-empty simplicial complex
    assert(sc.contains([]));

    // A newly constructed simplicial complex is empty
    SimplicialComplex scInit;
    assert(scInit.numFacets == 0);
    assert(scInit.facets == []);

    // star and link of an empty simplex give the entire (empty)
    // complex, just as before.
    assert(scInit.star([]) == []);
    assert(scInit.link([]) == []);


    // Note that an empty simplicial copmlex has no simplices
    // so none of them have the empty simplex as a face
    assert(!scInit.contains([]));

    // TO DO: Think about whether or not some of the operations above
    // that use the empty simplex should be errors instead?

    // --------------------------------------------------------------
    // Restrictions
    // --------------------------------------------------------------

    // vertices in an inserted facet must be sorted
    assertThrown!Error(sc.insertFacet([1,5,2,3]));

    // vertices may not be repeated
    assertThrown!Error(sc.insertFacet([1,3,3]));

    // it is an error to insert a simplex that is
    // already the face of an existing facet
    assertThrown!Error(sc.insertFacet([4,5]));
    assertThrown!Error(sc.insertFacet([2]));

    // since the empty facet is contained in any non-empty complex
    // it is an error to insert it too
    assertThrown!Error(sc.insertFacet([]));

    // For consistency and sanity we forbid inserting an empty
    // simplex into an empty complex as well
    assertThrown!Error(SimplicialComplex.init.insertFacet([]));
}