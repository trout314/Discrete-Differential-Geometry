struct Simplex
{
    this(R)(R vertices_)
    {
        import std.traits : isArray;
        import std.range : isInputRange, ElementType;
        import std.algorithm : isSorted, findAdjacent;
        import std.range : walkLength, empty;

        static assert(isInputRange!R || isArray!R,
            "Simplex must be constructed from a range or array");

        static assert(is(ElementType!R : int), "Element type of range or "
            ~ "array used to constuct a Simplex must be implicitly "
            ~ "convertible to int");

        assert(vertices_.isSorted, "vertices must be sorted");
        assert(vertices_.findAdjacent.walkLength == 0,
            "repeated vertices not allowed");
        assert(!vertices_.empty, "empty simplices not allowed");

        slice = vertices_.dup;
    }

    int[] vertices()
    {
        return slice[];
    }

    int opCmp(Simplex s) const
    {
        return this.slice < s.slice;
    }

    string toString() const
    {
        import std.conv : to;
        return slice.to!string;
    }

    private:
    int[] slice;
}


/// A simplicial complex type ... yada yada
struct SimplicialComplex
{
    // Inserts a facet into the simplicial complex.
    void insertFacet(Simplex simplex)
    in
    {
        import std.algorithm : any;
        assert(!facets.any!(f => f.hasFace(simplex)),
            "inserted simplex must not be a face of an existing facet");
    }
    body
    {
        import std.range : walkLength;
        immutable numVerts = simplex.vertices.walkLength; 
        facetLists[numVerts] ~= simplex;

        import std.algorithm : sort;
        sort(facetLists[numVerts]);
    }

    // Returns the facets of the simplicial complex.
    // These are simplicies that are not the face of another simplex.
    // They are returned in increasing order of dimension and in
    // lexicographic order within dimensions.
    // e.g. [8],[2,3],[2,4],[1,6,7]
    Simplex[] facets()
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
    Simplex[] link(Simplex simplex)
    {
        import std.range : array;
        import std.algorithm : map;
        return this.star(simplex)
            .map!(facet => facet.oppositeFace(simplex)).array;
    }

    // Returns the star of the given simplex as a list of facets
    // in same order as they appear in facets() 
    Simplex[] star(Simplex simplex)
    {
        import std.range : array;
        import std.algorithm : filter;
        return this.facets.filter!(f => f.hasFace(simplex)).array;
    }

    // Returns true if simplex is in this simplicial complex and false otherwise
    bool contains(Simplex simplex)
    {
        import std.algorithm : any;
        return this.facets.any!(f => f.hasFace(simplex));
    }

    string toString()
    {
        import std.conv : to;
        return "SComp" ~ this.facets.to!string;
    }

    private:
    // Lists of facets, indexed by number of vertices in the facet
    Simplex[][size_t] facetLists;

    // TO DO: Make the above into an array of (size_t, Simplex[])
    // pairs. Probably won't ever have enough distince facet
    // dimenstion to justify using an associative array
}

bool hasFace(Simplex s1, Simplex s2)
{
    import std.algorithm : all, canFind;
    return s2.vertices.all!(vertex => s1.vertices.canFind(vertex));
}

Simplex oppositeFace(Simplex simplex, Simplex face)
{
    import std.algorithm : filter, canFind;
    import std.range : array;
    import std.conv : to;

    return simplex.vertices.filter!(v => !face.vertices.canFind(v)).array.to!Simplex;
}

version(unittest)
{
    import std.exception : assertThrown;
    import std.stdio : writeln;
    import std.conv : to;
    import std.algorithm : map, each;
    import std.range : array;
}

unittest
{

    // To help with readability, we define a short aliases for Simplex and
    // Simplex[].
    // NOTE: Unfortunately we can't locally alias to!Simplex and
    // to!(Simplex[]) and also use UFCS.
    alias S = Simplex;
    alias SA = Simplex[];

    SimplicialComplex sc;
    assert(sc.numFacets == 0);
    assert(sc.facets == []);

    // insert some facets
    sc.insertFacet([1,2,3].to!S);
    sc.insertFacet([2,3,4].to!S);
    sc.insertFacet([4,5].to!S);
    sc.insertFacet([5,6].to!S);

    // get list of facets back, in order of increasing dimension
    // and dictionary order within a dimension 
    assert(sc.facets == [[4,5], [5,6], [1,2,3], [2,3,4]].to!SA);
    assert(sc.numFacets == 4);

    // check for the presence of simplices
    assert(sc.contains([4,5].to!S));
    assert(sc.contains([2,3,4].to!S));
    assert(sc.contains([6].to!S));
    assert(!sc.contains([1,2,5].to!S));
    assert(!sc.contains([4,6].to!S));

    // get the star of a simplex as list of facets
    assert(sc.star([5].to!S) == [[4,5], [5,6]].to!SA);
    assert(sc.star([4].to!S) == [[4,5], [2,3,4]].to!SA);
    assert(sc.star([2,3].to!S) == [[1,2,3], [2,3,4]].to!SA);

    // get link of a simplex as list of facets
    assert(sc.link([5].to!S) == [[4], [6]].to!SA);
    assert(sc.link([4].to!S) == [[5], [2,3]].to!SA);
    assert(sc.link([2,3].to!S) == [[1], [4]].to!SA);

    // --------------------------------------------------------------
    // Restrictions
    // --------------------------------------------------------------

    // vertices in an inserted facet must be sorted
    assertThrown!Error(sc.insertFacet([1,5,2,3].to!S));

    // vertices may not be repeated
    assertThrown!Error(sc.insertFacet([1,3,3].to!S));

    // it is an error to insert a simplex that is
    // already the face of an existing facet
    assertThrown!Error(sc.insertFacet([4,5].to!S));
    assertThrown!Error(sc.insertFacet([2].to!S));

    // since the empty facet is contained in any non-empty complex
    // it is an error to insert it too
//    auto s = Simplex([]);
//    assertThrown!Error(sc.insertFacet([].to!S));

    // For consistency and sanity we forbid inserting an empty
    // simplex into an empty complex as well
//    assertThrown!Error(SimplicialComplex.init.insertFacet([]));
}