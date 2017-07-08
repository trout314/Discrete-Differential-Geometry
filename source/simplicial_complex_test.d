auto test(alias SimpComp)()
{
    import std.range : array, walkLength;
    import utility : throwsWithMsg;
    import simplex : simplex;

    alias s = simplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimpComp!()();
    static assert(is(sc.VertexType == int)); // Default vertex type is int
    assert(sc.facets == []);
    assert(sc.facets!0.array == []);
    assert(sc.facets!1.array == []);
    assert(sc.facets(0) == []);
    assert(sc.facets(1) == []);
    assert(sc.numFacets == 0);
    
    throwsWithMsg(sc.link(s(1,2)),
        "expected a simplex in the simplicial complex");

    throwsWithMsg(sc.facets(-1), "facets expected a non-negative dimension but "
        ~ "got dimension -1"); 

    // ---------------- ADD SOME FACETS ---------------------
    sc.insertFacet(s(1,2));
    sc.insertFacet([2,3]);
    sc.insertFacet([3, 4, 5]);
    assert(sc.facets.walkLength == 3);
    
    auto sc1 = SimpComp!()([[1,2], [2,3], [2,3,4]]);
//    assert(sc == sc1);    

    return true;
}