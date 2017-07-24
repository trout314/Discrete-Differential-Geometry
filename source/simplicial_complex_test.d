import unit_threaded : Types;
import simplicial_complex : SimplicialComplex;

@Types!(SimplicialComplex)
void testSimpComp(alias SimpComp)()
{
    import unit_threaded : shouldEqual, shouldBeEmpty;
    import std.range : array, empty, walkLength;
    import utility : throwsWithMsg;
    import simplex : simplex, Simplex;

    ///
    alias s = simplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimpComp!()();

    static assert(is(sc.VertexType == int)); // Default vertex type is int

    assert(sc.facets.empty);

    import fluent.asserts;

    Simplex!(0, sc.VertexType)[] Empty;
    sc.facets!0.should.equal(Empty);
    // sc.facets!0.shouldBeEmpty;

    assert(sc.facets!0.empty);
    assert(sc.facets!1.empty);
    assert(sc.facets!2.empty);

    assert(sc.facets!0.empty);
    assert(sc.facets!1.empty);
    assert(sc.facets!2.empty);

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
    
    auto sc1 = SimpComp!()([[1,2], [2,3], [3,4,5]]);
    assert(sc == sc1);
}