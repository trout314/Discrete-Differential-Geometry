import unit_threaded : Types;
import simplicial_complex : SimplicialComplex, simplicialComplex;

@Types!(SimplicialComplex)
void testSimpComp(alias SimpComp)()
{
    import std.algorithm : equal;
    import std.range : array, empty, walkLength;
    import utility : throwsWithMsg, staticIota;
    import simplex : simplex, Simplex;
    import fluent.asserts;

    alias s = simplex;
    alias sComp = simplicialComplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimpComp!()();
    static assert(is(sc.VertexType == int)); // Default vertex type is int

    sc.facets.empty.should.equal(true);
    sc.numFacets.should.equal(0);
    foreach(d; staticIota!(0, 16))
    {
        sc.facets!d.empty.should.equal(true);
        sc.simplices!d.empty.should.equal(true);
    }
   
    throwsWithMsg(sc.link(s(1,2)),
        "expected a simplex in the simplicial complex");

    throwsWithMsg(sc.facets(-1), "facets expected a non-negative dimension but "
        ~ "got dimension -1"); 

    // ---------------------------- ADD SOME FACETS ----------------------------
    sc.insertFacet(s(1,2));
    sc.insertFacet([2,3]);
    sc.insertFacet([3, 4, 5]);
    
    sc.should.equal(sComp([[1,2], [2,3], [3,4,5]]));

    // -------------------------- TEST FUNCTIONALITY ---------------------------
    sc.facets.walkLength.should.equal(3);
    sc.numFacets.should.equal(3);

    sc.connectedComponents.walkLength.should.equal(1);
    sc.isConnected.should.equal(true);

    sc.facets!0.empty.should.equal(true);
    sc.facets!1.should.equal([s(1,2), s(2,3)]);
    sc.facets!2.should.equal([s(3,4,5)]);

    sc.simplices!0.should.containOnly([s(1), s(2), s(3), s(4), s(5)]);
    sc.simplices!1.should.containOnly([s(1,2), s(2,3), s(3,4), s(3,5), s(4,5)]);
    sc.simplices!2.should.containOnly([s(3,4,5)]);

    foreach(d; staticIota!(3, 16))
    {
        sc.facets!d.empty.should.equal(true);
        sc.facets(d).empty.should.equal(true);
        sc.simplices!d.empty.should.equal(true);
        sc.simplices(d).empty.should.equal(true);
    }

}