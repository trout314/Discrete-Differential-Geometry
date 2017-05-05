import simplex : Simplex, simplex;


/// A simplicial complex type ... yada yada
struct SimplicialComplex
{
    // Inserts a facet into the simplicial complex. T
    void insertFacet(int[] vertices) pure nothrow @safe @nogc
    in
    {
        import std.algorithm : isSorted, findAdjacent;
        assert(vertices.isSorted, "vertices must be sorted");
        assert(vertices.findAdjacent.length == 0, "repeated vertices not allowed");
    }
    body
    {
    }
    private:
    
}

pure nothrow @safe unittest
{
    SimplicialComplex sc;
    sc.insertFacet([1,3,5]);
    sc.insertFacet([1,3,3]);
}