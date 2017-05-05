void main()
{
    import std.stdio : writeln;
    writeln("Entering main...");

    import simplicial_complex : SimplicialComplex;
    SimplicialComplex sc;
    sc.insertFacet([1,2,3]);
}
