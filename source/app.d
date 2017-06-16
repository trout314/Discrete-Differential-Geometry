void main()
{
    import std.stdio : writeln;
    "Entering main...".writeln;
    scope(exit) {"Exiting main.".writeln;}

    import simplicial_complex : SimplicialComplex;
    import simplex : simplex;
    // SimplicialComplex!() sc;

    // sc.insertFacet(simplex(1,2,3));
    // writeln(sc);
}
