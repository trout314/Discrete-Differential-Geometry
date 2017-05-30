void main()
{
    import std.stdio : writeln;
    "Entering main...".writeln;
    scope(exit) {"Exiting main.".writeln;}

    import simplicial_complex : SimplicialComplex;
    import simplex : simplex;
    SimplicialComplex!() sc;

    sc.insertFacet(simplex(1,2,3));
    writeln(sc);

/*    import simplex : simplex;
    import std.algorithm : map;
    writeln(simplex(1,2));
    writeln([1,2].simplex);
    writeln(typeof([1,2,3,4].simplex).stringof);

    writeln(typeof([1,2].map!simplex[0]).stringof);
*/
}
