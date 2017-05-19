void main()
{
    import std.stdio : writeln;
    "Entering main...".writeln;
    scope(exit) {"Exiting main.".writeln;}

    import simplicial_complex : SimplicialComplex, Simplex;
    SimplicialComplex sc;

    import std.conv : to;
    auto s = Simplex([1,2,3]);
    sc.insertFacet(s);
    writeln(sc);

/*    import simplex : simplex;
    import std.algorithm : map;
    writeln(simplex(1,2));
    writeln([1,2].simplex);
    writeln(typeof([1,2,3,4].simplex).stringof);

    writeln(typeof([1,2].map!simplex[0]).stringof);
*/
}
