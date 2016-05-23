void main()
{
    import std.stdio : writeln;

    writeln("Entering main...");

    import simplex : simplex;

    simplex(1, 2, 4).writeln;

    import rational_extension_vector : simplexVecs, simplexRoots;
    import std.algorithm : each;

    simplexVecs!3().each!(v => v.writeln);
    simplexRoots(10).writeln;
}
