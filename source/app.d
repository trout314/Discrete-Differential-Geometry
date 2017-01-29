void main()
{
    import std.stdio : writeln;

    writeln("Entering main...");

    import simplex : simplex;

    // simplex(1, 2, 4).writeln;

    import rational_extension_vector : simplexPoints, simplexRoots, simplexCoefs;
    import factoring : primeFactorsRange;
    import std.algorithm : each;

    // simplexPoints!10.each!writeln;    
    // simplexRoots!10.each!writeln;
    // simplexCoefs!10.each!writeln;
    
    primeFactorsRange(360).writeln;
}
