void main()
{
    import std.stdio : writeln;

    writeln("Entering main...");

    import simplex : simplex;

    simplex(1, 2, 4).writeln;

    import rational_extension_vector : simplexVecs, simplexRoots, simplexCoefs;
    import std.algorithm : each;

    simplexVecs!10.each!(v => v.writeln);    
    simplexRoots!10.each!(r => r.writeln);
    simplexCoefs!10.each!(c => c.writeln);
   
}
