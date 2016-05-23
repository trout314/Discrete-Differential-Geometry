void main()
{
    import std.stdio : writeln;

    writeln("Entering main...");

    import simplex : simplex;

    simplex(1, 2, 4).writeln;

    import rational_extension_vector : simplexPoints, simplexRoots, simplexCoefs;
    import std.algorithm : each;

    simplexPoints!10.each!(v => v.writeln);    
    simplexRoots!10.each!(r => r.writeln);
    simplexCoefs!10.each!(c => c.writeln);
   
}
