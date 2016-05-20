void main()
{
    import std.stdio : writeln;
    writeln("Entering main...");
    
    import simplex : simplex; 
    simplex(1,2,4).writeln;
    
    import std.range : iota;
    iota(5).writeln;
    
    import rational_extension_vector : simplexVecs, simplexRoots;
    auto sVecs = simplexVecs!10(); 
    foreach(sVec; sVecs)
    {
        sVec.writeln;
    }
    
    import std.rational : rational;
    auto a = rational(2,4);
    writeln(a);
    writeln(typeof(sVecs).stringof);
    
    simplexRoots(19).writeln;
}
