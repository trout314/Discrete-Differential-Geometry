void main()
{
    import std.stdio : writeln;
    writeln("Entering main...");
    
    import simplex : test; 
    test();
    
    import std.range : iota;
    iota(5).writeln;
    
    import rational_extension_vector : simplexVecs;
    writeln(simplexVecs!10());
}
