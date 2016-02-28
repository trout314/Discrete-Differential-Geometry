void main()
{
    import std.stdio : writeln;
    writeln("Entering main...");
    
    import simplex : simplex; 
    simplex(1,2,4).writeln;
    
    import std.range : iota;
    iota(5).writeln;
    
    import rational_extension_vector : simplexVecs;
    foreach(sVec; simplexVecs!10())
    {
        writeln(sVec);
    }
}
