void main()
{
    import std.stdio : writeln;
    writeln("Entering main...");
    
    import simplex; 
    test();
    
    
    import std.range : iota;
    iota(5).writeln;
}
