version(unittest) {}
else
{
    void main()
    {
        // import sampler : sample;
        // import gperftools_d.profiler;
        // ProfilerStart();
        // Test comment
        // sample;
        // ProfilerStop();
        import std.stdio : writeln;

        import simplicial_complex;
        auto sc = simplicialComplex([[1,2], [2,3,4], [2,3,5]]);
        sc.toDetailedString.writeln;
        
    }
}

