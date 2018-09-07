version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();

        import std.stdio : File;
        import sampler;
        import manifold : standardSphere, saveTo;
        
        auto s = Sampler!(int, 3)(standardSphere!3);
        s.sample;

        auto fileName = "test.dat";
        s.manifold.saveTo(fileName);
        auto saveFile = File(fileName, "a");
        s.report(saveFile);
    
        // ProfilerStop();
    }
}
