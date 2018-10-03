version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();

        import std.stdio : File, write, writeln, stdout;
        import sampler;
        import manifold : findProblems, standardSphere, saveTo, saveEdgeGraphTo, loadManifold;
        import std.range : array, empty, front, back, popFront, popBack;
        import std.algorithm : each;

        // auto m = loadManifold!3("mfd_S3_16000_V.dat");
        // m.saveEdgeGraphTo("edge_graph_S3_16k_V.dat");

        // --------------------------------------------

        // "checking input manifold for problems...".write;
        // stdout.flush;

        // auto problems = m.findProblems.array;
        // if(!problems.empty)
        // {
        //     problems.each!writeln;
        //     assert(0);
        // }
        // "done.".writeln;
        // stdout.flush;


        // --------------------------------------------

        auto s = Sampler!(int, 3)(standardSphere!3);
        s.sample;

        auto ending = "test.dat";
        
        auto fileName = "mfd_" ~ ending;
        s.manifold.saveTo(fileName);
        auto saveFile = File(fileName, "a");
        saveFile.writeln;
        s.report(saveFile);
        s.manifold.saveEdgeGraphTo("edge_graph_" ~ ending);


        // ProfilerStop();
    }
}
