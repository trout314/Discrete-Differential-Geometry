version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();

        import std.stdio : File;
        import sampler;
        import manifold : standardSphere, saveTo, saveEdgeGraphTo, loadManifold;

        auto m = loadManifold!3("mfd_S3_128000_VGC.dat");
        m.saveEdgeGraphTo("edge_graph_S3_128000_VGC.dat");

        // auto s = Sampler!(int, 3)(standardSphere!3);
        // s.sample;

        // auto fileName = "edge_graph_test.dat";
        // s.manifold.saveEdgeGraphTo(fileName);
        // auto saveFile = File(fileName, "a");
   
        // ProfilerStop();
    }
}
