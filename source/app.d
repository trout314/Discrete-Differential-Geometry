version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();
        // sampleOld;
        // ProfilerStop();

        import sampler;

        import manifold : Manifold, standardSphereFacets;
        auto m = Manifold!3(standardSphereFacets(3));
        auto s = Sampler!(int, 3)(m);
        s.sample;
    }
}