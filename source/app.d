version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();

        import sampler;

        import manifold : standardSphere;
        auto s = Sampler!(int, 3)(standardSphere!3);
        s.sample;

        // ProfilerStop();
    }
}
