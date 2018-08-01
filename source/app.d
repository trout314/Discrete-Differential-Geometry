version(unittest) {}
else
{
    void main()
    {
        import sampler : sampleOld;
        // import gperftools_d.profiler;
        // ProfilerStart();
        sampleOld;
        // ProfilerStop();
    }
}