version(unittest) {}
else
{
    void main()
    {
        import sampler : sample;
        // import gperftools_d.profiler;
        // ProfilerStart();
        // Test comment
        sample;
        // ProfilerStop();
    }
}

