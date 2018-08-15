version(unittest) {}
else
{
    void main()
    {
        // import gperftools_d.profiler;
        // ProfilerStart();

        import sampler;

        import manifold : Manifold, standardSphereFacets;
        auto m = Manifold!3(standardSphereFacets(3));
        auto s = Sampler!(int, 3)(m);
        s.sample;

        // ProfilerStop();
    }
}

/* BUG REPORT?
import std.range : array, ElementType, isForwardRange, walkLength;
import std.algorithm  : copy, filter, map;

// return a newly constructed immutable static array
immutable(ElementType!R)[arrayLength] toImmutStaticArray(size_t arrayLength, R)(R range)
if (isForwardRange!R)
{
    assert(range.walkLength <= arrayLength,
        "range has too many elements for static array");

    ElementType!R[arrayLength] r;
    copy(range, r[]);
    return r;
}

void main()
{
    int[int[]] aaHeapKeys;
    int[int[]] aaStackKeys;
    int[int[2]] aaStaticKeys;
    
    int[] setA = [2,3,5];
    int[] setB = [2,4,3,9,5,25];
    
    // range of ranges: [[1,2,3], [1,2], [1]]
    auto sets = setA.map!(j => setB.filter!(i => i % j == 0));

    foreach(s; sets)
    {
        // insert using a slice of imuutable(int) on the heap as the key
        auto heapSlice = s.array.idup;
        static assert(is(typeof(heapSlice) == immutable(int)[]));        
        aaHeapKeys[heapSlice] = 0;

        // insert using an identical slice of immutable(int) on the stack
        // as the key
        auto buffer = s.toImmutStaticArray!2;
        auto stackSlice = buffer[0 .. s.walkLength];
        static assert(is(typeof(stackSlice) == immutable(int)[]));
        assert(stackSlice == heapSlice);
        aaStackKeys[stackSlice] = 0;

        aaStaticKeys[buffer] = 0;
    }

    import std.stdio : writeln;
    aaHeapKeys.writeln;
    aaStackKeys.writeln;
    aaStaticKeys.writeln;

    // assert(aaHeapKeys  == [[1]:0, [1, 2]:0, [1, 2, 3]:0]); // OK
    // assert(aaStackKeys == [[1]:0, [1, 2]:0, [1, 2, 3]:0]); // fails!

    // The associative array seems corrupted, since the following:
    // import std.stdio : writeln;
    // aaStackKeys.writeln;
    // prints: [[1, 0]:0, [1, 0, 0]:0, [1]:0] but still this fails:
    // assert(aaStackKeys == [[1, 0]:0, [1, 0, 0]:0, [1]:0]); // fails!
}

*/