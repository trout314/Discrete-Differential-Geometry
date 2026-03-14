/// Contains functionality for working with combinatorial n-manifolds
module manifold;

import std.algorithm, std.array, std.conv, std.exception, std.math, std.range,
    std.stdio, std.traits, std.typecons;
import algorithms, hashmap, manifold_examples, manifold_moves,
    simplicial_complex, utility;

alias isIRof = isInputRangeOf;
alias isIRofIRof = isInputRangeOfInputRangeOf;


//dfmt off
/*******************************************************************************
Manifold type... TO DO: More info here.

Manifold models a compact combinatorial n-manifold. Supports doing pachner moves
quickly
*/
struct Manifold(int dimension_, Vertex_ = int)
{
private:
    alias SimpComp = SimplicialComplex!Vertex_;
    // Dummy field: forces DMD to instantiate SimpComp postblit in this context,
    // which is needed for correct purity inference in algorithms.d tests.
    // Can be removed once SmallMap/HashMap get explicit pure @safe attributes.
    SimpComp _scDummy;

    // --- Per-dimension index maps ---
    // For each dimension d (0..dimension), a separate hash map keyed by
    // Vertex[d+1] (exact-size static array keys). For the ridge dimension
    // (d == dimension-1), the value bundles degree + link info. For all
    // other dimensions, the value is just the degree (size_t).

    alias RidgeLink = StackArray!(Vertex, 2);

    struct RidgeInfo
    {
        uint degree;
        RidgeLink link;
    }

    // Flat array for vertex degrees (dim 0): _vertexDegrees[v] = degree of vertex v.
    // Zero means vertex not present. Avoids hashing for the most frequent lookups.
    uint[] _vertexDegrees;

    // Per-dimension HashMaps for dimensions 1 .. dimension-1 only.
    // Dimension 0 uses _vertexDegrees; dimension `dimension` uses _facetArrayIdx.
    static foreach (_d_; 1 .. dimension_)
    {
        static if (_d_ == dimension_ - 1)
            mixin("HashMap!(Vertex_[" ~ to!string(_d_ + 1) ~ "], RidgeInfo) _dimMap" ~ to!string(_d_) ~ ";");
        else
            mixin("HashMap!(Vertex_[" ~ to!string(_d_ + 1) ~ "], uint) _dimMap" ~ to!string(_d_) ~ ";");
    }

    // Compile-time accessor for per-dimension maps (dimensions 1..dimension-1 only)
    ref auto dimMap(int d)() if (d >= 1 && d <= dimension_ - 1)
    {
        return mixin("_dimMap" ~ to!string(d));
    }

    ref auto dimMap(int d)() const if (d >= 1 && d <= dimension_ - 1)
    {
        return mixin("_dimMap" ~ to!string(d));
    }

    // Extract degree from a per-dimension map value (size_t or RidgeInfo)
    private static uint extractDegree(V)(V value)
    {
        static if (is(V == RidgeInfo) || is(V == const(RidgeInfo)))
            return value.degree;
        else
            return value;
    }

    // Legacy aliases kept for compatibility
    alias NSimplex = StackArray!(Vertex, dimension + 1);
    alias Ridge = Vertex[dimension_];
    alias Facet = Vertex[dimension_ + 1];

    uint[dimension_ + 1] numSimplices;
    ulong[dimension - 1] totSqrDegrees;

    // Facet array for O(1) random access
    Facet[] _facetArray;
    HashMap!(Facet, uint) _facetArrayIdx;

    version (TrackValidMoves)
    {
        Facet[][Vertex] _vertexFacets;

        // Per-dimension valid move arrays (swap-with-last) and reverse index maps
        static foreach (_d_vm; 0 .. dimension_)
        {
            mixin("Vertex_[" ~ to!string(_d_vm + 1) ~ "][] _validMoves" ~ to!string(_d_vm) ~ ";");
            mixin("HashMap!(Vertex_[" ~ to!string(_d_vm + 1) ~ "], size_t) _validMoveIdx" ~ to!string(_d_vm) ~ ";");
        }
    }

    //--------------- Helper Functions ---------------

    /// Find the vertex in `facet` that is not in `subset` (for ridge computations).
    private static Vertex_ oppositeVertex(const Vertex_[] facet, const Vertex_[dimension_] subset) pure nothrow @nogc @safe
    {
        foreach (v; facet)
        {
            bool found = false;
            foreach (s; subset)
                if (v == s) { found = true; break; }
            if (!found) return v;
        }
        assert(false, "no opposite vertex found");
    }

    static NSimplex toNSimp(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.toStackArray!(Vertex_, dimension + 1, R);
    }

    static Ridge toRidge(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.staticArray!dimension;
    }

    static Facet toFacet(R)(R range) if (isIRof!(R, Vertex_))
    {
        return range.staticArray!(dimension_ + 1);
    }
public:
    /// Dimension of the manifold
    static immutable dimension = dimension_;

    /// Vertex type used in the manifold
    alias Vertex = Vertex_;

    /// Move type used for the manifold
    // alias Move = BistellarMove!(dimension_, Vertex_);

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    /// We can initialize the manifold from an input range of input ranges
    /// of vertices
    this(F)(F initialFacets) if (isIRofIRof!(F, const(Vertex)))
    {
        initialFacets.each!(f => this.insertFacet(f));
        foreach (d; 0 .. dimension - 1)
        {
            totSqrDegrees[d] = simplices(d).map!(s => this.degree(s)^^2).sum;
        }
        version (TrackValidMoves)
        {
            // _vertexFacets and _facetArray already populated by insertFacet.
            // Build per-dimension valid move arrays.
            // Dimension 0: iterate _vertexDegrees
            {
                enum _targetDeg = dimension + 1;
                foreach (v; 0 .. cast(int) _vertexDegrees.length)
                {
                    if (_vertexDegrees[v] == _targetDeg)
                    {
                        Vertex[1] key = [cast(Vertex) v];
                        auto cc = localCoCenter!(0, Vertex, dimension)(this, key);
                        if (!this.contains(cc[]))
                            this.addValidMoveToArray!0(key);
                    }
                }
            }
            // Dimensions 1..dimension-1: HashMap-backed
            static foreach (_d_init; 1 .. dimension)
            {{
                foreach (kv; dimMap!_d_init.byKeyValue())
                {
                    auto deg = extractDegree(kv.value);
                    if (deg == dimension + 1 - _d_init)
                    {
                        auto key = kv.key;
                        auto cc = localCoCenter!(_d_init, Vertex, dimension)(this, key);
                        if (!this.contains(cc[]))
                            this.addValidMoveToArray!_d_init(key);
                    }
                }
            }}
        }
    }

    this(this)
    {
        _vertexDegrees = _vertexDegrees.dup;
        static foreach (d; 1 .. dimension_)
            mixin("_dimMap" ~ to!string(d) ~ " = _dimMap" ~ to!string(d) ~ ".dup;");

        _facetArray = _facetArray.dup;
        _facetArrayIdx = _facetArrayIdx.dup;

        version (TrackValidMoves)
        {
            Facet[][Vertex] newVF;
            foreach (k, v; _vertexFacets)
                newVF[k] = v.dup;
            _vertexFacets = newVF;

            static foreach (_d_dup; 0 .. dimension_)
            {{
                mixin("_validMoves" ~ to!string(_d_dup) ~ " = _validMoves" ~ to!string(_d_dup) ~ ".dup;");
                mixin("_validMoveIdx" ~ to!string(_d_dup) ~ " = _validMoveIdx" ~ to!string(_d_dup) ~ ".dup;");
            }}
        }
    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        if (simplex.empty) return false;
        return containsInDimMaps(simplex);
    }

    /// Runtime-dispatch lookup into per-dimension maps.
    private bool containsInDimMaps(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        Vertex[dimension + 1] buf;
        size_t len = 0;
        foreach (v; simplex) buf[len++] = v;

        final switch (cast(int)(len - 1))
        {
            case 0:
            {
                auto v = buf[0];
                return v < _vertexDegrees.length && _vertexDegrees[v] > 0;
            }
            static foreach (d; 1 .. dimension)
            {
                case d:
                {
                    Vertex[d + 1] key = buf[0 .. d + 1];
                    return (key in dimMap!d) !is null;
                }
            }
            case dimension:
            {
                Facet key = buf[0 .. dimension + 1];
                return (key in _facetArrayIdx) !is null;
            }
        }
    }

    /***************************************************************************
    Returns the degree of a simplex in the manifold.
    */
    size_t degree(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        assert(this.contains(simplex),
            "called degree on a simplex not in the manifold");
        return degreeFromDimMaps(simplex);
    }

    /// Runtime-dispatch degree lookup.
    private size_t degreeFromDimMaps(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        Vertex[dimension + 1] buf;
        size_t len = 0;
        foreach (v; simplex) buf[len++] = v;

        final switch (cast(int)(len - 1))
        {
            case 0:
            {
                auto v = buf[0];
                assert(v < _vertexDegrees.length && _vertexDegrees[v] > 0, "Vertex not in manifold");
                return _vertexDegrees[v];
            }
            static foreach (d; 1 .. dimension)
            {
                case d:
                {
                    Vertex[d + 1] key = buf[0 .. d + 1];
                    auto ptr = key in dimMap!d;
                    assert(ptr !is null, "Key not found in dimMap");
                    return extractDegree(*ptr);
                }
            }
            case dimension:
            {
                Facet key = buf[0 .. dimension + 1];
                assert((key in _facetArrayIdx) !is null, "Facet not in manifold");
                return 1;  // facets always have degree 1
            }
        }
    }

    /***************************************************************************
    Returns the degree of a simplex, or 0 if the simplex is not in the manifold.
    Avoids the double lookup of contains() + degree().
    */
    size_t degreeOrZero(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        Vertex[dimension + 1] buf;
        size_t len = 0;
        foreach (v; simplex) buf[len++] = v;

        final switch (cast(int)(len - 1))
        {
            case 0:
            {
                auto v = buf[0];
                if (v >= _vertexDegrees.length) return 0;
                return _vertexDegrees[v];
            }
            static foreach (d; 1 .. dimension)
            {
                case d:
                {
                    Vertex[d + 1] key = buf[0 .. d + 1];
                    auto ptr = key in dimMap!d;
                    if (ptr is null) return 0;
                    return extractDegree(*ptr);
                }
            }
            case dimension:
            {
                Facet key = buf[0 .. dimension + 1];
                return (key in _facetArrayIdx) !is null ? 1 : 0;
            }
        }
    }

    /// Fast compile-time dimension-aware degreeOrZero (avoids runtime dispatch).
    size_t degreeOrZero(int d, S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        static if (d == 0)
        {
            auto key = simplex.staticArray!1;
            auto v = key[0];
            if (v >= _vertexDegrees.length) return 0;
            return _vertexDegrees[v];
        }
        else static if (d == dimension)
        {
            auto key = simplex.staticArray!(dimension + 1);
            return (key in _facetArrayIdx) !is null ? 1 : 0;
        }
        else
        {
            auto key = simplex.staticArray!(d + 1);
            auto ptr = key in dimMap!d;
            if (ptr is null) return 0;
            return extractDegree(*ptr);
        }
    }

    /***************************************************************************
    Returns an array containing the number of simplices in each dimension.
    This is called the "fVector" of the manifold.
    */
    const(uint)[] fVector()() const
    {
        return numSimplices[];
    }

    version (TrackValidMoves)
    {
        /// Valid move count computed from array sizes.
        size_t validMoveCount()() const
        {
            size_t count = _facetArray.length;
            static foreach (_d_c; 0 .. dimension)
                count += mixin("_validMoves" ~ to!string(_d_c)).length;
            return count;
        }

        /// Accessor for valid move array in dimension d.
        ref auto validMoves(int d)()
        {
            return mixin("_validMoves" ~ to!string(d));
        }

        ref auto validMoves(int d)() const
        {
            return mixin("_validMoves" ~ to!string(d));
        }

        /// Accessor for valid move reverse index in dimension d.
        ref auto validMoveIdx(int d)()
        {
            return mixin("_validMoveIdx" ~ to!string(d));
        }

        ref auto validMoveIdx(int d)() const
        {
            return mixin("_validMoveIdx" ~ to!string(d));
        }

        /// Append a valid center to the swap-with-last array.
        private void addValidMoveToArray(int d)(Vertex[d + 1] key)
        {
            assert(key !in validMoveIdx!d, "adding duplicate valid move");
            validMoveIdx!d[key] = validMoves!d.length;
            validMoves!d ~= key;
        }

        /// Remove a valid center via swap-with-last.
        private void removeValidMoveFromArray(int d)(Vertex[d + 1] key)
        {
            auto idxPtr = key in validMoveIdx!d;
            assert(idxPtr !is null, "removing non-existent valid move");
            auto idx = *idxPtr;
            auto lastIdx = validMoves!d.length - 1;
            if (idx != lastIdx)
            {
                auto last = validMoves!d[lastIdx];
                validMoves!d[idx] = last;
                validMoveIdx!d[last] = idx;
            }
            validMoves!d = validMoves!d[0 .. lastIdx];
            validMoveIdx!d.remove(key);
        }
    }

    // Special version of insertFacet to update tracked info.
    // Self-contained: also updates numSimplices (f-vector).
    // Uses per-dimension maps with exact-size static array keys.
    private void insertFacet(F)(F facet_) if (isIRof!(F, Vertex))
    {
        assert(facet_.walkLength == dimension + 1,
            "facet has wrong dimension");

        auto facetBuffer = toFacet(facet_);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);
        assert(facetBuffer !in _facetArrayIdx,
            "tried to insert a facet already in the manifold");

        // Dimension 0: vertex degrees via flat array
        foreach (v; facet)
        {
            if (v >= _vertexDegrees.length)
                _vertexDegrees.length = v + 1;

            if (_vertexDegrees[v] == 0)
            {
                numSimplices[0]++;
                _vertexDegrees[v] = 1;
                static if (dimension > 1)
                    totSqrDegrees[0] += 1; // 2*1 - 1
            }
            else
            {
                ++_vertexDegrees[v];
                static if (dimension > 1)
                    totSqrDegrees[0] += 2 * _vertexDegrees[v] - 1;
            }
        }

        // Dimensions 1 .. dimension-1: HashMap-backed
        static foreach (d; 1 .. dimension)
        {{
            foreach (subset; facet[].subsetsOfSize(d + 1))
            {
                auto key = subset.staticArray!(d + 1);

                static if (d == dimension - 1)
                {
                    // Ridge dimension: bundled degree + link
                    auto ptr = key in dimMap!d;
                    if (ptr is null)
                    {
                        numSimplices[d]++;
                        auto oppVert = oppositeVertex(facet, key);
                        RidgeInfo ri;
                        ri.degree = 1;
                        ri.link ~= oppVert;
                        dimMap!d[key] = ri;
                    }
                    else
                    {
                        ptr.degree++;
                        assert(ptr.degree == 2);
                        auto oppVert = oppositeVertex(facet, key);
                        ptr.link ~= oppVert;
                        ptr.link[].sort;
                    }
                }
                else
                {
                    // Sub-ridge dimension: degree + totSqrDegrees
                    auto ptr = key in dimMap!d;
                    if (ptr is null)
                    {
                        dimMap!d[key] = 1;
                        numSimplices[d]++;
                        totSqrDegrees[d] += 1; // 2*1 - 1
                    }
                    else
                    {
                        ++(*ptr);
                        totSqrDegrees[d] += 2 * (*ptr) - 1;
                    }
                }
            }
        }}

        // Facet dimension: tracked by _facetArrayIdx, just count
        numSimplices[dimension]++;

        _facetArrayIdx[facetBuffer] = cast(uint) _facetArray.length;
        _facetArray ~= facetBuffer;

        version (TrackValidMoves)
        {
            foreach (v; facet)
                _vertexFacets[v] ~= facetBuffer;
        }
    }

    // Special version of removeFacet to update tracked info.
    // Self-contained: also updates numSimplices (f-vector).
    private void removeFacet(F)(F facet_) if (isIRof!(F, Vertex))
    {
        assert(facet_.walkLength == dimension + 1);

        auto facetBuffer = toFacet(facet_);
        auto facet = facetBuffer[];

        facet.assertValidSimplex(dimension);
        assert(facetBuffer in _facetArrayIdx);

        // Dimension 0: vertex degrees via flat array
        foreach (v; facet)
        {
            assert(v < _vertexDegrees.length && _vertexDegrees[v] > 0);
            --_vertexDegrees[v];
            static if (dimension > 1)
                totSqrDegrees[0] -= 2 * _vertexDegrees[v] + 1;

            if (_vertexDegrees[v] == 0)
                numSimplices[0]--;
        }

        // Dimensions 1 .. dimension-1: HashMap-backed
        static foreach (d; 1 .. dimension)
        {{
            foreach (subset; facet[].subsetsOfSize(d + 1))
            {
                auto key = subset.staticArray!(d + 1);

                static if (d == dimension - 1)
                {
                    // Ridge dimension
                    auto ptr = key in dimMap!d;
                    assert(ptr !is null);
                    ptr.degree--;

                    if (ptr.degree == 1)
                    {
                        auto oppVert = oppositeVertex(facet, key);
                        assert(ptr.link[].canFind(oppVert));
                        if (ptr.link.length == 2)
                        {
                            if (ptr.link[0] == oppVert)
                                ptr.link[0] = ptr.link[1];
                            else
                                assert(ptr.link[1] == oppVert);
                            ptr.link.length = 1;
                        }
                        else
                        {
                            assert(ptr.link.length == 1);
                        }
                    }

                    if (ptr.degree == 0)
                    {
                        numSimplices[d]--;
                        dimMap!d.remove(key);
                    }
                }
                else
                {
                    // Sub-ridge dimension
                    auto ptr = key in dimMap!d;
                    assert(ptr !is null);
                    (*ptr)--;
                    totSqrDegrees[d] -= 2 * (*ptr) + 1;

                    if (*ptr == 0)
                    {
                        numSimplices[d]--;
                        dimMap!d.remove(key);
                    }
                }
            }
        }}

        // Facet dimension: tracked by _facetArrayIdx, just count
        numSimplices[dimension]--;

        // Remove facet from _facetArray (swap with last)
        {
            auto fIdxPtr = facetBuffer in _facetArrayIdx;
            assert(fIdxPtr !is null);
            auto fIdx = *fIdxPtr;
            auto fLastIdx = _facetArray.length - 1;
            if (fIdx != fLastIdx)
            {
                auto fLast = _facetArray[fLastIdx];
                _facetArray[fIdx] = fLast;
                _facetArrayIdx[fLast] = fIdx;
            }
            _facetArray = _facetArray[0 .. fLastIdx];
            _facetArrayIdx.remove(facetBuffer);
        }

        version (TrackValidMoves)
        {
            // Remove facet from per-vertex adjacency
            foreach (v; facet)
            {
                auto arr = _vertexFacets[v];
                foreach (i; 0 .. arr.length)
                {
                    if (arr[i] == facetBuffer)
                    {
                        arr[i] = arr[$ - 1];
                        _vertexFacets[v] = arr[0 .. $ - 1];
                        break;
                    }
                }
                if (_vertexFacets[v].length == 0)
                    _vertexFacets.remove(v);
            }
        }

        assert(facetBuffer !in _facetArrayIdx);
    }

    /// Build a SimplicialComplex from the manifold's facets (on demand).
    SimpComp toSimplicialComplex()() const
    {
        return SimpComp(facets);
    }

    /// Return all simplices of the given dimension as arrays of vertices.
    const(Vertex)[][] simplices()(int dim) const
    {
        final switch (cast(int) dim)
        {
            case 0:
            {
                const(Vertex)[][] result;
                foreach (v; 0 .. cast(int) _vertexDegrees.length)
                    if (_vertexDegrees[v] > 0)
                        result ~= [cast(const(Vertex)) v];
                return result;
            }
            static foreach (d; 1 .. dimension)
            {
                case d:
                {
                    const(Vertex)[][] result;
                    foreach (kv; dimMap!d.byKeyValue())
                        result ~= cast(const(Vertex)[]) kv.key[].dup;
                    return result;
                }
            }
            case dimension:
            {
                const(Vertex)[][] result;
                foreach (ref f; _facetArray)
                    result ~= cast(const(Vertex)[]) f[].dup;
                return result;
            }
        }
    }

    /// Return all facets as an array of vertex arrays, sorted for determinism.
    const(Vertex)[][] facets()() const
    {
        const(Vertex)[][] result;
        foreach (ref f; _facetArray)
            result ~= cast(const(Vertex)[]) f[].dup;
        result.sort();
        return result;
    }

    /// Return the number of facets.
    size_t numFacets()() const
    {
        return numSimplices[dimension];
    }

    /// Return true if the given simplex is a facet.
    bool containsFacet(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        auto key = toFacet(simplex);
        return (key in _facetArrayIdx) !is null;
    }

    /// Return a randomly chosen facet of the given dimension.
    const(Vertex)[] randomFacetOfDim()(int dim) const
    {
        assert(dim == dimension, "manifolds only have facets of one dimension");
        import std.random : uniform;
        auto idx = uniform(0, _facetArray.length);
        return cast(const(Vertex)[]) _facetArray[idx][].dup;
    }

    /// Return the star of a simplex (facets containing it).
    auto star(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        return facets.filter!(f => simplex.isSubsetOf(f));
    }

    /// Return the link of a simplex.
    auto link(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        return star(simplex).map!(f => f.setDifference(simplex));
    }

    /// Save the manifold to a file.
    void saveTo()(string fileName) const
    {
        toSimplicialComplex.saveTo(fileName);
    }

    /// Save the edge graph to a file.
    void saveEdgeGraphTo()(string fileName) const
    {
        toSimplicialComplex.saveEdgeGraphTo(fileName);
    }

    /// Euler characteristic computed from the f-vector.
    int eulerCharacteristic()() const
    {
        int chi = 0;
        foreach (d; 0 .. dimension + 1)
            chi += (d % 2 == 0 ? 1 : -1) * cast(int) numSimplices[d];
        return chi;
    }

    /// Equality comparison.
    bool opEquals()(const ref Manifold rhs) const
    {
        return numSimplices == rhs.numSimplices
            && facets.sort.array == rhs.facets.sort.array;
    }

    /// String representation (list of facets).
    string toString()() const
    {
        return facets.to!string;
    }
}



///
pure @safe unittest
{
    static foreach (d; 2 .. 8)
    {
        {
            auto m = standardSphere!d;
            m.allBistellarMoves.shouldBeEmpty;
        }
    }

    // trigonal bipyramid
    auto tb = Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,4]]);
    tb.allBistellarMoves.shouldBeSameSetAs([
        BistellarMove!2([0],[1,2,3]),
        BistellarMove!2([4],[1,2,3]),
        BistellarMove!2([1,2],[0,4]),
        BistellarMove!2([1,3],[0,4]),
        BistellarMove!2([2,3],[0,4])
    ]);

    // octahedron
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    oct.allBistellarMoves.shouldBeSameSetAs([
        BistellarMove!2([0,1],[2,4]),
        BistellarMove!2([0,2],[1,3]),
        BistellarMove!2([0,3],[2,4]),
        BistellarMove!2([0,4],[1,3]),
        BistellarMove!2([1,5],[2,4]),
        BistellarMove!2([2,5],[1,3]),
        BistellarMove!2([3,5],[2,4]),
        BistellarMove!2([4,5],[1,3]),
        BistellarMove!2([1,2],[0,5]),
        BistellarMove!2([2,3],[0,5]),
        BistellarMove!2([3,4],[0,5]),
        BistellarMove!2([1,4],[0,5])
    ]);

    auto m = Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,5],
        [2,4,5],[3,4,5]]);    
    m.allBistellarMoves.shouldBeSameSetAs([
        BistellarMove!2([0],[1,2,3]),
        BistellarMove!2([5],[2,3,4]),
        BistellarMove!2([1,2],[0,4]),
        BistellarMove!2([1,3],[0,4]),
        BistellarMove!2([2,3],[0,5]),
        BistellarMove!2([2,4],[1,5]),
        BistellarMove!2([3,4],[1,5])
    ]);

    // two-point suspension over boundary of 3-simplex
    auto tps = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],
        [1,2,3,5],[1,2,4,5],[1,3,4,5],[2,3,4,5]]);
}


/*******************************************************************************
Returns a list of all the bistellar moves in this manifold (except for the
1->(dim+1) moves.
*/
BistellarMove!(dim, Vertex)[] allBistellarMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    BistellarMove!(dim, Vertex)[] result;

    // Dimension 0: vertices from flat array
    {
        enum targetDeg = dim + 1;
        foreach (v; 0 .. cast(int) mfd._vertexDegrees.length)
        {
            if (mfd._vertexDegrees[v] == targetDeg)
            {
                auto simp = [cast(const(Vertex)) v];
                auto coCenter_ = mfd.findCoCenter(simp);
                if (!mfd.contains(coCenter_))
                    result ~= BistellarMove!(dim, Vertex)(simp, coCenter_);
            }
        }
    }

    // Dimensions 1..dim-1: HashMap-backed
    static foreach (d; 1 .. dim)
    {{
        foreach (kv; mfd.dimMap!d.byKeyValue())
        {
            auto deg = mfd.extractDegree(kv.value);
            if (deg == dim + 1 - d)
            {
                auto simp = kv.key[];
                auto coCenter_ = mfd.findCoCenter(simp);
                if (!mfd.contains(coCenter_))
                {
                    result ~= BistellarMove!(dim, Vertex)(simp, coCenter_);
                }
            }
        }
    }}
    return result;
}

/*******************************************************************************
Returns the number of valid bistellar moves (including 1->(dim+1) stellar
subdivisions). Used for exact Hastings correction.
*/
size_t countValidBistellarMoves(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd)
{
    size_t count = 0;

    // Dimension 0: vertices from flat array
    {
        enum targetDeg = dim + 1;
        foreach (v; 0 .. cast(int) mfd._vertexDegrees.length)
        {
            if (mfd._vertexDegrees[v] == targetDeg)
            {
                auto simp = [cast(const(Vertex)) v];
                auto coCenter_ = mfd.findCoCenter(simp);
                if (!mfd.contains(coCenter_))
                    count++;
            }
        }
    }

    // Dimensions 1..dim-1: HashMap-backed
    static foreach (d; 1 .. dim)
    {{
        foreach (kv; mfd.dimMap!d.byKeyValue())
        {
            auto deg = mfd.extractDegree(kv.value);
            if (deg == dim + 1 - d)
            {
                auto simp = kv.key[];
                auto coCenter_ = mfd.findCoCenter(simp);
                if (!mfd.contains(coCenter_))
                    count++;
            }
        }
    }}
    // Every facet is a valid center for a 1->(dim+1) stellar subdivision
    count += mfd.fVector[dim];
    return count;
}

///
pure unittest
{
    // Standard sphere has no moves except stellar subdivisions
    auto sphere = standardSphere!2;
    // allBistellarMoves is empty for standard sphere, so count = nFacets
    sphere.countValidBistellarMoves.shouldEqual(sphere.fVector[2]);

    // Trigonal bipyramid
    auto tb = Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,4]]);
    // 5 moves from allBistellarMoves + 6 facets
    tb.countValidBistellarMoves.shouldEqual(tb.allBistellarMoves.length + tb.fVector[2]);

    // Octahedron
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    oct.countValidBistellarMoves.shouldEqual(oct.allBistellarMoves.length + oct.fVector[2]);
}

version (TrackValidMoves)
{

///
unittest
{
    import std.random : uniform, Mt19937;
    alias BM = BistellarMove!2;

    // Start with octahedron
    auto mfd = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    // Verify initial count
    mfd.validMoveCount.shouldEqual(mfd.countValidBistellarMoves);

    // Do a sequence of moves and verify after each
    auto rng = Mt19937(42);
    int nextV = 6;

    foreach (_; 0 .. 30)
    {
        auto moves = mfd.allBistellarMoves;
        if (moves.empty) break;

        // Pick a random valid move
        auto idx = uniform(0, cast(int) moves.length, rng);
        auto move = moves[idx];
        mfd.doMove(move);

        // Add stellar subdivisions too
        if (uniform(0, 3, rng) == 0)
        {
            auto facet = mfd.facets.front.array;
            mfd.doMove(BM(facet, [nextV]));
            nextV++;
        }

        assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
            "incremental count %d != full recount %d"
            .format(mfd.validMoveCount, mfd.countValidBistellarMoves));
    }
}

///
unittest
{
    import std.random : uniform, Mt19937;
    alias BM = BistellarMove!3;

    // Start with a 3-sphere
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);

    mfd.validMoveCount.shouldEqual(mfd.countValidBistellarMoves);

    auto rng = Mt19937(123);
    int nextV = 5;

    foreach (_; 0 .. 30)
    {
        // Do a stellar subdivision to grow
        auto facet = mfd.facets.front.array;
        mfd.doMove(BM(facet, [nextV]));
        nextV++;

        assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
            "incremental count %d != full recount %d after stellar subdiv"
            .format(mfd.validMoveCount, mfd.countValidBistellarMoves));

        // Try a random valid move if available
        auto moves = mfd.allBistellarMoves;
        if (!moves.empty)
        {
            auto idx = uniform(0, cast(int) moves.length, rng);
            mfd.doMove(moves[idx]);

            assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
                "incremental count %d != full recount %d after bistellar move"
                .format(mfd.validMoveCount, mfd.countValidBistellarMoves));
        }
    }
}

} // version (TrackValidMoves) — tests

version (TrackValidMoves)
{

/*******************************************************************************
Check whether a d-simplex (given as a static array key) is a valid bistellar
center. Requires: degree == dim+1-d AND coCenter not in manifold.
Does NOT count facet-dimension centers (stellar subdivisions).
*/
private bool isValidNonFacetCenter(int d, Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd,
    const ref Vertex[d + 1] key)
{
    static assert(d < dim, "use fVector[dim] for facet centers");

    uint deg;
    static if (d == 0)
    {
        auto v = key[0];
        if (v >= mfd._vertexDegrees.length) return false;
        deg = mfd._vertexDegrees[v];
    }
    else
    {
        auto ptr = key in mfd.dimMap!d;
        if (ptr is null) return false;
        deg = mfd.extractDegree(*ptr);
    }

    if (deg != dim + 1 - d) return false;
    auto coCenter_ = mfd.localCoCenter!d(key);
    return !mfd.contains(coCenter_[]);
}

/*******************************************************************************
Compute the coCenter of a d-simplex using the per-vertex adjacency list
(_vertexFacets) instead of the O(n) simplicial complex link computation.
Only valid for simplices with degree == dim+1-d.
Returns a StackArray of coCenter vertices, sorted.
*/
private auto localCoCenter(int d, Vertex, int dim)(
    const ref Manifold!(dim, Vertex) mfd,
    const ref Vertex[d + 1] center)
{
    StackArray!(Vertex, dim + 1 - d) result;

    static if (d == dim - 1)
    {
        // Ridge: coCenter is directly from ridgeLinks
        auto ptr = center in mfd.dimMap!(dim - 1);
        assert(ptr !is null);
        foreach (v; (*ptr).link[])
            result ~= v;
    }
    else
    {
        // Find star(center) via _vertexFacets[center[0]]
        auto ptr = center[0] in mfd._vertexFacets;
        assert(ptr !is null);
        foreach (ref facet; *ptr)
        {
            // Check facet contains all of center
            bool containsAll = true;
            foreach (v; center[1 .. $])
                if (!facet[].canFind(v)) { containsAll = false; break; }
            if (!containsAll) continue;

            // Extract vertices not in center
            foreach (v; facet)
                if (!center[].canFind(v) && !result[].canFind(v))
                    result ~= v;
        }
    }

    result[].sort();
    return result;
}

/*******************************************************************************
Update valid move arrays for all affected simplices in the neighborhood of
allVerts. When adding=true, adds newly valid moves. When adding=false,
removes moves that are about to become invalid.

Type A: subsets of allVerts (degree/link changes directly).
Type B: non-subsets of allVerts with coCenter ⊆ allVerts (coCenter
containment may change).
*/
private void updateValidMoveArrays(bool adding, Vertex, int dim)(
    ref Manifold!(dim, Vertex) mfd,
    const(Vertex)[] allVerts)
{
    // --- Type A: subsets of allVerts ---
    static foreach (d; 0 .. dim)
    {{
        foreach (subset; allVerts.subsetsOfSize(d + 1))
        {
            auto key = subset.staticArray!(d + 1);
            static if (adding)
            {
                if (mfd.isValidNonFacetCenter!d(key))
                    mfd.addValidMoveToArray!d(key);
            }
            else
            {
                if (key in mfd.validMoveIdx!d)
                    mfd.removeValidMoveFromArray!d(key);
            }
        }
    }}

    // --- Type B: non-subsets with coCenter ⊆ allVerts ---
    // Collect facets incident to any vertex in allVerts.
    // Use static thread-local buffers to avoid GC allocation per call.
    alias MFacet = Vertex[dim + 1];
    static MFacet[] incidentFacets;
    incidentFacets.length = 0;
    incidentFacets.assumeSafeAppend;
    foreach (v; allVerts)
    {
        auto ptr = v in mfd._vertexFacets;
        if (ptr is null) continue;
        foreach (ref f; *ptr)
            incidentFacets ~= f;
    }
    incidentFacets.sort();

    static foreach (d; 0 .. dim)
    {{
        enum targetDeg = dim + 1 - d;

        // Collect d-faces of unique incident facets
        static Vertex[d + 1][] candidates;
        candidates.length = 0;
        candidates.assumeSafeAppend;
        MFacet prevFacet;
        bool havePrevFacet = false;
        foreach (ref facet; incidentFacets)
        {
            if (havePrevFacet && facet == prevFacet) continue;
            prevFacet = facet;
            havePrevFacet = true;

            foreach (subset; facet[].subsetsOfSize(d + 1))
                candidates ~= subset.staticArray!(d + 1);
        }
        candidates.sort();

        // Check each unique candidate
        Vertex[d + 1] prevKey;
        bool havePrevKey = false;
        foreach (ref key; candidates)
        {
            if (havePrevKey && key == prevKey) continue;
            prevKey = key;
            havePrevKey = true;

            // Skip Type A (already handled above)
            bool isSubset = true;
            foreach (v; key)
                if (!allVerts.canFind(v)) { isSubset = false; break; }
            if (isSubset) continue;

            static if (adding)
            {
                // Check full validity: degree, coCenter ⊆ allVerts, coCenter ∉ manifold
                uint deg;
                static if (d == 0)
                {
                    auto v0 = key[0];
                    if (v0 >= mfd._vertexDegrees.length || mfd._vertexDegrees[v0] == 0) continue;
                    deg = mfd._vertexDegrees[v0];
                }
                else
                {
                    auto kptr = key in mfd.dimMap!d;
                    if (kptr is null) continue;
                    deg = mfd.extractDegree(*kptr);
                }
                if (deg != targetDeg) continue;

                auto cc = mfd.localCoCenter!d(key);
                bool ccInAllVerts = true;
                foreach (v; cc[])
                    if (!allVerts.canFind(v)) { ccInAllVerts = false; break; }
                if (!ccInAllVerts) continue;

                if (!mfd.contains(cc[]))
                    mfd.addValidMoveToArray!d(key);
            }
            else
            {
                // Remove if in index and affected (coCenter ⊆ allVerts)
                if (key !in mfd.validMoveIdx!d) continue;

                auto cc = mfd.localCoCenter!d(key);
                bool ccInAllVerts = true;
                foreach (v; cc[])
                    if (!allVerts.canFind(v)) { ccInAllVerts = false; break; }
                if (!ccInAllVerts) continue;

                mfd.removeValidMoveFromArray!d(key);
            }
        }
    }}
}

/*******************************************************************************
Sample a uniformly random valid bistellar move in O(1).
For non-facet centers, computes the coCenter via localCoCenter.
For stellar subdivisions (facet center), uses newVertex as the coCenter.
*/
BistellarMove!(dim, Vertex) sampleValidMove(Vertex, int dim, Rng)(
    const ref Manifold!(dim, Vertex) mfd,
    ref Rng rng,
    Vertex newVertex)
{
    import std.random : uniform;
    alias BM = BistellarMove!(dim, Vertex);

    auto total = mfd.validMoveCount;
    assert(total > 0, "no valid moves to sample");
    auto idx = uniform(0, total, rng);

    size_t cumSum = 0;
    static foreach (d; 0 .. dim)
    {{
        auto nMoves = mfd.validMoves!d.length;
        if (idx < cumSum + nMoves)
        {
            auto center = mfd.validMoves!d[idx - cumSum];
            auto cc = mfd.localCoCenter!d(center);
            return BM(center[], cc[]);
        }
        cumSum += nMoves;
    }}

    // Stellar subdivision: center is a random facet, coCenter is newVertex
    auto facetIdx = idx - cumSum;
    auto facetCenter = mfd._facetArray[facetIdx];
    return BM(facetCenter[], [newVertex]);
}

} // version (TrackValidMoves)

///
pure @safe unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));

    octahedron.fVector.shouldEqual([6UL,12,8]);
    octahedron.eulerCharacteristic.shouldEqual(2);

    auto tetrahedron = Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    octahedron.star([1,2]).shouldBeSameSetAs([[0,1,2], [1,2,5]]);    

    octahedron.allBistellarMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[0,1], [0,2], [0,3], [0,4], [1,2], [1,4],  // 1-simplices
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5]]);      

    tetrahedron.allBistellarMoves.shouldBeEmpty;
    
    // TO DO: FINISH CHECKS!
    // octahedron.doMove([1,2]);
    // octahedron.doMove([0,5]);
}

// NOTE: The following unittest cannot be @safe since throwsWithMsg 
// catches an Error
///
pure @system unittest
{
    Manifold!2([[1,2,3,4]]).throwsWithMsg("facet has wrong dimension");

    auto sphere3 = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];    
    Manifold!3(chain(sphere3, sphere3)).throwsWithMsg(
        "tried to insert a facet already in the manifold");
    
}

///
pure @safe unittest
{
    auto m = Manifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    m.allBistellarMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[1], [5], [2, 3], [2, 4], [3, 4]]);

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.allBistellarMoves.map!(mv => mv.center.array).shouldBeSameSetAs(
        [[2, 3], [0, 1], [1, 5], [4, 5], [0, 3], [1, 4],
        [1, 2], [0, 4], [0, 2], [2, 5], [3, 5], [3, 4]]);
}

// NOTE: The following unittest cannot be @safe since throwsWithMsg 
// catches an Error
///
pure @system unittest
{
    // Can't do 2->2 move on the boundary of a 3-simplex
    auto manifold = Manifold!2([[1,2,3],[1,2,4], [1,3,4], [2,3,4]]);
    auto move = BistellarMove!2([1,2], [3,4]); 
    manifold.doMove(move).throwsWithMsg("coCenter of move in manifold");
}

///
pure unittest
{
    alias BM = BistellarMove!2;
    auto manifold = standardSphere!2;
    manifold.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    manifold.doMove(BM([1,2,3], [4]));
    manifold.facets.shouldBeSameSetAs(
        [[0,1,2],[0,1,3], [0,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    manifold.doMove(BM([0,2,3], [7]));
    manifold.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,7],
        [0,3,7], [1,2,4], [1,3,4], [2,3,4], [2,3,7]]);

    manifold.doMove(BM([7],[0,2,3]));
    manifold.doMove(BM([4],[1,2,3]));
    manifold.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);

    manifold = standardSphere!2;
    manifold.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    
    manifold.doMove(BM([1,2,3], [4]));
    manifold.doMove(BM([1,2], [0,4]));
    manifold.facets.shouldBeSameSetAs(
        [[0,1,3], [0,1,4], [0,2,3], [0,2,4], [1,3,4], [2,3,4]]);
    manifold.doMove(BM([0,4], [1,2]));
    manifold.doMove(BM([4], [1,2,3]));
    manifold.facets.shouldBeSameSetAs([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);      

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.doMove(BM([1,2],[0,5]));
    octahedron.degree([0]).shouldEqual(5);
    octahedron.degree([5]).shouldEqual(5);
    octahedron.degree([1]).shouldEqual(3);
    octahedron.degree([2]).shouldEqual(3);

    // We can undo the 2->2 move
    octahedron.doMove(BM([0,5], [1,2]));
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.doMove(BM([0,1,2], [99]));
    octahedron.doMove(BM([99], [0,1,2]));
}

///
unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    auto move = BistellarMove!2([1,2],[0,5]);
    octahedron.doMove(move);
    octahedron.facets.shouldBeSameSetAs([[0, 1, 4], [0, 1, 5],
         [0, 2, 3], [0, 2, 5], [0, 3, 4], [1, 4, 5], [2, 3, 5], [3, 4, 5]]);

    // TO DO: More tests
}


/******************************************************************************
Do a bistellar move replacing the star(center) with star(coCenter)
TO DO: Better docs here!
*/
void doMove(Vertex, int dim)(
    ref Manifold!(dim, Vertex) manifold,
    BistellarMove!(dim, Vertex) move
)
{
    auto center = move.center;
    auto coCenter = move.coCenter;

    assert(manifold.contains(center), "center of move not in manifold");
    assert(!manifold.contains(coCenter), "coCenter of move in manifold");
    version (ExpensiveAsserts)
    {
        if (coCenter.walkLength > 1)
        {
            auto pm = BistellarMove!(dim, Vertex)(center, coCenter);
            assert(manifold.allBistellarMoves.canFind(pm), "not a valid pachner move");
        }
    }

    auto coCenDim = coCenter.length.to!int - 1;
    auto cenDim = center.length.to!int - 1;
    auto oldPiece = productUnion(coCenter.subsetsOfSize(coCenDim), center.only);
    auto newPiece = productUnion(center.subsetsOfSize(cenDim), coCenter.only);

    version (ExpensiveAsserts)
    {
        alias SC = SimplicialComplex!(Vertex, dim);
        alias MFD = Manifold!(dim, Vertex);

        assert(SC(oldPiece).isPureOfDim(dim));
        assert(SC(newPiece).isPureOfDim(dim));
        assert(MFD(chain(oldPiece, newPiece)).numFacets == dim + 2);
        assert(manifold.star(center).map!array.array.sort
            .equal!equal(oldPiece.map!array.array.sort));
    }

    version (TrackValidMoves)
    {
        // Build sorted allVerts = center ∪ coCenter
        Vertex[dim + 2] allVertsBuf;
        int avLen = 0;
        foreach (v; center) allVertsBuf[avLen++] = v;
        foreach (v; coCenter) allVertsBuf[avLen++] = v;
        auto allVerts = allVertsBuf[0 .. avLen];
        allVerts.sort();

        // Remove all valid moves in the neighborhood before the move
        manifold.updateValidMoveArrays!false(allVerts);
    }

    foreach (f; oldPiece)
        manifold.removeFacet(f);
    foreach (f; newPiece)
        manifold.insertFacet(f);

    version (TrackValidMoves)
    {
        // Add back valid moves in the neighborhood after the move
        manifold.updateValidMoveArrays!true(allVerts);
    }
}

auto findCoCenter(Vertex, int dim, C)(
    const ref Manifold!(dim, Vertex) manifold,
    C center
)
if (isIRof!(C, const(Vertex)))
{
    assert(manifold.contains(center));
    return manifold.link(center).joiner.array.dup.sort.uniq.array;
}

/*******************************************************************************
Returns the coCenter for the Pachner move with given center. For efficiency we
also need to know a facet with face center.
*/
auto coCenter(Vertex, int dim, C, F)(
    const ref Manifold!(dim, Vertex) mfd,
    C center,
    F facet
)
if (isIRof!(C, const(Vertex)) && isIRof!(F, const(Vertex)))
{
    assert(mfd.contains(facet));
    assert(mfd.contains(center));
    assert(center.isSubsetOf(facet));
    
    // The coCenter of a facet is a new vertex not in the manifold.
    assert(center.walkLength < dim + 1);

    // TO DO: Clean this up!
    auto ridges = facet.subsetsOfSize(dim)
        .filter!(r => center.isSubsetOf(r)).map!(r => mfd.toRidge(r));
    auto coCenterVerts = ridges.map!(r => mfd.dimMap!(dim - 1)[r].link[])
        .joiner.array.dup.sort.uniq.array;

    assert(coCenterVerts.equal(mfd.findCoCenter(center.array)));
    return coCenterVerts;
}

///
pure @safe unittest
{
    // Octahedron: 6 vertices, 8 triangular facets
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    // Co-center of vertex 0 is its link: the ring of neighbors [1,2,3,4]
    oct.findCoCenter([0]).sort.shouldEqual([1, 2, 3, 4]);

    // Co-center of vertex 5 is its link: [1,2,3,4]
    oct.findCoCenter([5]).sort.shouldEqual([1, 2, 3, 4]);

    // Co-center of vertex 1 is its link: [0,2,4,5]
    oct.findCoCenter([1]).sort.shouldEqual([0, 2, 4, 5]);
}

///
pure @safe unittest
{
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    // Edge [0,1] is shared by facets [0,1,2] and [0,1,4].
    // Link of [0,1] is the two opposite vertices: {2, 4}
    oct.findCoCenter([0, 1]).sort.shouldEqual([2, 4]);

    // Edge [1,2] is shared by facets [0,1,2] and [1,2,5].
    // Link of [1,2] is {0, 5}
    oct.findCoCenter([1, 2]).sort.shouldEqual([0, 5]);

    // Edge [2,3] is shared by facets [0,2,3] and [2,3,5].
    // Link of [2,3] is {0, 5}
    oct.findCoCenter([2, 3]).sort.shouldEqual([0, 5]);
}

///
pure @safe unittest
{
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    // Link of a facet (triangle) in a 2-manifold is a single vertex
    // (the one that would be added in a 3->1 move... but that vertex
    // must not already be in the manifold for the move to be valid).
    // For the octahedron, the link of every facet is empty since
    // all vertices are already present. findCoCenter returns the
    // flattened link vertices regardless of move validity.

    // Facet [0,1,2]: the link is the set of vertices completing each ridge.
    // Ridges [0,1]->{2,4}, [0,2]->{1,3}, [1,2]->{0,5} → flattened unique: {0,1,2,3,4,5} minus center {0,1,2} ...
    // Actually, link of a facet = vertices v such that facet ∪ {v} would be a (dim+1)-simplex in the complex.
    // In a 2-manifold there are no 3-simplices, so the link of a facet is empty.
    oct.findCoCenter([0, 1, 2]).shouldBeEmpty;
}

///
pure @safe unittest
{
    // Triangle as a 1-manifold (circle with 3 edges)
    auto circle = Manifold!1([[0,1], [1,2], [0,2]]);

    // Link of vertex 0 in a circle is its two neighbors: {1, 2}
    circle.findCoCenter([0]).sort.shouldEqual([1, 2]);
    circle.findCoCenter([1]).sort.shouldEqual([0, 2]);
    circle.findCoCenter([2]).sort.shouldEqual([0, 1]);

    // Link of an edge (facet) in a 1-manifold is empty
    circle.findCoCenter([0, 1]).shouldBeEmpty;
}

///
pure @safe unittest
{
    // Standard 1-sphere: boundary of a triangle (3 edges, 3 vertices)
    immutable s1 = standardSphere!1;
    // Each vertex links to the other two
    foreach (v; 0 .. 3)
        s1.findCoCenter([v]).length.shouldEqual(2);

    // Standard 2-sphere: boundary of a tetrahedron (4 triangles, 4 vertices)
    immutable s2 = standardSphere!2;
    // Each vertex links to the other three
    foreach (v; 0 .. 4)
        s2.findCoCenter([v]).length.shouldEqual(3);
    // Each edge links to the other two vertices
    s2.findCoCenter([0, 1]).sort.shouldEqual([2, 3]);
    s2.findCoCenter([0, 2]).sort.shouldEqual([1, 3]);
    s2.findCoCenter([1, 3]).sort.shouldEqual([0, 2]);
}

///
pure @safe unittest
{
    // For any center, |co-center| = degree(center) - |center| + 1
    // because the link of a k-simplex in a d-manifold consists of
    // (d-k)-simplices, and its vertex set has degree - |center| + 1 elements.
    // More precisely: |co-center| = degree(center) - |center| + 1...
    // Let's just verify the relationship empirically on the octahedron.
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    // Vertices: degree 4, co-center length 4
    foreach (v; 0 .. 6)
    {
        auto cc = oct.findCoCenter([v]);
        cc.length.shouldEqual(oct.degree([v]));
    }

    // Edges: degree 2 (each edge in 2 triangles), co-center length 2
    auto edges = oct.simplices(1);
    foreach (e; edges)
    {
        auto cc = oct.findCoCenter(e);
        cc.length.shouldEqual(oct.degree(e));
    }
}

///
pure @safe unittest
{
    auto mfd = Manifold!1([[0,1],[0,2],[1,2]]);
    mfd.facets.sort.shouldBeSameSetAs([[0,1],[0,2],[1,2]]);
}

///
pure @safe unittest
{
    auto mfd = Manifold!1([[0,1],[0,2],[1,2]]);
    mfd.star([1]).map!array.array.sort.shouldBeSameSetAs([[0,1],[1,2]]);
}

///
pure @safe unittest
{
    auto mfd = Manifold!1([[0,1],[0,2],[1,2]]);
    mfd.link([1]).map!(r => r.array).array.sort.shouldBeSameSetAs([[0],[2]]);

    assert(mfd.link([1]).empty == false);
}

///
pure @safe unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge07 = [0,7];

    () pure nothrow @nogc @safe {
        assert(sc.contains(edge01[]));
        assert(!sc.contains(edge07[]));
    }();  
}

/******************************************************************************
* Returns a manifold loaded from the file specified by fileName. If fileName
* is the empty string, the returned manifold is the standard sphere.
*/
Manifold!(dim, Vertex) loadManifold(int dim, Vertex = int)(string fileName)
{
    SimplicialComplex!Vertex sc = loadSimplicialComplex!Vertex(fileName);
    return Manifold!(dim, Vertex)(sc.facets);
}

///
@system unittest
{
    auto m = loadManifold!2(
        "data/manifold_sampler_unittest_load.dat");
    auto expected = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5],
        [0, 5, 6], [0, 1, 6], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
    assert(m == expected);

    assertThrown(loadManifold!2(
        "data/manifold_sampler_unittest_bad_load.dat"));
}

///
@system unittest
{
    auto fileName = "data/manifold_sampler_unittest_save.dat";
    auto sphere = standardSphere!4;
    sphere.saveTo(fileName);
    Manifold!4 loadedSphere = loadManifold!4(fileName);
    assert(loadedSphere == sphere);
}

/******************************************************************************
Returns the mean degree of simplices of the given dimension 'dim'
*/
real meanDegree(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfdDim);
    immutable nSimps = mfd.fVector[dim];
    immutable nFacets = mfd.fVector.back;
    immutable simpsPerFacet = binomial(mfdDim + 1, dim + 1);
    return simpsPerFacet * nFacets / real(nSimps);
}
///
pure @safe unittest
{
    foreach (dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        assert((dim + 1).iota.all!(d => m.meanDegree(d) == (dim - d + 1)));
    }

    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    // pyramid has two vertices with deg(v)=3 and three vertices with deg(v)=4
    assert(pyramid.meanDegree(0).isClose((2 * 3 + 3 * 4)/5.0));
    assert(pyramid.meanDegree(1).isClose(2.0));
    assert(pyramid.meanDegree(2).isClose(1.0));
}

/******************************************************************************
Returns total (degree(s))^^2 for all simplices s of given dimension 'dim'.
*/
ulong totalSquareDegree(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfd.dimension);
    if (dim == mfdDim - 1)
    {
        // All codimension-1 simplices (ridges) have degree 2
        return 4 * mfd.numSimplices[dim];
    }
    else if (dim == mfdDim)
    {
        // All facets have degree 1
        return mfd.numSimplices[dim];
    }
    return mfd.totSqrDegrees[dim];
}
///
pure @safe unittest
{
    foreach (dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        foreach (d; 0 .. dim + 1)
        {
            // (dim + 2) choose (d + 1) faces of dimension d
            // each with degree (dim - d + 1)
            m.totalSquareDegree(d).shouldEqual(
                binomial(dim + 2, d + 1) * (dim - d + 1)^^2);
        }
    }

    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    pyramid.totalSquareDegree(0).shouldEqual(2 * 3^^2 + 3 * 4^^2);
    pyramid.totalSquareDegree(1).shouldEqual(9 * 2^^2);
    pyramid.totalSquareDegree(2).shouldEqual(6 * 1^^2);
}

/******************************************************************************
Returns the variance in the degree of simplices of the given dimension 'dim'
*/
real degreeVariance(Vertex, int mfdDim)(
    const ref Manifold!(mfdDim, Vertex) mfd,
    int dim)
{
    assert(dim >= 0);
    assert(dim <= mfdDim);
    if (dim >= mfdDim - 1)
    {
        return 0;   // No variance in ridge and facet dimension
    }

    immutable meanDeg = mfd.meanDegree(dim);
    immutable meanSqrDeg = mfd.totalSquareDegree(dim) / real(mfd.fVector[dim]);
    return meanSqrDeg - meanDeg^^2;
}
///
pure @safe unittest
{
    foreach (dim; staticIota!(1, 8))
    {
        auto m = standardSphere!dim;
        assert((dim+1).iota.all!(d => m.degreeVariance(d) == 0));
    }

    // pyramid has two vertices with deg(v)=3 and three vertices with deg(v)=4
    auto pyramid = Manifold!2(
        [[0,1,2], [0,2,3], [0,1,3], [1,2,4], [2,3,4], [1,3,4]]);

    immutable real md0 = pyramid.meanDegree(0);
    assert(pyramid.degreeVariance(0)
        .isClose((2 * (3 - md0)^^2 + 3 * (4 - md0)^^2)/5.0));

    assert(pyramid.degreeVariance(1).isClose(0.0));
    assert(pyramid.degreeVariance(2).isClose(0.0));

}

void undoMove(int dim, Vertex)(
    ref Manifold!(dim, Vertex) manifold,
    ref const(BistellarMove!(dim, Vertex)) move)
{
    auto inverseMove = BistellarMove!(dim, Vertex)(move.coCenter, move.center);
    manifold.doMove(inverseMove);
}

///
pure @safe unittest
{
    alias BM = BistellarMove!2;

    // 1->3 move then undo on standard sphere
    auto sphere = standardSphere!2;
    auto origFacets = sphere.facets.sort;
    auto origFVector = sphere.fVector;

    auto move1to3 = BM([1,2,3], [4]);
    sphere.doMove(move1to3);
    assert(sphere.facets.sort != origFacets);

    sphere.undoMove(move1to3);
    sphere.facets.sort.shouldEqual(origFacets);
    sphere.fVector.shouldEqual(origFVector);
}

///
pure @safe unittest
{
    alias BM = BistellarMove!2;

    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    auto origFacets = oct.facets.sort;

    // 2->2 move swaps the diagonal of a quadrilateral
    auto move2to2 = BM([1,2], [0,5]);
    oct.doMove(move2to2);

    // Degrees changed
    oct.degree([1]).shouldEqual(3);
    oct.degree([2]).shouldEqual(3);
    oct.degree([0]).shouldEqual(5);
    oct.degree([5]).shouldEqual(5);

    oct.undoMove(move2to2);
    oct.facets.sort.shouldEqual(origFacets);

    // All degrees restored to 4
    assert(oct.simplices(0).all!(s => oct.degree(s) == 4));
}

///
pure @safe unittest
{
    alias BM = BistellarMove!2;

    auto sphere = standardSphere!2;
    auto origFacets = sphere.facets.sort;
    auto origFVector = sphere.fVector;

    // Apply two moves
    auto move1 = BM([1,2,3], [4]);
    sphere.doMove(move1);

    auto move2 = BM([0,2,3], [7]);
    sphere.doMove(move2);

    // Undo in reverse order
    sphere.undoMove(move2);
    sphere.undoMove(move1);

    sphere.facets.sort.shouldEqual(origFacets);
    sphere.fVector.shouldEqual(origFVector);
}

///
pure @safe unittest
{
    alias BM = BistellarMove!2;

    // Start with octahedron, do a 1->3 move, then undo it
    auto oct = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    auto origFacets = oct.facets.sort;

    auto move1to3 = BM([0,1,2], [99]);
    oct.doMove(move1to3);
    oct.numFacets.shouldEqual(10);

    // Now do a 3->1 to undo it
    oct.undoMove(move1to3);
    oct.facets.sort.shouldEqual(origFacets);
    oct.numFacets.shouldEqual(8);
}




@system unittest
{
    auto hyperbolicDodecahedral = loadManifold!3(
        "data/manifold_sampler_unittest_dodecahedral.dat");

    () pure @safe {
        // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
        assert(standardSphere!7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
                9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
                (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

        // TO DO: Get reference for this...
        assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);
    } ();
}

string[] findProblems(Vertex, int dim)(const ref Manifold!(dim, Vertex) mfd)
{
    alias SC = SimplicialComplex!Vertex;
    string[] problems;

    // Build a SimplicialComplex for topology checks
    auto sc = mfd.toSimplicialComplex;

    if (!sc.isPureOfDim(dim))
    {
        problems ~= "not all facets have the correct dimension";
    }

    if (!sc.isConnected)
    {
        problems ~= "facets do not define a connected simplicial complex";
    }

    static if (dim >= 1)
    {
        if (!mfd.simplices(dim - 1).all!(
            s => sc.star(s).walkLength == 2))
        {
            problems ~= "found a ridge with degree not equal to 2";
        }
    }

    static if (dim >= 2)
    {
        if (!mfd.simplices(dim - 2).all!(s => SC(sc.link(s)).isCircle))
        {
            problems ~= "found a hinge whose link is not a circle";
        }
    }

    static if (dim >= 3)
    {
        if (!mfd.simplices(dim - 3).all!(s => SC(sc.link(s)).is2Sphere))
        {
            problems ~= "found a codimension-3 simplex whose link is not a 2-sphere";
        }
    }

    if (mfd.numSimplices[] != sc.fVector[])
    {
        problems ~= "number of simplices incorrect";
    }

    foreach (d; 0 .. dim - 1)
    {
        auto correct = mfd.simplices(d).map!(s => mfd.degree(s)^^2).sum;
        if (mfd.totSqrDegrees[d] != correct)
        {
            problems ~= "found incorrect total squared degree";
            break;
        }
    }

    // Check SC simplices are in dimMaps
    foreach (d; 0 .. dim + 1)
    {
        foreach (s; sc.simplices(d))
        {
            if (!mfd.contains(s))
            {
                problems ~= "found a simplex in SC that is not in dimMaps";
                goto done;
            }
        }
    }
    done:

    // Check vertex degrees match SC
    foreach (v; 0 .. cast(int) mfd._vertexDegrees.length)
    {
        if (mfd._vertexDegrees[v] > 0)
        {
            if (!sc.contains([cast(const(Vertex)) v]))
            {
                problems ~= "found a vertex in _vertexDegrees that is not in SC";
                goto done2;
            }
            if (mfd._vertexDegrees[v] != sc.star([cast(const(Vertex)) v]).walkLength)
            {
                problems ~= "found a vertex with incorrect degree";
                goto done2;
            }
        }
    }
    // Check dimMap simplices (dims 1..dim-1) are in SC, and degrees match
    static foreach (d; 1 .. dim)
    {{
        foreach (kv; mfd.dimMap!d.byKeyValue)
        {
            if (!sc.contains(kv.key[]))
            {
                problems ~= "found a simplex in dimMaps that is not in SC";
                goto done2;
            }
            auto deg = mfd.extractDegree(kv.value);
            if (deg != sc.star(kv.key[]).walkLength)
            {
                problems ~= "found a simplex in dimMaps with incorrect degree";
                goto done2;
            }
        }
    }}
    // Check facets in _facetArrayIdx are in SC
    foreach (ref f; mfd._facetArray)
    {
        if (!sc.contains(f[]))
        {
            problems ~= "found a facet in _facetArray that is not in SC";
            goto done2;
        }
    }
    done2:

    // Check ridge links
    foreach (kv; mfd.dimMap!(dim - 1).byKeyValue)
    {
        auto ridge = kv.key[];
        auto link = kv.value.link[];

        if (link.walkLength != sc.star(ridge).walkLength)
        {
            problems ~= "found a ridge in dimMaps whose link has incorrect number of vertices";
            break;
        }

        if (link.array != sc.link(ridge).joiner.array.dup.sort.array)
        {
            problems ~= "found a ridge in dimMaps whose link has the wrong vertices";
            break;
        }
    }

    foreach (ridge; sc.simplices(dim - 1))
    {
        if (mfd.toRidge(ridge) !in mfd.dimMap!(dim - 1))
        {
            problems ~= "found a ridge in SC that is not in dimMaps";
            break;
        }
    }

    return problems;
}
///
pure @safe unittest
{
    auto m3 = standardSphere!3;
    assert(m3.findProblems.empty);

    m3 = Manifold!3([[1,2,3,4]]);
    m3.findProblems.shouldBeSameSetAs([
        "found a hinge whose link is not a circle",
        "found a ridge with degree not equal to 2",
        "found a codimension-3 simplex whose link is not a 2-sphere"
    ]);

    // These facets separately define 2-spheres, but share the vertex 1
    auto sphere2a = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]];
    auto sphere2b = [[1,5,6], [1,5,7], [1,6,7], [5,6,7]];
    auto m2 = Manifold!2(chain(sphere2a, sphere2b));
    m2.findProblems.shouldBeSameSetAs([
        "found a hinge whose link is not a circle"
    ]);

    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5],
        [3,4,5], [1,4,5]];
    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];
    m2 = Manifold!2(chain(octahedron, sphere));
    m2.findProblems.shouldBeSameSetAs([
        "facets do not define a connected simplicial complex"
    ]);

    // These facets separately define 2-spheres, but share the vertex 1
    auto sphere3a = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];
    auto sphere3b = [[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]];
    m3 = Manifold!3(chain(sphere3a, sphere3b));
    m3.findProblems.shouldBeSameSetAs([
        "found a codimension-3 simplex whose link is not a 2-sphere"
    ]);

    // TO DO: This test trips an internal assert. I don't want to remove
    // the internal assert. Decide what to do...

    // auto m1 = Manifold!1([[0,1],[1,2],[0,2],[0,3]]);
    // m3.findProblems.shouldBeSameSetAs([
    //     "found ...."
    // ]);
}

size_t[] degreeHistogram(Vertex, int mfdDim)(
    auto ref const(Manifold!(mfdDim, Vertex)) mfd,
    size_t dim)
{
    size_t[] histogram;

    void addToHistogram(size_t deg)
    {
        if (deg > histogram.length)
            histogram.length = deg;
        ++histogram[deg - 1];
    }

    final switch (cast(int) dim)
    {
        case 0:
        {
            foreach (deg; mfd._vertexDegrees)
                if (deg > 0) addToHistogram(deg);
            return histogram;
        }
        static foreach (d; 1 .. mfdDim)
        {
            case d:
            {
                foreach (kv; mfd.dimMap!d.byKeyValue)
                    addToHistogram(mfd.extractDegree(kv.value));
                return histogram;
            }
        }
        case mfdDim:
        {
            // All facets have degree 1
            if (mfd._facetArray.length > 0)
            {
                histogram.length = 1;
                histogram[0] = mfd._facetArray.length;
            }
            return histogram;
        }
    }
}
///
unittest
{
    // TO DO: More tests...

    auto m = standardSphere!3;
    assert(4.iota.map!(k => m.degreeHistogram(k)).equal!equal(
        [[0,0,0,5], [0,0,10], [0,10], [5]]));
}


void saveDualGraphTo(int dimension, Vertex = int)(
    const ref Manifold!(dimension, Vertex) mfd,
    string fileName)
{
    size_t[Vertex[]] facetIndex;
    foreach(index, facet; mfd.facets.enumerate)
    {
        facetIndex[facet] = index;
    }

    size_t[][] edges;
    foreach(index, facet; mfd.facets.enumerate)
    {
        foreach(ridge; facet.subsetsOfSize(dimension))
        {
            auto linkVertices = mfd.dimMap!(dimension - 1)[mfd.toRidge(ridge)].link;
            Vertex[] oppositeFacet;
            if (facet.canFind(linkVertices.front))
            {
                oppositeFacet = chain(ridge, linkVertices.back.only).array.dup.sort.array;
            }
            else
            {
                oppositeFacet = chain(ridge, linkVertices.front.only).array.dup.sort.array;
            }

            auto oppositeIndx = facetIndex[oppositeFacet];

            // To avoid duplicate edges, only write edge if second vertex
            // corresponds to facet further down the facet list
            if (oppositeIndx > index)
            {
                edges ~= [index, oppositeIndx];
            }
        }
    }

    auto saveFile = File(fileName, "w"); // Open in write-only mode
    foreach(edge; edges)
    {
        saveFile.writeln(edge.front, " ", edge.back);
    }
}
/// TO DO: Unittests for saveDualGraph


pure /* @safe */ unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    auto m1 = octahedron;
    auto m2 = m1;

    auto move = BistellarMove!2([0,1],[2,4]);
    m2.doMove(move);

    assert(m2.facets.array.sort != octahedron.facets.array.sort);
    assert(m2.facets.array.sort != m1.facets.array.sort);
    assert(m1.facets.array.sort == octahedron.facets.array.sort);
}

pure @safe unittest
{
    auto vertices = [0,1,2];
    auto edge1 = vertices[0..2];
    auto edge2 = vertices[1..3];
    // edge1 and edge2 slices now overlap at index 1

    auto edge3 = [0,2];

    auto m = Manifold!1([edge1, edge2, edge3]);
    m.facets.shouldBeSameSetAs([[0,1],[1,2],[0,2]]);
    vertices[1] = 42;
    assert(!m.simplices(0).canFind([42]));
    m.facets.shouldBeSameSetAs([[0,1],[1,2],[0,2]]);
}