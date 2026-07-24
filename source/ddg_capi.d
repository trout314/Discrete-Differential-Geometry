/// C API layer for the discrete differential geometry library.
/// Provides extern(C) functions for use via ctypes/FFI from Python.
module ddg_capi;

import std.algorithm, std.array, std.conv, std.range, std.stdio, std.string,
    std.math, std.random;
import core.memory : GC;
import algorithms, manifold, manifold_examples, manifold_moves,
    sampler, simplicial_complex, utility;

// ---------------------------------------------------------------------------
// Thread-local error string
// ---------------------------------------------------------------------------

private string lastErrorMsg;

private void setError(string msg) nothrow
{
    try { lastErrorMsg = msg; } catch (Exception) {}
}

private void clearError() nothrow
{
    try { lastErrorMsg = null; } catch (Exception) {}
}

extern(C) const(char)* ddg_last_error() nothrow
{
    if (lastErrorMsg.length == 0) return null;
    return lastErrorMsg.ptr;
}

// ---------------------------------------------------------------------------
// Opaque handle types
// ---------------------------------------------------------------------------

private struct ManifoldHandle
{
    int dim;
    void* ptr;
}

// GC-allocated wrappers so we can take stable pointers.
private struct ManifoldWrapper(int dim)
{
    Manifold!(dim, int) mfd;
}

// ---------------------------------------------------------------------------
// GC root management
// ---------------------------------------------------------------------------
// Objects allocated with `new` live on D's GC heap, but when the only
// reference is held by Python (via ctypes void*), the GC cannot see it.
// We pin objects with GC.addRoot when handing them to C and unpin with
// GC.removeRoot when the C side calls _free.

private void pinForC(void* ptr) nothrow
{
    if (ptr !is null) GC.addRoot(ptr);
}

private void unpinForC(void* ptr) nothrow
{
    if (ptr !is null) GC.removeRoot(ptr);
}

// ---------------------------------------------------------------------------
// GC control
// ---------------------------------------------------------------------------

/// Trigger a D garbage collection cycle. Call this periodically from Python
/// to reclaim temporaries created by chooseRandomMove, coCenter, etc.
extern(C) void ddg_gc_collect() nothrow
{
    try { GC.collect(); } catch (Exception) {}
}

/// Minimize the D GC heap, returning free pages to the OS.
extern(C) void ddg_gc_minimize() nothrow
{
    try { GC.minimize(); } catch (Exception) {}
}

/// Return D GC heap stats: used bytes, free bytes, total mapped bytes.
extern(C) void ddg_gc_stats(long* used_bytes, long* free_bytes) nothrow
{
    try
    {
        auto stats = GC.stats;
        if (used_bytes !is null) *used_bytes = cast(long) stats.usedSize;
        if (free_bytes !is null) *free_bytes = cast(long) stats.freeSize;
    }
    catch (Exception) {}
}

/// Seed the thread-local RNG driving all move proposals and Metropolis
/// accepts. Chains are reproducible given (initial state, params, seed);
/// without a call, the default unpredictable per-process seed applies.
extern(C) void ddg_set_random_seed(uint seed) nothrow
{
    import std.random : rndGen;
    try { rndGen.seed(seed); } catch (Exception) {}
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

private size_t strlen(const(char)* s) pure nothrow @nogc
{
    size_t len = 0;
    while (s[len] != '\0') len++;
    return len;
}

private string toStr(const(char)* s)
{
    return s[0 .. strlen(s)].idup;
}


// ---------------------------------------------------------------------------
// Manifold lifecycle
// ---------------------------------------------------------------------------

extern(C) void* ddg_manifold_standard_sphere(int dim) nothrow
{
    clearError();
    try
    {
        ManifoldHandle* h = new ManifoldHandle;
        h.dim = dim;
        switch (dim)
        {
            case 2: auto w = new ManifoldWrapper!2; w.mfd = standardSphere!2; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = standardSphere!3; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = standardSphere!4; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            default:
                setError("unsupported dimension: " ~ dim.to!string ~ " (must be 2, 3, or 4)");
                return null;
        }
        pinForC(cast(void*) h);
        return cast(void*) h;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_manifold_from_facets(int dim, const(int)* data, int num_facets) nothrow
{
    clearError();
    try
    {
        int verts_per_facet = dim + 1;
        const(int)[][] facets;
        foreach (i; 0 .. num_facets)
            facets ~= data[i * verts_per_facet .. (i + 1) * verts_per_facet];

        ManifoldHandle* h = new ManifoldHandle;
        h.dim = dim;
        switch (dim)
        {
            case 2: auto w = new ManifoldWrapper!2; w.mfd = Manifold!(2, int)(facets); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = Manifold!(3, int)(facets); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = Manifold!(4, int)(facets); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            default:
                setError("unsupported dimension: " ~ dim.to!string);
                return null;
        }
        pinForC(cast(void*) h);
        return cast(void*) h;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_manifold_load(const(char)* path, int dim) nothrow
{
    clearError();
    try
    {
        string fileName = toStr(path);
        ManifoldHandle* h = new ManifoldHandle;
        h.dim = dim;
        switch (dim)
        {
            case 2: auto w = new ManifoldWrapper!2; w.mfd = loadManifold!2(fileName); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = loadManifold!3(fileName); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = loadManifold!4(fileName); h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            default:
                setError("unsupported dimension: " ~ dim.to!string);
                return null;
        }
        pinForC(cast(void*) h);
        return cast(void*) h;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_manifold_copy(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return null; }
        ManifoldHandle* src = cast(ManifoldHandle*) handle;
        ManifoldHandle* h = new ManifoldHandle;
        h.dim = src.dim;
        switch (src.dim)
        {
            case 2: auto w = new ManifoldWrapper!2; w.mfd = (cast(ManifoldWrapper!2*) src.ptr).mfd; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = (cast(ManifoldWrapper!3*) src.ptr).mfd; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = (cast(ManifoldWrapper!4*) src.ptr).mfd; h.ptr = cast(void*) w; pinForC(cast(void*) w); break;
            default: setError("bad dimension"); return null;
        }
        pinForC(cast(void*) h);
        return cast(void*) h;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void ddg_manifold_free(void* handle) nothrow
{
    if (handle is null) return;
    auto h = cast(ManifoldHandle*) handle;
    unpinForC(h.ptr);
    unpinForC(handle);
}

// ---------------------------------------------------------------------------
// Manifold queries
// ---------------------------------------------------------------------------

extern(C) int ddg_manifold_dimension(void* handle) nothrow
{
    if (handle is null) return -1;
    return (cast(ManifoldHandle*) handle).dim;
}

extern(C) long ddg_manifold_num_facets(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return cast(long)(cast(ManifoldWrapper!2*) h.ptr).mfd.numFacets;
            case 3: return cast(long)(cast(ManifoldWrapper!3*) h.ptr).mfd.numFacets;
            case 4: return cast(long)(cast(ManifoldWrapper!4*) h.ptr).mfd.numFacets;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_euler_characteristic(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return 0; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return (cast(ManifoldWrapper!2*) h.ptr).mfd.eulerCharacteristic;
            case 3: return (cast(ManifoldWrapper!3*) h.ptr).mfd.eulerCharacteristic;
            case 4: return (cast(ManifoldWrapper!4*) h.ptr).mfd.eulerCharacteristic;
            default: setError("bad dimension"); return 0;
        }
    }
    catch (Exception e) { setError(e.msg); return 0; }
}

/// Freeze (frozen=1) or unfreeze (frozen=0) a list of vertices. The sampler
/// then rejects any move whose support contains a frozen vertex, preserving
/// the frozen set's entire closed star (facets and hinge degrees/curvature).
extern(C) int ddg_manifold_freeze_vertices(void* handle, const(int)* verts,
                                           long n, int frozen) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        if (verts is null && n > 0) { setError("null verts"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2:
                foreach (i; 0 .. n)
                    (cast(ManifoldWrapper!2*) h.ptr).mfd.setVertexFrozen(verts[i], frozen != 0);
                return 0;
            case 3:
                foreach (i; 0 .. n)
                    (cast(ManifoldWrapper!3*) h.ptr).mfd.setVertexFrozen(verts[i], frozen != 0);
                return 0;
            case 4:
                foreach (i; 0 .. n)
                    (cast(ManifoldWrapper!4*) h.ptr).mfd.setVertexFrozen(verts[i], frozen != 0);
                return 0;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Unfreeze all vertices (drops the flag storage; restores the fast path).
extern(C) int ddg_manifold_clear_frozen(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: (cast(ManifoldWrapper!2*) h.ptr).mfd.clearFrozenVertices(); return 0;
            case 3: (cast(ManifoldWrapper!3*) h.ptr).mfd.clearFrozenVertices(); return 0;
            case 4: (cast(ManifoldWrapper!4*) h.ptr).mfd.clearFrozenVertices(); return 0;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Query: 1 if vertex is frozen, 0 if not, -1 on error.
extern(C) int ddg_manifold_vertex_frozen(void* handle, int v) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return (cast(ManifoldWrapper!2*) h.ptr).mfd.vertexFrozen(v) ? 1 : 0;
            case 3: return (cast(ManifoldWrapper!3*) h.ptr).mfd.vertexFrozen(v) ? 1 : 0;
            case 4: return (cast(ManifoldWrapper!4*) h.ptr).mfd.vertexFrozen(v) ? 1 : 0;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Number of frozen vertices (-1 on error).
extern(C) long ddg_manifold_num_frozen(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return cast(long)(cast(ManifoldWrapper!2*) h.ptr).mfd.numFrozenVertices();
            case 3: return cast(long)(cast(ManifoldWrapper!3*) h.ptr).mfd.numFrozenVertices();
            case 4: return cast(long)(cast(ManifoldWrapper!4*) h.ptr).mfd.numFrozenVertices();
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_f_vector(void* handle, long* out_buf, int buf_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        const(uint)[] fv;
        switch (h.dim)
        {
            case 2: fv = (cast(ManifoldWrapper!2*) h.ptr).mfd.fVector; break;
            case 3: fv = (cast(ManifoldWrapper!3*) h.ptr).mfd.fVector; break;
            case 4: fv = (cast(ManifoldWrapper!4*) h.ptr).mfd.fVector; break;
            default: setError("bad dimension"); return -1;
        }

        int n = cast(int) fv.length;
        if (n > buf_len) n = buf_len;
        foreach (i; 0 .. n) out_buf[i] = fv[i];
        return n;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_is_orientable(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return (cast(ManifoldWrapper!2*) h.ptr).mfd.isOrientable ? 1 : 0;
            case 3: return (cast(ManifoldWrapper!3*) h.ptr).mfd.isOrientable ? 1 : 0;
            case 4: return (cast(ManifoldWrapper!4*) h.ptr).mfd.isOrientable ? 1 : 0;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_num_connected_components(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        int countCC(int dim)(ManifoldWrapper!dim* mw)
        {
            auto sc = mw.mfd.toSimplicialComplex;
            return cast(int) sc.connectedComponents.walkLength;
        }

        switch (h.dim)
        {
            case 2: return countCC!2(cast(ManifoldWrapper!2*) h.ptr);
            case 3: return countCC!3(cast(ManifoldWrapper!3*) h.ptr);
            case 4: return countCC!4(cast(ManifoldWrapper!4*) h.ptr);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// Manifold data — facets and simplices
// ---------------------------------------------------------------------------

extern(C) long ddg_manifold_facets(void* handle, int* out_data) nothrow
{
    clearError();
    if (handle is null) { setError("null handle"); return -1; }
    auto h = cast(ManifoldHandle*) handle;

    switch (h.dim)
    {
        case 2: return (cast(ManifoldWrapper!2*) h.ptr).mfd.writeFacetsToBuffer(out_data);
        case 3: return (cast(ManifoldWrapper!3*) h.ptr).mfd.writeFacetsToBuffer(out_data);
        case 4: return (cast(ManifoldWrapper!4*) h.ptr).mfd.writeFacetsToBuffer(out_data);
        default: setError("bad dimension"); return -1;
    }
}

extern(C) long ddg_manifold_simplices(void* handle, int simp_dim, int* out_data) nothrow
{
    clearError();
    if (handle is null) { setError("null handle"); return -1; }
    auto h = cast(ManifoldHandle*) handle;

    switch (h.dim)
    {
        case 2: return (cast(ManifoldWrapper!2*) h.ptr).mfd.writeSimplicesToBuffer(simp_dim, out_data);
        case 3: return (cast(ManifoldWrapper!3*) h.ptr).mfd.writeSimplicesToBuffer(simp_dim, out_data);
        case 4: return (cast(ManifoldWrapper!4*) h.ptr).mfd.writeSimplicesToBuffer(simp_dim, out_data);
        default: setError("bad dimension"); return -1;
    }
}

extern(C) long ddg_manifold_degree(void* handle, const(int)* simplex_data, int simplex_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        auto simplex = simplex_data[0 .. simplex_len];

        switch (h.dim)
        {
            case 2: return cast(long)(cast(ManifoldWrapper!2*) h.ptr).mfd.degree(simplex);
            case 3: return cast(long)(cast(ManifoldWrapper!3*) h.ptr).mfd.degree(simplex);
            case 4: return cast(long)(cast(ManifoldWrapper!4*) h.ptr).mfd.degree(simplex);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Count the number of valid Pachner moves (including stellar subdivisions).
extern(C) long ddg_manifold_count_valid_moves(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return cast(long)(cast(ManifoldWrapper!2*) h.ptr).mfd.countValidBistellarMoves;
            case 3: return cast(long)(cast(ManifoldWrapper!3*) h.ptr).mfd.countValidBistellarMoves;
            case 4: return cast(long)(cast(ManifoldWrapper!4*) h.ptr).mfd.countValidBistellarMoves;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Return 1/V(x): the importance weight that corrects the sampler's
/// stationary distribution (exp(-objective)*V) back to exp(-objective).
extern(C) double ddg_manifold_importance_weight(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        double compute(int dim)(ManifoldWrapper!dim* mw)
        {
            return 1.0 / cast(double) mw.mfd.countValidBistellarMoves;
        }

        switch (h.dim)
        {
            case 2: return compute!2(cast(ManifoldWrapper!2*) h.ptr);
            case 3: return compute!3(cast(ManifoldWrapper!3*) h.ptr);
            case 4: return compute!4(cast(ManifoldWrapper!4*) h.ptr);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) double ddg_manifold_mean_degree(void* handle, int simp_dim) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return cast(double)(cast(ManifoldWrapper!2*) h.ptr).mfd.meanDegree(simp_dim);
            case 3: return cast(double)(cast(ManifoldWrapper!3*) h.ptr).mfd.meanDegree(simp_dim);
            case 4: return cast(double)(cast(ManifoldWrapper!4*) h.ptr).mfd.meanDegree(simp_dim);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// Manifold I/O
// ---------------------------------------------------------------------------

extern(C) int ddg_manifold_save(void* handle, const(char)* path) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        string fileName = toStr(path);
        switch (h.dim)
        {
            case 2: (cast(ManifoldWrapper!2*) h.ptr).mfd.saveTo(fileName); break;
            case 3: (cast(ManifoldWrapper!3*) h.ptr).mfd.saveTo(fileName); break;
            case 4: (cast(ManifoldWrapper!4*) h.ptr).mfd.saveTo(fileName); break;
            default: setError("bad dimension"); return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_save_with_comments(
    void* handle, const(char)* path,
    const(char*)* comments, int num_comments) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        string fileName = toStr(path);
        string[] commentStrings;
        foreach (i; 0 .. num_comments)
            commentStrings ~= toStr(comments[i]);
        switch (h.dim)
        {
            case 2: (cast(ManifoldWrapper!2*) h.ptr).mfd.saveTo(fileName, commentStrings); break;
            case 3: (cast(ManifoldWrapper!3*) h.ptr).mfd.saveTo(fileName, commentStrings); break;
            case 4: (cast(ManifoldWrapper!4*) h.ptr).mfd.saveTo(fileName, commentStrings); break;
            default: setError("bad dimension"); return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_save_edge_graph(void* handle, const(char)* path) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        string fileName = toStr(path);
        switch (h.dim)
        {
            case 2: (cast(ManifoldWrapper!2*) h.ptr).mfd.saveEdgeGraphTo(fileName); break;
            case 3: (cast(ManifoldWrapper!3*) h.ptr).mfd.saveEdgeGraphTo(fileName); break;
            case 4: (cast(ManifoldWrapper!4*) h.ptr).mfd.saveEdgeGraphTo(fileName); break;
            default: setError("bad dimension"); return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_save_dual_graph(void* handle, const(char)* path) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        string fileName = toStr(path);
        switch (h.dim)
        {
            case 2: (cast(ManifoldWrapper!2*) h.ptr).mfd.saveDualGraphTo(fileName); break;
            case 3: (cast(ManifoldWrapper!3*) h.ptr).mfd.saveDualGraphTo(fileName); break;
            case 4: (cast(ManifoldWrapper!4*) h.ptr).mfd.saveDualGraphTo(fileName); break;
            default: setError("bad dimension"); return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// Manifold moves
// ---------------------------------------------------------------------------

extern(C) int ddg_manifold_do_pachner_move(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        int doRandomMove(int dim)(ManifoldWrapper!dim* mw)
        {
            auto moves = mw.mfd.allBistellarMoves;
            if (moves.length == 0)
            {
                auto facet = mw.mfd.randomFacetOfDim(dim);
                int newVertex = cast(int) mw.mfd.fVector[0];
                while (mw.mfd.contains(newVertex.only))
                    newVertex++;
                auto bm = BistellarMove!(dim, int)(facet, newVertex.only);
                mw.mfd.doMove(bm);
                return 1;
            }

            auto idx = uniform(0, moves.length);
            mw.mfd.doMove(moves[idx]);
            return 1;
        }

        switch (h.dim)
        {
            case 2: return doRandomMove!2(cast(ManifoldWrapper!2*) h.ptr);
            case 3: return doRandomMove!3(cast(ManifoldWrapper!3*) h.ptr);
            case 4: return doRandomMove!4(cast(ManifoldWrapper!4*) h.ptr);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_undo_pachner_move(void* handle) nothrow
{
    setError("undo_pachner_move requires tracking move history; not yet implemented");
    return -1;
}

// ---------------------------------------------------------------------------
// Targeted moves (worm program): apply a SPECIFIED bistellar or 4-4 hinge
// move.  All preconditions are checked explicitly here (doMove/doHingeMove
// rely on asserts, which are compiled out in release builds).
//
// WARNING: these act on the bare manifold.  If a ManifoldSampler with cocycle
// tracking wraps this manifold, targeted moves applied here bypass that
// tracking (a genuine detachment) -- use for analysis/catalog work, or
// re-enable the cocycle afterwards.
// ---------------------------------------------------------------------------

private string checkBistellar(int dim)(ManifoldWrapper!dim* mw,
    const(int)[] center, const(int)[] coCenter)
{
    import std.algorithm : canFind, sort;

    if (center.length < 1 || coCenter.length < 1
        || center.length + coCenter.length != dim + 2)
        return "center + coCenter must total dim+2 vertices, each nonempty";
    foreach (v; center)
        if (coCenter.canFind(v))
            return "center and coCenter share a vertex";
    auto cen = center.dup;
    cen.sort();
    auto coc = coCenter.dup;
    coc.sort();
    foreach (i; 1 .. cen.length)
        if (cen[i] == cen[i - 1]) return "repeated vertex in center";
    foreach (i; 1 .. coc.length)
        if (coc[i] == coc[i - 1]) return "repeated vertex in coCenter";
    if (!mw.mfd.contains(cen))
        return "center is not a simplex of the manifold";
    if (mw.mfd.contains(coc))
        return "coCenter is already a simplex of the manifold";
    if (mw.mfd.degree(cen) != coCenter.length)
        return "degree(center) != |coCenter|: star does not match this move";
    // every facet of star(center) must be center + (coCenter minus one vertex)
    foreach (skip; 0 .. coc.length)
    {
        int[] facet;
        facet ~= cen;
        foreach (i, v; coc)
            if (i != skip) facet ~= v;
        facet.sort();
        if (!mw.mfd.contains(facet))
            return "star of center does not match coCenter (facet missing)";
    }
    return null;
}

extern(C) int ddg_manifold_do_bistellar_move(void* handle,
    const(int)* center, int center_len,
    const(int)* cocenter, int cocenter_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        int doIt(int dim)(ManifoldWrapper!dim* mw)
        {
            auto cen = center[0 .. center_len];
            auto coc = cocenter[0 .. cocenter_len];
            auto err = checkBistellar!dim(mw, cen, coc);
            if (err !is null) { setError(err); return -1; }
            auto bm = BistellarMove!(dim, int)(cen, coc);
            mw.mfd.doMove(bm);
            return 0;
        }

        switch (h.dim)
        {
            case 2: return doIt!2(cast(ManifoldWrapper!2*) h.ptr);
            case 3: return doIt!3(cast(ManifoldWrapper!3*) h.ptr);
            case 4: return doIt!4(cast(ManifoldWrapper!4*) h.ptr);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_has_bistellar_move(void* handle,
    const(int)* center, int center_len,
    const(int)* cocenter, int cocenter_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        int checkIt(int dim)(ManifoldWrapper!dim* mw)
        {
            auto err = checkBistellar!dim(mw, center[0 .. center_len],
                                          cocenter[0 .. cocenter_len]);
            return err is null ? 1 : 0;
        }

        switch (h.dim)
        {
            case 2: return checkIt!2(cast(ManifoldWrapper!2*) h.ptr);
            case 3: return checkIt!3(cast(ManifoldWrapper!3*) h.ptr);
            case 4: return checkIt!4(cast(ManifoldWrapper!4*) h.ptr);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

private string checkHinge(ManifoldWrapper!3* mw,
    const(int)[] re, const(int)[] cyc, int diagonal, out HingeMove!int mv)
{
    import std.algorithm : canFind, sort;
    import std.algorithm.comparison : min, max;

    if (re.length != 2 || cyc.length != 4)
        return "hinge move needs a 2-vertex edge and a 4-vertex link cycle";
    if (diagonal != 0 && diagonal != 1)
        return "diagonal must be 0 (cycle[0]-cycle[2]) or 1 (cycle[1]-cycle[3])";
    foreach (i, v; cyc)
    {
        if (re.canFind(v)) return "link cycle vertex repeats an edge vertex";
        foreach (j; 0 .. i)
            if (cyc[j] == v) return "repeated vertex in link cycle";
    }
    if (re[0] == re[1]) return "degenerate edge";

    int[2] edge = [min(re[0], re[1]), max(re[0], re[1])];
    if (!mw.mfd.contains(edge[]))
        return "removed edge is not in the manifold";
    if (mw.mfd.degree(edge[]) != 4)
        return "removed edge does not have degree 4";
    foreach (i; 0 .. 4)
    {
        int[4] tet = [edge[0], edge[1], cast(int) cyc[i],
                      cast(int) cyc[(i + 1) % 4]];
        tet[].sort();
        if (!mw.mfd.contains(tet[]))
            return "link cycle does not match the star of the removed edge";
    }
    int a = diagonal == 0 ? cyc[0] : cyc[1];
    int b = diagonal == 0 ? cyc[2] : cyc[3];
    int[2] added = [min(a, b), max(a, b)];
    if (mw.mfd.contains(added[]))
        return "added edge (chosen diagonal) already exists";

    mv.removedEdge = edge;
    mv.addedEdge = added;
    foreach (i; 0 .. 4) mv.linkCycle[i] = cyc[i];
    return null;
}

extern(C) int ddg_manifold_do_hinge_move(void* handle,
    const(int)* removed_edge, const(int)* link_cycle, int diagonal) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("hinge moves are dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        HingeMove!int mv;
        auto err = checkHinge(mw, removed_edge[0 .. 2], link_cycle[0 .. 4],
                              diagonal, mv);
        if (err !is null) { setError(err); return -1; }
        mw.mfd.doHingeMove(mv);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_manifold_has_hinge_move(void* handle,
    const(int)* removed_edge, const(int)* link_cycle, int diagonal) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("hinge moves are dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        HingeMove!int mv;
        auto err = checkHinge(mw, removed_edge[0 .. 2], link_cycle[0 .. 4],
                              diagonal, mv);
        return err is null ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// SimplicialComplex lifecycle
// ---------------------------------------------------------------------------

extern(C) void* ddg_sc_create() nothrow
{
    clearError();
    try
    {
        auto w = new SimplicialComplex!int;
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_sc_from_facets(const(int)* data, int num_facets, int verts_per_facet) nothrow
{
    clearError();
    try
    {
        const(int)[][] facets;
        foreach (i; 0 .. num_facets)
            facets ~= data[i * verts_per_facet .. (i + 1) * verts_per_facet];

        auto w = new SimplicialComplex!int;
        *w = SimplicialComplex!int(facets);
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_sc_load(const(char)* path) nothrow
{
    clearError();
    try
    {
        string fileName = toStr(path);
        auto w = new SimplicialComplex!int;
        *w = loadSimplicialComplex!int(fileName);
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void* ddg_sc_copy(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return null; }
        auto src = cast(SimplicialComplex!int*) handle;
        auto w = new SimplicialComplex!int;
        *w = *src;
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void ddg_sc_free(void* handle) nothrow
{
    unpinForC(handle);
}

// ---------------------------------------------------------------------------
// SimplicialComplex queries
// ---------------------------------------------------------------------------

extern(C) long ddg_sc_num_facets(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return cast(long) sc.numFacets;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_f_vector(void* handle, long* out_buf, int buf_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto fv = (*sc).fVector;
        int n = cast(int) fv.length;
        if (n > buf_len) n = buf_len;
        foreach (i; 0 .. n) out_buf[i] = fv[i];
        return n;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_contains(void* handle, const(int)* data, int len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return sc.contains(data[0 .. len]) ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_contains_facet(void* handle, const(int)* data, int len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return sc.containsFacet(data[0 .. len]) ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// SimplicialComplex data — facets, simplices, star, link
// ---------------------------------------------------------------------------

extern(C) long ddg_sc_facets(void* handle, int* out_data, int* out_sizes) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto f = sc.facets;
        if (out_data !is null)
        {
            int idx = 0;
            int sIdx = 0;
            foreach (facet; f)
            {
                if (out_sizes !is null)
                    out_sizes[sIdx++] = cast(int) facet.walkLength;
                foreach (v; facet)
                    out_data[idx++] = v;
            }
        }
        return cast(long) f.walkLength;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sc_facets_of_dim(void* handle, int dim, int* out_data) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto f = sc.facets(dim);
        long count = 0;
        if (out_data !is null)
        {
            int idx = 0;
            foreach (facet; f)
            {
                foreach (v; facet)
                    out_data[idx++] = v;
                count++;
            }
        }
        else
        {
            count = f.walkLength;
        }
        return count;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sc_simplices(void* handle, int simp_dim, int* out_data) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto s = sc.simplices(simp_dim);
        if (out_data !is null)
        {
            int idx = 0;
            foreach (simp; s)
                foreach (v; simp)
                    out_data[idx++] = v;
        }
        return cast(long) s.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sc_star(void* handle, const(int)* simplex, int simplex_len,
    int* out_data, int* out_sizes) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto s = sc.star(simplex[0 .. simplex_len]).array;
        if (out_data !is null)
        {
            int idx = 0;
            int sIdx = 0;
            foreach (facet; s)
            {
                if (out_sizes !is null)
                    out_sizes[sIdx++] = cast(int) facet.walkLength;
                foreach (v; facet)
                    out_data[idx++] = v;
            }
        }
        return cast(long) s.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sc_link(void* handle, const(int)* simplex, int simplex_len,
    int* out_data, int* out_sizes) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto lk = sc.link(simplex[0 .. simplex_len]).map!array.array;
        if (out_data !is null)
        {
            int idx = 0;
            int sIdx = 0;
            foreach (simp; lk)
            {
                if (out_sizes !is null)
                    out_sizes[sIdx++] = cast(int) simp.length;
                foreach (v; simp)
                    out_data[idx++] = v;
            }
        }
        return cast(long) lk.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// SimplicialComplex mutation
// ---------------------------------------------------------------------------

extern(C) int ddg_sc_insert_facet(void* handle, const(int)* data, int len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        sc.insertFacet(data[0 .. len]);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_remove_facet(void* handle, const(int)* data, int len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        sc.removeFacet(data[0 .. len]);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// SimplicialComplex I/O
// ---------------------------------------------------------------------------

extern(C) int ddg_sc_save(void* handle, const(char)* path) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        string fileName = toStr(path);
        (*sc).saveTo(fileName);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_save_edge_graph(void* handle, const(char)* path) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        string fileName = toStr(path);
        (*sc).saveEdgeGraphTo(fileName);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// SimplicialComplex algorithms
// ---------------------------------------------------------------------------

extern(C) int ddg_sc_euler_characteristic(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return 0; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).eulerCharacteristic;
    }
    catch (Exception e) { setError(e.msg); return 0; }
}

extern(C) int ddg_sc_is_connected(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).isConnected ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_connected_components(void* handle, void** out_handles) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        auto components = (*sc).connectedComponents.array;
        if (out_handles !is null)
        {
            foreach (i, ref comp; components)
            {
                auto w = new SimplicialComplex!int;
                *w = comp;
                pinForC(cast(void*) w);
                out_handles[i] = cast(void*) w;
            }
        }
        return cast(int) components.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_pure(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).isPure ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_pure_of_dim(void* handle, int d) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).isPureOfDim(d) ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_circle(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).isCircle ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_2_sphere(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).is2Sphere ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_2_torus(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).is2Torus ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sc_is_orientable_surface_of_genus(void* handle, int g) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto sc = cast(SimplicialComplex!int*) handle;
        return (*sc).isOrientableSurfaceOfGenus(g) ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) void* ddg_sc_join(void* handle1, void* handle2) nothrow
{
    clearError();
    try
    {
        if (handle1 is null || handle2 is null) { setError("null handle"); return null; }
        auto sc1 = cast(SimplicialComplex!int*) handle1;
        auto sc2 = cast(SimplicialComplex!int*) handle2;
        auto result = join(*sc1, *sc2);
        auto w = new SimplicialComplex!int;
        *w = SimplicialComplex!int(result.facets);
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

// ---------------------------------------------------------------------------
// Sampler
// ---------------------------------------------------------------------------

alias CCallback = extern(C) int function(long, long, void*) nothrow;

private struct SamplerState
{
    int dim;
    void* manifoldHandle;
    int[] unusedVertices;

    int numFacetsTarget;
    double hingeDegreeTarget;
    double numFacetsCoef;
    double numHingesCoef;
    double hingeDegreeVarianceCoef;
    double coDim3DegreeVarianceCoef;
    // Fixed-target (strictly local) degree penalties; 0 = off (default).
    // Edge target reuses hingeDegreeTarget; the codim-3 target is its own knob.
    double hingeDegreeTargetCoef = 0.0;
    double coDim3DegreeTargetCoef = 0.0;
    double coDim3DegreeTarget = 0.0;
    double hingeMoveProb = 0.0;

    // MCMC state
    real currentObjective = real.nan;

    // Callback interval: how many moves between callback invocations
    long callbackInterval = 1000;

    // Statistics (cumulative across run() calls; reset with ddg_sampler_reset_stats)
    ulong hingeTries, hingeAccepts;
    ulong[5] bistellarTries;   // indexed by coCenter.length - 1; max dim=4 -> 5 slots
    ulong[5] bistellarAccepts;
    long totalAccepted, totalTried;

    // Per-vertex move-attribution counters (measured combinatorial lapse).
    // Opt-in (small per-proposal AA overhead); dim=3 only. See sampler.MoveCounters.
    bool trackMoveCounts = false;
    MoveCounters!int moveCounters;

    // Role-resolved geometry ledger + event log (opt-in; dim=3 only).
    // See sampler.GeometryLedger for the role taxonomy and record layout.
    GeometryLedger!int geoLedger;

    // Vertex 6-valence potential (Z-legality + chemical tilts + impurity
    // valence; opt-in, dim=3 only). See sampler.VertexPot.
    bool potEnabled = false;
    VertexPot vertexPot;
    VertexPotState!int vertexPotState;

    // Integer 1-cocycle tracking (T^3 winding forms; opt-in, dim=3 only).
    // See sampler.CocycleState.
    CocycleState!int cocycle;
}

extern(C) void* ddg_sampler_create(void* manifold_handle,
    int numFacetsTarget, double hingeDegreeTarget,
    double numFacetsCoef, double numHingesCoef,
    double hingeDegreeVarianceCoef, double coDim3DegreeVarianceCoef) nothrow
{
    return ddg_sampler_create_ext(manifold_handle, numFacetsTarget,
        hingeDegreeTarget, numFacetsCoef, numHingesCoef,
        hingeDegreeVarianceCoef, coDim3DegreeVarianceCoef, 0.0);
}

extern(C) void* ddg_sampler_create_ext(void* manifold_handle,
    int numFacetsTarget, double hingeDegreeTarget,
    double numFacetsCoef, double numHingesCoef,
    double hingeDegreeVarianceCoef, double coDim3DegreeVarianceCoef,
    double hingeMoveProb) nothrow
{
    clearError();
    try
    {
        if (manifold_handle is null) { setError("null manifold handle"); return null; }
        auto mh = cast(ManifoldHandle*) manifold_handle;

        auto state = new SamplerState;
        state.dim = mh.dim;
        state.manifoldHandle = ddg_manifold_copy(manifold_handle);
        if (state.manifoldHandle is null) return null;

        state.numFacetsTarget = numFacetsTarget;
        state.hingeDegreeTarget = hingeDegreeTarget;
        state.numFacetsCoef = numFacetsCoef;
        state.numHingesCoef = numHingesCoef;
        state.hingeDegreeVarianceCoef = hingeDegreeVarianceCoef;
        state.coDim3DegreeVarianceCoef = coDim3DegreeVarianceCoef;
        state.hingeMoveProb = hingeMoveProb;

        // Pre-allocate _facetArray so doMove doesn't trigger GC
        void reserveAndGetUnused(int dim)(SamplerState* s)
        {
            auto mw = cast(ManifoldWrapper!dim*)(cast(ManifoldHandle*) s.manifoldHandle).ptr;
            mw.mfd.reserveCapacity(s.numFacetsTarget);
            s.unusedVertices = getUnusedVertices!dim(mw.mfd);
        }

        switch (mh.dim)
        {
            case 2: reserveAndGetUnused!2(state); break;
            case 3: reserveAndGetUnused!3(state); break;
            case 4: reserveAndGetUnused!4(state); break;
            default: setError("bad dimension"); return null;
        }

        // Compute initial objective
        void initObjective(int dim)(SamplerState* s)
        {
            auto mw = cast(ManifoldWrapper!dim*)(cast(ManifoldHandle*) s.manifoldHandle).ptr;
            struct Params
            {
                int numFacetsTarget;
                real hingeDegreeTarget;
                real numFacetsCoef;
                real numHingesCoef;
                real hingeDegreeVarianceCoef;
                real coDim3DegreeVarianceCoef;
                real hingeDegreeTargetCoef;
                real coDim3DegreeTargetCoef;
                real coDim3DegreeTarget;
            }
            auto p = Params(s.numFacetsTarget, cast(real) s.hingeDegreeTarget,
                cast(real) s.numFacetsCoef, cast(real) s.numHingesCoef,
                cast(real) s.hingeDegreeVarianceCoef, cast(real) s.coDim3DegreeVarianceCoef,
                cast(real) s.hingeDegreeTargetCoef, cast(real) s.coDim3DegreeTargetCoef,
                cast(real) s.coDim3DegreeTarget);
            s.currentObjective = mw.mfd.objective(p);
        }

        switch (mh.dim)
        {
            case 2: initObjective!2(state); break;
            case 3: initObjective!3(state); break;
            case 4: initObjective!4(state); break;
            default: break;
        }

        pinForC(cast(void*) state);
        return cast(void*) state;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) long ddg_sampler_run(void* sampler_handle, long num_moves,
    CCallback callback, void* user_data) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;

        // dim=3 with hinge moves: use mcmcStep from sampler.d
        if (state.dim == 3)
            return runSamplerDim3(state, num_moves, callback, user_data);

        // dim=2 or dim=4: bistellar-only path
        long runBistellar(int dim)(SamplerState* s, long numMoves,
            CCallback cb, void* ud)
        {
            auto smh = cast(ManifoldHandle*) s.manifoldHandle;
            auto mw = cast(ManifoldWrapper!dim*) smh.ptr;

            struct Params
            {
                int numFacetsTarget;
                real hingeDegreeTarget;
                real numFacetsCoef;
                real numHingesCoef;
                real hingeDegreeVarianceCoef;
                real coDim3DegreeVarianceCoef;
                real hingeDegreeTargetCoef;
                real coDim3DegreeTargetCoef;
                real coDim3DegreeTarget;
            }

            auto params = Params(
                s.numFacetsTarget,
                cast(real) s.hingeDegreeTarget,
                cast(real) s.numFacetsCoef,
                cast(real) s.numHingesCoef,
                cast(real) s.hingeDegreeVarianceCoef,
                cast(real) s.coDim3DegreeVarianceCoef,
                cast(real) s.hingeDegreeTargetCoef,
                cast(real) s.coDim3DegreeTargetCoef,
                cast(real) s.coDim3DegreeTarget
            );

            if (s.currentObjective != s.currentObjective) // NaN check
                s.currentObjective = mw.mfd.objective(params);

            long accepted = 0;

            foreach (moveNum; 0 .. numMoves)
            {
                s.totalTried++;

                if (cb !is null && s.callbackInterval > 0 && moveNum % s.callbackInterval == 0 && moveNum > 0)
                {
                    if (cb(moveNum, numMoves, ud) != 0)
                        return accepted;
                }

                if (s.unusedVertices.length == 0)
                    s.unusedVertices ~= cast(int) mw.mfd.fVector[0];

                auto bm = mw.mfd.chooseRandomMove(s.unusedVertices[$ - 1], params);
                s.bistellarTries[bm.coCenter.length - 1]++;

                real deltaObj = mw.mfd.speculativeBistellarDelta(bm, s.currentObjective, params);
                real logAlpha = -deltaObj;

                if ((logAlpha >= 0) || (uniform01 <= exp(logAlpha)))
                {
                    mw.mfd.doMove(bm);
                    if (bm.coCenter.length == 1)
                    {
                        if (s.unusedVertices.length > 0)
                        {
                            s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                            s.unusedVertices.assumeSafeAppend;
                        }
                    }
                    if (bm.center.length == 1) s.unusedVertices ~= bm.center;
                    s.currentObjective += deltaObj;
                    s.bistellarAccepts[bm.coCenter.length - 1]++;
                    accepted++;
                    s.totalAccepted++;
                }
            }

            if (cb !is null)
                cb(numMoves, numMoves, ud);

            return accepted;
        }

        switch (state.dim)
        {
            case 2: return runBistellar!2(state, num_moves, callback, user_data);
            case 4: return runBistellar!4(state, num_moves, callback, user_data);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// dim=3 sampler using mcmcStep (supports hinge moves and full stats tracking).
private long runSamplerDim3(SamplerState* s, long numMoves,
    CCallback cb, void* ud) nothrow
{
    try
    {
        auto smh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) smh.ptr;

        struct Params
        {
            int numFacetsTarget;
            real hingeDegreeTarget;
            real numFacetsCoef;
            real numHingesCoef;
            real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef;
            real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef;
            real coDim3DegreeTarget;
        }

        auto params = Params(
            s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget,
            cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef,
            cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget
        );

        if (s.currentObjective != s.currentObjective) // NaN check
        {
            s.currentObjective = mw.mfd.objective(params);
            if (s.potEnabled)
            {
                mw.mfd.recomputeVertexPotState(s.vertexPotState, s.vertexPot);
                s.currentObjective += s.vertexPotState.total;
            }
        }

        long accepted = 0;
        long acceptedSinceWriteback = 0;
        // Narrow slices for mcmcStep ref params
        auto hT = s.hingeTries;
        auto hA = s.hingeAccepts;
        ulong[4] bT = s.bistellarTries[0 .. 4];
        ulong[4] bA = s.bistellarAccepts[0 .. 4];

        foreach (moveNum; 0 .. numMoves)
        {
            s.totalTried++;

            if (cb !is null && s.callbackInterval > 0 && moveNum % s.callbackInterval == 0 && moveNum > 0)
            {
                // Write back stats before callback
                s.hingeTries = hT;
                s.hingeAccepts = hA;
                s.bistellarTries[0 .. 4] = bT;
                s.bistellarAccepts[0 .. 4] = bA;
                s.totalAccepted += acceptedSinceWriteback;
                acceptedSinceWriteback = 0;

                if (cb(moveNum, numMoves, ud) != 0)
                    return accepted;
            }

            if (mw.mfd.mcmcStep(s.currentObjective, s.unusedVertices, params,
                    s.hingeMoveProb, hT, hA, bT, bA,
                    s.trackMoveCounts ? &s.moveCounters : null,
                    (s.geoLedger.trackRoles || s.geoLedger.logEvents
                        || s.geoLedger.logSixFlips)
                        ? &s.geoLedger : null,
                    s.potEnabled ? &s.vertexPotState : null,
                    s.potEnabled ? &s.vertexPot : null,
                    s.cocycle.enabled ? &s.cocycle : null))
            {
                accepted++;
                acceptedSinceWriteback++;
            }
        }

        // Write back stats
        s.hingeTries = hT;
        s.hingeAccepts = hA;
        s.bistellarTries[0 .. 4] = bT;
        s.bistellarAccepts[0 .. 4] = bA;
        s.totalAccepted += acceptedSinceWriteback;

        if (cb !is null)
            cb(numMoves, numMoves, ud);

        return accepted;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Run the sampler with exact Hastings correction using countValidBistellarMoves.
/// Execute-then-undo approach: slower per move but produces the exact target distribution.
extern(C) long ddg_sampler_run_exact(void* sampler_handle, long num_moves,
    CCallback callback, void* user_data) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (state.potEnabled)
        {
            setError("run_exact does not support the n6 potential");
            return -1;
        }
        if (state.cocycle.enabled)
        {
            setError("run_exact does not support cocycle tracking");
            return -1;
        }

        long runExact(int dim)(SamplerState* s, long numMoves,
            CCallback cb, void* ud)
        {
            auto smh = cast(ManifoldHandle*) s.manifoldHandle;
            auto mw = cast(ManifoldWrapper!dim*) smh.ptr;

            struct Params
            {
                int numFacetsTarget;
                real hingeDegreeTarget;
                real numFacetsCoef;
                real numHingesCoef;
                real hingeDegreeVarianceCoef;
                real coDim3DegreeVarianceCoef;
                real hingeDegreeTargetCoef;
                real coDim3DegreeTargetCoef;
                real coDim3DegreeTarget;
            }

            auto params = Params(
                s.numFacetsTarget,
                cast(real) s.hingeDegreeTarget,
                cast(real) s.numFacetsCoef,
                cast(real) s.numHingesCoef,
                cast(real) s.hingeDegreeVarianceCoef,
                cast(real) s.coDim3DegreeVarianceCoef,
                cast(real) s.hingeDegreeTargetCoef,
                cast(real) s.coDim3DegreeTargetCoef,
                cast(real) s.coDim3DegreeTarget
            );

            auto currentObjective = mw.mfd.objective(params);
            long accepted = 0;

            foreach (moveNum; 0 .. numMoves)
            {
                if (cb !is null && s.callbackInterval > 0 && moveNum % s.callbackInterval == 0 && moveNum > 0)
                {
                    if (cb(moveNum, numMoves, ud) != 0)
                        return accepted;
                }

                if (s.unusedVertices.length == 0)
                    s.unusedVertices ~= cast(int) mw.mfd.fVector[0];

                auto bm = mw.mfd.chooseRandomMove(s.unusedVertices[$ - 1], params);

                // Exact Hastings: execute, compute V_after, accept or undo
                immutable vBefore = cast(real) mw.mfd.countValidBistellarMoves;

                mw.mfd.doMove(bm);
                if (bm.coCenter.length == 1)
                {
                    if (s.unusedVertices.length > 0)
                    {
                        s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                        s.unusedVertices.assumeSafeAppend;
                    }
                }
                if (bm.center.length == 1) s.unusedVertices ~= bm.center;

                real newObjective = mw.mfd.objective(params);
                real deltaObj = newObjective - currentObjective;
                immutable vAfter = cast(real) mw.mfd.countValidBistellarMoves;

                real logAlpha = -deltaObj + log(vBefore) - log(vAfter);

                if ((logAlpha < 0) && (uniform01 > exp(logAlpha)))
                {
                    // Rejected — undo
                    mw.mfd.undoMove(bm);
                    if (bm.coCenter.length == 1) s.unusedVertices ~= bm.coCenter;
                    if (bm.center.length == 1)
                    {
                        assert(bm.center.front == s.unusedVertices[$ - 1]);
                        s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                        s.unusedVertices.assumeSafeAppend;
                    }
                }
                else
                {
                    currentObjective = newObjective;
                    accepted++;
                }
            }

            if (cb !is null)
                cb(numMoves, numMoves, ud);

            return accepted;
        }

        switch (state.dim)
        {
            case 2: return runExact!2(state, num_moves, callback, user_data);
            case 3: return runExact!3(state, num_moves, callback, user_data);
            case 4: return runExact!4(state, num_moves, callback, user_data);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) void* ddg_sampler_get_manifold(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return null; }
    return (cast(SamplerState*) sampler_handle).manifoldHandle;
}

extern(C) void ddg_sampler_free(void* handle) nothrow
{
    if (handle is null) return;
    auto state = cast(SamplerState*) handle;
    // Free the sampler's internal manifold copy
    ddg_manifold_free(state.manifoldHandle);
    unpinForC(handle);
}


// ---------------------------------------------------------------------------
// Sampler direct queries (avoid copying the manifold)
// ---------------------------------------------------------------------------

extern(C) int ddg_sampler_f_vector(void* sampler_handle, long* out_buf, int buf_len) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
        return ddg_manifold_f_vector(mfd_handle, out_buf, buf_len);
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) double ddg_sampler_importance_weight(void* sampler_handle) nothrow
{
    if (sampler_handle is null) return double.nan;
    auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
    return ddg_manifold_importance_weight(mfd_handle);
}

extern(C) long ddg_sampler_simplices(void* sampler_handle, int dim, int* out_data) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
        return ddg_manifold_simplices(mfd_handle, dim, out_data);
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sampler_degree(void* sampler_handle, const(int)* simplex, int len) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
        return ddg_manifold_degree(mfd_handle, simplex, len);
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sampler_set_num_facets_target(void* sampler_handle, int target) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.numFacetsTarget = target;
    state.currentObjective = real.nan; // force recompute on next run
    return 0;
}

extern(C) int ddg_sampler_set_callback_interval(void* sampler_handle, long interval) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    (cast(SamplerState*) sampler_handle).callbackInterval = interval;
    return 0;
}

extern(C) int ddg_sampler_set_hinge_move_prob(void* sampler_handle, double prob) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    (cast(SamplerState*) sampler_handle).hingeMoveProb = prob;
    return 0;
}

/// Enable/disable per-vertex move-attribution counters (dim=3 samplers only;
/// small per-proposal overhead, off by default). Enabling does not clear
/// previously accumulated counts; use ddg_sampler_reset_move_counts.
extern(C) int ddg_sampler_track_move_counts(void* sampler_handle, int enable) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    if (enable != 0 && state.dim != 3)
    {
        setError("move-count tracking is only supported for dim=3 samplers");
        return -1;
    }
    state.trackMoveCounts = enable != 0;
    return 0;
}

extern(C) int ddg_sampler_reset_move_counts(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    try
    {
        (cast(SamplerState*) sampler_handle).moveCounters.clear();
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Fetch the per-vertex move counters. Two-call pattern (like degree_histogram):
/// call with labels==null to get the entry count, then with all five buffers
/// (each of that length) to fill them. Entries are sorted by vertex label.
/// Every event contributes total weight 1 to its ledger, so ledger sums equal
/// event counts; a 1-4 move's created-vertex label appears like any other.
extern(C) long ddg_sampler_move_counts(void* sampler_handle, int* labels,
    double* proposed, double* valid,
    double* accepted_bistellar, double* accepted_hinge) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;

        bool[int] keySet;
        foreach (k; s.moveCounters.proposed.byKey) keySet[k] = true;
        foreach (k; s.moveCounters.valid.byKey) keySet[k] = true;
        foreach (k; s.moveCounters.acceptedBistellar.byKey) keySet[k] = true;
        foreach (k; s.moveCounters.acceptedHinge.byKey) keySet[k] = true;

        if (labels is null)
            return keySet.length;

        auto keys = keySet.keys;
        keys.sort();
        foreach (i, k; keys)
        {
            labels[i] = k;
            proposed[i] = s.moveCounters.proposed.get(k, 0.0);
            valid[i] = s.moveCounters.valid.get(k, 0.0);
            accepted_bistellar[i] = s.moveCounters.acceptedBistellar.get(k, 0.0);
            accepted_hinge[i] = s.moveCounters.acceptedHinge.get(k, 0.0);
        }
        return keys.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Enable/disable the role-resolved geometry ledger (dim=3 only; heavier than
/// move-count tracking). Does not clear accumulated data; see reset_geometry.
extern(C) int ddg_sampler_track_geometry(void* sampler_handle, int enable) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    if (enable != 0 && state.dim != 3)
    {
        setError("geometry tracking is only supported for dim=3 samplers");
        return -1;
    }
    state.geoLedger.trackRoles = enable != 0;
    return 0;
}

extern(C) int ddg_sampler_reset_geometry(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    try
    {
        (cast(SamplerState*) sampler_handle).geoLedger.clearRoles();
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Vertex role counts. Two-call pattern: labels==null returns the entry count;
/// otherwise fills labels[n] (sorted) and counts[n*11] (row-major, columns in
/// sampler.VRole order).
extern(C) long ddg_sampler_vertex_role_counts(void* sampler_handle,
    int* labels, double* counts) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
        bool[int] keySet;
        foreach (ref aa; g.vertexRoles)
            foreach (k; aa.byKey) keySet[k] = true;
        if (labels is null) return keySet.length;
        auto keys = keySet.keys;
        keys.sort();
        enum nc = VRole.max + 1;
        foreach (i, k; keys)
        {
            labels[i] = k;
            foreach (c; 0 .. nc)
                counts[i * nc + c] = g.vertexRoles[c].get(k, 0.0);
        }
        return keys.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Edge role counts. Two-call pattern: labels_a==null returns the entry count;
/// otherwise fills labels_a[n], labels_b[n] (sorted pairs, lexicographic) and
/// counts[n*15] (row-major, columns in sampler.ERole order).
extern(C) long ddg_sampler_edge_role_counts(void* sampler_handle,
    int* labels_a, int* labels_b, double* counts) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
        bool[int[2]] keySet;
        foreach (ref aa; g.edgeRoles)
            foreach (k; aa.byKey) keySet[k] = true;
        if (labels_a is null) return keySet.length;
        auto keys = keySet.keys;
        keys.sort();
        enum nc = ERole.max + 1;
        foreach (i, k; keys)
        {
            labels_a[i] = k[0];
            labels_b[i] = k[1];
            foreach (c; 0 .. nc)
                counts[i * nc + c] = g.edgeRoles[c].get(k, 0.0);
        }
        return keys.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Tet aggregates: created[5] / destroyed[5] by move type code, the log2
/// lifetime histogram [64] (age in attempted moves), the number of currently
/// living tracked tets, censored deaths, and the ledger clock.
extern(C) int ddg_sampler_tet_stats(void* sampler_handle,
    long* created, long* destroyed, long* lifetime_hist,
    long* living, long* censored, long* clock) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
        foreach (i; 0 .. 5)
        {
            created[i] = cast(long) g.tetsCreated[i];
            destroyed[i] = cast(long) g.tetsDestroyed[i];
        }
        foreach (i; 0 .. 64)
            lifetime_hist[i] = cast(long) g.tetLifetimeHist[i];
        *living = cast(long) g.tetBirth.length;
        *censored = cast(long) g.tetCensoredDeaths;
        *clock = cast(long) g.clock;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Enable the move event log with the given buffer capacity in bytes
/// (rounded down to whole records); capacity 0 disables. Enabling resets the
/// buffer and the overflow flag.
extern(C) int ddg_sampler_event_log_enable(void* sampler_handle,
    long capacity_bytes) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (capacity_bytes > 0 && state.dim != 3)
        {
            setError("event log is only supported for dim=3 samplers");
            return -1;
        }
        auto g = &state.geoLedger;
        if (capacity_bytes <= 0)
        {
            g.logEvents = false;
            g.eventBuf = null;
            g.eventUsed = 0;
            g.eventOverflow = false;
            return 0;
        }
        immutable cap = (cast(size_t) capacity_bytes / eventRecordBytes)
                        * eventRecordBytes;
        if (cap == 0) { setError("capacity smaller than one record"); return -1; }
        g.eventBuf = new ubyte[cap];
        g.eventUsed = 0;
        g.eventOverflow = false;
        g.logEvents = true;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Drain the event log. buf==null returns the number of buffered bytes.
/// Otherwise copies up to buf_len bytes (whole records), removes them from the
/// buffer, and returns the byte count copied. Check _event_log_overflowed to
/// detect records dropped between drains.
extern(C) long ddg_sampler_event_log_drain(void* sampler_handle,
    ubyte* buf, long buf_len) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
        if (buf is null) return g.eventUsed;
        import core.stdc.string : memcpy, memmove;
        immutable take = (min(cast(size_t) buf_len, g.eventUsed)
                          / eventRecordBytes) * eventRecordBytes;
        if (take > 0)
        {
            memcpy(buf, g.eventBuf.ptr, take);
            if (take < g.eventUsed)
                memmove(g.eventBuf.ptr, g.eventBuf.ptr + take, g.eventUsed - take);
            g.eventUsed -= take;
        }
        return take;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Returns 1 if records were dropped since the log was enabled/drained-clear,
/// else 0. Clears the flag.
extern(C) int ddg_sampler_event_log_overflowed(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
    immutable r = g.eventOverflow ? 1 : 0;
    g.eventOverflow = false;
    return r;
}

// ---------------------------------------------------------------------------
// Six-edge flip log (disclination-network rewiring stream; dim=3 only)
// ---------------------------------------------------------------------------

/// Enable the six-edge flip log: one fixed-size record (see
/// sampler.sixRecordBytes: u64 clock + i32 u + i32 v + i32 dir) per edge that
/// crosses the degree 5<->6 threshold during an accepted move. The stream is
/// the complete rewiring history of the disclination network. capacity <= 0
/// disables. Drain regularly; on overflow records are dropped (flagged).
extern(C) int ddg_sampler_six_flip_log_enable(void* sampler_handle,
    long capacity_bytes) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (capacity_bytes > 0 && state.dim != 3)
        {
            setError("six-flip log is only supported for dim=3 samplers");
            return -1;
        }
        auto g = &state.geoLedger;
        if (capacity_bytes <= 0)
        {
            g.logSixFlips = false;
            g.sixBuf = null;
            g.sixUsed = 0;
            g.sixOverflow = false;
            return 0;
        }
        immutable cap = (cast(size_t) capacity_bytes / sixRecordBytes)
                        * sixRecordBytes;
        if (cap == 0) { setError("capacity smaller than one record"); return -1; }
        g.sixBuf = new ubyte[cap];
        g.sixUsed = 0;
        g.sixOverflow = false;
        g.logSixFlips = true;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Drain the six-flip log. buf==null returns the number of buffered bytes.
/// Otherwise copies up to buf_len bytes (whole records), removes them from
/// the buffer, and returns the byte count copied.
extern(C) long ddg_sampler_six_flip_log_drain(void* sampler_handle,
    ubyte* buf, long buf_len) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null sampler handle"); return -1; }
        auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
        if (buf is null) return g.sixUsed;
        import core.stdc.string : memcpy, memmove;
        immutable take = (min(cast(size_t) buf_len, g.sixUsed)
                          / sixRecordBytes) * sixRecordBytes;
        if (take > 0)
        {
            memcpy(buf, g.sixBuf.ptr, take);
            if (take < g.sixUsed)
                memmove(g.sixBuf.ptr, g.sixBuf.ptr + take, g.sixUsed - take);
            g.sixUsed -= take;
        }
        return take;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Returns 1 if six-flip records were dropped, else 0. Clears the flag.
extern(C) int ddg_sampler_six_flip_log_overflowed(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto g = &(cast(SamplerState*) sampler_handle).geoLedger;
    immutable r = g.sixOverflow ? 1 : 0;
    g.sixOverflow = false;
    return r;
}

// ---------------------------------------------------------------------------
// Integer 1-cocycle tracking (T^3 winding forms; dim=3 only)
// ---------------------------------------------------------------------------

/// Enable cocycle tracking from an initial assignment: edges = 2*n_edges ints
/// (u, v pairs), values = 3*n_edges ints (omega(u->v) per edge). The edge set
/// must exactly cover the sampler's current manifold and the cochain must be
/// closed on every triangle — both are verified here (error + disabled state
/// on failure). n_edges <= 0 disables tracking and frees the state.
extern(C) int ddg_sampler_cocycle_enable(void* sampler_handle,
    const(int)* edges, const(int)* values, long n_edges) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (n_edges <= 0)
        {
            state.cocycle.enabled = false;
            state.cocycle.omega = null;
            return 0;
        }
        if (state.dim != 3)
        {
            setError("cocycle tracking is only supported for dim=3 samplers");
            return -1;
        }
        if (edges is null || values is null)
        {
            setError("null edges/values");
            return -1;
        }
        state.cocycle.omega = null;
        foreach (i; 0 .. cast(size_t) n_edges)
        {
            immutable u = edges[2 * i], v = edges[2 * i + 1];
            if (u == v) { setError("degenerate edge"); return -1; }
            int[2] key = u < v ? [u, v] : [v, u];
            immutable s = u < v ? 1 : -1;
            int[3] val = [s * values[3 * i], s * values[3 * i + 1],
                          s * values[3 * i + 2]];
            if (key in state.cocycle.omega)
            {
                state.cocycle.omega = null;
                setError("duplicate edge in cocycle init");
                return -1;
            }
            state.cocycle.omega[key] = val;
        }
        state.cocycle.enabled = true;
        auto smh = cast(ManifoldHandle*) state.manifoldHandle;
        auto prob = cocycleProblems(
            (cast(ManifoldWrapper!3*) smh.ptr).mfd, state.cocycle);
        if (prob !is null)
        {
            state.cocycle.enabled = false;
            state.cocycle.omega = null;
            setError(prob);
            return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Read the current cocycle. With null buffers, returns the edge count.
/// Otherwise fills edges_out (2 ints/edge, sorted u < v) and values_out
/// (3 ints/edge, omega(u->v)) for up to cap_edges edges and returns the
/// number written (order unspecified).
extern(C) long ddg_sampler_cocycle_read(void* sampler_handle,
    int* edges_out, int* values_out, long cap_edges) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (!state.cocycle.enabled) { setError("cocycle not enabled"); return -1; }
        if (edges_out is null || values_out is null)
            return cast(long) state.cocycle.omega.length;
        long i = 0;
        foreach (key, val; state.cocycle.omega)
        {
            if (i >= cap_edges) break;
            edges_out[2 * i] = key[0];
            edges_out[2 * i + 1] = key[1];
            values_out[3 * i] = val[0];
            values_out[3 * i + 1] = val[1];
            values_out[3 * i + 2] = val[2];
            ++i;
        }
        return i;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Audit the tracked cocycle against the current manifold (edge-set match +
/// closedness on every triangle). Returns 0 if clean, -1 with error set if
/// not — the production drift check.
extern(C) int ddg_sampler_cocycle_check(void* sampler_handle) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (!state.cocycle.enabled) { setError("cocycle not enabled"); return -1; }
        auto smh = cast(ManifoldHandle*) state.manifoldHandle;
        auto prob = cocycleProblems(
            (cast(ManifoldWrapper!3*) smh.ptr).mfd, state.cocycle);
        if (prob !is null) { setError(prob); return -1; }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sampler_set_num_facets_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.numFacetsCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_num_hinges_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.numHingesCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_hinge_degree_variance_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.hingeDegreeVarianceCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_codim3_degree_variance_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.coDim3DegreeVarianceCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_hinge_degree_target(void* sampler_handle, double target) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.hingeDegreeTarget = target;
    state.currentObjective = real.nan;
    return 0;
}

/// Fixed-target (strictly local) degree penalties. Coefficients default to 0
/// (off). The hinge term targets hingeDegreeTarget; the codim-3 term targets
/// coDim3DegreeTarget (set it consistently with the pinned f-vector, e.g. in
/// dim 3: t = 4/(6/hingeDegreeTarget - 1)).
extern(C) int ddg_sampler_set_hinge_degree_target_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.hingeDegreeTargetCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_codim3_degree_target_coef(void* sampler_handle, double coef) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.coDim3DegreeTargetCoef = coef;
    state.currentObjective = real.nan;
    return 0;
}

extern(C) int ddg_sampler_set_codim3_degree_target(void* sampler_handle, double target) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto state = cast(SamplerState*) sampler_handle;
    state.coDim3DegreeTarget = target;
    state.currentObjective = real.nan;
    return 0;
}

/// Reset cumulative statistics counters.
extern(C) int ddg_sampler_reset_stats(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    s.hingeTries = 0;
    s.hingeAccepts = 0;
    s.bistellarTries[] = 0;
    s.bistellarAccepts[] = 0;
    s.totalAccepted = 0;
    s.totalTried = 0;
    return 0;
}

/// Get sampler statistics. All output pointers are optional (may be null).
extern(C) int ddg_sampler_get_stats(void* sampler_handle,
    long* out_total_tried, long* out_total_accepted,
    long* out_hinge_tries, long* out_hinge_accepts,
    long* out_bistellar_tries, long* out_bistellar_accepts,
    int out_len) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;

    if (out_total_tried !is null) *out_total_tried = s.totalTried;
    if (out_total_accepted !is null) *out_total_accepted = s.totalAccepted;
    if (out_hinge_tries !is null) *out_hinge_tries = cast(long) s.hingeTries;
    if (out_hinge_accepts !is null) *out_hinge_accepts = cast(long) s.hingeAccepts;
    if (out_bistellar_tries !is null)
        foreach (i; 0 .. out_len)
            out_bistellar_tries[i] = cast(long) s.bistellarTries[i];
    if (out_bistellar_accepts !is null)
        foreach (i; 0 .. out_len)
            out_bistellar_accepts[i] = cast(long) s.bistellarAccepts[i];
    return 0;
}

extern(C) double ddg_sampler_current_objective(void* sampler_handle) nothrow
{
    if (sampler_handle is null) return double.nan;
    return cast(double)(cast(SamplerState*) sampler_handle).currentObjective;
}

/// Configure the vertex 6-valence potential (Z-legality + chemical tilts +
/// impurity valence; see sampler.VertexPot). dim=3 samplers only. tilt5 may be
/// null (all-zero tilts). Passing all-zero coefficients disables the term.
extern(C) int ddg_sampler_set_n6_potential(void* sampler_handle,
    double zleg_coef, double imp_coef, const(double)* tilt5) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (state.dim != 3)
        {
            setError("n6 potential requires a dim=3 sampler");
            return -1;
        }
        state.vertexPot.zlegCoef = cast(real) zleg_coef;
        state.vertexPot.impCoef = cast(real) imp_coef;
        if (tilt5 !is null)
            foreach (i; 0 .. 5)
                state.vertexPot.tilt[i] = cast(real) tilt5[i];
        else
            state.vertexPot.tilt[] = 0;
        state.potEnabled = state.vertexPot.enabled;
        state.currentObjective = real.nan;  // recompute (incl. state) next run
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// Degree variance
// ---------------------------------------------------------------------------

extern(C) double ddg_manifold_degree_variance(void* handle, int simp_dim) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return double.nan; }
        auto h = cast(ManifoldHandle*) handle;
        switch (h.dim)
        {
            case 2: return cast(double)(cast(ManifoldWrapper!2*) h.ptr).mfd.degreeVariance(simp_dim);
            case 3: return cast(double)(cast(ManifoldWrapper!3*) h.ptr).mfd.degreeVariance(simp_dim);
            case 4: return cast(double)(cast(ManifoldWrapper!4*) h.ptr).mfd.degreeVariance(simp_dim);
            default: setError("bad dimension"); return double.nan;
        }
    }
    catch (Exception e) { setError(e.msg); return double.nan; }
}

extern(C) double ddg_sampler_degree_variance(void* sampler_handle, int simp_dim) nothrow
{
    if (sampler_handle is null) return double.nan;
    auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
    return ddg_manifold_degree_variance(mfd_handle, simp_dim);
}

// ---------------------------------------------------------------------------
// Degree histogram
// ---------------------------------------------------------------------------

/// Returns degree histogram for simplices of given dimension.
/// If out_buf is null, returns the histogram length.
/// Otherwise fills out_buf with counts: out_buf[i] = number of simplices with degree i+1.
extern(C) long ddg_manifold_degree_histogram(void* handle, int simp_dim, long* out_buf) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        size_t[] getHist(int dim)()
        {
            return (cast(ManifoldWrapper!dim*) h.ptr).mfd.degreeHistogram(simp_dim);
        }

        size_t[] hist;
        switch (h.dim)
        {
            case 2: hist = getHist!2(); break;
            case 3: hist = getHist!3(); break;
            case 4: hist = getHist!4(); break;
            default: setError("bad dimension"); return -1;
        }

        if (out_buf is null)
            return cast(long) hist.length;

        foreach (i, v; hist)
            out_buf[i] = cast(long) v;

        return cast(long) hist.length;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sampler_degree_histogram(void* sampler_handle, int simp_dim, long* out_buf) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
    return ddg_manifold_degree_histogram(mfd_handle, simp_dim, out_buf);
}

// ---------------------------------------------------------------------------
// Disclination-network censuses (dim=3 only)
// ---------------------------------------------------------------------------

/// Joint (n6, m) vertex census (see sampler.valenceCensus). Fills out_buf
/// with (n6_cap+1)*(m_cap+1) longs, row-major [min(n6, n6_cap)][min(m, m_cap)].
extern(C) int ddg_manifold_valence_census(void* handle, long* out_buf,
    int n6_cap, int m_cap) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        if (out_buf is null) { setError("null out_buf"); return -1; }
        if (n6_cap < 1 || m_cap < 1) { setError("caps must be >= 1"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3)
        {
            setError("valence census requires a dim=3 manifold");
            return -1;
        }
        auto len = (n6_cap + 1) * (m_cap + 1);
        (cast(ManifoldWrapper!3*) h.ptr).mfd.valenceCensus(
            out_buf[0 .. len], n6_cap, m_cap);
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) int ddg_sampler_valence_census(void* sampler_handle, long* out_buf,
    int n6_cap, int m_cap) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
    return ddg_manifold_valence_census(mfd_handle, out_buf, n6_cap, m_cap);
}

/// Disclination-network census (see sampler.DisclinationCensus). host_mask
/// bit k marks n6-class k as a native host class (C15: (1<<0)|(1<<4)); 0 for
/// no host/dopant split. With out_buf null, returns the required slot count
/// (sampler.disclinationCensusSlots); otherwise fills the flattened layout
/// documented at sampler.flattenCensus and returns the slot count.
extern(C) long ddg_manifold_disclination_census(void* handle, int host_mask,
    long* out_buf, long buf_len) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        if (out_buf is null) return disclinationCensusSlots;
        if (buf_len < disclinationCensusSlots)
        {
            setError("buffer too small for census layout");
            return -1;
        }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3)
        {
            setError("disclination census requires a dim=3 manifold");
            return -1;
        }
        auto c = (cast(ManifoldWrapper!3*) h.ptr).mfd.disclinationCensus(host_mask);
        flattenCensus(c, out_buf[0 .. disclinationCensusSlots]);
        return disclinationCensusSlots;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_sampler_disclination_census(void* sampler_handle,
    int host_mask, long* out_buf, long buf_len) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto mfd_handle = (cast(SamplerState*) sampler_handle).manifoldHandle;
    return ddg_manifold_disclination_census(mfd_handle, host_mask, out_buf, buf_len);
}
