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
// Runtime initialization
// ---------------------------------------------------------------------------

/// Initialize the D runtime. The Python loader MUST call this once right
/// after dlopen, before any other ddg_* call.
///
/// A dlopen'd D shared library does not get rt_init() the way a D executable
/// does from its main() -- and without it the calling thread is never
/// registered with the runtime, so the stop-the-world collector silently
/// bails: GC.collect() frees NOTHING, ever, and every GC allocation grows
/// the heap permanently. This stayed invisible for years because the hot
/// paths are allocation-free by design; the worm channel's enumeration
/// garbage turned it into a multi-GB-per-hour leak (2026-07-30).
///
/// Returns 1 on success, 0 on failure. Safe to call more than once
/// (druntime refcounts rt_init). There is deliberately no matching
/// terminate: tearing down druntime under a live Python interpreter at
/// exit is riskier than letting the OS reclaim the process.
extern(C) int ddg_runtime_init() nothrow
{
    import core.runtime : Runtime;
    try { return Runtime.initialize() ? 1 : 0; }
    catch (Throwable) { return 0; }
}

// ---------------------------------------------------------------------------
// GC control
// ---------------------------------------------------------------------------

/// Trigger a D garbage collection cycle. Call this periodically from Python
/// to reclaim temporaries created by chooseRandomMove, coCenter, etc.
/// Only effective after ddg_runtime_init (see above).
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
    catch (Throwable e) { setError(e.msg); return -1; }
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
            import std.algorithm : sort;
            auto cen = center[0 .. center_len];
            auto coc = cocenter[0 .. cocenter_len];
            auto err = checkBistellar!dim(mw, cen, coc);
            if (err !is null) { setError(err); return -1; }
            // CANONICALIZE: the core's simplex invariant is ascending vertex
            // order (assertValidSimplex), and doMove's productUnion is a
            // sorted MERGE -- unsorted center/coCenter here would insert
            // phantom out-of-order keys into the dimension maps (the
            // worm_slide corruption bug, fixed 2026-07-25).
            auto cenS = cen.dup;
            cenS.sort();
            auto cocS = coc.dup;
            cocS.sort();
            auto bm = BistellarMove!(dim, int)(cenS, cocS);
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

/// Edges whose degree is NOT in {5, 6} (the FK-illegal set), dim=3 only.
/// Returns the count; if out_pairs/out_degs are non-null they receive 2 ints
/// per edge (sorted pair) and the degree.  Reads the manifold's own
/// incremental degree map -- no recomputation.
extern(C) long ddg_manifold_illegal_edges(void* handle,
    int* out_pairs, int* out_degs) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("illegal_edges is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        // One allocation-free scan of the edge map (degrees read in place);
        // callers size buffers with a first null-buffer call, so pass an
        // unbounded cap on the write pass.
        return mw.mfd.writeIllegalEdges(out_pairs, out_degs,
            out_pairs is null ? 0 : long.max);
    }
    catch (Throwable e) { setError(e.msg); return -1; }
}

/// Link of an edge (dim=3): one UNORDERED vertex pair per tetrahedron in the
/// edge's star, written flat (2 ints per pair).  Returns the pair count
/// (= edge degree), or -1 with an error if the edge is absent.
extern(C) long ddg_manifold_edge_link(void* handle, int a, int b,
    int* out_pairs) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("edge_link is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        import std.algorithm.comparison : max, min;
        int[2] eb = [min(a, b), max(a, b)];
        if (!mw.mfd.contains(eb[]))
        {
            setError("edge not in manifold");
            return -1;
        }
        long n = 0;
        foreach (pr; mw.mfd.link(eb[]))
        {
            if (out_pairs !is null)
            {
                int i = 0;
                foreach (v; pr) out_pairs[2 * n + i++] = v;
            }
            n++;
        }
        return n;
    }
    catch (Throwable e) { setError(e.msg); return -1; }
}

/// Ordered link cycle of edge (a,b) via the O(degree) ridge walk, dim=3.
/// hint_tet: 4 vertices of a CURRENT facet containing a and b (validated).
/// Writes the cycle vertices to out_cycle, returns the count (= degree).
extern(C) long ddg_manifold_edge_link_cycle(void* handle, int a, int b,
    const(int)* hint_tet, int* out_cycle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("edge_link_cycle is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        import std.algorithm : sort;
        int[4] tet = hint_tet[0 .. 4];
        tet[].sort();
        bool hasA, hasB;
        foreach (v; tet)
        {
            hasA |= v == a;
            hasB |= v == b;
        }
        if (!hasA || !hasB) { setError("hint tet lacks the edge"); return -1; }
        if (!mw.mfd.contains(tet[]))
        {
            setError("hint tet is not a current facet");
            return -1;
        }
        return mw.mfd.writeEdgeLinkCycle(a, b, tet[], out_cycle);
    }
    catch (Throwable e) { setError(e.msg); return -1; }
}

/// The two apexes of interior triangle (a,b,c): O(1) ridge-link lookup.
extern(C) long ddg_manifold_face_apexes(void* handle, int a, int b, int c,
    int* out_apexes) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("face_apexes is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        auto n = mw.mfd.writeFaceApexes(a, b, c, out_apexes);
        if (n < 0) { setError("triangle not in manifold"); return -1; }
        return n;
    }
    catch (Throwable e) { setError(e.msg); return -1; }
}

/// Census of slide-frame outcomes over all degree-3 chords x SLIDE_SLOTS
/// slots (dim = 3 only). Writes 8 longs to out8:
///   [0] ok   [1] degenerate-frame   [2] missing-face     (all chords)
///   [3] number of degree-3 chords   [4] number of degree-4 edges
///   [5] ok   [6] degenerate         [7] missing   (chords touching a deg-4 edge)
/// A `degenerate` outcome is a slide blocked by a folded apex-walk -- i.e. by a
/// nearby degree-4 disclination edge (a defect complex), NOT by a broken window
/// (`missing`). The [5..7] split isolates that near the complexes. O(1) per
/// frame; read-only. Returns 0 on success, -1 on error.
extern(C) long ddg_slide_frame_census(void* handle, long* out8) nothrow
{
    clearError();
    if (handle is null) { setError("null handle"); return -1; }
    if (out8 is null) { setError("null out buffer"); return -1; }
    auto h = cast(ManifoldHandle*) handle;
    if (h.dim != 3) { setError("slide_frame_census is dim=3 only"); return -1; }
    try
    {
        import std.algorithm : sort, map;
        import std.array : array;
        auto mfd = &(cast(ManifoldWrapper!3*) h.ptr).mfd;
        foreach (i; 0 .. 8) out8[i] = 0;
        static immutable int[2][6] picks =
            [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];

        // vertices incident to a degree-4 edge (a defect disclination line)
        bool[int] deg4vert;
        foreach (e; mfd.simplices(1).map!(x => [x[0], x[1]]).array)
            if (mfd.degree(e) == 4)
            {
                out8[4]++;
                deg4vert[e[0]] = true;
                deg4vert[e[1]] = true;
            }

        foreach (edge; mfd.simplices(1).map!(x => [x[0], x[1]]).array)
        {
            if (mfd.degree(edge) != 3) continue;
            out8[3]++;
            immutable bool near =
                (edge[0] in deg4vert) !is null || (edge[1] in deg4vert) !is null;

            // hint tet: any current facet containing the edge
            int[4] hint = 0;
            bool got = false;
            foreach (f; mfd.facets)
            {
                bool a0 = false, a1 = false;
                foreach (v; f)
                {
                    if (v == edge[0]) a0 = true;
                    if (v == edge[1]) a1 = true;
                }
                if (a0 && a1)
                {
                    int j = 0;
                    foreach (v; f) hint[j++] = v;
                    got = true;
                    break;
                }
            }
            if (!got) continue;

            int[8] lb = 0;
            if (mfd.writeEdgeLinkCycle(edge[0], edge[1], hint[], lb.ptr) != 3)
                continue;
            int[3] lk = [lb[0], lb[1], lb[2]];
            lk[].sort();

            foreach (slot; 0 .. SLIDE_SLOTS)
            {
                immutable int orient = slot / 6;
                immutable int pick = slot % 6;
                immutable int c0 = orient == 0 ? edge[0] : edge[1];
                immutable int c4 = orient == 0 ? edge[1] : edge[0];
                immutable int c2 = lk[picks[pick][0]];
                immutable int c3 = lk[picks[pick][1]];
                SlideFrame!int fr;
                immutable cause =
                    deriveSlideFrameCause(*mfd, c0, c4, c2, c3, fr);
                final switch (cause)
                {
                case SlideFrameCause.ok:
                    out8[0]++; if (near) out8[5]++; break;
                case SlideFrameCause.degenerate:
                    out8[1]++; if (near) out8[6]++; break;
                case SlideFrameCause.missingFace:
                    out8[2]++; if (near) out8[7]++; break;
                }
            }
        }
        return 0;
    }
    catch (Throwable ex) { setError(ex.msg); return -1; }
}

/// Audit the manifold's internal hash maps; returns 0 healthy, -1 with an
/// error message describing the first inconsistency.
extern(C) int ddg_manifold_validate_maps(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        string err;
        switch (h.dim)
        {
            case 2: err = (cast(ManifoldWrapper!2*) h.ptr).mfd.validateMaps; break;
            case 3: err = (cast(ManifoldWrapper!3*) h.ptr).mfd.validateMaps; break;
            case 4: err = (cast(ManifoldWrapper!4*) h.ptr).mfd.validateMaps; break;
            default: setError("bad dimension"); return -1;
        }
        if (err !is null) { setError(err); return -1; }
        return 0;
    }
    catch (Throwable e) { setError(e.msg); return -1; }
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

    // Knot-slide move class (dim=3 only): probability of proposing a slide
    // once the unified proposal lands on a degree-3 edge, plus its counters.
    SlideConfig slideCfg;

    // Non-local slide channel (dim=3 only): picks a degree-3 chord uniformly
    // from `deg3Chords` (kept live only while nlSlideCfg.prob > 0), plus its
    // config/counters. See sampler.NonlocalSlideConfig / tryNonlocalSlide.
    NonlocalSlideConfig nlSlideCfg;
    Deg3Set!int deg3Chords;

    // Deg-4 worm channel (dim=3 only): catalysed 2-move transport step, anchor
    // uniform over `deg4Edges` (kept live only while wormCfg.prob > 0).  See
    // sampler.tryWormMove.  Not cocycle-safe (gated off in mcmcStep).
    WormConfig wormCfg;
    Deg3Set!int deg4Edges;

    // Contract/split channel (dim=3 only): edge contraction + vertex split,
    // the f0-changing block-move pair. See sampler.ContractSplitConfig /
    // tryContractSplit. Not cocycle-safe (gated off in mcmcStep).
    ContractSplitConfig csCfg;

    // f0 worm channel (scheme C; dim=3 only): extended-ensemble vertex
    // removal/insertion with a frozen umbrella table. See sampler.wormF0Episode.
    // Not cocycle- or ledger-safe (gated in the episode entry point).
    WormF0Params wormF0;
    bool wormF0On = false;
    WF0Applied!int[] wormF0Undo;

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

    // Hard cap on the illegal-edge count (dim=3 only; opt-in). See
    // sampler.IllegalBudget. Requires the vertex potential (it maintains the
    // sum_v m counter the gate reads) and is incompatible with the channels
    // mcmcStep runs before the gated ones.
    IllegalBudget illegalBudget;
}

/// Recompute the tracked objective from scratch, including the n6 potential
/// state for dim=3. Parameter setters invalidate currentObjective to NaN;
/// this is the single re-derivation used at create, at the top of run(), and
/// lazily by ddg_sampler_current_objective.
private void recomputeObjective(SamplerState* s)
{
    void doIt(int dim)(SamplerState* s)
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
        static if (dim == 3)
        {
            if (s.potEnabled)
            {
                mw.mfd.recomputeVertexPotState(s.vertexPotState, s.vertexPot);
                s.currentObjective += s.vertexPotState.total;
            }
        }
    }
    switch (s.dim)
    {
        case 2: doIt!2(s); break;
        case 3: doIt!3(s); break;
        case 4: doIt!4(s); break;
        default: break;
    }
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

        recomputeObjective(state);

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
                recomputeObjective(s);

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
            recomputeObjective(s);

        // Keep the degree-3 chord set live only while the non-local channel is
        // on; rebuild it fresh at run start (O(E), once per run).
        if (s.nlSlideCfg.prob > 0)
            rebuildDeg3(mw.mfd, s.deg3Chords);
        if (s.wormCfg.prob > 0)
            rebuildDegSet(mw.mfd, s.deg4Edges, 4);

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
                    s.cocycle.enabled ? &s.cocycle : null,
                    (s.dim == 3 && s.slideCfg.prob > 0) ? &s.slideCfg : null,
                    (s.dim == 3 && s.nlSlideCfg.prob > 0) ? &s.nlSlideCfg : null,
                    (s.dim == 3 && s.nlSlideCfg.prob > 0) ? &s.deg3Chords : null,
                    (s.dim == 3 && s.wormCfg.prob > 0) ? &s.wormCfg : null,
                    (s.dim == 3 && s.wormCfg.prob > 0) ? &s.deg4Edges : null,
                    (s.dim == 3 && s.csCfg.prob > 0) ? &s.csCfg : null,
                    (s.dim == 3 && s.illegalBudget.enabled)
                        ? &s.illegalBudget : null))
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

/******************************************************************************
Enable the incrementally-maintained vertex lift (dim = 3, cocycle required).

Integrates the cocycle over a spanning forest ONCE to seed a position per
vertex, then keeps it exact move-by-move at O(1) -- see CocycleState.pos. This
replaces the pattern of re-deriving positions from the whole cochain on every
sample, which costs O(V+E) per call and, because it rebuilds the spanning tree
each time, can reassign a persisting vertex (the "gauge glitch" that consumers
otherwise have to detect and discard). Here the gauge is fixed for the run.

`box` is the torus period per axis, in the same units as the cocycle values.
*/
extern(C) int ddg_sampler_cocycle_enable_positions(void* sampler_handle,
    const(int)* box) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (!state.cocycle.enabled) { setError("cocycle not enabled"); return -1; }
        if (box is null) { setError("null box"); return -1; }
        int[3] b = [box[0], box[1], box[2]];
        auto prob = cocycleSeedPositions(state.cocycle, b);
        if (prob !is null)
        {
            state.cocycle.posEnabled = false;
            state.cocycle.pos = null;
            setError(prob);
            return -1;
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Positions of the given vertices, 3 ints each, written to out_pos in the
/// order requested. Returns 0, or -1 with the error set if any vertex has no
/// position (i.e. is not in the current triangulation).
extern(C) int ddg_sampler_cocycle_positions(void* sampler_handle,
    const(int)* verts, long n_verts, int* out_pos) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        if (!state.cocycle.posEnabled)
        { setError("cocycle positions not enabled"); return -1; }
        if (verts is null || out_pos is null)
        { setError("null buffer"); return -1; }
        foreach (i; 0 .. cast(size_t) n_verts)
        {
            auto p = verts[i] in state.cocycle.pos;
            if (p is null)
            {
                import std.conv : to;
                setError("vertex " ~ verts[i].to!string ~ " has no position");
                return -1;
            }
            out_pos[3 * i] = (*p)[0];
            out_pos[3 * i + 1] = (*p)[1];
            out_pos[3 * i + 2] = (*p)[2];
        }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Audit the lift: pos(v) - pos(u) == omega(u->v) mod the torus period, on
/// every edge. Returns 0 if clean, -1 with the error set otherwise.
extern(C) int ddg_sampler_cocycle_pos_check(void* sampler_handle) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto state = cast(SamplerState*) sampler_handle;
        auto prob = cocyclePosProblems(state.cocycle);
        if (prob !is null) { setError(prob); return -1; }
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/******************************************************************************
The subcomplex INDUCED by a vertex set: every simplex of the manifold all of
whose vertices lie in `verts`. Returns a NEW SimplicialComplex handle (free it
with ddg_sc_free), or null on error.

This is what a defect IS -- the subcomplex spanned by its own vertices -- as
opposed to ddg_manifold_closed_star, which is the region it OCCUPIES. Uses the
transient vertex->facet incidence, so only facets meeting `verts` are examined.
*/
extern(C) void* ddg_manifold_induced_subcomplex(void* handle,
    const(int)* verts, long n_verts) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return null; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3)
        { setError("induced_subcomplex is dim=3 only"); return null; }
        if (verts is null && n_verts > 0) { setError("null verts"); return null; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        auto vs = verts[0 .. cast(size_t) n_verts].dup;
        auto w = new SimplicialComplex!int;
        *w = mw.mfd.inducedSubcomplex(vs);
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

/******************************************************************************
The CLOSED STAR of a vertex set: every facet incident to any of `verts`, as a
new SimplicialComplex handle (free with ddg_sc_free), or null on error. The
region a defect occupies; its boundary is the surface flux would cross.
*/
extern(C) void* ddg_manifold_closed_star(void* handle,
    const(int)* verts, long n_verts) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return null; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3)
        { setError("closed_star is dim=3 only"); return null; }
        if (verts is null && n_verts > 0) { setError("null verts"); return null; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        auto vs = verts[0 .. cast(size_t) n_verts].dup;
        auto w = new SimplicialComplex!int;
        *w = mw.mfd.closedStar(vs);
        pinForC(cast(void*) w);
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
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
    s.slideCfg.tries = 0;
    s.slideCfg.accepts = 0;
    return 0;
}

/******************************************************************************
Enable the knot-slide move class (dim = 3 only).

`prob` is the probability of proposing a slide rather than the ordinary 3->2
bistellar move, given that the unified proposal has landed on a degree-3 edge.
0 (the default) disables slides entirely. Every valid slide is in the class
and acceptance is plain Metropolis on the exact action change; restrict to
CLEAN (species-preserving) slides with ddg_sampler_set_slide_clean_only.
See sampler.trySlideMove.
*/
extern(C) int ddg_sampler_set_slide_clean_only(void* sampler_handle,
    int clean_only) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (clean_only == 0 && s.slideCfg.prob > 0 && s.illegalBudget.enabled)
    {
        setError("cannot drop the clean-slide restriction while an "
                 ~ "illegality budget is active and the slide channel is on");
        return -1;
    }
    s.slideCfg.cleanOnly = (clean_only != 0);
    return 0;
}

extern(C) int ddg_sampler_set_slide_prob(void* sampler_handle, double prob) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (!(prob >= 0.0 && prob <= 1.0))
    { setError("slide probability must be in [0, 1]"); return -1; }
    if (s.dim != 3 && prob > 0)
    { setError("knot slides are dim=3 only"); return -1; }
    // A CLEAN slide preserves the multiset of illegal degrees over its changed
    // edges exactly (sampler.trySlideMove), hence the global illegal-edge
    // count, so it cannot breach the cap -- and it is precisely the transport
    // move a budgeted reservoir needs: birth/death alone only ever retraces.
    // A dirty slide carries no such guarantee.
    if (prob > 0 && s.illegalBudget.enabled && !s.slideCfg.cleanOnly)
    {
        setError("with an illegality budget the slide channel must be "
                 ~ "restricted to CLEAN slides (they preserve the illegal "
                 ~ "count exactly); call ddg_sampler_set_slide_clean_only(1)");
        return -1;
    }
    s.slideCfg.prob = cast(real) prob;
    return 0;
}

/******************************************************************************
Enable the non-local slide channel (dim = 3 only): each mcmcStep proposes it
with probability `prob`, picking a degree-3 chord uniformly (1/n_3) and moving
it up to `max_step` tets along its BC chain (annihilate + re-create, with the
1/n_3 Hastings factor). max_step = 4 recovers the local slide's displacement.
Set prob = 0 (default) to disable. See sampler.tryNonlocalSlide.
*/
extern(C) int ddg_sampler_set_nonlocal_slide_prob(void* sampler_handle,
    double prob, int max_step) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (!(prob >= 0.0 && prob <= 1.0))
    { setError("nonlocal slide probability must be in [0, 1]"); return -1; }
    if (s.dim != 3 && prob > 0)
    { setError("nonlocal slides are dim=3 only"); return -1; }
    if (prob > 0 && max_step < 1)
    { setError("nonlocal slide max_step must be >= 1"); return -1; }
    if (prob > 0 && s.illegalBudget.enabled)
    { setError("this channel bypasses the illegality budget; disable the budget first"); return -1; }
    s.nlSlideCfg.prob = cast(real) prob;
    if (max_step >= 1) s.nlSlideCfg.maxStep = max_step;
    return 0;
}

/// Non-local slide counters: proposals that formed a legal move (reached
/// Metropolis), and how many were accepted. Both pointers optional.
extern(C) int ddg_sampler_nonlocal_slide_stats(void* sampler_handle,
    long* out_tries, long* out_accepts) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (out_tries !is null) *out_tries = cast(long) s.nlSlideCfg.tries;
    if (out_accepts !is null) *out_accepts = cast(long) s.nlSlideCfg.accepts;
    return 0;
}

/// Size of the live degree-3 chord set (the 1/n_3 denominator). Only populated
/// while the non-local channel is enabled + a run has started; a diagnostic for
/// checking the incremental maintenance against a full degree-3 recount.
extern(C) long ddg_sampler_deg3_count(void* sampler_handle) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    return cast(long) s.deg3Chords.length;
}

/******************************************************************************
Enable the contract/split channel (dim = 3 only): each mcmcStep proposes it
with probability `prob`, choosing edge contraction or vertex split with a
fair coin. max_ring caps deg(uv) on the contract side and the splitting
cycle length |gamma| on the split side (the pair must be capped together for
detailed balance); pass max_ring <= 0 to keep the current value (default 6).
Cocycle-safe: the winding cochain and vertex lift are maintained exactly
through both directions. The channel is automatically inert while six-flip
logging is on. Set prob = 0 (default) to disable. See
sampler.ContractSplitConfig / tryContractSplit.
*/
extern(C) int ddg_sampler_set_contract_split(void* sampler_handle,
    double prob, int max_ring) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (!(prob >= 0.0 && prob <= 1.0))
    { setError("contract/split probability must be in [0, 1]"); return -1; }
    if (s.dim != 3 && prob > 0)
    { setError("the contract/split channel is dim=3 only"); return -1; }
    if (max_ring > 0 && (max_ring < 3 || max_ring > 8))
    { setError("contract/split max_ring must be in [3, 8]"); return -1; }
    s.csCfg.prob = cast(real) prob;
    if (max_ring > 0) s.csCfg.maxRing = max_ring;
    return 0;
}

/// Contract/split counters: per-direction proposals that reached Metropolis,
/// accepts, and proposals that failed validity/geometry gates. All pointers
/// optional.
/// (hits, misses) of the link-cycle memo (link_cycles.countCyclesCached).
/// Process-wide diagnostic for sizing the memo; not per-sampler.
extern(C) int ddg_cycle_memo_stats(ulong* hits, ulong* misses) nothrow
{
    clearError();
    import link_cycles : cycleMemoStats;
    try
    {
        immutable st = cycleMemoStats();
        if (hits !is null) *hits = st[0];
        if (misses !is null) *misses = st[1];
    }
    catch (Exception) { setError("cycle memo stats failed"); return -1; }
    return 0;
}

extern(C) int ddg_sampler_contract_split_stats(void* sampler_handle,
    long* out_contract_tries, long* out_contract_accepts,
    long* out_split_tries, long* out_split_accepts,
    long* out_no_valid) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (out_contract_tries !is null)
        *out_contract_tries = cast(long) s.csCfg.contractTries;
    if (out_contract_accepts !is null)
        *out_contract_accepts = cast(long) s.csCfg.contractAccepts;
    if (out_split_tries !is null)
        *out_split_tries = cast(long) s.csCfg.splitTries;
    if (out_split_accepts !is null)
        *out_split_accepts = cast(long) s.csCfg.splitAccepts;
    if (out_no_valid !is null)
        *out_no_valid = cast(long) s.csCfg.noValid;
    return 0;
}

/******************************************************************************
Cap the illegal-edge count at `cap` (dim = 3 only): any move whose post-move
count of edges with degree outside {5, 6} would exceed the cap is rejected
outright.  Pass a negative cap to disable (the default).

This is a BUDGET, not a price.  A fugacity on illegal edges suppresses them by
making them expensive, but they are the only mobile objects, so a large
fugacity freezes the chain before it reaches zero.  A cap instead keeps a fixed
small reservoir circulating, so the chain returns to n_ill = 0 repeatedly and
those returns are exact samples of the uniform measure on FK-legal states (the
action is exactly degenerate there at fixed f0, f3).

Requires the n6 potential to be enabled -- its per-vertex state maintains the
sum_v m counter the gate reads (sum_v m = 2 * #illegal edges exactly).  It is
also refused while the slide, non-local slide or deg-4 worm channels are on:
those run BEFORE the gated move types in mcmcStep and would tunnel through the
cap.  Enable them in either order; whichever comes second is refused.
*/
extern(C) int ddg_sampler_set_illegal_budget(void* sampler_handle,
    long cap) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (cap >= 0)
    {
        if (s.dim != 3)
        { setError("the illegality budget is dim=3 only"); return -1; }
        if (!s.potEnabled)
        {
            setError("the illegality budget needs the n6 potential enabled "
                     ~ "(it maintains the counter the gate reads); call "
                     ~ "ddg_sampler_set_n6_potential first");
            return -1;
        }
        if ((s.slideCfg.prob > 0 && !s.slideCfg.cleanOnly)
            || s.nlSlideCfg.prob > 0 || s.wormCfg.prob > 0)
        {
            setError("the illegality budget cannot be combined with dirty "
                     ~ "slides or with the non-local-slide / worm channels: "
                     ~ "they bypass the gate. CLEAN slides are allowed (they "
                     ~ "preserve the illegal count exactly).");
            return -1;
        }
    }
    s.illegalBudget.cap = cap;
    s.illegalBudget.blocked = 0;
    return 0;
}

/// Illegality-budget readout: (cap, current illegal-edge count, #proposals
/// rejected by the cap alone). The count is the maintained sum_v m / 2 when
/// the potential is on, else -1.
extern(C) int ddg_sampler_illegal_budget_stats(void* sampler_handle,
    long* out_cap, long* out_n_illegal, long* out_blocked) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (out_cap !is null) *out_cap = s.illegalBudget.cap;
    if (out_n_illegal !is null)
        *out_n_illegal = s.potEnabled ? s.vertexPotState.sumM / 2 : -1;
    if (out_blocked !is null)
        *out_blocked = cast(long) s.illegalBudget.blocked;
    return 0;
}

/******************************************************************************
Enable the deg-4 worm channel (dim = 3 only): each mcmcStep proposes it with
probability `prob`, drawing an anchor uniformly from the live deg-4 edge set
and attempting the catalysed 2-move transport step (one 2->3 + one 3->2,
content-neutral, escaped landing) with anchor-sum Hastings weights.  The
channel is automatically inert while a cocycle is attached.  Set prob = 0
(default) to disable.  See sampler.tryWormMove.
*/
extern(C) int ddg_sampler_set_worm_prob(void* sampler_handle,
    double prob) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (!(prob >= 0.0 && prob <= 1.0))
    { setError("worm probability must be in [0, 1]"); return -1; }
    if (s.dim != 3 && prob > 0)
    { setError("the worm channel is dim=3 only"); return -1; }
    if (prob > 0 && s.illegalBudget.enabled)
    { setError("this channel bypasses the illegality budget; disable the budget first"); return -1; }
    s.wormCfg.prob = cast(real) prob;
    return 0;
}

/// Worm counters: proposals with >= 1 candidate (reached Metropolis),
/// accepts, and proposals rejected for lack of candidates.
extern(C) int ddg_sampler_worm_stats(void* sampler_handle,
    long* out_tries, long* out_accepts, long* out_nocands) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (out_tries !is null) *out_tries = cast(long) s.wormCfg.tries;
    if (out_accepts !is null) *out_accepts = cast(long) s.wormCfg.accepts;
    if (out_nocands !is null) *out_nocands = cast(long) s.wormCfg.noCands;
    return 0;
}

/******************************************************************************
Crossval / scripted entry into the worm move at anchor (a, b) (dim = 3).

mode 0 = enumerate only: returns the candidate count; if the out arrays are
non-null (capacity `cap`), writes per-candidate landing edges and exact dS
values (dS computed via a trial application of that candidate).
mode 1 = trial of candidate `cand`: dS written to out_ds[0]; state restored.
mode 2 = commit candidate `cand` unconditionally (SlideAccept.force).
Returns the candidate count (mode 0) or 0/1 = rejected/committed, -1 error.
*/
extern(C) long ddg_sampler_worm_at(void* sampler_handle, int a, int b,
    int mode, int cand, int* out_la, int* out_lb, double* out_ds,
    long cap) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("worm is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*)(
            cast(ManifoldHandle*) s.manifoldHandle).ptr;
        import std.algorithm.comparison : max, min;
        int[2] anchor = [min(a, b), max(a, b)];
        if (!mw.mfd.contains(anchor[]))
        { setError("edge not in manifold"); return -1; }
        int[4] hint = -1;
        if (!findTetOfEdge(mw.mfd, anchor[0], anchor[1], hint))
        { setError("no facet on the edge"); return -1; }
        hint[].sort();
        if (s.currentObjective != s.currentObjective)   // NaN
            recomputeObjective(s);

        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);

        if (mode == 0)
        {
            auto cands = new WormCand!int[](WORM_MAX_CANDS);
            immutable long n = mw.mfd.wormEnumerate(anchor, hint[], cands);
            foreach (i; 0 .. n)
            {
                if (out_la !is null && i < cap) out_la[i] = cands[i].landing[0];
                if (out_lb !is null && i < cap) out_lb[i] = cands[i].landing[1];
                if (out_ds !is null && i < cap)
                {
                    bool v = false;
                    real ds = real.nan;
                    mw.mfd.tryWormMove(s.currentObjective, anchor, hint[],
                        params, v,
                        s.potEnabled ? &s.vertexPotState : null,
                        s.potEnabled ? &s.vertexPot : null,
                        null, null, SlideAccept.trialOnly, &ds, cast(int) i);
                    out_ds[i] = v ? cast(double) ds : double.nan;
                }
            }
            return n;
        }
        if (mode != 1 && mode != 2)
        { setError("mode must be 0 (enum), 1 (trial) or 2 (commit)"); return -1; }
        bool valid = false;
        real ds = real.nan;
        immutable ok = mw.mfd.tryWormMove(s.currentObjective, anchor, hint[],
            params, valid,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            (s.wormCfg.prob > 0) ? &s.deg4Edges : null,
            (s.nlSlideCfg.prob > 0) ? &s.deg3Chords : null,
            mode == 1 ? SlideAccept.trialOnly : SlideAccept.force,
            &ds, cand, null, null,
            (s.geoLedger.trackRoles || s.geoLedger.logEvents)
                ? &s.geoLedger : null);
        if (!valid) { setError("candidate invalid"); return 0; }
        if (out_ds !is null) out_ds[0] = cast(double) ds;
        return ok ? 1 : 0;
    }
    catch (Exception e)
    {
        setError(e.msg);
        return -1;
    }
}

/******************************************************************************
Attempt the knot slide at chord (a, b) in slot `slot` directly (dim = 3).

This is the scripted / crossval entry point into the same code path the
sampler's slide move uses -- not a sampling path. `mode` selects what happens
once a legal clean slide is in hand:

    0  trial only: measure dS and cleanliness, restore the state exactly
    1  force: always commit (objective, potential, cocycle and ledger all
       updated as for an accepted move)

Writes the exact action change to *out_dS when the slot yields a legal clean
slide. Returns 1 if it did (and, for mode 1, was committed), 0 if the slot is
not a clean slide of this state, or -1 on error.
*/
extern(C) int ddg_sampler_slide_at(void* sampler_handle,
    int a, int b, int slot, int mode, double* out_dS) nothrow
{
    return ddg_sampler_slide_at2(sampler_handle, a, b, slot, mode, out_dS,
                                 null, null);
}

/// slide_at with the slot's ARRIVAL chord (c4, c8) reported from the frame
/// decode (valid whenever the frame derives, even if the slide itself is
/// rejected later). Lets drivers identify the exact inverse of a slide
/// without guessing (FPKMC walk-back).
extern(C) int ddg_sampler_slide_at2(void* sampler_handle,
    int a, int b, int slot, int mode, double* out_dS,
    int* out_c4, int* out_c8) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("knot slides are dim=3 only"); return -1; }
        if (slot < 0 || slot >= SLIDE_SLOTS)
        { setError("slot out of range"); return -1; }
        if (mode != 0 && mode != 1)
        { setError("mode must be 0 (trial) or 1 (force)"); return -1; }

        auto mw = cast(ManifoldWrapper!3*)(
            cast(ManifoldHandle*) s.manifoldHandle).ptr;
        import std.algorithm.comparison : max, min;
        int[2] eb = [min(a, b), max(a, b)];
        if (!mw.mfd.contains(eb[])) { setError("edge not in manifold"); return -1; }

        // any facet on the edge serves as the link-walk hint
        int[4] hint = 0;
        {
            bool got = false;
            foreach (pr; mw.mfd.link(eb[]))
            {
                int i = 0;
                foreach (v; pr) hint[2 + i++] = v;
                got = true;
                break;
            }
            if (!got) { setError("edge has empty link"); return -1; }
            hint[0] = eb[0]; hint[1] = eb[1];
            hint[].sort();
        }

        if (s.currentObjective != s.currentObjective) // NaN check
            recomputeObjective(s);

        if (out_c4 !is null || out_c8 !is null)
        {
            SlideFrame!int fr;
            int[3] lk;
            if (slideDecode(mw.mfd, a, b, hint[], slot, fr, lk))
            {
                if (out_c4 !is null) *out_c4 = fr.c4;
                if (out_c8 !is null) *out_c8 = fr.c8;
            }
            else
            {
                if (out_c4 !is null) *out_c4 = -1;
                if (out_c8 !is null) *out_c8 = -1;
            }
        }

        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);

        bool valid = false;
        real dS = real.nan;
        mw.mfd.trySlideMove(s.currentObjective, eb[0], eb[1], hint[], slot,
            params, valid, s.trackMoveCounts ? &s.moveCounters : null,
            (s.geoLedger.trackRoles || s.geoLedger.logEvents
                || s.geoLedger.logSixFlips) ? &s.geoLedger : null,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            s.cocycle.enabled ? &s.cocycle : null,
            mode == 0 ? SlideAccept.trialOnly : SlideAccept.force,
            &dS, s.slideCfg.cleanOnly);
        if (!valid) return 0;
        if (out_dS !is null) *out_dS = cast(double) dS;
        return 1;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Slide-move counters: proposals that formed a legal clean slide, and how
/// many of those were accepted. Both pointers optional.
extern(C) int ddg_sampler_slide_stats(void* sampler_handle,
    long* out_tries, long* out_accepts) nothrow
{
    clearError();
    if (sampler_handle is null) { setError("null handle"); return -1; }
    auto s = cast(SamplerState*) sampler_handle;
    if (out_tries !is null) *out_tries = cast(long) s.slideCfg.tries;
    if (out_accepts !is null) *out_accepts = cast(long) s.slideCfg.accepts;
    return 0;
}

/******************************************************************************
Non-local slide (dim = 3): annihilate the degree-3 chord (a,b) with a 3->2 and
re-create it `steps` tets down the BC chain with a 2->3. `slot` (0..11) selects
the walk (orientation x link order). `mode`: 0 = trial (measure, restore),
1 = force (commit). Writes the exact action change to *out_dS, the exact change
in the degree-3 edge count to *out_dn3 (both optional), and the arrival chord to
*out_ta,*out_tb. Returns 1 if the slot yielded a legal move (and, for mode 1,
committed), 0 if not a legal move, -1 on error. The 1/n_3 Hastings factor for a
Metropolis channel is applied by the caller using out_dn3. See
sampler.tryNonlocalSlide.
*/
extern(C) int ddg_sampler_nonlocal_slide_at(void* sampler_handle,
    int a, int b, int slot, int steps, int mode, double* out_dS,
    long* out_dn3, int* out_ta, int* out_tb) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("nonlocal slide is dim=3 only"); return -1; }
        if (mode != 0 && mode != 1)
        { setError("mode must be 0 (trial) or 1 (force)"); return -1; }

        auto mw = cast(ManifoldWrapper!3*)(
            cast(ManifoldHandle*) s.manifoldHandle).ptr;
        import std.algorithm.comparison : max, min;
        int[2] eb = [min(a, b), max(a, b)];
        if (!mw.mfd.contains(eb[])) { setError("edge not in manifold"); return -1; }

        // any facet on the edge serves as the link-walk hint
        int[4] hint = 0;
        {
            bool got = false;
            foreach (pr; mw.mfd.link(eb[]))
            {
                int i = 0;
                foreach (v; pr) hint[2 + i++] = v;
                got = true;
                break;
            }
            if (!got) { setError("edge has empty link"); return -1; }
            hint[0] = eb[0]; hint[1] = eb[1];
            hint[].sort();
        }

        if (s.currentObjective != s.currentObjective)   // NaN
            recomputeObjective(s);

        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);

        bool valid = false;
        real dS = real.nan;
        long dn3 = 0;
        int ta = -1, tb = -1;
        mw.mfd.tryNonlocalSlide(s.currentObjective, eb[0], eb[1], hint[],
            slot, steps, params, valid,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            mode == 0 ? SlideAccept.trialOnly : SlideAccept.force, &dS,
            -1, &dn3, &ta, &tb);
        if (!valid) return 0;
        if (out_dS !is null) *out_dS = cast(double) dS;
        if (out_dn3 !is null) *out_dn3 = dn3;
        if (out_ta !is null) *out_ta = ta;
        if (out_tb !is null) *out_tb = tb;
        return 1;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

// ---------------------------------------------------------------------------
// FPKMC infrastructure (notes/FPKMC_DESIGN.md, M1)
// ---------------------------------------------------------------------------

/// Sliding-window BC-chain walk from `window` (4 vertices of a facet, in
/// walk order), n steps. Writes n+4 vertices; window k is
/// out_verts[k .. k+4]. Returns the number of vertices written, or -1 on
/// error (window not a facet -- the walk itself cannot break on a closed
/// 3-manifold).
extern(C) long ddg_manifold_chain_walk(void* handle, const(int)* window,
    long n, int* out_verts) nothrow
{
    clearError();
    try
    {
        if (handle is null || window is null || out_verts is null)
        { setError("null argument"); return -1; }
        auto h = cast(ManifoldHandle*) handle;
        if (h.dim != 3) { setError("chain_walk is dim=3 only"); return -1; }
        auto mw = cast(ManifoldWrapper!3*) h.ptr;
        int[4] w = [window[0], window[1], window[2], window[3]];
        int[4] chk = w;
        chk[].sort();
        if (!mw.mfd.contains(chk[]))
        { setError("window is not a facet"); return -1; }
        foreach (i; 0 .. 4) out_verts[i] = w[i];
        long len = 4;
        foreach (k; 0 .. n)
        {
            if (!mw.mfd.nextChainWindow(w))
            { setError("chain walk broke (window desynced?)"); return -1; }
            out_verts[len++] = w[3];
        }
        return len;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Stride of one site_survey output row: dS_create + SLIDE_SLOTS x
/// (dS, dest c4, dest c8, clean).
enum int SURVEY_STRIDE = 1 + 4 * SLIDE_SLOTS;

/******************************************************************************
Washboard site survey along a BC chain (notes/FPKMC_DESIGN.md).

For each window k (0 <= k <= chain_len-5): CREATE the knot by the 2->3 on
the window face -- through the same speculative-delta path every sampler
move uses, potential state updated -- measure all SLIDE_SLOTS slide slots
in trial-only mode, then undo the creation. Manifold, potential state and
sampler objective are restored exactly (audited; -1 on drift).

Output, per window (SURVEY_STRIDE doubles):
  [0]            dS_create, or NaN if the window is not a valid creation
                 site (face missing, apexes disagree with the chain, move
                 invalid, frozen vertex in support);
  [1+4s .. ]     slot s: (dS trial, dest chord c4, dest chord c8, clean);
                 dS NaN if the slot does not form a legal slide; clean is
                 1.0 / 0.0 (species-preserving or not; NaN when illegal);
                 destinations from the frame derivation whenever it
                 succeeds. Comparing destinations against chain arithmetic
                 classifies slots as chain / off-chain -- the slot census
                 (design invariant I1 and R3(iii)). MEASURED on R m4: every
                 site has exactly 1 clean chain-forward + 1 clean
                 chain-backward slot (I1 exact) PLUS 2-3 clean OFF-CHAIN
                 slots and ~7 dirty ones: the clean slide network is a
                 branching graph, not a union of chains.

Returns the number of windows surveyed, or -1 on error.
*/
extern(C) long ddg_sampler_site_survey(void* sampler_handle,
    const(int)* chain, long chain_len, double* out_data) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null || chain is null || out_data is null)
        { setError("null argument"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("site_survey is dim=3 only"); return -1; }
        if (chain_len < 5) { setError("chain too short"); return -1; }
        auto mw = cast(ManifoldWrapper!3*)(
            cast(ManifoldHandle*) s.manifoldHandle).ptr;

        if (s.currentObjective != s.currentObjective)   // NaN
            recomputeObjective(s);

        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef;
            real hingeDegreeVarianceCoef; real coDim3DegreeVarianceCoef;
            real hingeDegreeTargetCoef; real coDim3DegreeTargetCoef;
            real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);

        alias BM = BistellarMove!(3, int);
        immutable long nwin = chain_len - 4;
        immutable real objBefore = s.currentObjective;
        immutable real potBefore =
            s.potEnabled ? s.vertexPotState.total : 0.0L;
        static immutable int[2][6] picks =
            [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];

        foreach (k; 0 .. nwin)
        {
            double* row = out_data + k * SURVEY_STRIDE;
            foreach (i; 0 .. SURVEY_STRIDE) row[i] = double.nan;

            int[3] face = [chain[k + 1], chain[k + 2], chain[k + 3]];
            face[].sort();
            int[2] apx = [chain[k], chain[k + 4]];
            apx[].sort();
            {
                int[2] ap = 0;
                if (mw.mfd.writeFaceApexes(face[0], face[1], face[2],
                                           ap.ptr) != 2)
                    continue;
                if (!((ap[0] == apx[0] && ap[1] == apx[1])
                      || (ap[0] == apx[1] && ap[1] == apx[0])))
                    continue;
            }
            auto bm = BM(face[], apx[]);
            if (!mw.mfd.hasValidMove(bm)) continue;
            {
                int[5] sup = [face[0], face[1], face[2], apx[0], apx[1]];
                if (mw.mfd.anyFrozen(sup[])) continue;
            }

            immutable real potNow =
                s.potEnabled ? s.vertexPotState.total : 0.0L;
            immutable real baseRun = s.currentObjective - potNow;
            immutable real dBase =
                mw.mfd.speculativeBistellarDelta(bm, baseRun, params);
            real dPot = 0.0L;
            if (s.potEnabled)
                dPot = mw.mfd.potentialBistellarDelta(bm, s.vertexPotState,
                    s.vertexPot, true);
            mw.mfd.doMove(bm);
            immutable real dSc = dBase + dPot;
            row[0] = cast(double) dSc;

            // the chord's 3 post-creation tets are {chord, face_i, face_j}
            int[4] hint = [apx[0], apx[1], face[0], face[1]];
            hint[].sort();
            real obj = s.currentObjective + dSc;
            foreach (slot; 0 .. SLIDE_SLOTS)
            {
                immutable int orient = slot / 6;
                immutable int pick = slot % 6;
                immutable int c0 = orient == 0 ? apx[0] : apx[1];
                immutable int c4 = orient == 0 ? apx[1] : apx[0];
                immutable int c2 = face[picks[pick][0]];
                immutable int c3 = face[picks[pick][1]];
                SlideFrame!int fr;
                if (deriveSlideFrame(mw.mfd, c0, c4, c2, c3, fr))
                {
                    row[1 + 4 * slot + 1] = cast(double) fr.c4;
                    row[1 + 4 * slot + 2] = cast(double) fr.c8;
                }
                bool valid = false;
                real dS = real.nan;
                mw.mfd.trySlideMove(obj, apx[0], apx[1], hint[], slot,
                    params, valid, null, null,
                    s.potEnabled ? &s.vertexPotState : null,
                    s.potEnabled ? &s.vertexPot : null,
                    null, SlideAccept.trialOnly, &dS, false);
                if (valid)
                {
                    row[1 + 4 * slot] = cast(double) dS;
                    // cleanliness: re-trial under cleanOnly -- valid iff
                    // the slide preserves the illegal-degree multiset
                    bool vClean = false;
                    real dS2 = real.nan;
                    mw.mfd.trySlideMove(obj, apx[0], apx[1], hint[], slot,
                        params, vClean, null, null,
                        s.potEnabled ? &s.vertexPotState : null,
                        s.potEnabled ? &s.vertexPot : null,
                        null, SlideAccept.trialOnly, &dS2, true);
                    row[1 + 4 * slot + 3] = vClean ? 1.0 : 0.0;
                }
            }

            auto inv = BM(apx[], face[]);
            if (!mw.mfd.hasValidMove(inv))
            { setError("site_survey: creation undo invalid (bug)"); return -1; }
            if (s.potEnabled)
                mw.mfd.potentialBistellarDelta(inv, s.vertexPotState,
                    s.vertexPot, true);
            mw.mfd.doMove(inv);
        }

        // restoration audit
        if (s.potEnabled)
        {
            immutable real drift = s.vertexPotState.total - potBefore;
            if (drift > 1e-6L || drift < -1e-6L)
            { setError("site_survey: potState drift after restore (bug)");
              return -1; }
        }
        if (s.currentObjective != objBefore)
        { setError("site_survey: objective touched (bug)"); return -1; }
        return nwin;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/******************************************************************************
Slide-graph scan: bounded DFS over the defect states reachable by LEGAL
slides (dirty included -- the physical channel) from the CURRENT state,
with exact per-edge dS and exact rollback (notes/FPKMC_DESIGN.md M2).

The mobile defect is identified by `root_a, root_b` (a degree-3 chord of
the current state). Nodes are keyed EXACTLY by the overlay: the sorted
list of (edge, degree) deviations from the scan-start state -- no hashing,
no collisions. Expansion is pruned at cum_dS > dS_max (recorded as
boundary, not traversed), depth > max_depth, or the node/edge caps.

Outputs (caller-allocated):
  node_dS[n]      cumulative action vs the root state
  node_depth[n]   BFS depth (slides from root)
  node_nchords[n] number of degree-3 chords (the state's slide handles)
  node_sig[n]     packed illegal-degree histogram of the overlay
                  (4 bits per degree 1..16 -- the species signature)
  edge_src/dst, edge_dS, edge_chord_a/b, edge_slot   per directed edge
Returns the node count, sets *n_edges_out; -1 on error. The scan restores
the manifold and potential state exactly (audited; error on drift).
*/
extern(C) long ddg_sampler_slide_graph_scan(void* sampler_handle,
    int root_a, int root_b, double dS_max, int max_depth, long max_nodes,
    long max_edges,
    double* node_dS, int* node_depth, int* node_nchords, long* node_sig,
    int* node_chord_a, int* node_chord_b,
    int* edge_src, int* edge_dst, double* edge_dS,
    int* edge_chord_a, int* edge_chord_b, int* edge_slot,
    long* n_edges_out,
    const(int)* blocked_verts, long n_blocked, int* node_dock) nothrow
{
    // FP mode (blocked_verts non-null): node_dock[i] = 1 iff node i's
    // knot complex (every current deg-3 chord + its link vertices)
    // intersects the one-tet-layer neighborhood of the blocked vertex
    // set, computed ONCE from the root state. Dock nodes are absorbing:
    // enumerated with exact dS but never expanded. The static rule is a
    // state-independent boundary (exactness needs only that FP and any
    // reference dynamics absorb on the SAME rule) and is one tet-layer
    // conservative: interior states cannot be tet-adjacent to the
    // frozen defect, so their energies are additive (no-halo).
    // blocked_verts null (HB path) => behavior identical to before.
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("graph scan is dim=3 only"); return -1; }
        if (blocked_verts !is null && n_blocked > 0 && node_dock is null)
        { setError("blocked_verts given but node_dock is null"); return -1; }
        auto mw = cast(ManifoldWrapper!3*)(
            cast(ManifoldHandle*) s.manifoldHandle).ptr;

        if (s.currentObjective != s.currentObjective)
            recomputeObjective(s);

        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef;
            real hingeDegreeVarianceCoef; real coDim3DegreeVarianceCoef;
            real hingeDegreeTargetCoef; real coDim3DegreeTargetCoef;
            real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);

        auto potState = s.potEnabled ? &s.vertexPotState : null;
        auto pot = s.potEnabled ? &s.vertexPot : null;
        immutable real potBefore = s.potEnabled ? s.vertexPotState.total : 0.0L;

        int[2] root = [min(root_a, root_b), max(root_a, root_b)];
        if (mw.mfd.degreeOrZero!1(root[]) != 3)
        { setError("root chord is not a degree-3 edge"); return -1; }

        // one-tet-layer neighborhood of the blocked set (root state)
        bool[int] blockedNbhd;
        if (blocked_verts !is null && n_blocked > 0)
        {
            bool[int] blockedSet;
            foreach (i; 0 .. n_blocked)
                blockedSet[blocked_verts[cast(size_t) i]] = true;
            auto nf = cast(size_t) mw.mfd.numFacets;
            auto fbuf = new int[](4 * nf);
            fbuf[] = 0;
            mw.mfd.writeFacetsToBuffer(fbuf.ptr);
            foreach (fi; 0 .. nf)
            {
                auto f = fbuf[4 * fi .. 4 * fi + 4];
                bool hit = false;
                foreach (v; f) if (v in blockedSet) { hit = true; break; }
                if (hit) foreach (v; f) blockedNbhd[v] = true;
            }
        }

        // lazily captured baseline degrees + tracked degree-3 set
        size_t[int[2]] baseline;
        bool[int[2]] deg3;
        deg3[root] = true;

        void captureBaseline(scope int[2][] pairs)
        {
            foreach (p; pairs)
                if (p !in baseline)
                    baseline[p] = mw.mfd.degreeOrZero!1(p[]);
        }
        void refreshDeg3(scope int[2][] pairs)
        {
            foreach (p; pairs)
            {
                if (mw.mfd.degreeOrZero!1(p[]) == 3) deg3[p] = true;
                else deg3.remove(p);
            }
        }

        immutable(int)[] overlayKey()
        {
            int[] flat;
            int[2][] ks = baseline.keys;
            ks.sort();
            foreach (p; ks)
            {
                immutable cur = mw.mfd.degreeOrZero!1(p[]);
                if (cur != baseline[p])
                {
                    flat ~= p[0]; flat ~= p[1]; flat ~= cast(int) cur;
                }
            }
            return cast(immutable(int)[]) flat;
        }

        long sigOf()
        {
            long sig = 0;
            foreach (p, base; baseline)
            {
                immutable cur = mw.mfd.degreeOrZero!1(p[]);
                if (cur != 0 && cur != 5 && cur != 6 && cur >= 1 && cur <= 16)
                {
                    immutable shift = 4 * (cast(int) cur - 1);
                    immutable cnt = (sig >> shift) & 0xF;
                    if (cnt < 15) sig = (sig & ~(0xFL << shift))
                                      | ((cnt + 1) << shift);
                }
            }
            // root chord may predate baseline capture
            return sig;
        }

        long nNodes = 0, nEdges = 0;
        long[immutable(int)[]] seen;
        bool capped = false;

        // dock test against the CURRENT deg-3 chords (called at node
        // creation, when deg3 reflects the node's state)
        bool dockCheck()
        {
            if (blockedNbhd.length == 0) return false;
            foreach (p, _; deg3)
            {
                if ((p[0] in blockedNbhd) || (p[1] in blockedNbhd))
                    return true;
                foreach (pr; mw.mfd.link(p[]))
                    foreach (v; pr)
                        if (v in blockedNbhd) return true;
            }
            return false;
        }

        long addNode(immutable(int)[] key, real cumDS, int depth)
        {
            if (nNodes >= max_nodes) { capped = true; return -1; }
            immutable idx = nNodes++;
            seen[key] = idx;
            node_dS[cast(size_t) idx] = cast(double) cumDS;
            node_depth[cast(size_t) idx] = depth;
            int nc = 0;
            int[2] rep = [int.max, int.max];
            foreach (p, _; deg3)
            {
                nc++;
                if (p[0] < rep[0] || (p[0] == rep[0] && p[1] < rep[1]))
                    rep = p;
            }
            node_nchords[cast(size_t) idx] = nc;
            node_sig[cast(size_t) idx] = sigOf();
            node_chord_a[cast(size_t) idx] = nc ? rep[0] : -1;
            node_chord_b[cast(size_t) idx] = nc ? rep[1] : -1;
            if (node_dock !is null)
                node_dock[cast(size_t) idx] = dockCheck() ? 1 : 0;
            return idx;
        }

        // support pairs of a decoded frame (mirrors trySlideMove)
        int[2][] supportPairs(ref SlideFrame!int f, ref int[3] link)
        {
            int[9] sup = 0;
            int ns = 0;
            void add(int v)
            {
                foreach (i; 0 .. ns) if (sup[i] == v) return;
                sup[ns++] = v;
            }
            add(f.c0); add(f.c2); add(f.c3); add(f.c4);
            add(f.c5); add(f.c6); add(f.c7); add(f.c8);
            foreach (v; link) add(v);
            int[2][] pairs;
            foreach (i; 0 .. ns)
                foreach (j; i + 1 .. ns)
                {
                    int[2] e = [sup[i], sup[j]];
                    e[].sort();
                    pairs ~= e;
                }
            return pairs;
        }

        void dfs(long myIdx, real cumDS, int depth)
        {
            if (capped || depth >= max_depth) return;
            // chord menu snapshot (deg3 mutates during expansion)
            int[2][] menu = deg3.keys;
            menu.sort();
            foreach (chord; menu)
            {
                if (capped) break;
                if (mw.mfd.degreeOrZero!1(chord[]) != 3) continue;
                // hint tet: any facet on the chord
                int[4] hint = 0;
                {
                    bool got = false;
                    foreach (pr; mw.mfd.link(chord[]))
                    {
                        int i = 0;
                        foreach (v; pr) hint[2 + i++] = v;
                        got = true;
                        break;
                    }
                    if (!got) continue;
                    hint[0] = chord[0]; hint[1] = chord[1];
                    hint[].sort();
                }
                foreach (slot; 0 .. SLIDE_SLOTS)
                {
                    if (capped) break;
                    SlideFrame!int f;
                    int[3] link;
                    if (!slideDecode(mw.mfd, chord[0], chord[1],
                                     hint[], slot, f, link))
                        continue;
                    auto pairs = supportPairs(f, link);
                    captureBaseline(pairs);
                    SlideRec!int[4] recs;
                    real dS = real.nan;
                    if (!slideApplyKeep(mw.mfd,
                            s.currentObjective + cumDS,
                            chord[0], chord[1], hint[], slot, params,
                            potState, pot, recs, dS))
                        continue;
                    refreshDeg3(pairs);
                    auto key = overlayKey();
                    long child;
                    bool fresh = false;
                    if (auto p = key in seen)
                        child = *p;
                    else
                    {
                        child = addNode(key, cumDS + dS, depth + 1);
                        fresh = child >= 0;
                    }
                    if (child >= 0 && nEdges < max_edges)
                    {
                        edge_src[cast(size_t) nEdges] = cast(int) myIdx;
                        edge_dst[cast(size_t) nEdges] = cast(int) child;
                        edge_dS[cast(size_t) nEdges] = cast(double) dS;
                        edge_chord_a[cast(size_t) nEdges] = chord[0];
                        edge_chord_b[cast(size_t) nEdges] = chord[1];
                        edge_slot[cast(size_t) nEdges] = slot;
                        nEdges++;
                    }
                    else if (nEdges >= max_edges)
                        capped = true;
                    if (fresh && cumDS + dS <= cast(real) dS_max
                        && !(node_dock !is null
                             && node_dock[cast(size_t) child] == 1))
                        dfs(child, cumDS + dS, depth + 1);
                    slideRollback(mw.mfd, recs, potState, pot);
                    refreshDeg3(pairs);
                }
            }
        }

        auto rootKey = overlayKey();       // empty: no deviations yet
        immutable rootIdx = addNode(rootKey, 0.0L, 0);
        dfs(rootIdx, 0.0L, 0);

        // restoration audit
        foreach (p, base; baseline)
            if (mw.mfd.degreeOrZero!1(p[]) != base)
            { setError("graph scan: edge-degree drift after rollback (bug)");
              return -1; }
        if (s.potEnabled)
        {
            immutable real drift = s.vertexPotState.total - potBefore;
            if (drift > 1e-6L || drift < -1e-6L)
            { setError("graph scan: potState drift (bug)"); return -1; }
        }
        if (n_edges_out !is null) *n_edges_out = nEdges;
        if (capped)
            setError("graph scan CAPPED (node/edge limit) -- result is a "
                     ~ "truncation, not the full component");
        return nNodes;
    }
    catch (Exception e) { setError(e.msg); return -1; }
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
    auto s = cast(SamplerState*) sampler_handle;
    // Setters invalidate the tracked value to NaN; recompute lazily so a
    // read between configuration and the first run() sees the real objective.
    if (s.currentObjective != s.currentObjective)
    {
        try recomputeObjective(s);
        catch (Exception e) { setError(e.msg); return double.nan; }
    }
    return cast(double) s.currentObjective;
}

/// Targeted bistellar move applied THROUGH a sampler: validates and applies
/// the move on the sampler's internal manifold, applies the forced cocycle
/// update when cocycle tracking is enabled, and invalidates the tracked
/// objective (recomputed lazily, along with the n6-potential state). This is
/// the bookkeeping-safe entry point for externally proposed compound moves
/// (e.g. the knot slide) -- the bare ddg_manifold_* targeted moves applied to
/// a sampler's manifold would silently detach the cocycle. Vertex-changing
/// moves (1-4 / 4-1) maintain the unusedVertices label pool: a 4-1 pushes
/// the removed label; a 1-4 consumes its (caller-chosen) label from the
/// pool if present there (a label above the tracked pool is simply used --
/// the manifold-level validation guarantees it is fresh). The opt-in
/// geometry ledger / event logs do NOT see these moves.
extern(C) int ddg_sampler_do_bistellar_move(void* sampler_handle,
    const(int)* center, int center_len,
    const(int)* cocenter, int cocenter_len) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        auto rc = ddg_manifold_do_bistellar_move(s.manifoldHandle,
            center, center_len, cocenter, cocenter_len);
        if (rc != 0) return rc;   // validation failed; error already set
        if (center_len == 1)
        {
            // 4->1: the removed vertex label returns to the pool.
            s.unusedVertices ~= center[0];
        }
        else if (cocenter_len == 1)
        {
            // 1->4: consume the caller-chosen label from the pool. It may
            // legitimately be absent (a fresh label above the pool).
            foreach (i; 0 .. s.unusedVertices.length)
            {
                if (s.unusedVertices[i] == cocenter[0])
                {
                    s.unusedVertices[i] = s.unusedVertices[$ - 1];
                    s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                    s.unusedVertices.assumeSafeAppend;
                    break;
                }
            }
        }
        if (s.dim == 3 && s.cocycle.enabled)
        {
            import std.algorithm : sort;
            // cocycleBistellar reads only edges that persist through the
            // move, so running it just after doMove is valid. Sorted order
            // is fine: cocSet's sorted-key sign convention absorbs it.
            auto cen = center[0 .. center_len].dup;
            cen.sort();
            auto coc = cocenter[0 .. cocenter_len].dup;
            coc.sort();
            cocycleBistellar(s.cocycle, cen, coc);
        }
        s.currentObjective = real.nan;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Configure the f0 worm channel (scheme C; dim=3 samplers only). keys/vals
/// are the umbrella table (packed spoke-multiset bucket counts -> value;
/// see sampler.wf0Key). df1/df3 (optional, may be null) make the table
/// f-ADAPTIVE: vals become build-time cumulative dS, df1/df3 each state's
/// exact f-vector offset from the corridor start, and (tube_f1, tube_f3)
/// the f-vector the values were measured at; the engine then recompiles
/// the frozen scalar table at each episode open from the current f-vector
/// (sampler.wf0Compile -- exactly unbiased, episodes start/end closed).
/// Null df1/df3 = plain frozen table (zero deltas compile to themselves).
/// ufb6 the 6 linear-fallback coefficients (n3,n4,n5,n6,n7plus,Zdeficit),
/// ufbc its constant, z0 the Z reference. Mixture weights: aof (open-flag
/// share of openings), ph/pg/bcf/bc4 the open-state kernel shares (must
/// sum to 1). maxstep sizes the exact-undo buffer. Pass n = 0 to
/// reconfigure weights while keeping the table.
extern(C) int ddg_sampler_worm_f0_config(void* sampler_handle,
    const(ulong)* keys, const(double)* vals,
    const(double)* df1, const(double)* df3,
    double tube_f1, double tube_f3, long n,
    const(double)* ufb6, double ufbc, double z0, int lmax,
    double zeta, double aof, double ph, double pg, double bcf,
    double bc4, int maxstep, double ucap_hi, double ucap_lo,
    double mu) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (s.dim != 3) { setError("f0 worm requires dim=3"); return -1; }
        if (ph + pg + bcf + bc4 > 1.0 + 1e-9 || ph < 0 || pg < 0
            || bcf <= 0 || bc4 <= 0)
        { setError("bad kernel mixture"); return -1; }
        if (n > 0)
        {
            s.wormF0.utab = null;
            s.wormF0.skel = null;
            foreach (i; 0 .. n)
            {
                s.wormF0.utab[keys[i]] = vals[i];
                s.wormF0.skel[keys[i]] = WF0Skel(vals[i],
                    df1 !is null ? df1[i] : 0.0,
                    df3 !is null ? df3[i] : 0.0);
            }
            s.wormF0.tubeF1 = tube_f1;
            s.wormF0.tubeF3 = tube_f3;
        }
        if (ufb6 !is null)
            s.wormF0.ufb[] = ufb6[0 .. 6][];
        s.wormF0.ufbc = ufbc;
        s.wormF0.z0 = z0;
        s.wormF0.lmax = lmax;
        s.wormF0.zeta = zeta;
        s.wormF0.aof = aof;
        s.wormF0.ph = ph;
        s.wormF0.pg = pg;
        s.wormF0.bcf = bcf;
        s.wormF0.bc4 = bc4;
        s.wormF0.maxstep = maxstep;
        s.wormF0.ucapHi = ucap_hi;
        s.wormF0.ucapLo = ucap_lo;
        s.wormF0.mu = mu;
        if (s.wormF0Undo.length < maxstep + 4)
            s.wormF0Undo = new WF0Applied!int[](maxstep + 4);
        s.wormF0On = true;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Run one BILOCAL (two-ball) episode: a vertex is created at one ball
/// and destroyed at the other, so the net f-change vanishes and the
/// global pins cost exactly nothing. out12 has the same layout as
/// ddg_sampler_worm_f0_episode; closedHow = 4 means the pair closed.
/// Requires ddg_sampler_worm_f0_config first (shares the umbrella,
/// caps and mu; uses zeta2/bcp for the pair sector).
extern(C) int ddg_sampler_worm_pair_episode(void* sampler_handle,
    double* out12) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (!s.wormF0On) { setError("f0 worm not configured"); return -1; }
        if (s.dim != 3) { setError("pair worm requires dim=3"); return -1; }
        // The pair open prices its fugacity as log(zeta2), and the auto
        // path stores zeta2 = 0, so an auto request would silently make
        // every open log(0) = -inf and the channel would report a clean
        // zero-commit run forever. Auto is implemented for the CHORD
        // episode only -- the pair's proposal density carries a per-draw
        // log(W/wv) term, not just state counts, so cancelling it would
        // redefine the seed proposal and break the mirror at the close.
        // Fail loudly instead of running a dead channel.
        if (s.wormF0.zeta2Auto)
        { setError("pair worm: zeta2=NaN (auto) unsupported; pass a "
                   ~ "calibrated finite zeta2"); return -1; }
        if (s.cocycle.enabled)
        { setError("pair worm is not cocycle-safe"); return -1; }
        if (s.geoLedger.trackRoles || s.geoLedger.logEvents
            || s.geoLedger.logSixFlips)
        { setError("pair worm does not mirror the geometry ledger yet");
          return -1; }
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        if (s.currentObjective != s.currentObjective)
            recomputeObjective(s);
        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);
        wf0Compile(s.wormF0, params,
            mw.mfd.fVector[1], mw.mfd.fVector[3]);
        WormF0Result res;
        immutable changed = mw.mfd.wormPairEpisode(s.currentObjective,
            s.unusedVertices, params, s.wormF0,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            s.wormF0Undo, res);
        if (out12 !is null)
        {
            out12[0] = res.opened;   out12[1] = res.head;
            out12[2] = res.steps;    out12[3] = res.closedHow;
            out12[4] = res.dS;       out12[5] = res.umax;
            out12[6] = res.nH;       out12[7] = res.accH;
            out12[8] = res.nG;       out12[9] = res.accG;
            out12[10] = res.zmin;    out12[11] = res.nZ4;
            // best close log-alpha and its ingredients (res.df reused;
            // the pair episode has no f census of its own)
            foreach (k; 0 .. 4) out12[12 + k] = res.df[k];
        }
        return changed ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Run one CHORD (2<->3) bilocal episode -- the flicker carrier. A chord
/// is created at one ball and annihilated at the other, so the net
/// f-change is (0,+1,+2,+1)+(0,-1,-2,-1) = 0 and the pins cost nothing.
/// No umbrella is used (a chord's closure condition is what a 2->3
/// creates, so there is no barrier to flatten); zeta2 is read as a LOG
/// chemical potential here. out12 as for the other episodes;
/// closedHow = 5 means the chord pair closed.
extern(C) int ddg_sampler_worm_chord_episode(void* sampler_handle,
    double* out12) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (!s.wormF0On) { setError("worm not configured"); return -1; }
        if (s.dim != 3) { setError("chord worm requires dim=3"); return -1; }
        if (s.cocycle.enabled)
        { setError("chord worm is not cocycle-safe"); return -1; }
        if (s.geoLedger.trackRoles || s.geoLedger.logEvents
            || s.geoLedger.logSixFlips)
        { setError("chord worm does not mirror the geometry ledger yet");
          return -1; }
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        if (s.currentObjective != s.currentObjective)
            recomputeObjective(s);
        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);
        WormF0Result res;
        immutable changed = mw.mfd.wormChordPairEpisode(s.currentObjective,
            params, s.wormF0,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            s.wormF0Undo, res);
        if (out12 !is null)
        {
            out12[0] = res.opened;   out12[1] = res.head;
            out12[2] = res.steps;    out12[3] = res.closedHow;
            out12[4] = res.dS;       out12[5] = res.umax;
            out12[6] = res.nH;       out12[7] = res.accH;
            out12[8] = res.nG;       out12[9] = res.accG;
            out12[10] = res.zmin;    out12[11] = res.nZ4;
        }
        return changed ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Debug probe: count chain-targeted 2->3 creation sites (isolation
/// harness for wf0ChainSites; no state change).
extern(C) long ddg_sampler_chain_sites(void* sampler_handle,
    int kmax) nothrow
{
    clearError();
    try
    {
        auto s = cast(SamplerState*) sampler_handle;
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        return wf0ChainSitesProbe(mw.mfd, kmax);
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Run one STRICT-CLOSURE chord episode. Both marks are pure flags, so
/// the sector boundary makes no move; the close fires only when one mark
/// is a degree-3 chord, the other is ABSENT, and the absent one's region
/// carries no degree-3 edge -- i.e. the flicker relocated AND its old
/// neighbourhood came out clean. out12 as for the other episodes;
/// closedHow = 7 means the strict close fired.
extern(C) int ddg_sampler_worm_chord_strict(void* sampler_handle,
    double* out12) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (!s.wormF0On) { setError("worm not configured"); return -1; }
        if (s.dim != 3) { setError("requires dim=3"); return -1; }
        if (s.cocycle.enabled)
        { setError("not cocycle-safe"); return -1; }
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        if (s.currentObjective != s.currentObjective)
            recomputeObjective(s);
        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);
        WormF0Result res;
        immutable changed = mw.mfd.wormChordStrictEpisode(
            s.currentObjective, params, s.wormF0,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null, s.wormF0Undo, res);
        if (out12 !is null)
        {
            out12[0] = res.opened;   out12[1] = res.head;
            out12[2] = res.steps;    out12[3] = res.closedHow;
            out12[4] = res.dS;       out12[5] = res.umax;
            out12[6] = res.nH;       out12[7] = res.accH;
            out12[8] = res.nG;       out12[9] = res.accG;
            out12[10] = res.zmin;    out12[11] = res.nZ4;
            foreach (k; 0 .. 4) out12[12 + k] = res.df[k];
            // catalysis audit: [head, global] max accepted dS and
            // accepted-uphill counts (the global arm is the control)
            out12[16] = res.dsArm[0];  out12[17] = res.nUpArm[0];
            out12[18] = res.dsArm[1];  out12[19] = res.nUpArm[1];
        }
        return changed ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Upload the CHORD carrier's umbrella: keys are packed endpoint-spoke
/// multisets (sampler.wf0Key over the degrees of every edge at either
/// chord endpoint), vals the replayed cumulative dS of a measured
/// catalysed path. offpen is the flat off-tube value. n = 0 clears it.
extern(C) int ddg_sampler_worm_chord_config(void* sampler_handle,
    const(ulong)* keys, const(double)* vals, long n,
    double offpen) nothrow
{
    clearError();
    try
    {
        auto s = cast(SamplerState*) sampler_handle;
        if (s is null) { setError("null handle"); return -1; }
        s.wormF0.ctab = null;
        foreach (i; 0 .. n) s.wormF0.ctab[keys[i]] = vals[i];
        s.wormF0.cfb = offpen;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Set the chord channel's AGGREGATION knobs. `region_max` is how many
/// other degree-3 edges an EMPTY mark's region may hold (0 = the
/// original region-clean test, under which genuine aggregation is
/// impossible: the only tolerated neighbour is the adopted chord, i.e.
/// the one leaving). `agg_beta` weights the destination draw by
/// exp(beta * n) with n the flickers whose support meets the site's.
/// Both default to the certified behaviour (0, 0.0).
extern(C) int ddg_sampler_worm_chord_agg(void* sampler_handle,
    int region_max, double agg_beta) nothrow
{
    clearError();
    try
    {
        auto s = cast(SamplerState*) sampler_handle;
        if (s is null) { setError("null handle"); return -1; }
        if (region_max < 0) { setError("bad region_max"); return -1; }
        if (agg_beta != agg_beta) { setError("bad agg_beta"); return -1; }
        s.wormF0.regionMax = region_max;
        s.wormF0.aggBeta = agg_beta;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Set the PAIR sector knobs (zeta2 = pair fugacity, bcp = close-pair
/// share of open steps). Everything else comes from the f0 config.
extern(C) int ddg_sampler_worm_pair_config(void* sampler_handle,
    double zeta2, double bcp, int chain_k) nothrow
{
    clearError();
    try
    {
        auto s = cast(SamplerState*) sampler_handle;
        if (s is null) { setError("null handle"); return -1; }
        if (bcp <= 0 || bcp >= 1) { setError("bad bcp"); return -1; }
        // NaN requests auto-calibration from the proposal density
        s.wormF0.zeta2Auto = (zeta2 != zeta2);
        s.wormF0.zeta2 = s.wormF0.zeta2Auto ? 0.0 : zeta2;
        if (s.wormF0.zeta2Auto && bcp > 0) s.wormF0.pclTarget = 1.0 / bcp;
        s.wormF0.bcp = bcp;
        if (chain_k >= 1) s.wormF0.chainK = chain_k;
        return 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Debug probe: D-side umbrella value at a vertex's current star.
extern(C) double ddg_sampler_worm_f0_u(void* sampler_handle, int v) nothrow
{
    clearError();
    try
    {
        auto s = cast(SamplerState*) sampler_handle;
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef);
        wf0Compile(s.wormF0, params,
            mw.mfd.fVector[1], mw.mfd.fVector[3]);
        return mw.mfd.wormF0DebugU(s.wormF0, v);
    }
    catch (Exception e) { setError(e.msg); return double.nan; }
}

/// Run one f0-worm episode. out12 receives [opened, head, steps,
/// closedHow, dS, umax, nH, accH, nG, accG, zmin, nZ4]. Returns 1 if the closed
/// state changed, 0 otherwise, -1 on error. Gated: no cocycle, no
/// geometry ledger (v2 adds CHANNEL_F0 brackets).
extern(C) int ddg_sampler_worm_f0_episode(void* sampler_handle,
    double* out10) nothrow
{
    clearError();
    try
    {
        if (sampler_handle is null) { setError("null handle"); return -1; }
        auto s = cast(SamplerState*) sampler_handle;
        if (!s.wormF0On) { setError("f0 worm not configured"); return -1; }
        if (s.dim != 3) { setError("f0 worm requires dim=3"); return -1; }
        if (s.cocycle.enabled)
        { setError("f0 worm is not cocycle-safe"); return -1; }
        if (s.geoLedger.trackRoles || s.geoLedger.logEvents
            || s.geoLedger.logSixFlips)
        { setError("f0 worm does not mirror the geometry ledger yet");
          return -1; }
        auto mh = cast(ManifoldHandle*) s.manifoldHandle;
        auto mw = cast(ManifoldWrapper!3*) mh.ptr;
        if (s.currentObjective != s.currentObjective)
            recomputeObjective(s);
        struct Params { int numFacetsTarget; real hingeDegreeTarget;
            real numFacetsCoef; real numHingesCoef; real hingeDegreeVarianceCoef;
            real coDim3DegreeVarianceCoef; real hingeDegreeTargetCoef;
            real coDim3DegreeTargetCoef; real coDim3DegreeTarget; }
        auto params = Params(s.numFacetsTarget,
            cast(real) s.hingeDegreeTarget, cast(real) s.numFacetsCoef,
            cast(real) s.numHingesCoef, cast(real) s.hingeDegreeVarianceCoef,
            cast(real) s.coDim3DegreeVarianceCoef,
            cast(real) s.hingeDegreeTargetCoef,
            cast(real) s.coDim3DegreeTargetCoef,
            cast(real) s.coDim3DegreeTarget);
        // Recompile the frozen table from the f-adaptive skeleton at the
        // episode-start f-vector (no-op for plain zero-delta tables).
        wf0Compile(s.wormF0, params,
            mw.mfd.fVector[1], mw.mfd.fVector[3]);
        WormF0Result res;
        immutable changed = mw.mfd.wormF0Episode(s.currentObjective,
            s.unusedVertices, params, s.wormF0,
            s.potEnabled ? &s.vertexPotState : null,
            s.potEnabled ? &s.vertexPot : null,
            s.wormF0Undo, res);
        if (out10 !is null)
        {
            out10[0] = res.opened;
            out10[1] = res.head;
            out10[2] = res.steps;
            out10[3] = res.closedHow;
            out10[4] = res.dS;
            out10[5] = res.umax;
            out10[6] = res.nH;
            out10[7] = res.accH;
            out10[8] = res.nG;
            out10[9] = res.accG;
            out10[10] = res.zmin;
            out10[11] = res.nZ4;
        }
        return changed ? 1 : 0;
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Configure the vertex 6-valence potential (Z-legality + chemical tilts +
/// impurity valence; see sampler.VertexPot). dim=3 samplers only. tilt5 may be
/// null (all-zero tilts). Passing all-zero coefficients disables the term.
/// imp_offset shifts the impurity quadratic's foot: V(m) = imp_coef *
/// max(0, m - imp_offset)^2 + imp_lin_coef * m. 0 reproduces the bare
/// quadratic. imp_lin_coef is a pure chemical potential on impure edges
/// (sum_v m = 2 * #impure edges), with no arrangement dependence.
extern(C) int ddg_sampler_set_n6_potential(void* sampler_handle,
    double zleg_coef, double imp_coef, const(double)* tilt5,
    long imp_offset, double imp_lin_coef) nothrow
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
        if (imp_offset < 0) { setError("imp_offset must be >= 0"); return -1; }
        state.vertexPot.impOffset = imp_offset;
        state.vertexPot.impLinCoef = cast(real) imp_lin_coef;
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
