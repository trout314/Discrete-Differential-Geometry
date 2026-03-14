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
            case 2: auto w = new ManifoldWrapper!2; w.mfd = standardSphere!2; h.ptr = cast(void*) w; break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = standardSphere!3; h.ptr = cast(void*) w; break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = standardSphere!4; h.ptr = cast(void*) w; break;
            default:
                setError("unsupported dimension: " ~ dim.to!string ~ " (must be 2, 3, or 4)");
                return null;
        }
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
            case 2: auto w = new ManifoldWrapper!2; w.mfd = Manifold!(2, int)(facets); h.ptr = cast(void*) w; break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = Manifold!(3, int)(facets); h.ptr = cast(void*) w; break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = Manifold!(4, int)(facets); h.ptr = cast(void*) w; break;
            default:
                setError("unsupported dimension: " ~ dim.to!string);
                return null;
        }
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
            case 2: auto w = new ManifoldWrapper!2; w.mfd = loadManifold!2(fileName); h.ptr = cast(void*) w; break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = loadManifold!3(fileName); h.ptr = cast(void*) w; break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = loadManifold!4(fileName); h.ptr = cast(void*) w; break;
            default:
                setError("unsupported dimension: " ~ dim.to!string);
                return null;
        }
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
            case 2: auto w = new ManifoldWrapper!2; w.mfd = (cast(ManifoldWrapper!2*) src.ptr).mfd; h.ptr = cast(void*) w; break;
            case 3: auto w = new ManifoldWrapper!3; w.mfd = (cast(ManifoldWrapper!3*) src.ptr).mfd; h.ptr = cast(void*) w; break;
            case 4: auto w = new ManifoldWrapper!4; w.mfd = (cast(ManifoldWrapper!4*) src.ptr).mfd; h.ptr = cast(void*) w; break;
            default: setError("bad dimension"); return null;
        }
        return cast(void*) h;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void ddg_manifold_free(void* handle) nothrow
{
    // No-op: D's GC will collect.
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
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        long writeFacets(int dim)(ManifoldWrapper!dim* mw, int* buf)
        {
            auto f = mw.mfd.facets;
            if (buf !is null)
            {
                int idx = 0;
                foreach (facet; f)
                    foreach (v; facet)
                        buf[idx++] = v;
            }
            return cast(long) f.length;
        }

        switch (h.dim)
        {
            case 2: return writeFacets!2(cast(ManifoldWrapper!2*) h.ptr, out_data);
            case 3: return writeFacets!3(cast(ManifoldWrapper!3*) h.ptr, out_data);
            case 4: return writeFacets!4(cast(ManifoldWrapper!4*) h.ptr, out_data);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

extern(C) long ddg_manifold_simplices(void* handle, int simp_dim, int* out_data) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        long writeSimplices(int dim)(ManifoldWrapper!dim* mw, int sdim, int* buf)
        {
            auto s = mw.mfd.simplices(sdim);
            if (buf !is null)
            {
                int idx = 0;
                foreach (simp; s)
                    foreach (v; simp)
                        buf[idx++] = v;
            }
            return cast(long) s.length;
        }

        switch (h.dim)
        {
            case 2: return writeSimplices!2(cast(ManifoldWrapper!2*) h.ptr, simp_dim, out_data);
            case 3: return writeSimplices!3(cast(ManifoldWrapper!3*) h.ptr, simp_dim, out_data);
            case 4: return writeSimplices!4(cast(ManifoldWrapper!4*) h.ptr, simp_dim, out_data);
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
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
            case 2: return cast(long)(cast(ManifoldWrapper!2*) h.ptr).mfd.validMoveCount;
            case 3: return cast(long)(cast(ManifoldWrapper!3*) h.ptr).mfd.validMoveCount;
            case 4: return cast(long)(cast(ManifoldWrapper!4*) h.ptr).mfd.validMoveCount;
            default: setError("bad dimension"); return -1;
        }
    }
    catch (Exception e) { setError(e.msg); return -1; }
}

/// Return F(x)/V(x): the importance weight that corrects the approximate
/// sampler's stationary distribution back to exp(-objective(x)).
extern(C) double ddg_manifold_importance_weight(void* handle) nothrow
{
    clearError();
    try
    {
        if (handle is null) { setError("null handle"); return -1; }
        auto h = cast(ManifoldHandle*) handle;

        double compute(int dim)(ManifoldWrapper!dim* mw)
        {
            auto V = cast(double) mw.mfd.validMoveCount;
            auto F = cast(double) mw.mfd.numFacets;
            return F / V;
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
// SimplicialComplex lifecycle
// ---------------------------------------------------------------------------

extern(C) void* ddg_sc_create() nothrow
{
    clearError();
    try
    {
        auto w = new SimplicialComplex!int;
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
        return cast(void*) w;
    }
    catch (Exception e) { setError(e.msg); return null; }
}

extern(C) void ddg_sc_free(void* handle) nothrow
{
    // No-op: D's GC will collect.
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
}

extern(C) void* ddg_sampler_create(void* manifold_handle,
    int numFacetsTarget, double hingeDegreeTarget,
    double numFacetsCoef, double numHingesCoef,
    double hingeDegreeVarianceCoef, double coDim3DegreeVarianceCoef) nothrow
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

        switch (mh.dim)
        {
            case 2: state.unusedVertices = getUnusedVertices!2((cast(ManifoldWrapper!2*)(cast(ManifoldHandle*) state.manifoldHandle).ptr).mfd); break;
            case 3: state.unusedVertices = getUnusedVertices!3((cast(ManifoldWrapper!3*)(cast(ManifoldHandle*) state.manifoldHandle).ptr).mfd); break;
            case 4: state.unusedVertices = getUnusedVertices!4((cast(ManifoldWrapper!4*)(cast(ManifoldHandle*) state.manifoldHandle).ptr).mfd); break;
            default: setError("bad dimension"); return null;
        }

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

        long runSampler(int dim)(SamplerState* s, long numMoves,
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
            }

            auto params = Params(
                s.numFacetsTarget,
                cast(real) s.hingeDegreeTarget,
                cast(real) s.numFacetsCoef,
                cast(real) s.numHingesCoef,
                cast(real) s.hingeDegreeVarianceCoef,
                cast(real) s.coDim3DegreeVarianceCoef
            );

            auto currentObjective = mw.mfd.objective(params);
            long accepted = 0;

            foreach (moveNum; 0 .. numMoves)
            {
                if (cb !is null && moveNum % 1000 == 0 && moveNum > 0)
                {
                    if (cb(moveNum, numMoves, ud) != 0)
                        return accepted;
                }

                if (s.unusedVertices.length == 0)
                    s.unusedVertices ~= cast(int) mw.mfd.fVector[0];

                auto bm = mw.mfd.chooseRandomMove(s.unusedVertices[$ - 1], params);

                real deltaObj = mw.mfd.speculativeBistellarDelta(bm, currentObjective, params);
                immutable nFacetsBefore = cast(real) mw.mfd.fVector[dim];
                immutable deltaFacets = dim + 2 - 2 * cast(int) bm.center.length;
                immutable nFacetsAfter = nFacetsBefore + deltaFacets;
                real logAlpha = -deltaObj + log(nFacetsBefore) - log(nFacetsAfter);

                if ((logAlpha >= 0) || (uniform01 <= exp(logAlpha)))
                {
                    mw.mfd.doMove(bm);
                    if (bm.coCenter.length == 1)
                    {
                        if (s.unusedVertices.length > 0)
                            s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                    }
                    if (bm.center.length == 1) s.unusedVertices ~= bm.center;
                    currentObjective += deltaObj;
                    accepted++;
                }
            }

            if (cb !is null)
                cb(numMoves, numMoves, ud);

            return accepted;
        }

        switch (state.dim)
        {
            case 2: return runSampler!2(state, num_moves, callback, user_data);
            case 3: return runSampler!3(state, num_moves, callback, user_data);
            case 4: return runSampler!4(state, num_moves, callback, user_data);
            default: setError("bad dimension"); return -1;
        }
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
            }

            auto params = Params(
                s.numFacetsTarget,
                cast(real) s.hingeDegreeTarget,
                cast(real) s.numFacetsCoef,
                cast(real) s.numHingesCoef,
                cast(real) s.hingeDegreeVarianceCoef,
                cast(real) s.coDim3DegreeVarianceCoef
            );

            auto currentObjective = mw.mfd.objective(params);
            long accepted = 0;

            foreach (moveNum; 0 .. numMoves)
            {
                if (cb !is null && moveNum % 1000 == 0 && moveNum > 0)
                {
                    if (cb(moveNum, numMoves, ud) != 0)
                        return accepted;
                }

                if (s.unusedVertices.length == 0)
                    s.unusedVertices ~= cast(int) mw.mfd.fVector[0];

                auto bm = mw.mfd.chooseRandomMove(s.unusedVertices[$ - 1], params);

                // Exact Hastings: execute, compute V_after, accept or undo
                immutable vBefore = cast(real) mw.mfd.validMoveCount;

                mw.mfd.doMove(bm);
                if (bm.coCenter.length == 1)
                {
                    if (s.unusedVertices.length > 0)
                        s.unusedVertices = s.unusedVertices[0 .. $ - 1];
                }
                if (bm.center.length == 1) s.unusedVertices ~= bm.center;

                real newObjective = mw.mfd.objective(params);
                real deltaObj = newObjective - currentObjective;
                immutable vAfter = cast(real) mw.mfd.validMoveCount;

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
    // No-op: D's GC will collect.
}
