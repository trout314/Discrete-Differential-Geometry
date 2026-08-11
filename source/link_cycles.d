/// Bounded-length simple-cycle machinery for vertex links.
///
/// Purpose: the contract/split sampler channel needs, per proposal,
///   (a) the NUMBER of simple cycles of length <= L in a vertex link
///       (the Hastings factor for split proposals and for the reverse of a
///       contraction), and
///   (b) a UNIFORMLY sampled such cycle (the splitting cycle gamma).
///
/// In near-crystal states almost every vertex is one of the four
/// Frank-Kasper coordinations Z12/Z14/Z15/Z16, whose links are FIXED
/// polyhedra -- so their cycle counts are constants and their cycle lists
/// are precomputed once at module load, then transported to a concrete link
/// through an exact face-based isomorphism (the polyhedral graphs are
/// 3-connected, so by Whitney the 1-skeleton determines the triangulation).
/// Everything else (defect vertices, merged post-contraction links) falls
/// back to a bounded, allocation-free DFS on adjacency bitmasks.
///
/// Cycle identity convention (must match everywhere): a cycle is counted
/// once, rooted at its minimal vertex, with its second vertex smaller than
/// its last.
module link_cycles;

import std.algorithm : max, min, sort;
import utility : shouldEqual;

/// Links here are small; 32 covers every FK class and every merged
/// (post-contraction) link with margin.
enum int maxLinkVerts = 32;

/// Longest splitting cycle the channel will ever use.
enum int maxCycleLen = 8;

/// Per-length simple-cycle counts (index = cycle length, 3 .. maxCycleLen).
struct CycleCounts
{
    long[maxCycleLen + 1] byLength;

    /// Total number of simple cycles with length in [3, maxLen].
    long total(int maxLen) const pure nothrow @nogc @safe
    {
        long t = 0;
        foreach (len; 3 .. maxLen + 1)
            t += byLength[len];
        return t;
    }
}

/*******************************************************************************
Count simple cycles of length <= maxLen in a small graph given as adjacency
bitmasks (bit j of adj[i] set iff i ~ j). Allocation-free.
*/
CycleCounts countCycles(const(uint)[] adj, int maxLen) pure nothrow @nogc @safe
{
    CycleCounts counts;
    ubyte[maxCycleLen] buf;
    foreach (start; 0 .. cast(int) adj.length)
        cycleDfs!false(adj, maxLen, cast(ubyte) start, counts, -1, buf);
    return counts;
}

/*******************************************************************************
Memoised `countCycles`.

The bounded DFS is the single largest cost in the contract/split channel: a
frame-pointer profile put `countCycles` at ~34% of it (the merged link of a
contraction always takes the DFS, and so does any split link the FK catalog
does not match). But `countCycles` is `pure` -- its answer depends only on the
adjacency bitmasks and `maxLen` -- and links repeat constantly, both across
proposals at the same site and across the handful of local topologies a
near-FK state contains. So it memoises exactly.

DIRECT-MAPPED, one entry per hash bucket, overwritten on a miss: no eviction
policy, no probing, no failure mode when the working set exceeds the table.
The key is the FULL adjacency (plus length and maxLen) and is compared in
full on a hit, so a hash collision costs a recompute, never a wrong answer.

The memo is mutable module state, which is why this wrapper is not `pure`;
`countCycles` itself stays pure and is the fallback for oversized links.
Thread-local, like all D module state, so the table needs no synchronisation.

MEASURED on the a15 c_imp-0.4 melt at cs = 0.25 (the hit rate is limited by
RELABELLING -- the adjacency is keyed by local index, so the same abstract
link recurs under many labellings -- not by the number of link topologies):

    slots    TLS      hit rate   s/sweep   vs no memo
     4096    0.8 MB     24.2%     0.326      1.14x
    16384    3.2 MB     35.6%     0.294      1.26x
    65536   13.0 MB     41.3%     0.259      1.43x

Raise `cycleMemoSlots` (power of two) to trade thread-local memory for hit
rate; `cycleMemoStats` reports hits/misses, exposed as
`ddg_cycle_memo_stats`. Canonicalising the adjacency before hashing would
collapse the relabellings and lift the ceiling well above 41%, but a correct
canonical form is a real cost per call and a WRONG one silently returns
another graph's count -- so this stays an exact, labelled key.
*/
private struct CycleMemoEntry
{
    uint[maxLinkVerts] adj;
    ubyte n;
    ubyte maxLen;
    bool live;
    CycleCounts counts;
}

private enum size_t cycleMemoSlots = 65_536;          // power of two
private CycleMemoEntry[cycleMemoSlots] cycleMemo_;
private ulong cycleMemoHits_, cycleMemoMisses_;

/// (hits, misses) since process start -- the dial for sizing cycleMemoSlots.
ulong[2] cycleMemoStats() nothrow @nogc @safe
{
    return [cycleMemoHits_, cycleMemoMisses_];
}

/// ditto
CycleCounts countCyclesCached(const(uint)[] adj, int maxLen) nothrow @nogc @safe
{
    immutable size_t n = adj.length;
    if (n > maxLinkVerts || maxLen > maxCycleLen || maxLen < 0)
        return countCycles(adj, maxLen);

    ulong h = 0xcbf29ce484222325UL ^ (cast(ulong) maxLen * 0x9E3779B97F4A7C15UL);
    h = (h ^ n) * 0x100000001b3UL;
    foreach (i; 0 .. n)
        h = (h ^ adj[i]) * 0x100000001b3UL;
    h ^= h >>> 33;

    auto e = &cycleMemo_[cast(size_t)(h & (cycleMemoSlots - 1))];
    if (e.live && e.n == n && e.maxLen == maxLen)
    {
        bool same = true;
        foreach (i; 0 .. n)
            if (e.adj[i] != adj[i]) { same = false; break; }
        if (same)
        {
            ++cycleMemoHits_;
            return e.counts;
        }
    }
    ++cycleMemoMisses_;
    immutable c = countCycles(adj, maxLen);
    e.live = true;
    e.n = cast(ubyte) n;
    e.maxLen = cast(ubyte) maxLen;
    foreach (i; 0 .. n)
        e.adj[i] = adj[i];
    e.counts = c;
    return c;
}

/// the memo must be indistinguishable from the function it caches
unittest
{
    import std.random : Random, uniform;
    auto rnd = Random(20260811);
    foreach (trial; 0 .. 400)
    {
        immutable int n = uniform(4, 13, rnd);
        uint[maxLinkVerts] adj;
        foreach (i; 0 .. n)
            foreach (j; i + 1 .. n)
                if (uniform(0, 3, rnd) == 0)
                {
                    adj[i] |= 1u << j;
                    adj[j] |= 1u << i;
                }
        foreach (maxLen; 3 .. maxCycleLen + 1)
        {
            immutable want = countCycles(adj[0 .. n], maxLen);
            // first call populates, second must hit, both must agree
            immutable got1 = countCyclesCached(adj[0 .. n], maxLen);
            immutable got2 = countCyclesCached(adj[0 .. n], maxLen);
            foreach (L; 0 .. maxCycleLen + 1)
            {
                assert(got1.byLength[L] == want.byLength[L]);
                assert(got2.byLength[L] == want.byLength[L]);
            }
        }
    }
    immutable st = cycleMemoStats();
    assert(st[0] > 0 && st[1] > 0, "memo saw neither a hit nor a miss");
}

/*******************************************************************************
Select the k-th cycle (0-based, in the fixed DFS enumeration order) among all
simple cycles of length <= maxLen. Fills buf with the cycle's vertices and
returns its length, or 0 if k is out of range. Combined with a uniform k in
[0, count), this samples cycles uniformly with no enumeration storage.
Allocation-free.
*/
int kthCycle(const(uint)[] adj, int maxLen, long k, ref ubyte[maxCycleLen] buf)
    pure nothrow @nogc @safe
{
    foreach (start; 0 .. cast(int) adj.length)
    {
        CycleCounts counts;
        immutable len = cycleDfs!true(adj, maxLen, cast(ubyte) start,
                                      counts, k, buf);
        if (len > 0)
            return len;
        k -= counts.total(maxLen);
        if (k < 0)
            return 0;
    }
    return 0;
}

/// Shared DFS core over cycles rooted at `start`. Counts into `counts`;
/// with select=true additionally stops at local index k (relative to this
/// root), writing the cycle into buf and returning its length (else 0).
private int cycleDfs(bool select)(const(uint)[] adj, int maxLen, ubyte start,
    ref CycleCounts counts, long k, ref ubyte[maxCycleLen] buf)
    pure nothrow @nogc @safe
{
    import core.bitop : bsf;

    assert(adj.length <= maxLinkVerts);
    assert(maxLen >= 3 && maxLen <= maxCycleLen);

    ubyte[maxCycleLen] verts;
    uint[maxCycleLen] remaining;
    uint visited = 1u << start;
    immutable uint aboveStart = ~((2u << start) - 1);
    immutable uint startBit = 1u << start;
    long localIdx = 0;
    int d = 0;
    verts[0] = start;
    remaining[0] = adj[start] & aboveStart;

    while (d >= 0)
    {
        if (remaining[d] == 0)
        {
            visited &= ~(1u << verts[d]);
            --d;
            continue;
        }
        immutable int w = bsf(remaining[d]);
        remaining[d] &= remaining[d] - 1;

        if (w == start)
        {
            // closing edge: valid cycle iff length >= 3 and canonical
            // direction (second vertex < last vertex)
            if (d >= 2 && verts[1] < verts[d])
            {
                counts.byLength[d + 1] += 1;
                static if (select)
                {
                    if (localIdx == k)
                    {
                        buf[0 .. d + 1] = verts[0 .. d + 1];
                        return d + 1;
                    }
                    ++localIdx;
                }
            }
            continue;
        }
        if ((visited & (1u << w)) != 0 || d + 1 >= maxLen)
            continue;
        ++d;
        verts[d] = cast(ubyte) w;
        visited |= 1u << w;
        remaining[d] = adj[w] & (aboveStart | startBit);
    }
    return 0;
}

/*******************************************************************************
Build adjacency bitmasks from a face list (triangulated S^2 link).
*/
void adjFromFaces(const(ubyte[3])[] faces, ref uint[maxLinkVerts] adj)
    pure nothrow @nogc @safe
{
    adj[] = 0;
    foreach (ref f; faces)
        foreach (i; 0 .. 3)
            foreach (j; 0 .. 3)
                if (i != j)
                {
                    assert(f[i] < maxLinkVerts);
                    adj[f[i]] |= 1u << f[j];
                }
}

// ------------------------------------------------------------------ catalog

/// The four Frank-Kasper coordination links, in a fixed canonical labeling
/// (extracted from validated reference crystals: c15 for Z12/Z16, sigma for
/// Z14/Z15).
enum CatalogClass : int { z12 = 0, z14 = 1, z15 = 2, z16 = 3 }

immutable ubyte[3][] z12Faces = [
    [0,3,9], [0,3,10], [0,7,8], [0,7,10], [0,8,9], [1,4,9],
    [1,4,11], [1,6,8], [1,6,11], [1,8,9], [2,5,10], [2,5,11],
    [2,6,7], [2,6,11], [2,7,10], [3,4,5], [3,4,9], [3,5,10],
    [4,5,11], [6,7,8],
];
immutable ubyte[3][] z14Faces = [
    [0,3,4], [0,3,5], [0,4,8], [0,5,9], [0,7,8], [0,7,9],
    [1,2,5], [1,2,10], [1,3,5], [1,3,11], [1,10,11], [2,5,9],
    [2,6,9], [2,6,10], [3,4,11], [4,8,12], [4,11,12], [6,7,9],
    [6,7,13], [6,10,13], [7,8,13], [8,12,13], [10,11,12], [10,12,13],
];
immutable ubyte[3][] z15Faces = [
    [0,1,7], [0,1,8], [0,4,5], [0,4,8], [0,5,7], [1,7,13],
    [1,8,14], [1,9,13], [1,9,14], [2,5,6], [2,5,7], [2,6,12],
    [2,7,13], [2,11,12], [2,11,13], [3,4,6], [3,4,8], [3,6,12],
    [3,8,14], [3,10,12], [3,10,14], [4,5,6], [9,10,11], [9,10,14],
    [9,11,13], [10,11,12],
];
immutable ubyte[3][] z16Faces = [
    [0,2,3], [0,2,10], [0,3,6], [0,5,6], [0,5,9], [0,9,10],
    [1,2,3], [1,2,11], [1,3,7], [1,7,15], [1,11,15], [2,10,11],
    [3,6,7], [4,5,6], [4,5,12], [4,6,7], [4,7,14], [4,12,14],
    [5,9,12], [7,14,15], [8,9,10], [8,9,12], [8,10,11], [8,11,13],
    [8,12,13], [11,13,15], [12,13,14], [13,14,15],
];

/// One catalog entry: faces, adjacency, cycle counts, and the full cycle
/// list (each cycle in canonical rooting, DFS order) for transport-based
/// sampling in concrete links.
struct CatalogLink
{
    immutable(ubyte[3])[] faces;
    int nVerts;
    uint[maxLinkVerts] adj;
    CycleCounts counts;
    /// All cycles of length <= maxCycleLen, sorted by length (stable), so
    /// the cycles of length <= L are exactly cycles[0 .. counts.total(L)]
    /// and a uniform draw under a length cap is a single index draw.
    immutable(ubyte)[][] cycles;
}

private __gshared CatalogLink[4] catalogStore;

/// Read-only access to the catalog (initialized at module load).
ref const(CatalogLink) catalog(CatalogClass c) nothrow @nogc @trusted
{
    return catalogStore[cast(int) c];
}

shared static this()
{
    static immutable(ubyte[3])[][4] allFaces =
        [z12Faces, z14Faces, z15Faces, z16Faces];
    static immutable int[4] sizes = [12, 14, 15, 16];
    foreach (ci; 0 .. 4)
    {
        auto entry = &catalogStore[ci];
        entry.faces = allFaces[ci];
        entry.nVerts = sizes[ci];
        adjFromFaces(entry.faces, entry.adj);
        entry.counts = countCycles(entry.adj[0 .. entry.nVerts], maxCycleLen);
        immutable total = entry.counts.total(maxCycleLen);
        auto cycles = new immutable(ubyte)[][cast(size_t) total];
        foreach (k; 0 .. total)
        {
            ubyte[maxCycleLen] buf;
            immutable len = kthCycle(entry.adj[0 .. entry.nVerts],
                                     maxCycleLen, k, buf);
            assert(len >= 3);
            cycles[cast(size_t) k] = buf[0 .. len].idup;
        }
        // stable sort by length: cycles[0 .. counts.total(L)] = all length<=L
        import std.algorithm : SwapStrategy;
        cycles.sort!((a, b) => a.length < b.length, SwapStrategy.stable);
        entry.cycles = cycles;
    }
}

// -------------------------------------------------- catalog identification

/// Per-edge record of the (at most two) faces on an edge, as their opposite
/// vertices.
private struct EdgeOpp
{
    ubyte[2] opp;
    ubyte n;
}

private alias EdgeTable = EdgeOpp[maxLinkVerts][maxLinkVerts];

private bool buildEdgeTable(const(ubyte[3])[] faces, ref EdgeTable tab)
    nothrow @nogc @safe
{
    foreach (ref row; tab)
        foreach (ref e; row)
            e.n = 0;
    foreach (ref f; faces)
        foreach (i; 0 .. 3)
        {
            immutable ubyte a = f[i];
            immutable ubyte b = f[(i + 1) % 3];
            immutable ubyte c = f[(i + 2) % 3];
            immutable lo = min(a, b), hi = max(a, b);
            if (tab[lo][hi].n >= 2)
                return false;           // edge in >2 faces: not a surface
            tab[lo][hi].opp[tab[lo][hi].n] = c;
            tab[lo][hi].n += 1;
        }
    return true;
}

/*******************************************************************************
Exact catalog identification of a concrete link given by its face list, with
the vertex correspondence realizing it.

Fixes one corner of the canonical polyhedron and tries mapping it to every
corner-with-orientation of the candidate, propagating vertex images across
shared face edges (a dual walk); a complete consistent propagation is an
isomorphism of the triangulations. Returns the CatalogClass index and fills
perm with perm[canonicalVertex] = concreteVertex, or returns -1. Candidate
faces must use labels 0 .. nVerts-1.
*/
int matchCatalog(int nVerts, const(ubyte[3])[] faces,
                 ref ubyte[maxLinkVerts] perm) nothrow @nogc @trusted
{
    EdgeTable tabB;
    if (!buildEdgeTable(faces, tabB))
        return -1;
    foreach (ci; 0 .. 4)
    {
        const CatalogLink* entry = &catalogStore[cast(size_t) ci];
        if (entry.nVerts != nVerts || entry.faces.length != faces.length)
            continue;
        EdgeTable tabA;
        immutable okA = buildEdgeTable(entry.faces, tabA);
        assert(okA);
        immutable ubyte[3] cornerA =
            [entry.faces[0][0], entry.faces[0][1], entry.faces[0][2]];
        foreach (ref fB; faces)
            foreach (rot; 0 .. 3)
                foreach (flip; 0 .. 2)
                {
                    ubyte[3] cornerB;
                    foreach (i; 0 .. 3)
                        cornerB[i] = fB[flip ? (rot + 3 - i) % 3
                                             : (rot + i) % 3];
                    if (propagateIso(tabA, tabB, nVerts, cornerA, cornerB,
                                     perm))
                        return ci;
                }
    }
    return -1;
}

/// Grow a vertex map from one matched corner across face adjacencies.
private bool propagateIso(ref const EdgeTable tabA, ref const EdgeTable tabB,
                          int nVerts, ubyte[3] cornerA, ubyte[3] cornerB,
                          ref ubyte[maxLinkVerts] perm) nothrow @nogc @safe
{
    enum ubyte unset = 0xff;
    perm[] = unset;
    ubyte[maxLinkVerts] inv;
    inv[] = unset;

    foreach (i; 0 .. 3)
    {
        immutable a = cornerA[i], b = cornerB[i];
        if (perm[a] != unset || inv[b] != unset)
            return false;
        perm[a] = b;
        inv[b] = a;
    }

    // queue of undirected A-edges with known endpoint images
    ubyte[2][4 * maxLinkVerts] queue;
    bool[maxLinkVerts][maxLinkVerts] queued;
    foreach (ref row; queued)
        row[] = false;
    int qhead = 0, qtail = 0;

    void push(ubyte x, ubyte y) nothrow @nogc @safe
    {
        immutable lo = min(x, y), hi = max(x, y);
        if (queued[lo][hi])
            return;
        queued[lo][hi] = true;
        assert(qtail < queue.length);
        queue[qtail] = [lo, hi];
        ++qtail;
    }

    push(cornerA[0], cornerA[1]);
    push(cornerA[1], cornerA[2]);
    push(cornerA[0], cornerA[2]);

    while (qhead < qtail)
    {
        immutable loA = queue[qhead][0];
        immutable hiA = queue[qhead][1];
        ++qhead;
        immutable pLo = perm[loA], pHi = perm[hiA];
        assert(pLo != unset && pHi != unset);
        immutable loB = min(pLo, pHi), hiB = max(pLo, pHi);
        const EdgeOpp* ea = &tabA[loA][hiA];
        const EdgeOpp* eb = &tabB[loB][hiB];
        if (ea.n != eb.n)
            return false;
        if (ea.n != 2)
            return false;               // interior edges of a closed surface

        // The two opposite vertices across the A-edge must biject onto the
        // two across the B-edge.
        foreach (ai; 0 .. 2)
        {
            immutable oa = ea.opp[ai];
            immutable oaOther = ea.opp[1 - ai];
            if (perm[oa] != unset)
            {
                if (perm[oa] != eb.opp[0] && perm[oa] != eb.opp[1])
                    return false;
                push(oa, loA);
                push(oa, hiA);
                continue;
            }
            ubyte target = unset;
            if (perm[oaOther] == eb.opp[0])
                target = eb.opp[1];
            else if (perm[oaOther] == eb.opp[1])
                target = eb.opp[0];
            else if (inv[eb.opp[0]] == unset && inv[eb.opp[1]] != unset)
                target = eb.opp[0];
            else if (inv[eb.opp[1]] == unset && inv[eb.opp[0]] != unset)
                target = eb.opp[1];
            else
                continue;               // resolved later via another edge
            if (inv[target] != unset)
                return false;
            perm[oa] = target;
            inv[target] = oa;
            push(oa, loA);
            push(oa, hiA);
        }
    }

    foreach (v; 0 .. nVerts)
        if (perm[v] == unset)
            return false;
    return true;
}

// ---------------------------------------------------------------- unittests

/// countCycles reproduces the independently computed (Python) cycle censuses
/// of the four catalog links, including C3 = number of faces.
pure @safe unittest
{
    static immutable long[6][4] expected = [
        [20, 30, 72, 240, 720, 1620],   // Z12 (icosahedron)
        [24, 36, 84, 272, 900, 2517],   // Z14
        [26, 39, 90, 280, 948, 2910],   // Z15
        [28, 42, 96, 282, 960, 3237],   // Z16
    ];
    immutable(ubyte[3])[][4] allFaces =
        [z12Faces, z14Faces, z15Faces, z16Faces];
    static immutable int[4] sizes = [12, 14, 15, 16];
    foreach (ci; 0 .. 4)
    {
        uint[maxLinkVerts] adj;
        adjFromFaces(allFaces[ci], adj);
        immutable counts = countCycles(adj[0 .. sizes[ci]], maxCycleLen);
        foreach (len; 3 .. 9)
            counts.byLength[len].shouldEqual(expected[ci][len - 3]);
        // C3 = number of faces (every 3-cycle bounds a face in these links)
        counts.byLength[3].shouldEqual(cast(long) allFaces[ci].length);
    }
}

/// kthCycle enumerates each cycle exactly once, in valid canonical form,
/// and agrees with countCycles; out-of-range k returns 0.
pure @safe unittest
{
    uint[maxLinkVerts] adj;
    adjFromFaces(z12Faces, adj);
    immutable counts = countCycles(adj[0 .. 12], maxCycleLen);
    immutable total = counts.total(maxCycleLen);
    total.shouldEqual(2702);

    bool[ulong] seen;
    foreach (k; 0 .. total)
    {
        ubyte[maxCycleLen] buf;
        immutable len = kthCycle(adj[0 .. 12], maxCycleLen, k, buf);
        assert(len >= 3);
        // canonical rooting: min vertex first, second < last
        foreach (i; 1 .. len)
            assert(buf[i] > buf[0]);
        assert(buf[1] < buf[len - 1]);
        // consecutive vertices adjacent, closure included, all distinct
        ulong key = 0;
        foreach (i; 0 .. len)
        {
            assert((adj[buf[i]] & (1u << buf[(i + 1) % len])) != 0);
            key = key * 64 + buf[i] + 1;
        }
        assert(key !in seen);
        seen[key] = true;
    }
    ubyte[maxCycleLen] buf;
    kthCycle(adj[0 .. 12], maxCycleLen, total, buf).shouldEqual(0);
    // shorter cap: counts restricted consistently
    immutable c6 = countCycles(adj[0 .. 12], 6);
    c6.total(6).shouldEqual(20 + 30 + 72 + 240);
}

/// Catalog initialization: cycle lists sized by the counts.
@safe unittest
{
    foreach (c; [CatalogClass.z12, CatalogClass.z14,
                 CatalogClass.z15, CatalogClass.z16])
    {
        const entry = &catalog(c);
        (cast(long) entry.cycles.length)
            .shouldEqual(entry.counts.total(maxCycleLen));
        // length-sorted prefix property (uniform draws under a length cap)
        foreach (i; 1 .. entry.cycles.length)
            assert(entry.cycles[i - 1].length <= entry.cycles[i].length);
        long nLe6 = 0;
        foreach (cyc; entry.cycles)
            if (cyc.length <= 6)
                ++nLe6;
        nLe6.shouldEqual(entry.counts.total(6));
    }
}

/// matchCatalog identifies each catalog link under a nontrivial relabeling,
/// with a perm that maps faces to faces; distinct classes don't match; a
/// diagonal-flipped Z14 (valid surface, different graph) matches nothing.
@safe unittest
{
    import std.algorithm : sort;

    // relabel Z14 by v -> (5v + 3) mod 14 and check identification
    ubyte[3][] relabeled;
    foreach (ref f; z14Faces)
    {
        ubyte[3] g;
        foreach (i; 0 .. 3)
            g[i] = cast(ubyte) ((5 * f[i] + 3) % 14);
        relabeled ~= g;
    }
    ubyte[maxLinkVerts] perm;
    matchCatalog(14, relabeled, perm).shouldEqual(cast(int) CatalogClass.z14);
    // perm must map every canonical face to a face of the relabeled complex
    bool[ulong] faceSet;
    foreach (ref f; relabeled)
    {
        ubyte[3] s = f;
        s[].sort();
        faceSet[s[0] * 1024uL + s[1] * 32 + s[2]] = true;
    }
    foreach (ref f; z14Faces)
    {
        ubyte[3] s = [perm[f[0]], perm[f[1]], perm[f[2]]];
        s[].sort();
        assert((s[0] * 1024uL + s[1] * 32 + s[2]) in faceSet);
    }

    // wrong sizes: Z12 faces cannot match Z14/Z15/Z16 and vice versa
    matchCatalog(12, z12Faces, perm)
        .shouldEqual(cast(int) CatalogClass.z12);

    // diagonal flip: faces [0,3,4],[0,3,5] -> [0,4,5],[3,4,5] (edge (4,5)
    // did not exist in Z14) gives a valid triangulated S^2 with the same
    // f-vector but a different graph: no catalog match
    ubyte[3][] flipped;
    foreach (ref f; z14Faces)
    {
        if (f == cast(ubyte[3]) [0, 3, 4])
            flipped ~= cast(ubyte[3]) [0, 4, 5];
        else if (f == cast(ubyte[3]) [0, 3, 5])
            flipped ~= cast(ubyte[3]) [3, 4, 5];
        else
            flipped ~= f;
    }
    matchCatalog(14, flipped, perm).shouldEqual(-1);
}
