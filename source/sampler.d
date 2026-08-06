/// Shared MCMC sampler core: objective function, move selection, and helpers.
/// Used by manifold_sampler.d, bench_sampler.d, and the C API (ddg_capi.d).
module sampler;

import std.algorithm, std.array, std.conv, std.format, std.math, std.range, std.typecons;
import std.random : uniform, uniform01, rndChoice = choice;
import manifold, manifold_moves, simplicial_complex, utility;

alias isIRof = isInputRangeOf;

// ---------------------------------------------------------------------------
// Penalty computation
// ---------------------------------------------------------------------------

struct Penalty
{
    real volumePenalty;
    real globalCurvPenalty;
    real localCurvPenalty;
    real localSolidAngleCurvPenalty;
    real hingeDegTargetPenalty;
    real coDim3DegTargetPenalty;
}

/// Compute penalties from raw values (no manifold needed).
Penalty penaltiesFromValues(int dim_, P)(
    long nFacets, long nHinges, ulong hingeTotSqDeg,
    long nCoDim3, ulong coDim3TotSqDeg, P params)
{
    enum hingesPerFacet = dim_ * (dim_ + 1) / 2;
    enum coDim3PerFacet = binomial(dim_ + 1, dim_ - 2);

    Penalty penalty;
    penalty.volumePenalty = (nFacets - params.numFacetsTarget) ^^ 2;

    immutable nHingesTarget = hingesPerFacet * nFacets / params.hingeDegreeTarget;
    penalty.globalCurvPenalty = (nHinges - nHingesTarget) ^^ 2;

    immutable degTarget = hingesPerFacet * nFacets / cast(real) nHinges;
    real _;
    real x = modf(degTarget, _);
    real minPenalty = (x - x ^^ 2) * nHinges;

    penalty.localCurvPenalty = (
        degTarget ^^ 2 * nHinges - 2 * degTarget * hingesPerFacet * nFacets + hingeTotSqDeg) - minPenalty;
    // Intensive: divide by count to get mean variance per hinge
    if (nHinges > 0) penalty.localCurvPenalty /= nHinges;

    // Fixed-target hinge penalty: sum_e (deg_e - t)^2 with t a CONSTANT target
    // (params.hingeDegreeTarget), minus the integer-lattice floor x(1-x) per
    // hinge. Extensive and linear in (nFacets, nHinges, totSqDeg), so the
    // action term is strictly local: a move's delta depends only on degrees
    // inside its bistellar ball. Identity vs. the variance form:
    //   sum(d - t)^2 = sum(d - dbar)^2 + nHinges*(dbar - t)^2.
    immutable tHinge = cast(real) params.hingeDegreeTarget;
    x = modf(tHinge, _);
    penalty.hingeDegTargetPenalty =
        tHinge ^^ 2 * nHinges - 2 * tHinge * hingesPerFacet * nFacets + hingeTotSqDeg
        - (x - x ^^ 2) * nHinges;

    static if (dim_ > 2)
    {
        immutable coDim3DegTarget = coDim3PerFacet * nFacets / cast(real) nCoDim3;
        // Codim-3 face degrees are always EVEN: a codim-3 face's link is a
        // 2-sphere, and a triangulated 2-sphere always has an even number of
        // triangles — which equals the face's facet-degree. So the minimum
        // achievable variance splits across the nearest *even* integers, giving
        // floor 4 y(1-y) with y = frac(degTarget/2), not the integer-lattice
        // x(1-x). (The hinge term above keeps x(1-x): edge links are cycles of
        // any length, so edge degrees take any integer.)
        x = modf(coDim3DegTarget / 2.0, _);
        minPenalty = 4.0 * (x - x ^^ 2) * nCoDim3;

        penalty.localSolidAngleCurvPenalty = (
            coDim3DegTarget ^^ 2 * nCoDim3 - 2 * coDim3DegTarget * coDim3PerFacet * nFacets + coDim3TotSqDeg) - minPenalty;
        // Intensive: divide by count to get mean variance per codim-3 face
        if (nCoDim3 > 0) penalty.localSolidAngleCurvPenalty /= nCoDim3;

        // Fixed-target codim-3 penalty: sum_v (deg_v - t)^2 with t a CONSTANT
        // target (params.coDim3DegreeTarget), minus the even-lattice floor
        // 4y(1-y) per face (codim-3 degrees are even; see comment above).
        // Extensive and linear in the counters, hence strictly local.
        immutable tCoDim3 = cast(real) params.coDim3DegreeTarget;
        x = modf(tCoDim3 / 2.0, _);
        penalty.coDim3DegTargetPenalty =
            tCoDim3 ^^ 2 * nCoDim3 - 2 * tCoDim3 * coDim3PerFacet * nFacets + coDim3TotSqDeg
            - 4.0 * (x - x ^^ 2) * nCoDim3;
    }
    else
    {
        penalty.localSolidAngleCurvPenalty = 0;
        penalty.coDim3DegTargetPenalty = 0;
    }

    return penalty;
}

Penalty penalties(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    immutable nFacets = mfd.fVector[dim];
    immutable nHinges = mfd.fVector[dim - 2];
    immutable totSqDeg = mfd.totalSquareDegree(dim - 2);

    static if (dim > 2)
    {
        immutable nCoDim3 = mfd.fVector[dim - 3];
        immutable totSAsqDeg = mfd.totalSquareDegree(dim - 3);
    }
    else
    {
        enum nCoDim3 = 0;
        enum totSAsqDeg = 0;
    }

    return penaltiesFromValues!dim(
        nFacets, nHinges, totSqDeg, nCoDim3, totSAsqDeg, params);
}

real objectiveFromPenalty(P)(Penalty pen, P params)
{
    return params.numFacetsCoef * pen.volumePenalty
        + params.numHingesCoef * pen.globalCurvPenalty
        + params.hingeDegreeVarianceCoef * pen.localCurvPenalty
        + params.coDim3DegreeVarianceCoef * pen.localSolidAngleCurvPenalty
        + params.hingeDegreeTargetCoef * pen.hingeDegTargetPenalty
        + params.coDim3DegreeTargetCoef * pen.coDim3DegTargetPenalty;
}

real objective(int dim, Vertex, P)(const ref Manifold!(dim, Vertex) mfd, P params)
{
    auto pen = mfd.penalties(params);
    return objectiveFromPenalty(pen, params);
}

/// Fixed-target vs variance-about-mean identity:
///   sum(d - t)^2 = sum(d - dbar)^2 + n*(dbar - t)^2
/// checked through both floor conventions on a live manifold.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.7;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    alias BM = BistellarMove!3;
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    auto params = TestParams();
    auto pen = mfd.penalties(params);

    immutable nH = cast(real) mfd.fVector[1];
    immutable nV = cast(real) mfd.fVector[0];

    // Reconstruct raw (floor-free) sums of squared deviations from each form
    real _;
    immutable dbarH = 6.0L * mfd.fVector[3] / nH;
    real xm = modf(dbarH, _);
    real xt = modf(cast(real) params.hingeDegreeTarget, _);
    immutable rawVarH = pen.localCurvPenalty * nH + (xm - xm ^^ 2) * nH;
    immutable rawTgtH = pen.hingeDegTargetPenalty + (xt - xt ^^ 2) * nH;
    assert(isClose(rawTgtH,
        rawVarH + nH * (dbarH - params.hingeDegreeTarget) ^^ 2, 1e-9L));

    immutable dbarV = 4.0L * mfd.fVector[3] / nV;
    real ym = modf(dbarV / 2, _);
    real yt = modf(cast(real) params.coDim3DegreeTarget / 2, _);
    immutable rawVarV = pen.localSolidAngleCurvPenalty * nV + 4 * (ym - ym ^^ 2) * nV;
    immutable rawTgtV = pen.coDim3DegTargetPenalty + 4 * (yt - yt ^^ 2) * nV;
    assert(isClose(rawTgtV,
        rawVarV + nV * (dbarV - params.coDim3DegreeTarget) ^^ 2, 1e-9L));
}

// ---------------------------------------------------------------------------
// Speculative delta
// ---------------------------------------------------------------------------

/******************************************************************************
Compute the objective delta for a bistellar move without executing it.
Enumerates affected sub-simplices and looks up their degrees to compute
speculative totSqDeg deltas.
*/
real speculativeBistellarDelta(int dim, Vertex, P)(
    const ref Manifold!(dim, Vertex) mfd,
    const ref BistellarMove!(dim, Vertex) move,
    real currentObjective,
    P params)
{
    auto center = move.center;
    auto coCenter = move.coCenter;
    immutable cenLen = cast(int) center.length;
    immutable coCenLen = cast(int) coCenter.length;

    // Combined vertex set, sorted (needed for degreeMap lookup via subsetsOfSize)
    Vertex[dim + 2] allVertsBuf;
    allVertsBuf[0 .. cenLen] = center[];
    allVertsBuf[cenLen .. cenLen + coCenLen] = coCenter[];
    auto allVerts = allVertsBuf[0 .. cenLen + coCenLen];
    allVerts.sort();

    // Compute speculative f-vector
    uint[dim + 1] newFVector = mfd.fVector[0 .. dim + 1];
    newFVector[].modifyFVector(move);

    // Compute speculative totSqDeg for dimensions 0 through dim-2.
    // For each sub-simplex s of dimension d in the combined vertex set:
    //   delta(s) = (|C| - |s∩C|) - (|CC| - |s∩CC|)
    //   ΔtotSqDeg[d] += 2*deg(s)*delta + delta²
    long[dim - 1] newTotSqDeg;
    foreach (d; 0 .. dim - 1)
        newTotSqDeg[d] = cast(long) mfd.totalSquareDegree(d);

    // Enumerate subsets of each relevant dimension
    static foreach (d; 0 .. dim - 1)
    {{
        foreach (subset; allVerts[].subsetsOfSize(d + 1))
        {
            // Count how many vertices in this subset are from the center
            int s_C = 0;
            foreach (v; subset)
            {
                if (center.canFind(v)) s_C++;
            }
            int s_CC = d + 1 - s_C;
            int delta = (cenLen - s_C) - (coCenLen - s_CC);

            if (delta == 0) continue;

            long deg = cast(long) mfd.degreeOrZero!d(subset);
            newTotSqDeg[d] += 2 * deg * delta + cast(long) delta * delta;
        }
    }}

    // Compute new objective from speculative values
    static if (dim > 2)
    {
        auto newPen = penaltiesFromValues!dim(
            cast(long) newFVector[dim], cast(long) newFVector[dim - 2],
            cast(ulong) newTotSqDeg[dim - 2],
            cast(long) newFVector[dim - 3],
            cast(ulong) newTotSqDeg[dim - 3],
            params);
    }
    else
    {
        auto newPen = penaltiesFromValues!dim(
            cast(long) newFVector[dim], cast(long) newFVector[dim - 2],
            cast(ulong) newTotSqDeg[dim - 2],
            0, 0, params);
    }

    return objectiveFromPenalty(newPen, params) - currentObjective;
}

///
unittest
{
    import std.random : Mt19937;

    // Test that speculative delta matches actual delta for all move types
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    import manifold_examples : standardSphere;

    // Start from a sphere and do some 1→4 moves to get a nontrivial triangulation
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    // Grow the manifold a bit to create diverse degree distributions
    alias BM = BistellarMove!3;
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Now test speculative delta on all available bistellar moves
    foreach (move; mfd.allBistellarMoves)
    {
        auto currentObj = mfd.objective(params);
        auto specDelta = mfd.speculativeBistellarDelta(move, currentObj, params);

        // Actually execute and compute
        auto mfdCopy = mfd;
        mfdCopy.doMove(move);
        auto actualNewObj = mfdCopy.objective(params);
        auto actualDelta = actualNewObj - currentObj;

        assert(isClose(specDelta, actualDelta, 1e-6),
            "speculative delta mismatch: spec=%s actual=%s for move %s"
            .format(specDelta, actualDelta, move));
    }
}

// ---------------------------------------------------------------------------
// Speculative delta for hinge moves
// ---------------------------------------------------------------------------

/******************************************************************************
Compute the objective delta for a 4-4 hinge move without executing it.
A hinge move preserves the f-vector; only totSqDeg changes for dims 0..dim-2.
Enumerates all sub-simplices of the 6 involved vertices to compute degree deltas.
*/
real speculativeHingeDelta(Vertex, P)(
    const ref Manifold!(3, Vertex) mfd,
    const ref HingeMove!Vertex move,
    real currentObjective,
    P params)
{
    enum dim = 3;

    // Collect the 6 involved vertices (sorted for subset enumeration)
    Vertex[6] allVertsBuf;
    allVertsBuf[0] = move.removedEdge[0];
    allVertsBuf[1] = move.removedEdge[1];
    allVertsBuf[2 .. 6] = move.linkCycle;
    allVertsBuf[].sort();
    auto allVerts = allVertsBuf[];

    // Get old and new facets
    auto oldFacets = move.oldFacets;
    auto newFacets = move.newFacets;

    // f-vector is unchanged by a 4-4 move
    long[dim - 1] newTotSqDeg;
    foreach (d; 0 .. dim - 1)
        newTotSqDeg[d] = cast(long) mfd.totalSquareDegree(d);

    // For each sub-simplex dimension d (0 and 1 for dim=3),
    // enumerate all subsets of the 6 vertices and compute degree deltas.
    static foreach (d; 0 .. dim - 1)
    {{
        foreach (subset; allVerts[].subsetsOfSize(d + 1))
        {
            // Count how many old facets contain this subset
            int oldCount = 0;
            foreach (ref f; oldFacets)
                if (subset.isSubsetOf(f[]))
                    oldCount++;

            // Count how many new facets contain this subset
            int newCount = 0;
            foreach (ref f; newFacets)
                if (subset.isSubsetOf(f[]))
                    newCount++;

            int delta = newCount - oldCount;
            if (delta == 0) continue;

            long deg = cast(long) mfd.degreeOrZero!d(subset);
            newTotSqDeg[d] += 2 * deg * delta + cast(long) delta * delta;
        }
    }}

    // Compute new objective from unchanged f-vector and updated totSqDeg
    auto newPen = penaltiesFromValues!dim(
        cast(long) mfd.fVector[dim], cast(long) mfd.fVector[dim - 2],
        cast(ulong) newTotSqDeg[dim - 2],
        cast(long) mfd.fVector[dim - 3],
        cast(ulong) newTotSqDeg[dim - 3],
        params);

    return objectiveFromPenalty(newPen, params) - currentObjective;
}

///
unittest
{
    // Test that speculativeHingeDelta matches actual delta
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    alias BM = BistellarMove!3;

    // Build a 3-sphere and do 1-4 moves to create degree-4 edges
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Find all degree-4 edges and test speculative delta on valid hinge moves
    int tested = 0;
    foreach (edge; mfd.simplices(1))
    {
        if (mfd.degree(edge) != 4) continue;

        int[2] edgeArr = [edge[0], edge[1]];
        // Pick a start vertex from the link
        // Use any facet containing this edge to find a link vertex
        foreach (facet; mfd.facets)
        {
            if (!edgeArr[].isSubsetOf(facet[]))
                continue;

            // Find the two vertices not in the edge
            int[2] others;
            int oi = 0;
            foreach (v; facet)
                if (v != edgeArr[0] && v != edgeArr[1])
                    others[oi++] = v;

            foreach (diag; 0 .. 2)
            {
                auto hm = mfd.hingeMove(edgeArr, others[0], diag);
                if (!mfd.hasValidHingeMove(hm)) continue;

                auto currentObj = mfd.objective(params);
                auto specDelta = mfd.speculativeHingeDelta(hm, currentObj, params);

                // Actually execute and compute
                auto mfdCopy = mfd;
                mfdCopy.doHingeMove(hm);
                auto actualNewObj = mfdCopy.objective(params);
                auto actualDelta = actualNewObj - currentObj;

                assert(isClose(specDelta, actualDelta, 1e-6),
                    "speculative hinge delta mismatch: spec=%s actual=%s for move %s"
                    .format(specDelta, actualDelta, hm));
                tested++;
            }
            break; // only need one facet per edge
        }
    }
    assert(tested > 0, "should have tested at least one hinge move");
}

// ---------------------------------------------------------------------------
// Vertex 6-valence potential (Z-legality + chemical tilts + impurity valence)
// ---------------------------------------------------------------------------
//
// A strictly-local vertex potential on two per-vertex counters (dim = 3):
//   n6(v) = # incident edges with degree >= 6
//   m(v)  = # incident edges with degree not in {5, 6} ("impurity valence")
// Energy per vertex:
//   U(n6) = zlegCoef * dist^2(n6, {0,2,3,4}) + tilt[n6] (tilt only for n6 <= 4)
//   V(m)  = impLinCoef * m + impCoef * max(0, m - impOffset)^2
// Physics: by the link sum rule (sum over incident edges of (6 - deg) = 12,
// link = S^2), a vertex with m = 0 and n6 in {0,2,3,4} is EXACTLY a
// Frank-Kasper coordination (Z12/Z14/Z15/Z16); n6 = 1 is combinatorially
// impossible at m = 0. So U penalizes illegal 6-valences (with quadratic
// analytic tail — n6 is unbounded, a clamped table would leave hubs with no
// restoring force), the tilts are chemical potentials selecting AMONG the
// legal classes (phase selection: A15-type vs Laves-type stoichiometry), and
// The two impurity terms do DIFFERENT jobs, and separating them is the point.
//
// impLinCoef is a pure CHEMICAL POTENTIAL on impure edges: sum_v m(v) =
// 2 * #impure-edges exactly, so a linear term depends on how many impure
// edges exist and not at all on where they are. Zero arrangement dependence
// is a feature here -- it controls concentration without charging anything
// for two defects touching, so a move that leaves the impure-edge COUNT
// unchanged is exactly free under it.
//
// impCoef is the CLUSTERING term, and must be strictly convex to be one at
// all. It makes an octahedral vertex, m = 6, cost 36 rather than 6.
//
// Measured caveat on the clustering term (notes/memory/): m saturates. Over
// 73 dispersed-gas snapshots AND 6 fully percolated ones, m is 1 or 2 at
// ~97% of defect vertices with max 4, in EVERY class from an isolated
// flicker to a 3,600-vertex giant component -- the complexes are curve-like
// and extend rather than thicken. So the m-distribution barely differs
// between the dispersed and percolated phases, no function of m alone
// separates them, and sum_v m^2 = 2*n_impure + 2*P acts mostly through its
// concentration part. Prefer impLinCoef for concentration and reach for
// impCoef only when genuine convexity is wanted.
//
// impOffset shifts the quadratic's foot: V(m) = impCoef * max(0, m - off)^2 is
// FLAT for m <= off and steep beyond. The point is that sum_v m^2 =
// 2*n_impure + 2*P (P = pairs of impure edges sharing a vertex), so a bare
// quadratic charges one coefficient to BOTH the concentration of defects and
// their adjacency -- and much of P is INTERNAL to a single defect (an
// isolated (3,4,4) flicker's own impure edges share vertices). With off set
// to the largest m an isolated defect already carries, a defect pays nothing
// for its own structure and two defects touch cheaply, while deep burial
// (large m) stays forbidden. off = 0 reproduces the bare quadratic exactly.
//
// Zero coefficients = disabled, zero overhead.

struct VertexPot
{
    real zlegCoef = 0;
    real impCoef = 0;
    real impLinCoef = 0;
    long impOffset = 0;
    real[5] tilt = [0, 0, 0, 0, 0];

    bool enabled() const pure nothrow @nogc @safe
    {
        if (zlegCoef != 0 || impCoef != 0 || impLinCoef != 0) return true;
        foreach (t; tilt) if (t != 0) return true;
        return false;
    }

    real U(long k) const pure nothrow @nogc @safe
    {
        real d2;
        if (k == 0 || (k >= 2 && k <= 4)) d2 = 0;
        else if (k == 1 || k == 5) d2 = 1;
        else d2 = cast(real)(k - 4) ^^ 2;   // analytic quadratic tail
        immutable t = (k >= 0 && k <= 4) ? tilt[cast(size_t) k] : 0.0L;
        return zlegCoef * d2 + t;
    }

    real V(long m) const pure nothrow @nogc @safe
    {
        real e = impLinCoef * cast(real) m;          // V(0) = 0 still
        immutable long d = m - impOffset;
        if (d > 0) e += impCoef * cast(real)(d * d); // flat foot below off
        return e;
    }
}

unittest
{
    // impOffset = 0 must reproduce the bare quadratic EXACTLY: every existing
    // run and published number used that form.
    VertexPot p;
    p.impCoef = 0.5;
    foreach (m; 0 .. 8) assert(p.V(m) == 0.5 * m * m);

    // A flat foot up to the offset, quadratic beyond, and V(0) = 0 always --
    // the state machinery assumes a vertex absent from the counter maps
    // contributes V(0) (see VertexPotState.recompute).
    p.impOffset = 2;
    assert(p.V(0) == 0 && p.V(1) == 0 && p.V(2) == 0);
    assert(p.V(3) == 0.5 && p.V(4) == 2.0 && p.V(6) == 8.0);

    // The point of the offset: the FIRST unit past the foot is cheap, while
    // deep burial stays expensive. Bare quadratic charges 2m+1 to go m -> m+1
    // from m = 2 (i.e. 5); with the foot at 2 the same step costs 1.
    VertexPot bare;
    bare.impCoef = 0.5;
    assert(bare.V(3) - bare.V(2) == 2.5);
    assert(p.V(3) - p.V(2) == 0.5);
    assert(p.V(6) - p.V(5) == 3.5);       // burial still steep

    // The linear term is a pure chemical potential: sum_v m = 2*#impure
    // edges, so it charges per impure edge and nothing for arrangement. It
    // superposes on the quadratic, and V(0) = 0 must survive (the counter
    // state assumes an absent vertex contributes V(0)).
    VertexPot lin;
    lin.impLinCoef = 0.75;
    assert(lin.V(0) == 0);
    foreach (m; 0 .. 8) assert(lin.V(m) == 0.75 * m);
    lin.impCoef = 0.5;
    lin.impOffset = 2;
    assert(lin.V(0) == 0 && lin.V(2) == 1.5);          // linear only, foot flat
    assert(lin.V(4) == 0.75 * 4 + 0.5 * 4);            // both terms beyond
}

/// Per-vertex counter state for the vertex potential. Only nonzero counters
/// are stored; a vertex absent from both maps contributes U(0) + V(0)
/// (= tilt[0]) so `total` must account for the vertex count (see recompute).
struct VertexPotState(Vertex)
{
    int[Vertex] n6;
    int[Vertex] mImp;
    real total = 0;
}

private bool ind6(long d) pure nothrow @nogc @safe { return d >= 6; }
private bool indImp(long d) pure nothrow @nogc @safe
{
    return d > 0 && (d < 5 || d > 6);
}

/// Rebuild the counter state from scratch (init, coefficient changes, tests).
void recomputeVertexPotState(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot)
{
    st.n6.clear;
    st.mImp.clear;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d))
        {
            st.n6.require(e[0], 0)++;
            st.n6.require(e[1], 0)++;
        }
        if (indImp(d))
        {
            st.mImp.require(e[0], 0)++;
            st.mImp.require(e[1], 0)++;
        }
    }
    // total = f0 * U(0)  +  corrections for vertices with nonzero counters.
    st.total = cast(real) mfd.fVector[0] * pot.U(0);
    foreach (v, k; st.n6) st.total += pot.U(k) - pot.U(0);
    foreach (v, m; st.mImp) st.total += pot.V(m);
}

/// Potential delta for a bistellar move (dim 3). Must be called BEFORE the
/// move is applied (reads pre-move degrees). With commit = true, updates the
/// counter state and running total (still pre-move: all reads are speculative).
real potentialBistellarDelta(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    const ref BistellarMove!(3, Vertex) move,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot,
    bool commit)
{
    enum dim = 3;
    auto center = move.center;
    auto coCenter = move.coCenter;
    immutable cenLen = cast(int) center.length;
    immutable coCenLen = cast(int) coCenter.length;

    Vertex[dim + 2] allVertsBuf;
    allVertsBuf[0 .. cenLen] = center[];
    allVertsBuf[cenLen .. cenLen + coCenLen] = coCenter[];
    auto allVerts = allVertsBuf[0 .. cenLen + coCenLen];
    allVerts.sort();
    immutable nv = cast(int) allVerts.length;

    int[dim + 2] dn6;
    int[dim + 2] dm;
    foreach (subset; allVerts.subsetsOfSize(2))
    {
        int s_C = 0;
        foreach (v; subset)
            if (center.canFind(v)) s_C++;
        immutable s_CC = 2 - s_C;
        immutable delta = (cenLen - s_C) - (coCenLen - s_CC);
        if (delta == 0) continue;

        immutable oldDeg = cast(long) mfd.degreeOrZero!1(subset);
        immutable newDeg = oldDeg + delta;
        immutable d6 = (ind6(newDeg) ? 1 : 0) - (ind6(oldDeg) ? 1 : 0);
        immutable dI = (indImp(newDeg) ? 1 : 0) - (indImp(oldDeg) ? 1 : 0);
        if (d6 == 0 && dI == 0) continue;

        foreach (v; subset)
        {
            foreach (i; 0 .. nv)
            {
                if (allVerts[i] == v)
                {
                    dn6[i] += d6;
                    dm[i] += dI;
                    break;
                }
            }
        }
    }

    // 1->4 creates coCenter[0]; 4->1 destroys center[0].
    immutable hasCreated = (coCenLen == 1);
    immutable hasRemoved = (cenLen == 1);

    real dS = 0;
    foreach (i; 0 .. nv)
    {
        immutable v = allVerts[i];
        immutable isCreated = hasCreated && v == coCenter[0];
        immutable isRemoved = hasRemoved && v == center[0];
        if (dn6[i] == 0 && dm[i] == 0 && !isCreated && !isRemoved) continue;

        immutable long oldN6 = st.n6.get(v, 0);
        immutable long oldM = st.mImp.get(v, 0);
        immutable real oldE = isCreated ? 0 : pot.U(oldN6) + pot.V(oldM);
        immutable long newN6 = oldN6 + dn6[i];
        immutable long newM = oldM + dm[i];
        immutable real newE = isRemoved ? 0 : pot.U(newN6) + pot.V(newM);
        dS += newE - oldE;

        if (commit)
        {
            if (isRemoved)
            {
                st.n6.remove(v);
                st.mImp.remove(v);
            }
            else
            {
                if (newN6 != 0) st.n6[v] = cast(int) newN6;
                else st.n6.remove(v);
                if (newM != 0) st.mImp[v] = cast(int) newM;
                else st.mImp.remove(v);
            }
        }
    }
    if (commit) st.total += dS;
    return dS;
}

/// Potential delta for a 4-4 hinge move (dim 3). Same contract as the
/// bistellar version (pre-move reads; commit updates state).
real potentialHingeDelta(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    const ref HingeMove!Vertex move,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot,
    bool commit)
{
    Vertex[6] allVertsBuf;
    allVertsBuf[0] = move.removedEdge[0];
    allVertsBuf[1] = move.removedEdge[1];
    allVertsBuf[2 .. 6] = move.linkCycle;
    allVertsBuf[].sort();
    auto allVerts = allVertsBuf[];

    auto oldFacets = move.oldFacets;
    auto newFacets = move.newFacets;

    int[6] dn6;
    int[6] dm;
    foreach (subset; allVerts.subsetsOfSize(2))
    {
        int oldCount = 0;
        foreach (ref f; oldFacets)
            if (subset.isSubsetOf(f[])) oldCount++;
        int newCount = 0;
        foreach (ref f; newFacets)
            if (subset.isSubsetOf(f[])) newCount++;
        immutable delta = newCount - oldCount;
        if (delta == 0) continue;

        immutable oldDeg = cast(long) mfd.degreeOrZero!1(subset);
        immutable newDeg = oldDeg + delta;
        immutable d6 = (ind6(newDeg) ? 1 : 0) - (ind6(oldDeg) ? 1 : 0);
        immutable dI = (indImp(newDeg) ? 1 : 0) - (indImp(oldDeg) ? 1 : 0);
        if (d6 == 0 && dI == 0) continue;

        foreach (v; subset)
        {
            foreach (i; 0 .. 6)
            {
                if (allVerts[i] == v)
                {
                    dn6[i] += d6;
                    dm[i] += dI;
                    break;
                }
            }
        }
    }

    real dS = 0;
    foreach (i; 0 .. 6)
    {
        if (dn6[i] == 0 && dm[i] == 0) continue;
        immutable v = allVerts[i];
        immutable long oldN6 = st.n6.get(v, 0);
        immutable long oldM = st.mImp.get(v, 0);
        immutable long newN6 = oldN6 + dn6[i];
        immutable long newM = oldM + dm[i];
        dS += (pot.U(newN6) + pot.V(newM)) - (pot.U(oldN6) + pot.V(oldM));

        if (commit)
        {
            if (newN6 != 0) st.n6[v] = cast(int) newN6;
            else st.n6.remove(v);
            if (newM != 0) st.mImp[v] = cast(int) newM;
            else st.mImp.remove(v);
        }
    }
    if (commit) st.total += dS;
    return dS;
}

// ---------------------------------------------------------------------------
// Generic block-move deltas (dim = 3)
// ---------------------------------------------------------------------------
//
// The contract/split channel (and any future block move) replaces an
// arbitrary set of tets with another; its penalty deltas are driven by the
// PLANNED facet lists rather than a per-move-type formula. Per-proposal AA
// allocations here are amortized by the channel's small proposal
// probability; revisit if a profile ever says otherwise.

/// Edge/vertex degree deltas induced by (removedFacets, addedFacets).
private void blockDegreeDeltas(Vertex)(
    const(Vertex[4])[] removedFacets,
    const(Vertex[4])[] addedFacets,
    out int[Vertex[2]] eDelta,
    out int[Vertex] vDelta)
{
    static immutable int[2][6] pairIdx =
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
    void tally(const(Vertex[4])[] facets, int sign)
    {
        foreach (ref t; facets)
        {
            foreach (ref pi; pairIdx)
            {
                Vertex[2] e = [t[pi[0]], t[pi[1]]];
                if (e[0] > e[1])
                {
                    immutable tmp = e[0];
                    e[0] = e[1];
                    e[1] = tmp;
                }
                eDelta[e] = eDelta.get(e, 0) + sign;
            }
            foreach (v; t)
                vDelta[v] = vDelta.get(v, 0) + sign;
        }
    }
    tally(removedFacets, -1);
    tally(addedFacets, +1);
}

/// Base-objective delta of a planned block move (speculative: reads pre-move
/// degrees, no mutation). Same contract as speculativeBistellarDelta.
real speculativeBlockDelta(Vertex, P)(
    const ref Manifold!(3, Vertex) mfd,
    const(Vertex[4])[] removedFacets,
    const(Vertex[4])[] addedFacets,
    real currentObjective,
    P params)
{
    enum dim = 3;
    int[Vertex[2]] eDelta;
    int[Vertex] vDelta;
    blockDegreeDeltas!Vertex(removedFacets, addedFacets, eDelta, vDelta);

    long f3 = cast(long) mfd.fVector[dim]
        + cast(long) addedFacets.length - cast(long) removedFacets.length;
    long f1 = cast(long) mfd.fVector[1];
    long f0 = cast(long) mfd.fVector[0];
    long sqE = cast(long) mfd.totalSquareDegree(1);
    long sqV = cast(long) mfd.totalSquareDegree(0);

    foreach (e, d; eDelta)
    {
        if (d == 0)
            continue;
        immutable long old = cast(long) mfd.degreeOrZero!1(e[]);
        immutable long nw = old + d;
        assert(nw >= 0, "block move drives an edge degree negative");
        if (old == 0 && nw > 0) ++f1;
        if (old > 0 && nw == 0) --f1;
        sqE += nw * nw - old * old;
    }
    foreach (v, d; vDelta)
    {
        if (d == 0)
            continue;
        Vertex[1] vs = [v];
        immutable long old = cast(long) mfd.degreeOrZero!0(vs[]);
        immutable long nw = old + d;
        assert(nw >= 0, "block move drives a vertex degree negative");
        if (old == 0 && nw > 0) ++f0;
        if (old > 0 && nw == 0) --f0;
        sqV += nw * nw - old * old;
    }

    auto newPen = penaltiesFromValues!dim(f3, f1, cast(ulong) sqE,
                                          f0, cast(ulong) sqV, params);
    return objectiveFromPenalty(newPen, params) - currentObjective;
}

/// Vertex-potential delta of a planned block move. Same contract as
/// potentialBistellarDelta (pre-move reads; commit updates the counter
/// state). createdVerts/removedVerts get the created/destroyed-vertex energy
/// semantics (old resp. new energy is zero).
real potentialBlockDelta(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    const(Vertex[4])[] removedFacets,
    const(Vertex[4])[] addedFacets,
    const(Vertex)[] createdVerts,
    const(Vertex)[] removedVerts,
    ref VertexPotState!Vertex st,
    const ref VertexPot pot,
    bool commit)
{
    import std.algorithm : canFind;

    int[Vertex[2]] eDelta;
    int[Vertex] vDelta;
    blockDegreeDeltas!Vertex(removedFacets, addedFacets, eDelta, vDelta);

    int[Vertex] dn6;
    int[Vertex] dm;
    foreach (e, d; eDelta)
    {
        if (d == 0)
            continue;
        immutable long old = cast(long) mfd.degreeOrZero!1(e[]);
        immutable long nw = old + d;
        immutable d6 = (ind6(nw) ? 1 : 0) - (ind6(old) ? 1 : 0);
        immutable dI = (indImp(nw) ? 1 : 0) - (indImp(old) ? 1 : 0);
        if (d6 == 0 && dI == 0)
            continue;
        foreach (v; e)
        {
            dn6[v] = dn6.get(v, 0) + d6;
            dm[v] = dm.get(v, 0) + dI;
        }
    }

    bool[Vertex] affected;
    foreach (v, d; dn6) affected[v] = true;
    foreach (v, d; dm) affected[v] = true;
    foreach (v; createdVerts) affected[v] = true;
    foreach (v; removedVerts) affected[v] = true;

    real dS = 0;
    foreach (v, _; affected)
    {
        immutable isCreated = createdVerts.canFind(v);
        immutable isRemoved = removedVerts.canFind(v);
        immutable long oldN6 = st.n6.get(v, 0);
        immutable long oldM = st.mImp.get(v, 0);
        immutable real oldE = isCreated ? 0 : pot.U(oldN6) + pot.V(oldM);
        immutable long newN6 = oldN6 + dn6.get(v, 0);
        immutable long newM = oldM + dm.get(v, 0);
        immutable real newE = isRemoved ? 0 : pot.U(newN6) + pot.V(newM);
        dS += newE - oldE;

        if (commit)
        {
            if (isRemoved)
            {
                st.n6.remove(v);
                st.mImp.remove(v);
            }
            else
            {
                if (newN6 != 0) st.n6[v] = cast(int) newN6;
                else st.n6.remove(v);
                if (newM != 0) st.mImp[v] = cast(int) newM;
                else st.mImp.remove(v);
            }
        }
    }
    if (commit) st.total += dS;
    return dS;
}

// ---------------------------------------------------------------------------
// Disclination-network observables (dim = 3)
// ---------------------------------------------------------------------------
//
// Near the flat mean edge degree the curvature-carrying structure is the
// DISCLINATION NETWORK: the graph of degree>=6 edges ("six-edges"). By the
// link sum rule a legal (m = 0) Z14 vertex carries exactly 2 six-edge ends
// (a disclination line passing through), Z15 exactly 3 (a branch node), Z16
// exactly 4 — so the six-edge graph IS the Frank-Kasper skeleton, and its
// components / chains / loops / grafts are the physical line objects. This
// block provides snapshot-cadence censuses of that graph plus the joint
// (n6, m) vertex census. Everything is O(#edges) over a const manifold; no
// sampler state is required, so censuses work on loaded .mfd files too.

struct DisclinationCensus
{
    long nNetVerts;      // vertices with >= 1 incident six-edge
    long nSixEdges;      // edges of degree >= 6
    long nComponents;    // connected components of the six-edge graph
    long giantSize;      // largest component (in vertices)
    long secondSize;     // second-largest component
    long giantDiameter;  // double-BFS pseudo-diameter of the giant component
    long cycleRank;      // first Betti number E - V + C of the network
    long nSegments;      // maximal chains between non-degree-2 network vertices
    long sumSegLen;      // six-edges lying in such chains
    long nPureLoops;     // components that are cycles of degree-2 vertices
    long sumLoopLen;     // six-edges lying in pure loops
    long nEndpoints;     // network vertices with exactly 1 six-edge
    long nFrayVerts;     // network vertices with impurity valence m > 0
    long nImpEndEdges;   // six-edges with >= 1 impure endpoint
    long nZ14, nZ15, nZ16;  // legal (m == 0) network-vertex class populations
    // hostMask-resolved edge categories (bit k of hostMask = "n6-class k is a
    // native host class"). An "FK" endpoint is legal with n6 in {2, 3, 4}.
    long eDopDop;        // both endpoints dopant-FK (legal, class not in mask)
    long eDopHost;       // dopant-FK -- host-FK (the graft edges)
    long eHostHost;      // both endpoints host-FK (the host's own skeleton)
    long eImpAny;        // >= 1 endpoint impure or non-FK (n6 >= 5)
    long[8] netDegCensus;   // network degree 1..7; index 7 clamps >= 8
    long[64] segLenHist;    // chain length; index clamps at 63
    long[32] compSizeHist;  // log2-binned component size (index = bsr(size))
}

/// Slot count of the flattened C-API census layout (see flattenCensus).
enum disclinationCensusSlots = 24 + 8 + 64 + 32;

/// Flatten to the C-API layout: slots 0..20 = the 21 scalar fields in
/// declaration order, 21..23 reserved (zero), 24..31 netDegCensus,
/// 32..95 segLenHist, 96..127 compSizeHist.
void flattenCensus(const ref DisclinationCensus c, long[] outBuf)
in (outBuf.length >= disclinationCensusSlots)
{
    import std.traits : Unqual;
    outBuf[0 .. disclinationCensusSlots] = 0;
    size_t i = 0;
    foreach (field; c.tupleof)
        static if (is(Unqual!(typeof(field)) == long))
            outBuf[i++] = field;
    static assert(DisclinationCensus.tupleof.length == 24);
    assert(i == 21);
    outBuf[24 .. 32] = c.netDegCensus[];
    outBuf[32 .. 96] = c.segLenHist[];
    outBuf[96 .. 128] = c.compSizeHist[];
}

/// Joint (n6, m) vertex census over ALL vertices: outBuf is row-major
/// [min(n6, n6Cap)][min(m, mCap)] with (n6Cap+1)*(mCap+1) slots. Bin (0, 0)
/// is the FK-Z12 bulk; rows n6 >= 1 sum to the disclination-network vertex
/// count. Uncapped, sum_k k * row_k = 2 * #six-edges.
void valenceCensus(Vertex)(const ref Manifold!(3, Vertex) mfd,
    long[] outBuf, int n6Cap, int mCap)
in (outBuf.length >= (n6Cap + 1) * (mCap + 1))
{
    outBuf[0 .. (n6Cap + 1) * (mCap + 1)] = 0;
    int[Vertex] n6, mImp;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d)) { n6.require(e[0], 0)++; n6.require(e[1], 0)++; }
        if (indImp(d)) { mImp.require(e[0], 0)++; mImp.require(e[1], 0)++; }
    }
    long touched = 0;
    foreach (v, k; n6)
    {
        outBuf[min(k, n6Cap) * (mCap + 1) + min(mImp.get(v, 0), mCap)]++;
        ++touched;
    }
    foreach (v, m; mImp)
    {
        if (v in n6) continue;
        outBuf[min(m, mCap)]++;   // n6 == 0 row
        ++touched;
    }
    outBuf[0] += cast(long) mfd.fVector[0] - touched;
}

/// Full disclination-network census (see DisclinationCensus). hostMask marks
/// the host's native n6 classes (e.g. C15: bit 0 | bit 4); pass 0 for no
/// host/dopant split (all legal-legal edges then count as eDopDop).
DisclinationCensus disclinationCensus(Vertex)(
    const ref Manifold!(3, Vertex) mfd, int hostMask = 0)
{
    import core.bitop : bsr;

    DisclinationCensus c;

    // Pass 1: six-edges and per-vertex (n6, m) counters.
    Vertex[2][] six;
    int[Vertex] n6, mImp;
    foreach (e; mfd.simplices(1))
    {
        immutable d = cast(long) mfd.degree(e);
        if (ind6(d))
        {
            Vertex[2] pr = [e[0], e[1]];
            six ~= pr;
            n6.require(e[0], 0)++;
            n6.require(e[1], 0)++;
        }
        if (indImp(d))
        {
            mImp.require(e[0], 0)++;
            mImp.require(e[1], 0)++;
        }
    }
    immutable nv = cast(int) n6.length;
    c.nNetVerts = nv;
    c.nSixEdges = cast(long) six.length;
    if (nv == 0) return c;

    // Dense reindex + CSR adjacency.
    int[Vertex] idx;
    auto vlab = new Vertex[nv];
    {
        int k = 0;
        foreach (v, _; n6) { idx[v] = k; vlab[k] = v; ++k; }
    }
    auto ndeg = new int[nv];
    foreach (ref e; six) { ndeg[idx[e[0]]]++; ndeg[idx[e[1]]]++; }
    auto off = new int[nv + 1];
    foreach (i; 0 .. nv) off[i + 1] = off[i] + ndeg[i];
    auto adj = new int[2 * six.length];
    {
        auto fill = off[0 .. nv].dup;
        foreach (ref e; six)
        {
            immutable a = idx[e[0]], b = idx[e[1]];
            adj[fill[a]++] = b;
            adj[fill[b]++] = a;
        }
    }

    // Vertex classification. Note ndeg[i] == n6(vlab[i]) by construction.
    auto imp = new bool[nv];
    foreach (i; 0 .. nv) imp[i] = (vlab[i] in mImp) !is null;
    bool isFK(int i) { return !imp[i] && ndeg[i] >= 2 && ndeg[i] <= 4; }
    bool isHost(int i) { return isFK(i) && ((hostMask >> ndeg[i]) & 1); }

    foreach (i; 0 .. nv)
    {
        c.netDegCensus[min(ndeg[i], 8) - 1]++;
        if (ndeg[i] == 1) c.nEndpoints++;
        if (imp[i]) c.nFrayVerts++;
        else if (ndeg[i] == 2) c.nZ14++;
        else if (ndeg[i] == 3) c.nZ15++;
        else if (ndeg[i] == 4) c.nZ16++;
    }
    foreach (ref e; six)
    {
        immutable a = idx[e[0]], b = idx[e[1]];
        if (imp[a] || imp[b]) c.nImpEndEdges++;
        if (!isFK(a) || !isFK(b)) c.eImpAny++;
        else
        {
            immutable nHost = (isHost(a) ? 1 : 0) + (isHost(b) ? 1 : 0);
            if (nHost == 2) c.eHostHost++;
            else if (nHost == 1) c.eDopHost++;
            else c.eDopDop++;
        }
    }

    // Connected components (iterative DFS with a preallocated stack).
    auto comp = new int[nv];
    comp[] = -1;
    auto work = new int[nv];
    long[] compSizes;
    int giantId = -1;
    foreach (s; 0 .. nv)
    {
        if (comp[s] >= 0) continue;
        immutable cid = cast(int) compSizes.length;
        size_t sp = 0;
        work[sp++] = s;
        comp[s] = cid;
        long size = 0;
        while (sp)
        {
            immutable u = work[--sp];
            ++size;
            foreach (w; adj[off[u] .. off[u + 1]])
                if (comp[w] < 0) { comp[w] = cid; work[sp++] = w; }
        }
        compSizes ~= size;
        if (size > c.giantSize)
        {
            c.secondSize = c.giantSize;
            c.giantSize = size;
            giantId = cid;
        }
        else if (size > c.secondSize)
            c.secondSize = size;
    }
    c.nComponents = cast(long) compSizes.length;
    c.cycleRank = c.nSixEdges - c.nNetVerts + c.nComponents;
    foreach (sz; compSizes)
        c.compSizeHist[min(bsr(cast(ulong) sz), 31)]++;

    // Pseudo-diameter of the giant component: double BFS.
    auto dist = new int[nv];
    int bfsFar(int src)
    {
        dist[] = -1;
        size_t head = 0, tail = 0;
        work[tail++] = src;
        dist[src] = 0;
        int far = src;
        while (head < tail)
        {
            immutable u = work[head++];
            if (dist[u] > dist[far]) far = u;
            foreach (w; adj[off[u] .. off[u + 1]])
                if (dist[w] < 0) { dist[w] = dist[u] + 1; work[tail++] = w; }
        }
        return far;
    }
    foreach (i; 0 .. nv)
        if (comp[i] == giantId)
        {
            immutable a = bfsFar(i);
            immutable b = bfsFar(a);
            c.giantDiameter = dist[b];
            break;
        }

    // Segment decomposition: maximal chains of degree-2 vertices between
    // nodes (network degree != 2), then leftover pure loops (all-degree-2
    // cycles = free disclination loops). The six-graph is simple (edges of a
    // simplicial complex are unique), so chain walking needs no multi-edge
    // guards. Every six-edge lands in exactly one segment or one pure loop.
    bool[ulong] seen;
    ulong ekey(int a, int b)
    {
        return a < b ? (cast(ulong) a << 32) | cast(uint) b
                     : (cast(ulong) b << 32) | cast(uint) a;
    }
    int chainNext(int prev, int cur)
    {
        foreach (w; adj[off[cur] .. off[cur + 1]])
            if (w != prev) return w;
        assert(0, "degree-2 vertex with no forward neighbor");
    }
    foreach (i; 0 .. nv)
    {
        if (ndeg[i] == 2) continue;
        foreach (w0; adj[off[i] .. off[i + 1]])
        {
            if (ekey(i, w0) in seen) continue;
            seen[ekey(i, w0)] = true;
            int prev = i, cur = w0;
            long len = 1;
            while (ndeg[cur] == 2)
            {
                immutable nxt = chainNext(prev, cur);
                seen[ekey(cur, nxt)] = true;
                prev = cur;
                cur = nxt;
                ++len;
            }
            c.nSegments++;
            c.sumSegLen += len;
            c.segLenHist[min(len, 63)]++;
        }
    }
    foreach (i; 0 .. nv)
    {
        if (ndeg[i] != 2) continue;
        foreach (w0; adj[off[i] .. off[i + 1]])
        {
            if (ekey(i, w0) in seen) continue;
            seen[ekey(i, w0)] = true;
            int prev = i, cur = w0;
            long len = 1;
            while (cur != i)
            {
                immutable nxt = chainNext(prev, cur);
                seen[ekey(cur, nxt)] = true;
                prev = cur;
                cur = nxt;
                ++len;
            }
            c.nPureLoops++;
            c.sumLoopLen += len;
        }
    }
    return c;
}

///
unittest
{
    // Boundary of the 4-simplex: every edge has degree 3 — no six-edges, and
    // every vertex has impurity valence 4 (all 4 incident edges are degree 3).
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto c = mfd.disclinationCensus(0);
    c.nNetVerts.shouldEqual(0);
    c.nSixEdges.shouldEqual(0);
    c.nComponents.shouldEqual(0);

    long[5 * 5] vc;
    mfd.valenceCensus(vc[], 4, 4);
    vc[0 * 5 + 4].shouldEqual(5);   // all 5 vertices in bin (n6=0, m=4)
    long tot = 0;
    foreach (x; vc) tot += x;
    tot.shouldEqual(5);
}

// ---------------------------------------------------------------------------
// Integer 1-cocycle tracking (T^3 winding forms / emergent coordinates)
// ---------------------------------------------------------------------------
//
// Maintains representatives of a basis of H^1(M; Z^3) as integer 1-cochains
// omega: oriented edges -> Z^3, CLOSED on every triangle (signed sum = 0).
// For M = T^3 these are the three winding forms; the pairing sum_gamma omega
// over a closed edge loop gamma is its winding vector — gauge-invariant, so
// ANY representative serves for topology (winding, spanning, class checks).
// Geometry (edge direction vectors, Fourier space, nematic order) comes from
// the harmonic representative, computed OFFLINE in Python (graph-Laplacian
// solve); the D side only maintains an exact integer representative.
//
// The update under Pachner moves is FORCED by closedness (the move ball is
// simply connected, so path sums within it are path-independent):
//   1->4  new vertex w in tet {a,b,c,d}: gauge choice places w at a
//         (omega(a->w) = 0, omega(w->v) = omega(a->v)) — pure gauge.
//   2->3  new pole-pole edge: omega(p->q) = omega(p->a) + omega(a->q) for
//         any equator vertex a (choices agree by closedness of the old
//         cocycle on the two old tets).
//   3->2, 4->1  deletions only.
//   4-4   new diagonal: omega(x->z) = omega(x->p) + omega(p->z) via a pole p;
//         then the removed edge is deleted.
// Values are Z^3 so the class is preserved EXACTLY forever; no float drift.

struct CocycleState(Vertex)
{
    int[3][Vertex[2]] omega;   // key: sorted edge (u < v); value = omega(u->v)
    bool enabled;

    // --- Vertex lift (opt-in; see cocycleSeedPositions).
    //
    // omega is a min-imaged real displacement, not merely a winding indicator
    // (build_from_positions sets omega(u->v) = p(v) - p(u) - M k(u,v)), so its
    // integral is a genuine position on the cover. `pos` maintains that
    // integral INCREMENTALLY, which is possible because every move preserves
    //
    //     pos(v) - pos(u) == omega(u->v)   (mod boxM)   on every edge
    //
    // with at most one assignment: 2->3 forces the new pole-pole value through
    // an equator vertex, so the difference telescopes and nothing needs
    // updating; 3->2 and 4->1 only delete; the 4-4 diagonal telescopes through
    // a pole the same way; and 1->4 is the sole vertex creation, where the
    // gauge choice omega(a->w) = 0 made below puts the new vertex exactly at
    // its base vertex, so pos[w] = pos[a].
    //
    // Maintaining the lift rather than re-integrating it is not only cheaper
    // (O(1) per move and per query, against an O(V+E) spanning-tree rebuild
    // plus a full-cochain marshal) -- it FIXES THE GAUGE for the lifetime of
    // the run. A rebuilt spanning tree can reassign a persisting vertex, which
    // is the origin of the "gauge glitch" steps downstream trackers have to
    // detect and discard; here they cannot arise at all.
    int[3][Vertex] pos;
    int[3] boxM;               // torus period per axis, same units as omega
    bool posEnabled;
}

/// Seed the vertex lift by integrating omega over a BFS spanning forest.
/// O(V + E), run ONCE when positions are enabled; the lift is maintained
/// incrementally thereafter. Returns an error message, or null if clean.
string cocycleSeedPositions(Vertex)(ref CocycleState!Vertex st,
    const int[3] boxM)
{
    if (!st.enabled)
        return "cocycle not enabled";
    foreach (c; 0 .. 3)
        if (boxM[c] <= 0)
            return "torus period must be positive in every direction";

    Vertex[][Vertex] adj;
    foreach (key; st.omega.byKey)
    {
        adj[key[0]] ~= key[1];
        adj[key[1]] ~= key[0];
    }
    st.pos = null;
    st.posEnabled = true;          // cocGet/pos writes below are plain AA ops
    foreach (root; adj.byKey)
    {
        if (root in st.pos) continue;
        st.pos[root] = [0, 0, 0];
        Vertex[] stack = [root];
        while (stack.length)
        {
            immutable x = stack[$ - 1];
            stack = stack[0 .. $ - 1];
            foreach (w; adj[x])
            {
                if (w in st.pos) continue;
                int[3] pw = st.pos[x];
                pw[] += cocGet(st, x, w)[];
                st.pos[w] = pw;
                stack ~= w;
            }
        }
    }
    st.boxM = boxM;
    return null;
}

/// Audit the lift: every edge must satisfy pos(v) - pos(u) == omega(u->v)
/// modulo the torus period. Returns null if clean, else a message -- the
/// exact analogue of cocycleProblems for the position channel.
string cocyclePosProblems(Vertex)(const ref CocycleState!Vertex st)
{
    import std.conv : to;
    if (!st.posEnabled)
        return "cocycle positions not enabled";
    foreach (key, val; st.omega)
    {
        auto pu = key[0] in st.pos;
        auto pv = key[1] in st.pos;
        if (pu is null || pv is null)
            return "edge (" ~ key[0].to!string ~ "," ~ key[1].to!string
                 ~ ") has an endpoint with no position";
        foreach (c; 0 .. 3)
        {
            immutable long d = cast(long)(*pv)[c] - cast(long)(*pu)[c]
                             - cast(long) val[c];
            if (d % st.boxM[c] != 0)
                return "lift broken on edge (" ~ key[0].to!string ~ ","
                     ~ key[1].to!string ~ ") axis " ~ c.to!string
                     ~ ": residual " ~ (d % st.boxM[c]).to!string;
        }
    }
    return null;
}

/// omega(a->b) with sign handling for the sorted-key convention.
private int[3] cocGet(Vertex)(const ref CocycleState!Vertex st,
    const Vertex a, const Vertex b)
{
    auto p = mkEdge(a, b) in st.omega;
    assert(p !is null, "cocycle missing edge");
    if (a < b) return *p;
    int[3] r = [-(*p)[0], -(*p)[1], -(*p)[2]];
    return r;
}

/// Store omega(a->b) = val under the sorted-key convention.
private void cocSet(Vertex)(ref CocycleState!Vertex st,
    const Vertex a, const Vertex b, const int[3] val)
{
    if (a < b)
        st.omega[mkEdge(a, b)] = val;
    else
    {
        int[3] r = [-val[0], -val[1], -val[2]];
        st.omega[mkEdge(a, b)] = r;
    }
}

/// Apply the forced cocycle update for an ACCEPTED bistellar move. Reads only
/// edges that persist through the move, so it may run before or after doMove.
void cocycleBistellar(Vertex)(ref CocycleState!Vertex st,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    final switch (cast(int) coCenter.length - 1)
    {
    case 0: // 1->4: place the new vertex at center[0] (gauge choice)
        immutable w = coCenter[0];
        immutable a = center[0];
        immutable int[3] zero = [0, 0, 0];
        cocSet(st, a, w, zero);
        foreach (v; center[1 .. $])
            cocSet(st, w, v, cocGet(st, a, v));
        // The sole vertex creation. omega(a->w) = 0 above puts w exactly at a,
        // and each spoke then satisfies pos(v) - pos(w) = pos(v) - pos(a)
        // = omega(a->v) = omega(w->v), so one assignment keeps the lift exact.
        if (st.posEnabled)
            st.pos[w] = st.pos[a];
        break;
    case 1: // 2->3: new pole-pole edge, value forced via an equator vertex
        immutable p = coCenter[0], q = coCenter[1], a = center[0];
        int[3] val = cocGet(st, p, a);
        val[] += cocGet(st, a, q)[];
        cocSet(st, p, q, val);
        break;
    case 2: // 3->2: the center edge is destroyed
        st.omega.remove(mkEdge(center[0], center[1]));
        break;
    case 3: // 4->1: the 4 spokes at the removed vertex are destroyed
        immutable w = center[0];
        foreach (v; coCenter)
            st.omega.remove(mkEdge(w, v));
        if (st.posEnabled)
            st.pos.remove(w);
        break;
    }
}

/// Hinge (4-4) analog: diagonal forced via a pole, removed edge deleted.
void cocycleHinge(Vertex)(ref CocycleState!Vertex st,
    const ref HingeMove!Vertex hm)
{
    immutable p = hm.removedEdge[0];
    int[3] val = cocGet(st, hm.addedEdge[0], p);
    val[] += cocGet(st, p, hm.addedEdge[1])[];
    cocSet(st, hm.addedEdge[0], hm.addedEdge[1], val);
    st.omega.remove(mkEdge(hm.removedEdge[0], hm.removedEdge[1]));
}

/// Full audit: key set must equal the manifold's edge set and the cochain
/// must be closed on every triangle. Returns null if clean, else a message.
/// This is the drift check — cheap enough for test/production cadence.
string cocycleProblems(Vertex)(const ref Manifold!(3, Vertex) mfd,
    const ref CocycleState!Vertex st)
{
    long ne = 0;
    foreach (e; mfd.simplices(1))
    {
        ++ne;
        if (mkEdge(e[0], e[1]) !in st.omega)
            return format("cocycle missing edge [%s, %s]", e[0], e[1]);
    }
    if (ne != st.omega.length)
        return format("cocycle has %s edges, manifold has %s",
                      st.omega.length, ne);
    foreach (t; mfd.simplices(2))
    {
        // closedness: omega(a->b) + omega(b->c) - omega(a->c) == 0
        immutable ab = cocGet(st, t[0], t[1]);
        immutable bc = cocGet(st, t[1], t[2]);
        immutable ac = cocGet(st, t[0], t[2]);
        foreach (i; 0 .. 3)
            if (ab[i] + bc[i] - ac[i] != 0)
                return format("cocycle not closed on triangle [%s, %s, %s]",
                              t[0], t[1], t[2]);
    }
    return null;
}

///
unittest
{
    // Closedness is preserved under arbitrary churn. Start from an exact
    // cocycle omega = delta(phi) with random integer phi (closed by
    // construction; on S^3 every closed cochain is exact, so this loses no
    // generality), run mixed MCMC with the cocycle attached, and audit.
    struct TestParams
    {
        int numFacetsTarget = 48;
        real hingeDegreeTarget = 5.1;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    int[int] phi;
    foreach (v; 0 .. 5) phi[v] = uniform(-50, 50);

    CocycleState!int coc;
    coc.enabled = true;
    foreach (e; mfd.simplices(1))
    {
        int[3] val = [phi[e[1]] - phi[e[0]], 2 * (phi[e[1]] - phi[e[0]]), 0];
        coc.omega[mkEdge(e[0], e[1])] = val;
    }
    assert(cocycleProblems(mfd, coc) is null);

    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    foreach (step; 0 .. 3000)
        mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
            hingeTries, hingeAccepts, bTries, bAccepts, null, null, null, null,
            &coc);

    auto prob = cocycleProblems(mfd, coc);
    assert(prob is null, "cocycle drift after churn: " ~ prob);
    assert(hingeAccepts + bAccepts[].sum > 100, "churn too weak to test");
}

/******************************************************************************
Try to propose a hinge move on a 3-manifold. Picks a random facet, a random
edge from that facet, checks for degree 4, picks a random diagonal, and
checks validity. Returns null if no valid move found in this single attempt.

Proposal is symmetric: the forward and reverse proposal probabilities are
equal (same number of containing facets, same edge count per facet, same
f-vector), so no Hastings correction is needed beyond the objective delta.
*/
Nullable!(HingeMove!Vertex) tryProposeHingeMove(Vertex)(
    ref Manifold!(3, Vertex) mfd)
{
    // Pick a random facet
    auto facet = mfd.randomFacetOfDim(3);

    // Pick a random edge from the facet (one of C(4,2)=6 edges)
    static immutable int[2][6] edgePairs = [
        [0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
    auto pair = edgePairs[uniform(0, 6)];

    Vertex[2] edge = [facet[pair[0]], facet[pair[1]]];
    edge[].sort();

    // Check degree 4
    if (mfd.degree(edge[]) != 4)
        return typeof(return).init;

    // Find a start vertex (any facet vertex not in the edge)
    Vertex startVertex;
    foreach (v; facet)
        if (v != edge[0] && v != edge[1]) { startVertex = v; break; }

    // Pick a random diagonal and construct the move
    auto hm = mfd.hingeMove(edge, startVertex, uniform(0, 2));

    if (!mfd.hasValidHingeMove(hm))
        return typeof(return).init;

    return nullable(hm);
}

///
unittest
{
    alias BM = BistellarMove!3;

    // Build a manifold with some degree-4 edges
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    // Try many proposals; at least some should succeed
    int successes = 0;
    foreach (_; 0 .. 200)
    {
        auto result = mfd.tryProposeHingeMove();
        if (!result.isNull) successes++;
    }
    assert(successes > 0, "should have found at least one valid hinge move");
}

/******************************************************************************
Per-vertex move-attribution counters (the "measured combinatorial lapse").

Attribution rule ("Option A"): every event distributes TOTAL WEIGHT 1 uniformly
over the move's support vertices. For a bistellar move the support is the 5
vertices of its bistellar ball (center ∪ coCenter) — the vertex set of the one
4-simplex the move glues onto the triangulation, hence 1/5 per vertex. For a
4-4 hinge move the support is its 6 vertices (removedEdge ∪ linkCycle); its
pentachoron-stack 4-volume is 2, applied in ANALYSIS — the bistellar and hinge
ledgers are kept separate so any weighting convention is a linear combination
downstream. For a 1-4 move the coCenter label is the vertex the move CREATES
(it is attributed like the others; intersect with surviving vertices in
analysis).

Ladder: proposed (concrete move formed, post proposal-thinning, pre validity)
→ valid (passed hasValid*, i.e. counted as a "try") → accepted. Per vertex,
valid/proposed is the kinematic-availability field and accepted/valid the
energetic (Metropolis) field.
*/
struct MoveCounters(Vertex)
{
    double[Vertex] proposed;
    double[Vertex] valid;
    double[Vertex] acceptedBistellar;
    double[Vertex] acceptedHinge;

    void clear()()
    {
        proposed = null;
        valid = null;
        acceptedBistellar = null;
        acceptedHinge = null;
    }
}

/// Distribute total weight 1 uniformly over `support` in `ledger`.
private void addSupport(Vertex)(ref double[Vertex] ledger,
    scope const(Vertex)[] support)
{
    immutable w = 1.0 / support.length;
    foreach (v; support)
        ledger[v] = ledger.get(v, 0.0) + w;
}

///
unittest
{
    double[int] ledger;
    ledger.addSupport([0, 1, 2, 3, 4]);       // one bistellar-ball event
    ledger.addSupport([0, 1, 2, 3, 4]);
    ledger.addSupport([3, 4, 5, 6, 7, 8]);    // one hinge-support event
    ledger[0].shouldEqual(0.4);
    ledger[3].isClose(0.4 + 1.0 / 6).shouldEqual(true);
    // each event contributes total weight 1
    double total = 0;
    foreach (v; ledger.byValue) total += v;
    total.isClose(3.0).shouldEqual(true);
}

/******************************************************************************
Role-resolved geometry ledger (the "maximalist" move-participation record).

Every accepted move partitions its support simplices into orbits ("roles")
under the move's symmetry, and the degree change of a simplex is a FIXED
integer determined by (move type, role) — so the role-resolved ledger is a
lossless generating set for every linear geometry-change observable (volume
flux / trace-K, degree velocity, per-edge deficit flux, the lapse, channel
decompositions). Move type codes follow coCenter.length-1 for bistellar moves,
with 4 for the 4-4 hinge move: 0:1→4, 1:2→3, 2:3→2, 3:4→1, 4:4-4.

Vertex roles and their degree changes (degree = # incident facets):
  2→3: pole(+2)×2 = coCenter, equator(0)×3 = center
  3→2: pole(−2)×2 = center,   equator(0)×3 = coCenter
  1→4: base(+2)×4 = center,   created ×1 = coCenter (born at degree 4)
  4→1: base(−2)×4 = coCenter, destroyed ×1 = center (dies at degree 4)
  4-4: pole(−2)×2 = removedEdge, diag(+2)×2 = addedEdge, passive(0)×2

Edge roles (degree = # incident facets = deficit-angle carrier):
  2→3: equator(−1)×3, spoke(+1)×6, created (pole–pole, born at 3)
  3→2: triangle(+1)×3, spoke(−1)×6, destroyed (center edge, dies at 3)
  1→4: base(+1)×6, created spokes ×4 (born at 3)
  4→1: base(−1)×6, destroyed spokes ×4 (die at 3)
  4-4: destroyed (removedEdge, dies at 4), created (diagonal, born at 4),
       pole–diag(0)×4, pole–passive(−1)×4, equator(+1)×4

Tets have no partial roles (wholesale birth/death only) and are tracked as
AGGREGATES: created/destroyed counts by move type + a log2 lifetime histogram
(age in attempted moves; tets alive when tracking started count as censored).

The optional EVENT LOG appends one fixed-size record per accepted move
(clock: u64, type: u32, labels: 6×i32 = 36 bytes packed): bistellar labels are
center-then-coCenter (support size implied by type, unused slots = -1); 4-4
labels are removedEdge then linkCycle ROTATED so the added diagonal is
(labels[2], labels[4]). The log is the full 4D cobordism, one 4-simplex per
record; roles/ledgers/incarnations are reconstructible offline by replay.
*/
enum VRole
{
    v23Pole, v23Equator,
    v32Pole, v32Equator,
    v14Base, v14Created,
    v41Base, v41Destroyed,
    v44Pole, v44Diag, v44Passive
}

enum ERole
{
    e23Equator, e23Spoke, e23Created,
    e32Triangle, e32Spoke, e32Destroyed,
    e14Base, e14Created,
    e41Base, e41Destroyed,
    e44Destroyed, e44Created, e44PoleDiag, e44PolePassive, e44Equator
}

/// Degree change per role (birth/death roles listed as 0; they are separate
/// channels, not degree increments on a persisting simplex).
immutable int[VRole.max + 1] vRoleDegreeDelta =
    [+2, 0, -2, 0, +2, 0, -2, 0, -2, +2, 0];
immutable int[ERole.max + 1] eRoleDegreeDelta =
    [-1, +1, 0, +1, -1, 0, +1, 0, -1, 0, 0, 0, 0, -1, +1];

/// True for roles that create/destroy the simplex itself.
immutable bool[VRole.max + 1] vRoleIsBirthDeath =
    [false, false, false, false, false, true, false, true, false, false, false];
immutable bool[ERole.max + 1] eRoleIsBirthDeath =
    [false, false, true, false, false, true, false, true, false, true,
     true, true, false, false, false];

enum eventRecordBytes = 36;   // u64 clock + u32 type + 6 x i32 labels, packed
enum sixRecordBytes = 20;     // u64 clock + i32 u + i32 v + i32 dir, packed

struct GeometryLedger(Vertex)
{
    bool trackRoles;    // role-resolved AAs + tet aggregates
    bool logEvents;     // fixed-size event records
    bool logSixFlips;   // six-edge (degree 5<->6 crossing) flip records

    double[Vertex][VRole.max + 1] vertexRoles;
    double[Vertex[2]][ERole.max + 1] edgeRoles;

    // Tets: aggregates only (identities churn too fast to ledger usefully).
    ulong[5] tetsCreated;        // by move type code
    ulong[5] tetsDestroyed;
    ulong[Vertex[4]] tetBirth;   // living tets -> birth clock
    ulong[64] tetLifetimeHist;   // log2-binned age (in attempted moves)
    ulong tetCensoredDeaths;     // destroyed tets born before tracking began

    ulong clock;                 // attempted moves since tracking enabled

    ubyte[] eventBuf;
    size_t eventUsed;
    bool eventOverflow;

    ubyte[] sixBuf;
    size_t sixUsed;
    bool sixOverflow;

    void clearRoles()()
    {
        foreach (ref aa; vertexRoles) aa = null;
        foreach (ref aa; edgeRoles) aa = null;
        tetsCreated[] = 0; tetsDestroyed[] = 0;
        tetBirth = null; tetLifetimeHist[] = 0;
        tetCensoredDeaths = 0;
        clock = 0;
    }
}

private void bump(K)(ref double[K] aa, const K key)
{
    aa[key] = aa.get(key, 0.0) + 1.0;
}

private Vertex[2] mkEdge(Vertex)(const Vertex a, const Vertex b)
{
    Vertex[2] e = a < b ? [a, b] : [b, a];
    return e;
}

private void tetCreate(Vertex)(ref GeometryLedger!Vertex g, Vertex[4] key)
{
    key[].sort();
    g.tetBirth[key] = g.clock;
}

private void tetDestroy(Vertex)(ref GeometryLedger!Vertex g, Vertex[4] key)
{
    import core.bitop : bsr;
    key[].sort();
    if (auto p = key in g.tetBirth)
    {
        immutable age = g.clock - *p;
        g.tetLifetimeHist[age == 0 ? 0 : bsr(age) + 1]++;
        g.tetBirth.remove(key);
    }
    else
        g.tetCensoredDeaths++;   // born before tracking started
}

/// Record an accepted bistellar move. center/coCenter as in BistellarMove.
void recordBistellar(Vertex)(ref GeometryLedger!Vertex g,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    immutable typeCode = cast(int) coCenter.length - 1;
    final switch (typeCode)
    {
    case 0: // 1->4: center = tet, coCenter = created vertex
        foreach (v; center) bump(g.vertexRoles[VRole.v14Base], v);
        bump(g.vertexRoles[VRole.v14Created], coCenter[0]);
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                bump(g.edgeRoles[ERole.e14Base], mkEdge(center[i], center[j]));
        foreach (v; center)
            bump(g.edgeRoles[ERole.e14Created], mkEdge(v, coCenter[0]));
        g.tetsDestroyed[0]++;
        {
            Vertex[4] t;
            foreach (i; 0 .. 4) t[i] = center[i];
            tetDestroy(g, t);
        }
        g.tetsCreated[0] += 4;
        foreach (skip; 0 .. 4)
        {
            Vertex[4] t; size_t n = 0;
            foreach (i, v; center) if (i != skip) t[n++] = v;
            t[3] = coCenter[0];
            tetCreate(g, t);
        }
        break;
    case 1: // 2->3: center = triangle (equator), coCenter = poles
        foreach (v; center) bump(g.vertexRoles[VRole.v23Equator], v);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v23Pole], v);
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                bump(g.edgeRoles[ERole.e23Equator], mkEdge(center[i], center[j]));
        foreach (c; center)
            foreach (p; coCenter)
                bump(g.edgeRoles[ERole.e23Spoke], mkEdge(c, p));
        bump(g.edgeRoles[ERole.e23Created], mkEdge(coCenter[0], coCenter[1]));
        g.tetsDestroyed[1] += 2;
        foreach (p; coCenter)
        {
            Vertex[4] t = [center[0], center[1], center[2], p];
            tetDestroy(g, t);
        }
        g.tetsCreated[1] += 3;
        foreach (skip; 0 .. 3)
        {
            Vertex[4] t; size_t n = 0;
            foreach (i, v; center) if (i != skip) t[n++] = v;
            t[2] = coCenter[0]; t[3] = coCenter[1];
            tetCreate(g, t);
        }
        break;
    case 2: // 3->2: center = edge (poles), coCenter = triangle (equator)
        foreach (v; center) bump(g.vertexRoles[VRole.v32Pole], v);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v32Equator], v);
        bump(g.edgeRoles[ERole.e32Destroyed], mkEdge(center[0], center[1]));
        foreach (c; center)
            foreach (q; coCenter)
                bump(g.edgeRoles[ERole.e32Spoke], mkEdge(c, q));
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                bump(g.edgeRoles[ERole.e32Triangle], mkEdge(coCenter[i], coCenter[j]));
        g.tetsDestroyed[2] += 3;
        foreach (skip; 0 .. 3)
        {
            Vertex[4] t; size_t n = 0;
            foreach (i, v; coCenter) if (i != skip) t[n++] = v;
            t[2] = center[0]; t[3] = center[1];
            tetDestroy(g, t);
        }
        g.tetsCreated[2] += 2;
        foreach (p; center)
        {
            Vertex[4] t = [coCenter[0], coCenter[1], coCenter[2], p];
            tetCreate(g, t);
        }
        break;
    case 3: // 4->1: center = destroyed vertex, coCenter = base tet
        bump(g.vertexRoles[VRole.v41Destroyed], center[0]);
        foreach (v; coCenter) bump(g.vertexRoles[VRole.v41Base], v);
        foreach (v; coCenter)
            bump(g.edgeRoles[ERole.e41Destroyed], mkEdge(center[0], v));
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                bump(g.edgeRoles[ERole.e41Base], mkEdge(coCenter[i], coCenter[j]));
        g.tetsDestroyed[3] += 4;
        foreach (skip; 0 .. 4)
        {
            Vertex[4] t; size_t n = 0;
            foreach (i, v; coCenter) if (i != skip) t[n++] = v;
            t[3] = center[0];
            tetDestroy(g, t);
        }
        g.tetsCreated[3]++;
        {
            Vertex[4] t;
            foreach (i; 0 .. 4) t[i] = coCenter[i];
            tetCreate(g, t);
        }
        break;
    }
}

/// Record an accepted 4-4 hinge move.
void recordHinge(Vertex)(ref GeometryLedger!Vertex g,
    const Vertex[2] removedEdge, const Vertex[2] addedEdge,
    const Vertex[4] linkCycleIn)
{
    // Rotate the cycle so the added diagonal is (cycle[0], cycle[2]).
    Vertex[4] lc = linkCycleIn;
    immutable d0 = mkEdge(lc[0], lc[2]);
    if (!(d0[0] == min(addedEdge[0], addedEdge[1])
          && d0[1] == max(addedEdge[0], addedEdge[1])))
        lc = [linkCycleIn[1], linkCycleIn[2], linkCycleIn[3], linkCycleIn[0]];

    foreach (v; removedEdge) bump(g.vertexRoles[VRole.v44Pole], v);
    bump(g.vertexRoles[VRole.v44Diag], lc[0]);
    bump(g.vertexRoles[VRole.v44Diag], lc[2]);
    bump(g.vertexRoles[VRole.v44Passive], lc[1]);
    bump(g.vertexRoles[VRole.v44Passive], lc[3]);

    bump(g.edgeRoles[ERole.e44Destroyed], mkEdge(removedEdge[0], removedEdge[1]));
    bump(g.edgeRoles[ERole.e44Created], mkEdge(lc[0], lc[2]));
    foreach (p; removedEdge)
    {
        bump(g.edgeRoles[ERole.e44PoleDiag], mkEdge(p, lc[0]));
        bump(g.edgeRoles[ERole.e44PoleDiag], mkEdge(p, lc[2]));
        bump(g.edgeRoles[ERole.e44PolePassive], mkEdge(p, lc[1]));
        bump(g.edgeRoles[ERole.e44PolePassive], mkEdge(p, lc[3]));
    }
    foreach (i; 0 .. 4)
        bump(g.edgeRoles[ERole.e44Equator], mkEdge(lc[i], lc[(i + 1) % 4]));

    g.tetsDestroyed[4] += 4;
    foreach (i; 0 .. 4)
    {
        Vertex[4] t = [removedEdge[0], removedEdge[1], linkCycleIn[i],
                       linkCycleIn[(i + 1) % 4]];
        tetDestroy(g, t);
    }
    g.tetsCreated[4] += 4;
    foreach (p; removedEdge)
    {
        Vertex[4] t1 = [p, lc[0], lc[1], lc[2]];
        Vertex[4] t2 = [p, lc[0], lc[2], lc[3]];
        tetCreate(g, t1);
        tetCreate(g, t2);
    }
}

/// Append one fixed-size event record (see eventRecordBytes layout).
/// Event-log ANNOTATION records (state-neutral -- consumers that mirror
/// geometry must skip type >= EVT_CHANNEL_BEGIN). Composite move channels
/// (worm / slide / nonlocal slide) bracket their elementary move records:
///   BEGIN labels: [channel, k, anchor0, anchor1, landing0, landing1]
///   EDGE  labels: [role, u, v, -1, -1, -1]   (identity sets, see ROLE_*)
///   END   labels: [channel, nMoves, dS * 1e6 (int-clamped), -1, -1, -1]
/// so trackers can treat the enclosed moves as one atomic transaction and
/// carry worldline identity across discontinuous hops instead of recording
/// spurious death + birth pairs.
enum int EVT_CHANNEL_BEGIN = 5;
enum int EVT_CHANNEL_END = 6;
enum int EVT_CHANNEL_EDGE = 7;
enum int CHANNEL_WORM = 1;
enum int CHANNEL_SLIDE = 2;
enum int CHANNEL_NONLOCAL = 3;
enum int ROLE_GONE4 = 0;
enum int ROLE_NEW4 = 1;
enum int ROLE_GONE3 = 2;
enum int ROLE_NEW3 = 3;

/// dS packed into an int label slot (units of 1e-6, clamped).
private int dSQuantized(real dS) nothrow @nogc
{
    import std.math : lround;
    immutable real q = dS * 1e6L;
    if (q >= cast(real) int.max) return int.max;
    if (q <= cast(real) int.min) return int.min;
    return cast(int) lround(q);
}

void logEvent(Vertex)(ref GeometryLedger!Vertex g, int typeCode,
    scope const(Vertex)[] labelsA, scope const(Vertex)[] labelsB)
{
    if (g.eventUsed + eventRecordBytes > g.eventBuf.length)
    {
        g.eventOverflow = true;
        return;
    }
    import core.stdc.string : memcpy;
    auto p = g.eventBuf.ptr + g.eventUsed;
    immutable ulong clk = g.clock;
    immutable uint tc = cast(uint) typeCode;
    memcpy(p, &clk, 8);
    memcpy(p + 8, &tc, 4);
    int[6] lab = -1;
    size_t n = 0;
    foreach (v; labelsA) lab[n++] = cast(int) v;
    foreach (v; labelsB) lab[n++] = cast(int) v;
    memcpy(p + 12, lab.ptr, 24);
    g.eventUsed += eventRecordBytes;
}

/// Append one six-edge flip record: (clock, u < v, dir). dir = +1 when edge
/// (u, v) crosses degree 5 -> 6 (a disclination-line edge is born), -1 for
/// 6 -> 5 (one dies). The stream is the complete rewiring history of the
/// disclination network: E6(t) = E6(0) + sum of dir over records.
private void logSixFlip(Vertex)(ref GeometryLedger!Vertex g,
    const Vertex a, const Vertex b, int dir)
{
    if (g.sixUsed + sixRecordBytes > g.sixBuf.length)
    {
        g.sixOverflow = true;
        return;
    }
    import core.stdc.string : memcpy;
    auto p = g.sixBuf.ptr + g.sixUsed;
    immutable ulong clk = g.clock;
    immutable int u = cast(int) min(a, b);
    immutable int v = cast(int) max(a, b);
    memcpy(p, &clk, 8);
    memcpy(p + 8, &u, 4);
    memcpy(p + 12, &v, 4);
    memcpy(p + 16, &dir, 4);
    g.sixUsed += sixRecordBytes;
}

/// Emit six-edge flip records for an ACCEPTED bistellar move. Must be called
/// BEFORE the move is applied (reads pre-move degrees). Only persisting edges
/// can cross the 5/6 threshold: per the role tables, created edges are born
/// at degree 3 and destroyed edges die at degree 3, so the six-edge graph
/// changes exactly at the +-1 degree crossings enumerated here.
void sixFlipsBistellar(Vertex)(ref GeometryLedger!Vertex g,
    const ref Manifold!(3, Vertex) mfd,
    scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
{
    void crossing(const Vertex a, const Vertex b, int delta)
    {
        Vertex[2] e = a < b ? [a, b] : [b, a];
        immutable d = cast(long) mfd.degree(e[]);
        if (delta > 0 && d == 5) logSixFlip(g, a, b, +1);
        else if (delta < 0 && d == 6) logSixFlip(g, a, b, -1);
    }
    final switch (cast(int) coCenter.length - 1)
    {
    case 0: // 1->4: the 6 base edges gain a facet
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                crossing(center[i], center[j], +1);
        break;
    case 1: // 2->3: equator edges -1, spokes +1
        foreach (i; 0 .. center.length)
            foreach (j; i + 1 .. center.length)
                crossing(center[i], center[j], -1);
        foreach (cv; center)
            foreach (p; coCenter)
                crossing(cv, p, +1);
        break;
    case 2: // 3->2: coCenter triangle edges +1, spokes -1
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                crossing(coCenter[i], coCenter[j], +1);
        foreach (cv; center)
            foreach (q; coCenter)
                crossing(cv, q, -1);
        break;
    case 3: // 4->1: the 6 base edges lose a facet
        foreach (i; 0 .. coCenter.length)
            foreach (j; i + 1 .. coCenter.length)
                crossing(coCenter[i], coCenter[j], -1);
        break;
    }
}

/// Hinge (4-4) analog: equator edges +1, pole-passive edges -1. The removed
/// edge dies at degree 4 and the diagonal is born at 4 — no crossings there.
void sixFlipsHinge(Vertex)(ref GeometryLedger!Vertex g,
    const ref Manifold!(3, Vertex) mfd, const ref HingeMove!Vertex hm)
{
    void crossing(const Vertex a, const Vertex b, int delta)
    {
        Vertex[2] e = a < b ? [a, b] : [b, a];
        immutable d = cast(long) mfd.degree(e[]);
        if (delta > 0 && d == 5) logSixFlip(g, a, b, +1);
        else if (delta < 0 && d == 6) logSixFlip(g, a, b, -1);
    }
    foreach (i; 0 .. 4)
        crossing(hm.linkCycle[i], hm.linkCycle[(i + 1) % 4], +1);
    foreach (p; hm.removedEdge)
        foreach (v; hm.linkCycle)
            if (v != hm.addedEdge[0] && v != hm.addedEdge[1])
                crossing(p, v, -1);
}

///
unittest
{
    // Single controlled moves against a live manifold: measured degree changes
    // must reproduce the (type, role) tables exactly.
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    GeometryLedger!int g;

    long deg(int[] s) { return mfd.degree(s); }

    // --- 1->4 on facet [0,1,2,3], new vertex 5 ---
    int[] c14 = [0,1,2,3];
    auto degBefore = [deg([0]), deg([1]), deg([2]), deg([3])];
    auto edgeBefore = deg([0,1]);
    recordBistellar(g, c14, [5]);
    mfd.doMove(BistellarMove!3([0,1,2,3], [5]));
    foreach (i, v; [0,1,2,3])
        assert(deg([v]) - degBefore[i] == vRoleDegreeDelta[VRole.v14Base]);
    assert(deg([5]) == 4);                                   // born at 4
    assert(deg([0,1]) - edgeBefore == eRoleDegreeDelta[ERole.e14Base]);
    assert(deg([0,5]) == 3);                                 // spoke born at 3
    assert(g.vertexRoles[VRole.v14Created][5] == 1.0);
    assert(g.tetsCreated[0] == 4 && g.tetsDestroyed[0] == 1);

    // --- 2->3 on triangle [0,1,2] (in tets [0,1,2,4],[0,1,2,5]? use valid) ---
    // triangle [0,1,2] now has degree 2 (tets [0,1,2,4] and [0,1,2,5]).
    assert(deg([0,1,2]) == 2);
    auto dp4 = deg([4]); auto dp5 = deg([5]); auto de01 = deg([0,1]);
    recordBistellar(g, [0,1,2], [4,5]);
    mfd.doMove(BistellarMove!3([0,1,2], [4,5]));
    assert(deg([4]) - dp4 == vRoleDegreeDelta[VRole.v23Pole]);
    assert(deg([5]) - dp5 == vRoleDegreeDelta[VRole.v23Pole]);
    assert(deg([0,1]) - de01 == eRoleDegreeDelta[ERole.e23Equator]);
    assert(deg([4,5]) == 3);                                 // pole-pole born at 3
    assert(g.edgeRoles[ERole.e23Created][mkEdge(4,5)] == 1.0);

    // --- role totals: each move contributes its exact multiplicities ---
    double tot(double[int] aa) { double s=0; foreach(v; aa.byValue) s+=v; return s; }
    assert(tot(g.vertexRoles[VRole.v14Base]) == 4.0);
    assert(tot(g.vertexRoles[VRole.v23Pole]) == 2.0);
    assert(tot(g.vertexRoles[VRole.v23Equator]) == 3.0);
    double etot(double[int[2]] aa) { double s=0; foreach(v; aa.byValue) s+=v; return s; }
    assert(etot(g.edgeRoles[ERole.e23Spoke]) == 6.0);
    assert(etot(g.edgeRoles[ERole.e14Created]) == 4.0);

    assert(mfd.findProblems.length == 0);
}

/******************************************************************************
                            KNOT SLIDE (worm move)
*******************************************************************************

A "slide" translates a (3,4,4) illegal-degree knot one step along its local
Boerdijk-Coxeter chain. It is a COMPOSITE of four Pachner moves that, taken
together, destroy the degree-3 chord (c0,c4) and create a degree-3 chord
(c4,c8) four steps down the chain, translating the knot's entire local degree
pattern:

    M1  3->2  center (c0,c4), coCenter = its 3-vertex link
    M2  2->3  center face (c3,c4,c5), coCenter apexes (c2,c6)
    M3  2->3  center face (c5,c6,c7), coCenter apexes (c4,c8)
    M4  3->2  center (c2,c6), coCenter = its (then-current) 3-vertex link

Net effect on the f-vector is ZERO in every dimension (two 3->2 and two 2->3),
so N_3 and E are exactly preserved and every global/extensive term of the
action cancels identically. The frame c5..c8 is derived from (c0,c4,c2,c3) by
the sliding-window rule, each step an O(1) ridge-link lookup:

    c5 = apex(c2,c3,c4) != c0     c6 = apex(c3,c4,c5) != c2
    c7 = apex(c4,c5,c6) != c3     c8 = apex(c5,c6,c7) != c4

MOVE CLASS: by default EVERY valid slide is in the class, accepted by plain
Metropolis on the exact dS. Detailed balance needs only a symmetric proposal
and k_f = k_r, both of which hold regardless of whether the slide preserves
the illegal-degree multiset. Setting SlideConfig.cleanOnly restricts the class
to CLEAN (species-preserving) slides -- those that leave that multiset
unchanged. Cleanliness is tested locally on the O(1) set of edges the
composite can touch (an edge whose degree changes must be an edge of an added
or removed tet, hence has both endpoints in the move support); unchanged edges
contribute identically to the global multiset, so local == global.

ACCEPTANCE is plain Metropolis on the exact action change,

    alpha = min(1, exp(-dS)),

with no Hastings correction, because k_f = k_r = 1: exactly one (chord, slot)
at the departure end produces a given end state, and exactly one at the
arrival end comes back. Verified exhaustively on 124 transitions across four
states (a crystal + one knot and three thermal lam=0.40 snapshots) -- 58 clean
and 66 dirty, every one with k_f = k_r = 1 and none with k_r = 0 -- plus zero
local-vs-global cleanliness disagreements. That dirty slides close just as
clean ones do is why the clean restriction is not the default. See
scripts/defect_dynamics/worm_slide.py --closure-test (the Python oracle).

The reverse always has a handle: M3 creates (c4,c8) as a pole-pole edge born
at degree 3, and M4's three tets (on edge (c2,c6), whose link is then
{c3,c4,c5}) contain no c8 -- so the arrival chord is a degree-3 edge whatever
else the slide did. It is checked at runtime rather than assumed, since the
whole closure argument rests on it.

PROPOSAL: an independent channel, NOT part of the unified facet proposal (see
mcmcStep for why that distinction is load-bearing). Draw a facet uniformly and
one of its 6 vertex pairs uniformly; a degree-3 edge lies in exactly 3 facets,
so every chord is proposed with the SAME probability 3/(6*N_3) -- independent
of the state's defect content, with no global chord list to maintain. Since
N_3 is slide-invariant, that factor cancels between forward and reverse. A
slot j in [0, SLIDE_SLOTS) then fixes the frame: 2 chord orientations x 6
ordered (c2,c3) picks from the sorted 3-vertex link. The slot count is
deliberately constant, not a count of valid frames: an invalid slot is a
rejected proposal, so the denominator is identical in both states.
*/

/// 2 chord orientations x 6 ordered (c2,c3) picks from the 3-vertex link.
enum int SLIDE_SLOTS = 12;

/// Slide-move configuration + counters, passed to mcmcStep as an optional
/// trailing pointer (null = slides disabled, the default for every existing
/// call site). `prob` is the probability of proposing a slide rather than the
/// ordinary 3->2 bistellar move once the unified proposal has landed on a
/// degree-3 edge.
struct SlideConfig
{
    real prob = 0.0L;
    /// Restrict the move class to CLEAN (species-preserving) slides. OFF by
    /// default: cleanliness is not needed for detailed balance, and excluding
    /// dirty slides throws away ~2/3 of the class. See trySlideMove.
    bool cleanOnly = false;
    ulong tries;      // proposals that formed a legal slide (reached Metropolis)
    ulong accepts;
}

/// Non-local slide channel (dim=3): pick a degree-3 chord uniformly (1/n_3),
/// annihilate it and re-create it a random number of tets down its BC chain.
/// prob = probability of proposing this channel per mcmcStep; maxStep = the
/// upper bound on |k| (drawn uniformly in 1..maxStep; direction is the slot's
/// orientation, so the effective step distribution is symmetric about 0).
struct NonlocalSlideConfig
{
    real prob = 0.0L;
    int maxStep = 8;
    ulong tries;
    ulong accepts;
}

/// Contract/split channel (dim=3): edge contraction (u,v) -> u paired with
/// its inverse, vertex split along a cycle in the link -- the only channel
/// that changes f0 by +-1 in one accepted move (1<->4 needs a degree-4
/// vertex; this needs only the link condition). maxRing caps BOTH deg(uv)
/// on the contract side and |gamma| on the split side -- the pair must be
/// capped together or detailed balance breaks. Proposals: contract = random
/// facet + random edge (prob deg(uv)/(6 f3)); split = random facet + random
/// vertex (prob deg3(w)/(4 f3)), a uniform splitting cycle of length <=
/// maxRing in link(w) (FK-catalog table + transported cycle list when the
/// link is a Z12/Z14/Z15/Z16 polyhedron, bounded DFS otherwise), and a fair
/// coin for which side keeps w. Hastings factors follow ../notes: the
/// reverse of a contraction needs the cycle count of the MERGED link,
/// counted by DFS on the planned (not yet applied) state.
struct ContractSplitConfig
{
    real prob = 0.0L;
    int maxRing = 6;
    ulong contractTries;      // proposals reaching Metropolis
    ulong contractAccepts;
    ulong splitTries;
    ulong splitAccepts;
    ulong noValid;            // failed validity/geometry gates
}

/// A live set of the degree-3 edges (defect chords) for the 1/n_3 proposal:
/// O(1) uniform draw / add / remove, and -- crucially -- a WITNESS facet per
/// chord (a tet containing it) so the channel gets its slide hint in O(1)
/// instead of the O(N) `mfd.link` (which scans all facets via StarRange). A
/// witness can go stale when its tet is destroyed; the draw validates it in
/// O(1) with writeFaceApexes and only recomputes (once, then caches) on a miss.
/// Only touched while the non-local channel is enabled. Sentinel witness =
/// first component < 0 (vertex labels are non-negative).
struct Deg3Set(Vertex)
{
    Vertex[2][] edges;
    Vertex[4][] wit;                 // wit[i] contains edges[i]; [<0,..] = unknown
    size_t[Vertex[2]] index;

    size_t length() const @safe pure nothrow { return edges.length; }

    static bool edgeInFacet(Vertex[2] e, Vertex[4] f)
    {
        bool a0 = false, a1 = false;
        foreach (v; f) { if (v == e[0]) a0 = true; if (v == e[1]) a1 = true; }
        return a0 && a1;
    }

    void add(Vertex[2] e, Vertex[4] w)
    {
        immutable hasW = (w[0] >= 0 && edgeInFacet(e, w));
        auto p = e in index;
        if (p !is null)
        {
            if (hasW && wit[*p][0] < 0) wit[*p] = w;   // upgrade a sentinel
            return;
        }
        index[e] = edges.length;
        edges ~= e;
        if (hasW) wit ~= w;
        else { Vertex[4] s = -1; wit ~= s; }
    }
    void remove(Vertex[2] e)
    {
        auto p = e in index;
        if (p is null) return;
        immutable i = *p;
        immutable last = edges.length - 1;
        if (i != last)
        {
            edges[i] = edges[last]; wit[i] = wit[last];
            index[edges[i]] = i;
        }
        edges.length = last; wit.length = last;
        index.remove(e);
    }
    void reconcile(Vertex[2] e, bool isDeg3, Vertex[4] w)
    {
        if (isDeg3) add(e, w); else remove(e);
    }
    void clear() { edges.length = 0; wit.length = 0; index.clear(); }
}

/// Reconcile the degree-3 set over every edge among `verts` (all edges a move
/// could have changed lie within its support), using only the CURRENT degrees.
/// `witnessFacet` (a length-4 tet, optional) supplies the witness for any
/// newly-degree-3 edge it contains -- pass the facet the move just created so
/// the moved chord starts with a valid O(1) hint.
void reconcileDeg3(Vertex)(ref Manifold!(3, Vertex) mfd,
    ref Deg3Set!Vertex set, scope const(Vertex)[] verts,
    scope const(Vertex)[] witnessFacet = null)
{
    Vertex[4] w = -1;
    if (witnessFacet.length == 4)
        foreach (k; 0 .. 4) w[k] = witnessFacet[k];
    foreach (i; 0 .. verts.length)
        foreach (j; i + 1 .. verts.length)
        {
            Vertex[2] e = [verts[i], verts[j]]; e[].sort();
            set.reconcile(e, mfd.degreeOrZero!1(e[]) == 3, w);
        }
}

/// Populate a degree-3 set from scratch (O(N)); call once when the non-local
/// channel is enabled. Iterates FACETS (not edges) so every initial chord gets
/// a valid witness facet.
void rebuildDeg3(Vertex)(ref Manifold!(3, Vertex) mfd, ref Deg3Set!Vertex set)
{
    set.clear();
    foreach (f; mfd.facets)
    {
        Vertex[4] ft = [f[0], f[1], f[2], f[3]];
        foreach (a; 0 .. 4)
            foreach (b; a + 1 .. 4)
            {
                Vertex[2] e = [ft[a], ft[b]]; e[].sort();
                if (mfd.degreeOrZero!1(e[]) == 3) set.add(e, ft);
            }
    }
}

/// Degree-parameterised twins of reconcileDeg3 / rebuildDeg3.  The Deg3Set
/// struct is degree-agnostic (keyed edges + witness facets); the worm channel
/// keeps a second instance holding the DEGREE-4 edges for its 1/n_4 proposal.
void reconcileDegSet(Vertex)(ref Manifold!(3, Vertex) mfd,
    ref Deg3Set!Vertex set, int wantDeg, scope const(Vertex)[] verts,
    scope const(Vertex)[] witnessFacet = null)
{
    Vertex[4] w = -1;
    if (witnessFacet.length == 4)
        foreach (k; 0 .. 4) w[k] = witnessFacet[k];
    foreach (i; 0 .. verts.length)
        foreach (j; i + 1 .. verts.length)
        {
            Vertex[2] e = [verts[i], verts[j]]; e[].sort();
            set.reconcile(e, mfd.degreeOrZero!1(e[]) == wantDeg, w);
        }
}

/// ditto
void rebuildDegSet(Vertex)(ref Manifold!(3, Vertex) mfd,
    ref Deg3Set!Vertex set, int wantDeg)
{
    set.clear();
    foreach (f; mfd.facets)
    {
        Vertex[4] ft = [f[0], f[1], f[2], f[3]];
        foreach (a; 0 .. 4)
            foreach (b; a + 1 .. 4)
            {
                Vertex[2] e = [ft[a], ft[b]]; e[].sort();
                if (mfd.degreeOrZero!1(e[]) == wantDeg) set.add(e, ft);
            }
    }
}

/// What trySlideMove does once it has a legal clean slide in hand.
enum SlideAccept
{
    metropolis,   ///< the sampler's rule: accept with min(1, exp(-dS))
    trialOnly,    ///< never commit -- measure dS and cleanliness, then restore
    force,        ///< always commit (crossval / scripted transport only)
}

/// One applied Pachner move of a slide, kept for rollback. Public: the
/// FPKMC graph scan (ddg_capi) holds these across recursion for exact
/// rollback via slideRollback.
struct SlideRec(Vertex)
{
    Vertex[3] center;
    Vertex[3] coCenter;
    int centerLen;
    int coCenterLen;
}

/// The apex of triangle (a,b,c) that is not `excl`. Returns false if the
/// triangle is absent or `excl` is not one of its apexes (frame broken).
private bool apexExcluding(Vertex)(const ref Manifold!(3, Vertex) mfd,
    Vertex a, Vertex b, Vertex c, Vertex excl, out Vertex result)
{
    int[2] lk = 0;
    if (mfd.writeFaceApexes(a, b, c, lk.ptr) != 2) return false;
    if (lk[0] == excl) { result = cast(Vertex) lk[1]; return true; }
    if (lk[1] == excl) { result = cast(Vertex) lk[0]; return true; }
    return false;
}

/// The eight frame vertices of a slide, in template order.
struct SlideFrame(Vertex)
{
    Vertex c0, c2, c3, c4, c5, c6, c7, c8;
}

/// Derive c5..c8 from the chord (c0,c4) and the ordered link pick (c2,c3) by
/// the sliding-window rule. Returns false if any window face is missing, the
/// frame assumption breaks, or the eight vertices are not distinct.
/// Why a slide frame could not be derived. `degenerate` is the physically
/// interesting one: the eight frame vertices are not distinct, i.e. the
/// apex-walk folded back on itself. This happens exactly when the chain
/// runs into a degree-4 edge (whose link 4-cycle makes two triangles share
/// the same apex pair) -- the disclination lines of a defect complex. So a
/// `degenerate` count is a count of slides blocked by nearby defect
/// structure, not by a missing/broken window (`missingFace`).
enum SlideFrameCause : int { ok = 0, missingFace = 1, degenerate = 2 }

/// Frame derivation that reports WHY it failed (see SlideFrameCause). The
/// bool `deriveSlideFrame` below is the unchanged behaviour, kept so every
/// existing caller is untouched; this core exists for the slide-frame
/// census (ddg_slide_frame_census).
SlideFrameCause deriveSlideFrameCause(Vertex)(
    const ref Manifold!(3, Vertex) mfd,
    Vertex c0, Vertex c4, Vertex c2, Vertex c3, out SlideFrame!Vertex f)
{
    Vertex c5, c6, c7, c8;
    if (!apexExcluding(mfd, c2, c3, c4, c0, c5))
        return SlideFrameCause.missingFace;
    if (!apexExcluding(mfd, c3, c4, c5, c2, c6))
        return SlideFrameCause.missingFace;
    if (!apexExcluding(mfd, c4, c5, c6, c3, c7))
        return SlideFrameCause.missingFace;
    if (!apexExcluding(mfd, c5, c6, c7, c4, c8))
        return SlideFrameCause.missingFace;

    Vertex[8] all = [c0, c2, c3, c4, c5, c6, c7, c8];
    Vertex[8] sorted_ = all;
    sorted_[].sort();
    foreach (i; 0 .. 7)
        if (sorted_[i] == sorted_[i + 1])
            return SlideFrameCause.degenerate;   // apex-walk folded

    f = SlideFrame!Vertex(c0, c2, c3, c4, c5, c6, c7, c8);
    return SlideFrameCause.ok;
}

bool deriveSlideFrame(Vertex)(const ref Manifold!(3, Vertex) mfd,
    Vertex c0, Vertex c4, Vertex c2, Vertex c3, out SlideFrame!Vertex f)
{
    return deriveSlideFrameCause(mfd, c0, c4, c2, c3, f)
        == SlideFrameCause.ok;
}

/// Is this edge degree illegal (an FK/TCP defect)? Legal degrees are 5 and 6.
private bool isIllegalDegree(size_t d) { return d != 0 && d != 5 && d != 6; }

/******************************************************************************
Attempt one knot slide at the degree-3 chord (a,b) using slot `slot`.

Runs in two passes. The TRIAL pass applies the four Pachner moves for real,
accumulating the action change through the SAME speculative-delta path every
other move uses (there is no parallel reimplementation of the action to drift
out of sync), tests cleanliness, then rolls the composite back exactly. If the
Metropolis test passes, the COMMIT pass replays the recorded moves with the
full instrumentation in lockstep, so an accepted slide touches the potential
state, cocycle, six-flip ledger and event log through exactly the same calls,
in the same order, as four ordinary accepted bistellar moves.

Returns true if the slide was accepted (manifold, currentObjective, potState
and cocycle all advanced); false otherwise, with every structure restored.
`valid` is set when the proposal formed a legal CLEAN slide, i.e. when the
attempt counted as a Metropolis try rather than a malformed proposal.
*/
bool trySlideMove(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    Vertex a, Vertex b, const(Vertex)[] hintTet, int slot,
    P params,
    out bool valid,
    MoveCounters!Vertex* counters = null,
    GeometryLedger!Vertex* ledger = null,
    VertexPotState!Vertex* potState = null,
    const(VertexPot)* pot = null,
    CocycleState!Vertex* cocycle = null,
    SlideAccept policy = SlideAccept.metropolis,
    real* dSOut = null,
    bool cleanOnly = false,
    Deg3Set!Vertex* deg3Set = null)
{
    alias BM = BistellarMove!(3, Vertex);
    valid = false;
    if (dSOut !is null) *dSOut = real.nan;

    // --- the chord's link: 3 vertices, sorted so slot -> (c2,c3) is a
    // state-independent map (this is what makes k_f = k_r = 1 hold).
    int[8] linkBuf = 0;
    auto nl = mfd.writeEdgeLinkCycle(a, b, hintTet, linkBuf.ptr);
    if (nl != 3) return false;
    Vertex[3] link = [cast(Vertex) linkBuf[0], cast(Vertex) linkBuf[1],
                      cast(Vertex) linkBuf[2]];
    link[].sort();

    // --- decode the slot: orientation x ordered (c2,c3) pick.
    immutable int orient = slot / 6;
    immutable int pick = slot % 6;
    immutable Vertex c0 = orient == 0 ? a : b;
    immutable Vertex c4 = orient == 0 ? b : a;
    static immutable int[2][6] picks =
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];
    immutable Vertex c2 = link[picks[pick][0]];
    immutable Vertex c3 = link[picks[pick][1]];

    SlideFrame!Vertex f;
    if (!deriveSlideFrame(mfd, c0, c4, c2, c3, f)) return false;

    // --- support: the frame plus the chord's link (the third link vertex is
    // generally not a frame vertex). Every tet the composite adds or removes
    // has all four vertices in this set, so every edge whose degree can
    // change has both endpoints in it -- that is what makes the local
    // cleanliness test exact.
    Vertex[9] supBuf = 0;
    int nsup = 0;
    void addSup(Vertex v)
    {
        foreach (i; 0 .. nsup) if (supBuf[i] == v) return;
        supBuf[nsup++] = v;
    }
    addSup(f.c0); addSup(f.c2); addSup(f.c3); addSup(f.c4);
    addSup(f.c5); addSup(f.c6); addSup(f.c7); addSup(f.c8);
    foreach (v; link) addSup(v);
    auto support = supBuf[0 .. nsup];

    if (counters !is null)
        addSupport(counters.proposed, support);

    // Frozen-region rejection: as for the other move types, every facet the
    // composite adds or removes has all its vertices in the support.
    if (mfd.anyFrozen(support)) return false;

    // --- degrees of every support edge BEFORE the composite.
    size_t[36] degBefore = 0;
    int npair = 0;
    foreach (i; 0 .. nsup)
        foreach (j; i + 1 .. nsup)
        {
            Vertex[2] e = [supBuf[i], supBuf[j]];
            e[].sort();
            degBefore[npair++] = mfd.degreeOrZero!1(e[]);
        }

    // --- TRIAL PASS: apply the four moves, accumulating the action change
    // through the ordinary speculative-delta path.
    SlideRec!Vertex[4] recs;
    int nApplied = 0;
    real baseRun = currentObjective - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;

    void rollback()
    {
        foreach_reverse (k; 0 .. nApplied)
        {
            auto inv = BM(recs[k].coCenter[0 .. recs[k].coCenterLen],
                          recs[k].center[0 .. recs[k].centerLen]);
            assert(mfd.hasValidMove(inv), "slide rollback: inverse invalid");
            if (potState !is null)
                mfd.potentialBistellarDelta(inv, *potState, *pot, true);
            mfd.doMove(inv);
        }
    }

    /// Form, validate and apply one step. Returns false (after rolling the
    /// whole composite back) if the step is not a legal Pachner move here.
    bool step(scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
    {
        auto bm = BM(center, coCenter);
        if (!mfd.hasValidMove(bm)) { rollback(); return false; }
        real dBase = mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        recs[nApplied].centerLen = cast(int) center.length;
        recs[nApplied].coCenterLen = cast(int) coCenter.length;
        recs[nApplied].center[0 .. center.length] = center[];
        recs[nApplied].coCenter[0 .. coCenter.length] = coCenter[];
        nApplied++;
        return true;
    }

    // M1: 3->2 destroying the chord itself.
    Vertex[2] m1c = [c0, c4]; m1c[].sort();
    if (!step(m1c[], link[])) return false;

    // M2: 2->3 on face (c3,c4,c5) with apexes (c2,c6).
    Vertex[3] m2c = [f.c3, f.c4, f.c5]; m2c[].sort();
    Vertex[2] m2cc = [f.c2, f.c6]; m2cc[].sort();
    if (!step(m2c[], m2cc[])) return false;

    // M3: 2->3 on face (c5,c6,c7) with apexes (c4,c8) -- creates the arrival
    // chord (c4,c8).
    Vertex[3] m3c = [f.c5, f.c6, f.c7]; m3c[].sort();
    Vertex[2] m3cc = [f.c4, f.c8]; m3cc[].sort();
    if (!step(m3c[], m3cc[])) return false;

    // M4: 3->2 on (c2,c6); its link is resolved against the current state.
    Vertex[2] m4c = [f.c2, f.c6]; m4c[].sort();
    {
        auto degc = mfd.degreeOrZero!1(m4c[]);
        if (degc != 3) { rollback(); return false; }
        // any current facet on the edge serves as the walk hint
        Vertex[4] hint = 0;
        {
            int[2] ap = 0;
            if (mfd.writeFaceApexes(f.c2, f.c6, f.c3, ap.ptr) != 2)
            { rollback(); return false; }
            hint = [f.c2, f.c6, f.c3, cast(Vertex) ap[0]];
            hint[].sort();
        }
        int[8] lb4 = 0;
        auto n4 = mfd.writeEdgeLinkCycle(m4c[0], m4c[1], hint[], lb4.ptr);
        if (n4 != 3) { rollback(); return false; }
        Vertex[3] m4cc = [cast(Vertex) lb4[0], cast(Vertex) lb4[1],
                          cast(Vertex) lb4[2]];
        m4cc[].sort();
        if (!step(m4c[], m4cc[])) return false;
    }

    // --- LANDED: the arrival chord must be a degree-3 edge, i.e. the knot
    // came to rest with a usable handle. M3 creates (c4,c8) as a pole-pole
    // edge born at degree 3 and M4's three tets (on the edge (c2,c6), whose
    // link is then {c3,c4,c5}) contain no c8, so this always holds -- it is
    // checked rather than assumed because the whole inverse-closure argument
    // rests on it.
    {
        Vertex[2] arrival = [f.c4, f.c8]; arrival[].sort();
        if (mfd.degreeOrZero!1(arrival[]) != 3) { rollback(); return false; }
    }

    // --- CLEANLINESS (only gates the move when `cleanOnly` is set).
    //
    // A slide is "clean" when the multiset of illegal degrees over CHANGED
    // support edges is identical before and after. Unchanged edges contribute
    // identically to the global multiset, so this local test is exact.
    //
    // Cleanliness is NOT required for detailed balance, and is off by default.
    // The proposal is symmetric without it: a degree-3 chord lies in exactly 3
    // facets, so P(propose chord c) = (3/N_3)(1/6) regardless of the state's
    // defect content, N_3 is slide-invariant, and k_f = k_r = 1 holds for
    // dirty slides just as for clean ones (verified on 124 transitions across
    // four states: 66 dirty, all k_f = k_r = 1, none with k_r = 0). Dirty
    // slides simply cost energy and are turned away by the Metropolis test
    // rather than excluded by fiat, which is both richer and more honest --
    // species preservation becomes statistical instead of imposed.
    //
    // `cleanOnly` is kept because the clean subclass has a special analytic
    // status: it preserves the ENTIRE degree multiset (cleanliness fixes the
    // illegal part; E and sum(deg) = 6*N_3 being fixed then force n5 and n6
    // individually), so under a volume pin plus an edge-degree term alone its
    // dS vanishes identically -- an exactly zero-energy orbit, useful for
    // studying pure transport at fixed species.
    enum int MAXDEG = 64;
    int[MAXDEG + 1] histo = 0;
    bool overflow = false;
    {
        int k = 0;
        foreach (i; 0 .. nsup)
            foreach (j; i + 1 .. nsup)
            {
                Vertex[2] e = [supBuf[i], supBuf[j]];
                e[].sort();
                immutable d0 = degBefore[k++];
                immutable d1 = mfd.degreeOrZero!1(e[]);
                if (d0 == d1) continue;                 // unchanged: cancels
                if (d0 > MAXDEG || d1 > MAXDEG) { overflow = true; break; }
                if (isIllegalDegree(d0)) histo[d0]++;
                if (isIllegalDegree(d1)) histo[d1]--;
            }
    }
    bool clean = !overflow;
    if (clean)
        foreach (h; histo)
            if (h != 0) { clean = false; break; }

    // The trial is over either way: put the manifold back exactly as it was.
    rollback();

    if (cleanOnly && !clean) return false;

    valid = true;
    if (dSOut !is null) *dSOut = deltaTotal;
    if (counters !is null)
        addSupport(counters.valid, support);

    // --- Metropolis. k_f = k_r = 1 and n_3 is preserved on the clean class,
    // so the Hastings ratio is exactly exp(-dS) with no correction.
    final switch (policy)
    {
    case SlideAccept.trialOnly:
        return false;
    case SlideAccept.force:
        break;
    case SlideAccept.metropolis:
        immutable real logAlpha = -deltaTotal;
        if (logAlpha < 0 && uniform01 > exp(logAlpha))
            return false;
        break;
    }

    // --- COMMIT PASS: replay the recorded moves with full instrumentation,
    // exactly as four ordinary accepted bistellar moves.
    if (ledger !is null && ledger.logEvents)
    {
        Vertex[4] lb = [cast(Vertex) CHANNEL_SLIDE, cast(Vertex) 1, a, b];
        logEvent(*ledger, EVT_CHANNEL_BEGIN, lb[], null);
    }
    real committed = 0.0L;
    real baseCommit = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    foreach (k; 0 .. nApplied)
    {
        auto bm = BM(recs[k].center[0 .. recs[k].centerLen],
                     recs[k].coCenter[0 .. recs[k].coCenterLen]);
        assert(mfd.hasValidMove(bm), "slide commit: replayed move invalid");
        real dBase = mfd.speculativeBistellarDelta(bm, baseCommit, params);
        committed += dBase;
        baseCommit += dBase;
        if (potState !is null)
            committed += mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        if (ledger !is null && ledger.logSixFlips)
            sixFlipsBistellar(*ledger, mfd, bm.center, bm.coCenter);
        if (cocycle !is null && cocycle.enabled)
            cocycleBistellar(*cocycle, bm.center, bm.coCenter);
        mfd.doMove(bm);
        if (ledger !is null)
        {
            if (ledger.trackRoles)
                recordBistellar(*ledger, bm.center, bm.coCenter);
            if (ledger.logEvents)
                logEvent(*ledger, cast(int) bm.coCenter.length - 1,
                         bm.center, bm.coCenter);
        }
    }
    assert(abs(committed - deltaTotal) < 1e-9L,
        "slide commit delta disagrees with trial delta");
    if (ledger !is null && ledger.logEvents)
    {
        Vertex[3] le = [cast(Vertex) CHANNEL_SLIDE, cast(Vertex) nApplied,
                        cast(Vertex) dSQuantized(committed)];
        logEvent(*ledger, EVT_CHANNEL_END, le[], null);
    }

    currentObjective += committed;
    if (counters !is null)
        addSupport(counters.acceptedBistellar, support);
    if (deg3Set !is null)
        reconcileDeg3(mfd, *deg3Set, support);
    return true;
}

/******************************************************************************
                    NON-LOCAL SLIDE (undo / re-do worm move)
*******************************************************************************

Move a degree-3 chord by ANNIHILATING it with a 3->2 (undoing its creating
2->3) and RE-CREATING it `steps` tets down the Boerdijk-Coxeter chain with a
fresh 2->3. Just two Pachner moves, passing through the pristine intermediate,
so the exact action change is

    dS = dS(3->2 at source) + dS(2->3 at target)

with no frame derivation and no degree-4 collision (the local 4-move slide's
machinery existed only to cross the endpoint OVERLAP; going through pristine
sidesteps all of it). `steps` = 4 reproduces the local slide's displacement (a
chord that shares one vertex with the original); larger `steps` are the
non-local sampling move -- an excitation that hops directly to a distant
same-rung site, defeating the washboard caging.

USE FOR EQUILIBRIUM SAMPLING. It leaves pi ∝ exp(-objective) invariant, but the
3->2 annihilates the defect to the vacuum, so a distant re-do is re-nucleation,
not motion -- keep `steps` small (=4) for physical kinetics.

`slot` (0 .. SLIDE_SLOTS-1) selects the walk: orient = slot/6 picks which apex
leads the window, pick = slot%6 orders the chord's 3-vertex link. The walk is
nextChainWindow, a reversible permutation, so a matching reverse slot walks
back. Proposal-symmetry / detailed-balance wiring lives in the caller.
*/
bool tryNonlocalSlide(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    Vertex a, Vertex b, const(Vertex)[] hintTet, int slot, int steps,
    P params,
    out bool valid,
    VertexPotState!Vertex* potState = null,
    const(VertexPot)* pot = null,
    SlideAccept policy = SlideAccept.metropolis,
    real* dSOut = null,
    long n3Before = -1, long* dn3Out = null,
    Vertex* outTa = null, Vertex* outTb = null,
    Deg3Set!Vertex* deg3Set = null,
    GeometryLedger!Vertex* ledger = null)
{
    // `n3Before` (the current count of degree-3 edges) is used ONLY in
    // metropolis mode, where the acceptance carries the Hastings factor
    // n3Before/(n3Before + dn3) from the 1/n_3 proposal (pass -1 to omit it).
    // `dn3Out` reports the exact change in the degree-3 edge count (counted
    // over the two moves' O(1) supports); `outTa,outTb` report the arrival
    // chord. Both moves are always applied (to count dn3 exactly) and rolled
    // back on trial/reject.
    import std.math : exp, log;
    alias BM = BistellarMove!(3, Vertex);
    valid = false;
    if (dSOut !is null) *dSOut = real.nan;
    if (dn3Out !is null) *dn3Out = 0;
    if (steps <= 0 || slot < 0 || slot >= SLIDE_SLOTS) return false;

    // count degree-3 edges among all pairs of a small vertex set
    int countDeg3(scope const(Vertex)[] vs)
    {
        int c = 0;
        foreach (i; 0 .. vs.length)
            foreach (j; i + 1 .. vs.length)
            {
                Vertex[2] e = [vs[i], vs[j]]; e[].sort();
                if (mfd.degreeOrZero!1(e[]) == 3) c++;
            }
        return c;
    }

    // link of the degree-3 chord (a,b): 3 vertices in cycle order
    int[8] linkBuf = 0;
    auto nl = mfd.writeEdgeLinkCycle(a, b, hintTet, linkBuf.ptr);
    if (nl != 3) return false;
    Vertex[3] link = [cast(Vertex) linkBuf[0], cast(Vertex) linkBuf[1],
                      cast(Vertex) linkBuf[2]];

    // decode slot: orientation (leading apex) x ordered link pick
    immutable int orient = slot / 6;
    immutable int pick = slot % 6;
    static immutable int[2][6] picks =
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];
    immutable int i0 = picks[pick][0], i1 = picks[pick][1], i2 = 3 - i0 - i1;
    immutable Vertex w0 = orient == 0 ? a : b;
    Vertex[4] window = [w0, link[i0], link[i1], link[i2]];

    Vertex[2] chord = [a, b]; chord[].sort();
    Vertex[3] linkS = link; linkS[].sort();

    {
        Vertex[5] sup = [chord[0], chord[1], link[0], link[1], link[2]];
        if (mfd.anyFrozen(sup[])) return false;
    }

    // --- annihilate: 3->2 on (a,b), coCenter = its link ---
    BM annM = BM(chord[], linkS[]);
    if (!mfd.hasValidMove(annM)) return false;
    Vertex[5] annSup = [chord[0], chord[1], link[0], link[1], link[2]];
    immutable int cAnnBefore = countDeg3(annSup[]);
    real baseRun = currentObjective - (potState !is null ? potState.total : 0.0L);
    immutable real dAnnBase = mfd.speculativeBistellarDelta(annM, baseRun, params);
    real dAnnPot = 0.0L;
    if (potState !is null)
        dAnnPot = mfd.potentialBistellarDelta(annM, *potState, *pot, true);
    mfd.doMove(annM);                 // manifold now pristine locally
    baseRun += dAnnBase;
    immutable int dn3Ann = countDeg3(annSup[]) - cAnnBefore;

    // undo the annihilation (recreate the source chord) and return false
    bool undoAnnAndFail()
    {
        BM inv = BM(linkS[], chord[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(inv, *potState, *pot, true);
        mfd.doMove(inv);
        return false;
    }

    // the start window is now a facet (tet_n or tet_{n+1})
    {
        Vertex[4] ws = window; ws[].sort();
        if (!mfd.contains(ws[])) return undoAnnAndFail();
    }

    // --- walk `steps` tets down the chain ---
    Vertex[4] w = window;
    foreach (_; 0 .. steps)
        if (!mfd.nextChainWindow(w)) return undoAnnAndFail();

    // target chord = (w[0], apex of face (w1,w2,w3) opposite w[0])
    int[2] ap = 0;
    if (mfd.writeFaceApexes(w[1], w[2], w[3], ap.ptr) != 2)
        return undoAnnAndFail();
    Vertex nxt;
    if (ap[0] == w[0]) nxt = cast(Vertex) ap[1];
    else if (ap[1] == w[0]) nxt = cast(Vertex) ap[0];
    else return undoAnnAndFail();
    Vertex[3] tFace = [w[1], w[2], w[3]]; tFace[].sort();
    Vertex[2] tChord = [w[0], nxt]; tChord[].sort();

    {
        Vertex[5] sup = [tFace[0], tFace[1], tFace[2], w[0], nxt];
        if (mfd.anyFrozen(sup[])) return undoAnnAndFail();
    }

    // --- re-do: 2->3 on the target face, coCenter = target chord ---
    BM redoM = BM(tFace[], tChord[]);
    if (!mfd.hasValidMove(redoM)) return undoAnnAndFail();
    Vertex[5] redoSup = [tFace[0], tFace[1], tFace[2], w[0], nxt];
    immutable int cRedoBefore = countDeg3(redoSup[]);
    immutable real dRedoBase = mfd.speculativeBistellarDelta(redoM, baseRun, params);
    real dRedoPot = 0.0L;
    if (potState !is null)
        dRedoPot = mfd.potentialBistellarDelta(redoM, *potState, *pot, true);
    mfd.doMove(redoM);                // defect now at the target
    immutable int dn3Redo = countDeg3(redoSup[]) - cRedoBefore;

    immutable real dS = dAnnBase + dAnnPot + dRedoBase + dRedoPot;
    immutable long dn3 = cast(long)(dn3Ann + dn3Redo);
    valid = true;
    if (dSOut !is null) *dSOut = dS;
    if (dn3Out !is null) *dn3Out = dn3;
    if (outTa !is null) *outTa = tChord[0];
    if (outTb !is null) *outTb = tChord[1];

    // undo BOTH moves (both are applied by now): 3->2 the target chord, then
    // 2->3 recreate the source. Restores the state exactly.
    bool undoBothAndFail()
    {
        BM invRedo = BM(tChord[], tFace[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(invRedo, *potState, *pot, true);
        mfd.doMove(invRedo);
        BM invAnn = BM(linkS[], chord[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(invAnn, *potState, *pot, true);
        mfd.doMove(invAnn);
        return false;
    }

    final switch (policy)
    {
    case SlideAccept.trialOnly:
        return undoBothAndFail();
    case SlideAccept.force:
        break;
    case SlideAccept.metropolis:
        real logAlpha = -dS;
        if (n3Before > 0 && n3Before + dn3 > 0)   // 1/n_3 Hastings factor
            logAlpha += log(cast(real) n3Before / cast(real)(n3Before + dn3));
        if (logAlpha < 0 && uniform01 > exp(logAlpha))
            return undoBothAndFail();
        break;
    }

    // --- accepted: leave the defect at the target ---
    currentObjective += dS;
    // Mirror both Pachner moves into the ledger (post-move-safe records
    // only, in applied order), bracketed as one composite so trackers keep
    // worldline identity across the hop. sixFlipsBistellar needs pre-move
    // state this apply-then-decide structure cannot provide, so the channel
    // is gated off under logSixFlips in mcmcStep (as for the worm).
    if (ledger !is null)
    {
        if (ledger.logEvents)
        {
            Vertex[6] lb = [cast(Vertex) CHANNEL_NONLOCAL, cast(Vertex) 1,
                            a, b, tChord[0], tChord[1]];
            logEvent(*ledger, EVT_CHANNEL_BEGIN, lb[0 .. 4], lb[4 .. 6]);
        }
        BM annRec = BM(chord[], linkS[]);
        BM redoRec = BM(tFace[], tChord[]);
        void mirror(ref BM bm)
        {
            if (ledger.trackRoles)
                recordBistellar(*ledger, bm.center, bm.coCenter);
            if (ledger.logEvents)
                logEvent(*ledger, cast(int) bm.coCenter.length - 1,
                         bm.center, bm.coCenter);
        }
        mirror(annRec);
        mirror(redoRec);
        if (ledger.logEvents)
        {
            Vertex[3] ge = [cast(Vertex) ROLE_GONE3, chord[0], chord[1]];
            logEvent(*ledger, EVT_CHANNEL_EDGE, ge[], null);
            Vertex[3] ne = [cast(Vertex) ROLE_NEW3, tChord[0], tChord[1]];
            logEvent(*ledger, EVT_CHANNEL_EDGE, ne[], null);
            Vertex[3] le = [cast(Vertex) CHANNEL_NONLOCAL, cast(Vertex) 2,
                            cast(Vertex) dSQuantized(dS)];
            logEvent(*ledger, EVT_CHANNEL_END, le[], null);
        }
    }
    if (deg3Set !is null)
    {
        reconcileDeg3(mfd, *deg3Set, annSup[]);    // source healed to pristine
        // The 2->3 created three tets on the target chord; one of them is a
        // valid O(1) witness for the freshly-placed chord.
        Vertex[4] wt = [tChord[0], tChord[1], tFace[0], tFace[1]];
        reconcileDeg3(mfd, *deg3Set, redoSup[], wt[]);
    }
    return true;
}

/*******************************************************************************
                      DEG-4 WORM MOVE (catalysed transport step)
*******************************************************************************

One proposal = one 2->3 + one 3->2 (either order), f-vector neutral, moving
DEG-4 content out of its hinge-flip cavity.  Ported from the validated Python
oracle scripts/defect_dynamics/worm_deg4_slide.py (move class, inverse
closure 69/69, MH integration exact to 1e-14).

MOVE CLASS at an anchor deg-4 edge e:
  * e's final degree is legal (5/6) or the edge is removed;
  * k in {1,2} deg-4 edges lose deg-4 status and k gain it, everything
    balanced; the anchor is among the losers;
  * at most one deg-3 edge removed and one created (catalyst relocation);
  * every OTHER edge's illegal status is exactly unchanged (degree map,
    not set);
  * some new deg-4 f ("landing") lies OUTSIDE e's octahedron (e + its
    4-cycle) in the start state, with e outside f's octahedron in the end
    state (escape symmetry: the reverse is in the class).

PROPOSAL: anchor uniform over the live deg-4 set (1/n_4; the class preserves
n_4 so this factor cancels), then uniform over the enumerated candidates.
Because k = 2 transitions are proposable from either lost anchor, the
Hastings weight is the anchor sum q = (1/n_4) sum_a k(a)/n(a), computed by
enumeration at every lost (forward) / gained (reverse) anchor.  Acceptance:

    alpha = min(1, exp(-dS) * q_r / q_f).

dS comes from speculativeBistellarDelta + potentialBistellarDelta applied in
lockstep on the chosen pair only (enumeration is manifold-only).  Cocycle
tracking is NOT supported: the channel must be disabled when a cocycle is
attached.
*/

/// Worm channel configuration + counters (see mcmcStep).
struct WormConfig
{
    real prob = 0.0L;         ///< probability of proposing a worm per step
    ulong tries;              ///< proposals with >= 1 candidate
    ulong accepts;
    ulong noCands;            ///< proposals rejected for lack of candidates
}

enum int WORM_MAX_CANDS = 512;
private enum int WORM_MAX_NET = 10;

/// One enumerated worm candidate.
struct WormCand(Vertex)
{
    bool m1is23;              ///< m1 kind (m2 is the opposite kind)
    Vertex[3] f1;             ///< m1 face (2->3 center / 3->2 coCenter)
    Vertex[2] p1;             ///< m1 pair (2->3 coCenter / 3->2 center)
    Vertex[3] f2;
    Vertex[2] p2;
    int nPair;                ///< 1 or 2
    Vertex[2][2] gone4;       ///< deg-4 edges losing deg-4 status (anchor incl.)
    Vertex[2][2] new4;        ///< deg-4 edges gaining it
    Vertex[2] landing;        ///< the escaped new deg-4
    int nNet;                 ///< canonical net facet diff (end-state key)
    Vertex[4][WORM_MAX_NET] netT;
    byte[WORM_MAX_NET] netSg;
}

/// Two candidates lead to the same end state iff their net diffs agree.
private bool sameKey(Vertex)(ref const WormCand!Vertex a,
                             ref const WormCand!Vertex b)
{
    if (a.nNet != b.nNet) return false;
    foreach (i; 0 .. a.nNet)
        if (a.netT[i] != b.netT[i] || a.netSg[i] != b.netSg[i]) return false;
    return true;
}

/// Does candidate r's net diff equal the INVERSE of c's?
private bool inverseKey(Vertex)(ref const WormCand!Vertex c,
                                ref const WormCand!Vertex r)
{
    if (c.nNet != r.nNet) return false;
    foreach (i; 0 .. c.nNet)
    {
        bool found = false;
        foreach (j; 0 .. r.nNet)
            if (c.netT[i] == r.netT[j] && c.netSg[i] == -r.netSg[j])
            { found = true; break; }
        if (!found) return false;
    }
    return true;
}

/// Accumulate one move's facet changes into a (tet, sign) list with
/// cancellation; canonicalised by insertion (tets sorted internally).
private void accNet(Vertex)(ref Vertex[4][WORM_MAX_NET] ts,
    ref byte[WORM_MAX_NET] sg, ref int n, Vertex[4] t, byte s)
{
    t[].sort();
    foreach (i; 0 .. n)
        if (ts[i] == t)
        {
            sg[i] += s;
            if (sg[i] == 0)          // cancelled: remove entry
            {
                foreach (j; i .. n - 1) { ts[j] = ts[j + 1]; sg[j] = sg[j + 1]; }
                n--;
            }
            return;
        }
    ts[n] = t; sg[n] = s; n++;
}

/// Sort the net diff so equal end states compare equal structurally.
private void sortNet(Vertex)(ref Vertex[4][WORM_MAX_NET] ts,
    ref byte[WORM_MAX_NET] sg, int n)
{
    foreach (i; 0 .. n)
        foreach (j; i + 1 .. n)
            if (ts[j] < ts[i] || (ts[j] == ts[i] && sg[j] < sg[i]))
            {
                auto tt = ts[i]; ts[i] = ts[j]; ts[j] = tt;
                auto ss = sg[i]; sg[i] = sg[j]; sg[j] = ss;
            }
}

/// Find some facet containing edge (a,b): the O(N) fallback for hint
/// derivation off the hot path (stale witnesses, partner/reverse anchors).
bool findTetOfEdge(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex a, Vertex b, ref Vertex[4] outTet)
{
    foreach (f; mfd.facets)
    {
        bool ha = false, hb = false;
        foreach (v; f) { if (v == a) ha = true; if (v == b) hb = true; }
        if (ha && hb)
        {
            foreach (k; 0 .. 4) outTet[k] = f[k];
            return true;
        }
    }
    return false;
}

/// Append (dedup) every facet containing `v`, walking face-neighbours from
/// `seed` (a facet containing v).  The star of a vertex is a ball, so the
/// face-neighbour walk restricted to facets containing v covers it.  `arr`
/// may already hold other vertices' stars: the walk uses its own per-call
/// visited set (indices into arr), so shared facets do not block it.
/// collectStar must return the SAME SET for every valid seed, and that
/// set must be star(v).
///
/// This is what makes Manifold.someFacetContaining safe to substitute for
/// `star(v).front`: facets() sorts, so star(v).front was the
/// lexicographically least tet on v, while the witness is whichever was
/// inserted most recently. The seed feeds collectStar, so the star comes
/// back in a different ORDER, and the head kernel's uniform draw over the
/// candidate list therefore lands on a different candidate for the same
/// RNG stream -- observed as 96 -> 116 commits over 3000 episodes. That
/// is a relabelling of the proposal, not a change to it: the candidate
/// SET and its COUNT are what the uniform draw and the Hastings ratio
/// (nH0/nH) depend on, and this test pins both.
unittest
{
    import std.algorithm : sort, canFind;
    import manifold_examples : standardSphere;
    auto m = standardSphere!3;
    // churn a little so vertices have non-trivial stars
    int fresh = 0;
    foreach (sv; m.simplices(0))
        if (sv[0] >= fresh) fresh = sv[0] + 1;
    auto f0 = m.facets.front.dup;
    m.doMove(BistellarMove!(3, int)(f0, [fresh]));

    foreach (sv; m.simplices(0))
    {
        immutable v = sv[0];
        // every tet on v, as sorted quadruples, from star()
        int[4][] want;
        foreach (f; m.star(v.only))
        {
            int[4] t = [f[0], f[1], f[2], f[3]];
            t[].sort();
            want ~= t;
        }
        want.sort();
        assert(want.length > 0);
        // collectStar from EVERY valid seed must reproduce exactly that
        foreach (seed; want)
        {
            int[4][64] arr;
            immutable n = collectStar(m, v, seed, arr[], 0);
            assert(n == cast(int) want.length,
                "collectStar size depends on the seed");
            int[4][] got;
            foreach (i; 0 .. n) { auto t = arr[i]; t[].sort(); got ~= t; }
            got.sort();
            assert(got == want, "collectStar set depends on the seed");
        }
    }
}

private int collectStar(Vertex)(ref Manifold!(3, Vertex) mfd, Vertex v,
    Vertex[4] seed, Vertex[4][] arr, int n)
{
    int findOrAppend(Vertex[4] t)
    {
        foreach (i; 0 .. n) if (arr[i] == t) return i;
        if (n < arr.length) { arr[n] = t; return n++; }
        return -1;
    }
    seed[].sort();
    immutable int s0 = findOrAppend(seed);
    if (s0 < 0) return n;
    int[512] queue;
    int[512] visited;
    int qh = 0, qt = 0, nv = 0;
    bool seen(int idx)
    {
        foreach (i; 0 .. nv) if (visited[i] == idx) return true;
        return false;
    }
    void visit(int idx)
    {
        if (idx < 0 || seen(idx) || nv >= visited.length
            || qt >= queue.length) return;
        visited[nv++] = idx;
        queue[qt++] = idx;
    }
    visit(s0);
    while (qh < qt)
    {
        auto T = arr[queue[qh++]];
        // 3 faces of T containing v
        Vertex[3] oth;
        int no = 0;
        foreach (x; T) if (x != v) oth[no++] = x;
        if (no != 3) continue;                    // T does not contain v
        foreach (skip; 0 .. 3)
        {
            Vertex[3] face;
            int nf = 0;
            face[nf++] = v;
            foreach (k; 0 .. 3) if (k != skip) face[nf++] = oth[k];
            int[2] ap = 0;
            if (mfd.writeFaceApexes(face[0], face[1], face[2], ap.ptr) != 2)
                continue;
            Vertex fourth = -1;
            foreach (x; T)
            { if (x != face[0] && x != face[1] && x != face[2]) fourth = x; }
            Vertex nb = (cast(Vertex) ap[0] == fourth) ? cast(Vertex) ap[1]
                                                       : cast(Vertex) ap[0];
            Vertex[4] N = [face[0], face[1], face[2], nb];
            N[].sort();
            visit(findOrAppend(N));
        }
    }
    return n;
}

/// Enumerate every worm-class candidate at `anchor` (a deg-4 edge; `hintTet`
/// = a facet containing it).  Trial moves are applied to the manifold and
/// rolled back exactly; no objective/potential state is touched.  Returns
/// the candidate count written into `out_`.
int wormEnumerate(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] anchor, const(Vertex)[] hintTet, WormCand!Vertex[] out_)
{
    alias BM = BistellarMove!(3, Vertex);
    int nOut = 0;
    if (mfd.degreeOrZero!1(anchor[]) != 4) return 0;

    // anchor octahedron: the edge + its 4-cycle
    int[8] cycBuf = 0;
    if (mfd.writeEdgeLinkCycle(anchor[0], anchor[1], hintTet, cycBuf.ptr) != 4)
        return 0;
    Vertex[6] octaA = [anchor[0], anchor[1], cast(Vertex) cycBuf[0],
        cast(Vertex) cycBuf[1], cast(Vertex) cycBuf[2], cast(Vertex) cycBuf[3]];
    bool inOctaA(Vertex x)
    {
        foreach (v; octaA) if (v == x) return true;
        return false;
    }

    // stage-1 star: tets around anchor[0] and anchor[1]
    Vertex[4][160] star1;
    int n1 = 0;
    foreach (i; 0 .. 4)
    {
        Vertex[4] seed = [anchor[0], anchor[1], cast(Vertex) cycBuf[i],
                          cast(Vertex) cycBuf[(i + 1) % 4]];
        n1 = collectStar(mfd, anchor[0], seed, star1[], n1);
    }
    {
        Vertex[4] seed = [anchor[0], anchor[1], cast(Vertex) cycBuf[0],
                          cast(Vertex) cycBuf[1]];
        n1 = collectStar(mfd, anchor[1], seed, star1[], n1);
    }

    // ---- candidate first moves --------------------------------------------
    static struct Mv
    {
        bool is23;
        Vertex[3] f;      // face (2->3 center / 3->2 link)
        Vertex[2] p;      // pair (2->3 apexes / 3->2 edge)
    }
    Mv[384] m1s;
    int nM1 = 0;

    void gatherMoves(scope Vertex[4][] tets, int nT,
                     scope bool delegate(scope const(Vertex)[]) touches,
                     ref Mv[384] outM, ref int nM)
    {
        // 2->3 on faces of the tets
        Vertex[3][2048] seenF;
        int nSF = 0;
        foreach (ti; 0 .. nT)
        {
            auto T = tets[ti];
            foreach (skip; 0 .. 4)
            {
                Vertex[3] face;
                int k = 0;
                foreach (i; 0 .. 4) if (i != skip) face[k++] = T[i];
                bool dup = false;
                foreach (i; 0 .. nSF) if (seenF[i] == face) { dup = true; break; }
                if (dup) continue;
                if (nSF < seenF.length) seenF[nSF++] = face;
                int[2] ap = 0;
                if (mfd.writeFaceApexes(face[0], face[1], face[2], ap.ptr) != 2)
                    continue;
                Vertex[2] axis = [cast(Vertex) ap[0], cast(Vertex) ap[1]];
                axis[].sort();
                Vertex[5] sup = [face[0], face[1], face[2], axis[0], axis[1]];
                if (!touches(sup[])) continue;
                if (mfd.degreeOrZero!1(axis[]) != 0) continue;
                if (mfd.anyFrozen(sup[])) continue;
                if (nM < outM.length)
                {
                    outM[nM].is23 = true;
                    outM[nM].f = face;
                    outM[nM].p = axis;
                    nM++;
                }
            }
        }
        // 3->2 on degree-3 edges among the tets' pairs
        Vertex[2][1024] seenE;
        int nSE = 0;
        foreach (ti; 0 .. nT)
        {
            auto T = tets[ti];
            foreach (a; 0 .. 4)
                foreach (b; a + 1 .. 4)
                {
                    Vertex[2] e = [T[a], T[b]]; e[].sort();
                    bool dup = false;
                    foreach (i; 0 .. nSE) if (seenE[i] == e) { dup = true; break; }
                    if (dup) continue;
                    if (nSE < seenE.length) seenE[nSE++] = e;
                    if (mfd.degreeOrZero!1(e[]) != 3) continue;
                    int[8] lk = 0;
                    if (mfd.writeEdgeLinkCycle(e[0], e[1], T[], lk.ptr) != 3)
                        continue;
                    Vertex[3] link = [cast(Vertex) lk[0], cast(Vertex) lk[1],
                                      cast(Vertex) lk[2]];
                    link[].sort();
                    Vertex[5] sup = [e[0], e[1], link[0], link[1], link[2]];
                    if (!touches(sup[])) continue;
                    if (mfd.anyFrozen(sup[])) continue;
                    BM probe = BM(e[], link[]);
                    if (!mfd.hasValidMove(probe)) continue;
                    if (nM < outM.length)
                    {
                        outM[nM].is23 = false;
                        outM[nM].f = link;
                        outM[nM].p = e;
                        nM++;
                    }
                }
        }
    }

    bool touchesAnchor(scope const(Vertex)[] sup)
    {
        foreach (v; sup) if (v == anchor[0] || v == anchor[1]) return true;
        return false;
    }
    gatherMoves(star1[], n1, &touchesAnchor, m1s, nM1);

    // ---- try each (m1, m2) pair -------------------------------------------
    foreach (i1; 0 .. nM1)
    {
        auto m1 = m1s[i1];
        BM bm1 = m1.is23 ? BM(m1.f[], m1.p[]) : BM(m1.p[], m1.f[]);
        if (!mfd.hasValidMove(bm1)) continue;

        // snapshot original degrees over m1's support pairs.  Fixed-size
        // parallel arrays, NOT an AA: this loop runs ~70x per enumeration
        // and per-iteration AA allocation churn was leaking the heap out
        // from under long campaigns.  Capacity bound is exact: 10 m1 pairs
        // + anchor + 10 m2 extension pairs = 21.
        Vertex[2][21] origE;
        size_t[21] origD;
        int nOrig = 0;
        int origFind(ref const Vertex[2] e)
        {
            foreach (i; 0 .. nOrig) if (origE[i] == e) return i;
            return -1;
        }
        void origAdd(Vertex[2] e)
        {
            if (origFind(e) >= 0) return;
            assert(nOrig < origE.length);
            origE[nOrig] = e;
            origD[nOrig] = mfd.degreeOrZero!1(e[]);
            nOrig++;
        }
        Vertex[5] sup1 = [m1.f[0], m1.f[1], m1.f[2], m1.p[0], m1.p[1]];
        foreach (i; 0 .. 5)
            foreach (j; i + 1 .. 5)
            {
                Vertex[2] e = [sup1[i], sup1[j]]; e[].sort();
                origAdd(e);
            }
        {
            Vertex[2] e = anchor; e[].sort();
            origAdd(e);
        }
        mfd.doMove(bm1);

        // disturbed vertices: anchor + endpoints of changed pairs
        Vertex[16] dist;
        int nD = 0;
        void addD(Vertex v)
        {
            foreach (i; 0 .. nD) if (dist[i] == v) return;
            if (nD < dist.length) dist[nD++] = v;
        }
        addD(anchor[0]); addD(anchor[1]);
        foreach (oi; 0 .. nOrig)
            if (mfd.degreeOrZero!1(origE[oi][]) != origD[oi])
            { addD(origE[oi][0]); addD(origE[oi][1]); }

        // stage-2 stars around the disturbed vertices
        Vertex[4][768] star2;
        int n2 = 0;
        // seed tets: added tets of m1
        Vertex[4][3] added1;
        int nAdd1 = 0;
        if (m1.is23)
        {
            added1[nAdd1++] = [m1.p[0], m1.p[1], m1.f[0], m1.f[1]];
            added1[nAdd1++] = [m1.p[0], m1.p[1], m1.f[1], m1.f[2]];
            added1[nAdd1++] = [m1.p[0], m1.p[1], m1.f[0], m1.f[2]];
        }
        else
        {
            added1[nAdd1++] = [m1.f[0], m1.f[1], m1.f[2], m1.p[0]];
            added1[nAdd1++] = [m1.f[0], m1.f[1], m1.f[2], m1.p[1]];
        }
        foreach (di; 0 .. nD)
        {
            Vertex v = dist[di];
            Vertex[4] seed = -1;
            bool got = false;
            foreach (k; 0 .. nAdd1)
            {
                foreach (x; added1[k]) if (x == v) { seed = added1[k]; got = true; }
                if (got) break;
            }
            if (!got)
            {
                foreach (ti; 0 .. n1)
                {
                    bool hasV = false;
                    foreach (x; star1[ti]) if (x == v) hasV = true;
                    if (!hasV) continue;
                    Vertex[4] c = star1[ti]; c[].sort();
                    if (mfd.contains(c[])) { seed = star1[ti]; got = true; break; }
                }
            }
            if (got)
                n2 = collectStar(mfd, v, seed, star2[], n2);
        }

        bool touchesDist(scope const(Vertex)[] sup)
        {
            foreach (v; sup)
                foreach (i; 0 .. nD) if (dist[i] == v) return true;
            return false;
        }
        Mv[384] m2s;
        int nM2 = 0;
        gatherMoves(star2[], n2, &touchesDist, m2s, nM2);

        foreach (i2; 0 .. nM2)
        {
            auto m2 = m2s[i2];
            if (m2.is23 == m1.is23) continue;        // one 2->3 + one 3->2
            BM bm2 = m2.is23 ? BM(m2.f[], m2.p[]) : BM(m2.p[], m2.f[]);
            if (!mfd.hasValidMove(bm2)) continue;

            // extend the snapshot with m2's support pairs (their current
            // degree equals the original unless already recorded)
            immutable int nOrig1 = nOrig;
            Vertex[5] sup2 = [m2.f[0], m2.f[1], m2.f[2], m2.p[0], m2.p[1]];
            foreach (i; 0 .. 5)
                foreach (j; i + 1 .. 5)
                {
                    Vertex[2] e = [sup2[i], sup2[j]]; e[].sort();
                    origAdd(e);
                }
            mfd.doMove(bm2);

            // ---- class check ----------------------------------------------
            Vertex[2][2] g4, w4;
            int nG4 = 0, nW4 = 0;
            Vertex[2] g3 = -1, w3 = -1;
            int nG3 = 0, nW3 = 0;
            bool ok = true;
            foreach (oi; 0 .. nOrig)
            {
                Vertex[2] e = origE[oi];
                immutable size_t d0 = origD[oi];
                immutable dN = mfd.degreeOrZero!1(e[]);
                if (dN == d0) continue;
                immutable bool was4 = d0 == 4, is4 = dN == 4;
                immutable bool was3 = d0 == 3, is3 = dN == 3;
                if (was4 && !is4)
                {
                    if (nG4 >= 2) { ok = false; break; }
                    g4[nG4++] = e;
                    if (!(dN == 0 || dN == 5 || dN == 6) && !is3)
                    { ok = false; break; }
                    if (is3)                     // deg-4 dropping to 3 stays
                    {                            // illegal: allowed only as
                        if (nW3 >= 1) { ok = false; break; }   // the catalyst
                        w3 = e; nW3++;           // birth (counts as new deg-3)
                    }
                }
                else if (!was4 && is4)
                {
                    if (nW4 >= 2) { ok = false; break; }
                    w4[nW4++] = e;
                    if (was3)                    // deg-3 promoted to deg-4 =
                    {                            // the catalyst was consumed
                        if (nG3 >= 1) { ok = false; break; }
                        g3 = e; nG3++;
                    }
                    else if (!(d0 == 0 || d0 == 5 || d0 == 6))
                    { ok = false; break; }
                }
                else if (was3 && !is3)
                {
                    if (nG3 >= 1) { ok = false; break; }
                    g3 = e; nG3++;
                    if (!(dN == 0 || dN == 5 || dN == 6)) { ok = false; break; }
                }
                else if (!was3 && is3)
                {
                    if (nW3 >= 1) { ok = false; break; }
                    w3 = e; nW3++;
                    if (!(d0 == 0 || d0 == 5 || d0 == 6)) { ok = false; break; }
                }
                else
                {
                    // legal-to-legal shifts are free; anything else fails
                    immutable bool leg0 = d0 == 5 || d0 == 6 || d0 == 0;
                    immutable bool leg1 = dN == 5 || dN == 6 || dN == 0;
                    if (!(leg0 && leg1)) { ok = false; break; }
                }
            }
            if (ok)
            {
                // anchor must be among the losers, balances must hold
                immutable aDeg = mfd.degreeOrZero!1(anchor[]);
                bool anchorGone = false;
                foreach (i; 0 .. nG4)
                    if (g4[i] == anchor) anchorGone = true;
                ok = anchorGone && (aDeg == 0 || aDeg == 5 || aDeg == 6)
                    && nG4 == nW4 && nG4 >= 1 && nG3 == nW3;
            }
            Vertex[2] landing = -1;
            if (ok)
            {
                ok = false;
                foreach (i; 0 .. nW4)
                {
                    auto f = w4[i];
                    if (inOctaA(f[0]) && inOctaA(f[1])) continue;   // in-cavity
                    // anchor must be outside f's octahedron in the end state
                    int inside = 0;
                    foreach (av; anchor)
                    {
                        if (av == f[0] || av == f[1]) { inside++; continue; }
                        Vertex[3] fc = [f[0], f[1], av]; fc[].sort();
                        if (mfd.degreeOrZero!2(fc[]) > 0) inside++;
                    }
                    if (inside < 2) { landing = f; ok = true; break; }
                }
            }
            if (ok && nOut < out_.length)
            {
                auto c = &out_[nOut];
                c.m1is23 = m1.is23;
                c.f1 = m1.f; c.p1 = m1.p;
                c.f2 = m2.f; c.p2 = m2.p;
                c.nPair = nG4;
                foreach (i; 0 .. nG4) c.gone4[i] = g4[i];
                foreach (i; 0 .. nW4) c.new4[i] = w4[i];
                c.landing = landing;
                c.nNet = 0;
                void accMove(bool is23, Vertex[3] f, Vertex[2] p)
                {
                    if (is23)
                    {
                        accNet(c.netT, c.netSg, c.nNet,
                               [f[0], f[1], f[2], p[0]], cast(byte) -1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [f[0], f[1], f[2], p[1]], cast(byte) -1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[0], f[1]], cast(byte) 1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[1], f[2]], cast(byte) 1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[0], f[2]], cast(byte) 1);
                    }
                    else
                    {
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[0], f[1]], cast(byte) -1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[1], f[2]], cast(byte) -1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [p[0], p[1], f[0], f[2]], cast(byte) -1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [f[0], f[1], f[2], p[0]], cast(byte) 1);
                        accNet(c.netT, c.netSg, c.nNet,
                               [f[0], f[1], f[2], p[1]], cast(byte) 1);
                    }
                }
                accMove(m1.is23, m1.f, m1.p);
                accMove(m2.is23, m2.f, m2.p);
                sortNet(c.netT, c.netSg, c.nNet);
                nOut++;
            }

            // rollback m2 and forget its snapshot extensions
            BM inv2 = m2.is23 ? BM(m2.p[], m2.f[]) : BM(m2.f[], m2.p[]);
            mfd.doMove(inv2);
            nOrig = nOrig1;
        }
        BM inv1 = m1.is23 ? BM(m1.p[], m1.f[]) : BM(m1.f[], m1.p[]);
        mfd.doMove(inv1);
    }
    return nOut;
}

/*
Attempt one deg-4 worm move at `anchor`.  Enumerates the class candidates,
draws one uniformly (or `forceCand` for crossval), computes the anchor-sum
Hastings weights by enumeration at every forward (gone4) and reverse (new4)
anchor, applies the pair through the speculative + potential machinery for
the exact dS, and Metropolis-accepts on exp(-dS) * q_r/q_f.  A proposal that
does not commit restores the manifold exactly.

`valid` = a candidate existed and the Hastings machinery completed (whether
or not the move was accepted).  Cocycle tracking is not supported -- the
caller must gate the channel off when a cocycle is attached.
*/
bool tryWormMove(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    Vertex[2] anchor, const(Vertex)[] hintTet,
    P params,
    out bool valid,
    VertexPotState!Vertex* potState = null,
    const(VertexPot)* pot = null,
    Deg3Set!Vertex* deg4Set = null,
    Deg3Set!Vertex* deg3Set = null,
    SlideAccept policy = SlideAccept.metropolis,
    real* dSOut = null,
    int forceCand = -1,
    Vertex* outLa = null, Vertex* outLb = null,
    GeometryLedger!Vertex* ledger = null)
{
    import std.math : exp, log;
    alias BM = BistellarMove!(3, Vertex);
    valid = false;
    if (dSOut !is null) *dSOut = real.nan;

    // Scratch candidate buffers: allocated once per thread and reused, so a
    // long campaign's proposal stream produces no per-call GC garbage (the
    // ~150 KB/call churn here was leaking the heap over multi-hour runs).
    static WormCand!Vertex[] candsBuf, pcBuf, rcBuf;
    if (candsBuf.length == 0)
    {
        candsBuf = new WormCand!Vertex[](WORM_MAX_CANDS);
        pcBuf = new WormCand!Vertex[](WORM_MAX_CANDS);
        rcBuf = new WormCand!Vertex[](WORM_MAX_CANDS);
    }
    auto cands = candsBuf;
    immutable int nf = wormEnumerate(mfd, anchor, hintTet, cands);
    if (nf == 0) return false;
    immutable int ci = forceCand >= 0 ? forceCand : cast(int) uniform(0, nf);
    if (ci >= nf) return false;
    auto chosen = cands[ci];

    // ---- forward proposal weight: sum over the gone4 anchors -------------
    real qf = 0.0L;
    {
        int kf = 0;
        foreach (i; 0 .. nf)
            if (sameKey(chosen, cands[i])) kf++;
        qf += cast(real) kf / cast(real) nf;
    }
    if (chosen.nPair == 2)
    {
        Vertex[2] partner = chosen.gone4[0] == anchor ? chosen.gone4[1]
                                                      : chosen.gone4[0];
        Vertex[4] pw = -1;
        bool gotW = false;
        if (deg4Set !is null)
        {
            auto p = partner in deg4Set.index;
            if (p !is null)
            {
                pw = deg4Set.wit[*p];
                if (pw[0] >= 0 && Deg3Set!Vertex.edgeInFacet(partner, pw))
                {
                    Vertex[4] ps = pw; ps[].sort();
                    gotW = mfd.contains(ps[]);
                }
            }
        }
        if (!gotW)
            gotW = findTetOfEdge(mfd, partner[0], partner[1], pw);
        if (gotW)
        {
            auto pc = pcBuf;
            immutable int np = wormEnumerate(mfd, partner, pw[], pc);
            int kp = 0;
            foreach (i; 0 .. np)
                if (sameKey(chosen, pc[i])) kp++;
            if (np > 0)
                qf += cast(real) kp / cast(real) np;
        }
    }

    // ---- apply the chosen pair with objective/potential lockstep ---------
    BM bm1 = chosen.m1is23 ? BM(chosen.f1[], chosen.p1[])
                           : BM(chosen.p1[], chosen.f1[]);
    BM bm2 = chosen.m1is23 ? BM(chosen.p2[], chosen.f2[])
                           : BM(chosen.f2[], chosen.p2[]);
    // (m2 has the opposite kind of m1: 2->3 center is a face, 3->2 an edge)
    real baseRun = currentObjective - (potState !is null ? potState.total : 0.0L);
    if (!mfd.hasValidMove(bm1)) return false;
    immutable real dB1 = mfd.speculativeBistellarDelta(bm1, baseRun, params);
    real dP1 = 0.0L;
    if (potState !is null)
        dP1 = mfd.potentialBistellarDelta(bm1, *potState, *pot, true);
    mfd.doMove(bm1);
    baseRun += dB1;
    bool undo1AndFail()
    {
        BM inv = chosen.m1is23 ? BM(chosen.p1[], chosen.f1[])
                               : BM(chosen.f1[], chosen.p1[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(inv, *potState, *pot, true);
        mfd.doMove(inv);
        return false;
    }
    if (!mfd.hasValidMove(bm2)) return undo1AndFail();
    immutable real dB2 = mfd.speculativeBistellarDelta(bm2, baseRun, params);
    real dP2 = 0.0L;
    if (potState !is null)
        dP2 = mfd.potentialBistellarDelta(bm2, *potState, *pot, true);
    mfd.doMove(bm2);
    immutable real dS = dB1 + dP1 + dB2 + dP2;

    bool undoBoth()
    {
        BM inv2 = chosen.m1is23 ? BM(chosen.f2[], chosen.p2[])
                                : BM(chosen.p2[], chosen.f2[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(inv2, *potState, *pot, true);
        mfd.doMove(inv2);
        BM inv1 = chosen.m1is23 ? BM(chosen.p1[], chosen.f1[])
                                : BM(chosen.f1[], chosen.p1[]);
        if (potState !is null)
            mfd.potentialBistellarDelta(inv1, *potState, *pot, true);
        mfd.doMove(inv1);
        return false;
    }

    // ---- reverse proposal weight: sum over the new4 anchors --------------
    real qr = 0.0L;
    foreach (bi; 0 .. chosen.nPair)
    {
        Vertex[2] b = chosen.new4[bi];
        Vertex[4] bw = -1;
        bool gotW = false;
        // an added tet of the composite containing b is a valid witness
        foreach (i; 0 .. chosen.nNet)
            if (chosen.netSg[i] > 0
                && Deg3Set!Vertex.edgeInFacet(b, chosen.netT[i]))
            { bw = chosen.netT[i]; gotW = true; break; }
        if (!gotW)
            gotW = findTetOfEdge(mfd, b[0], b[1], bw);
        if (!gotW) continue;
        auto rc = rcBuf;
        immutable int nr = wormEnumerate(mfd, b, bw[], rc);
        int kr = 0;
        foreach (i; 0 .. nr)
            if (inverseKey(chosen, rc[i])) kr++;
        if (nr > 0)
            qr += cast(real) kr / cast(real) nr;
    }
    valid = true;
    if (dSOut !is null) *dSOut = dS;
    if (outLa !is null) *outLa = chosen.landing[0];
    if (outLb !is null) *outLb = chosen.landing[1];

    final switch (policy)
    {
    case SlideAccept.trialOnly:
        return undoBoth();
    case SlideAccept.force:
        break;
    case SlideAccept.metropolis:
        if (qr <= 0)                         // closure says impossible; reject
            return undoBoth();
        real logAlpha = -dS + log(qr) - log(qf);
        if (logAlpha < 0 && uniform01 > exp(logAlpha))
            return undoBoth();
        break;
    }

    // ---- accepted ---------------------------------------------------------
    currentObjective += dS;
    // Mirror the two Pachner moves into the ledger so event-log consumers
    // (e.g. reaction_census's incremental state tracker) stay in sync.
    // recordBistellar/logEvent are post-move-safe (the unified commit calls
    // them after doMove too); sixFlipsBistellar is NOT (it reads pre-move
    // state), so the worm channel is gated off when logSixFlips is on --
    // see the channel gate in mcmcStep.
    if (ledger !is null)
    {
        if (ledger.logEvents)
        {
            Vertex[6] lb = [cast(Vertex) CHANNEL_WORM,
                            cast(Vertex) chosen.nPair,
                            anchor[0], anchor[1],
                            chosen.landing[0], chosen.landing[1]];
            logEvent(*ledger, EVT_CHANNEL_BEGIN, lb[0 .. 4], lb[4 .. 6]);
        }
        void mirror(ref BM bm)
        {
            if (ledger.trackRoles)
                recordBistellar(*ledger, bm.center, bm.coCenter);
            if (ledger.logEvents)
                logEvent(*ledger, cast(int) bm.coCenter.length - 1,
                         bm.center, bm.coCenter);
        }
        mirror(bm1);
        mirror(bm2);
        if (ledger.logEvents)
        {
            foreach (i; 0 .. chosen.nPair)
            {
                Vertex[3] g = [cast(Vertex) ROLE_GONE4,
                               chosen.gone4[i][0], chosen.gone4[i][1]];
                logEvent(*ledger, EVT_CHANNEL_EDGE, g[], null);
                Vertex[3] w = [cast(Vertex) ROLE_NEW4,
                               chosen.new4[i][0], chosen.new4[i][1]];
                logEvent(*ledger, EVT_CHANNEL_EDGE, w[], null);
            }
            Vertex[3] le = [cast(Vertex) CHANNEL_WORM, cast(Vertex) 2,
                            cast(Vertex) dSQuantized(dS)];
            logEvent(*ledger, EVT_CHANNEL_END, le[], null);
        }
    }
    {
        // reconcile the live degree sets over everything the pair touched
        Vertex[12] sup;
        int nS = 0;
        void addS(Vertex v)
        {
            foreach (i; 0 .. nS) if (sup[i] == v) return;
            if (nS < sup.length) sup[nS++] = v;
        }
        foreach (v; chosen.f1) addS(v);
        foreach (v; chosen.p1) addS(v);
        foreach (v; chosen.f2) addS(v);
        foreach (v; chosen.p2) addS(v);
        addS(anchor[0]); addS(anchor[1]);
        Vertex[4] wt = -1;
        foreach (i; 0 .. chosen.nNet)
            if (chosen.netSg[i] > 0) { wt = chosen.netT[i]; break; }
        if (deg3Set !is null)
            reconcileDegSet(mfd, *deg3Set, 3, sup[0 .. nS], wt[]);
        if (deg4Set !is null)
            reconcileDegSet(mfd, *deg4Set, 4, sup[0 .. nS], wt[]);
    }
    return true;
}

/// Rollback exactness: a slide that does not commit must leave the manifold
/// bitwise as it found it, whatever the reason it declined -- malformed
/// frame, a mid-composite Pachner move that fails to validate, an unclean
/// (species-changing) end state, or a Metropolis rejection. The trial pass
/// applies all four moves for real, so this is the invariant that keeps a
/// rejected proposal from corrupting the state.
///
/// COVERAGE: a thermalized sphere has degree-3 chords whose frames derive, so
/// the composite really is applied and rolled back here (the test asserts as
/// much, so it cannot quietly decay into checking nothing). It does NOT
/// produce CLEAN slides -- those need a (3,4,4) knot embedded in a
/// Boerdijk-Coxeter chain, i.e. a crystal-derived state, not a random one.
/// The clean/commit half of the move is covered exhaustively against the
/// Python oracle by scripts/defect_dynamics/worm_slide.py --dcross.
unittest
{
    import std.random : rndGen;

    struct TestParams
    {
        int numFacetsTarget = 150;
        real hingeDegreeTarget = 5.105025;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.85;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 9.5;
    }

    static int[][] snapshot(ref Manifold!3 mfd)
    {
        int[][] fs;
        foreach (f; mfd.facets)
            fs ~= f.dup;
        fs.sort();
        return fs;
    }

    rndGen.seed(20260724);          // deterministic thermalization
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    {
        real obj = mfd.objective(params);
        int[] unused;
        ulong hT, hA;
        ulong[4] bT, bA;
        foreach (_; 0 .. 2000)
            mfd.mcmcStep(obj, unused, params, 0.0, hT, hA, bT, bA);
    }

    auto before = snapshot(mfd);
    auto fvec0 = mfd.fVector.dup;
    static immutable int[2][6] picks =
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];

    int nProbed = 0, nFramed = 0, nClean = 0;
    foreach (edge; mfd.simplices(1).map!(e => [e[0], e[1]]).array)
    {
        if (mfd.degree(edge) != 3) continue;

        // any facet on the edge serves as the link-walk hint
        int[] hint;
        foreach (f; mfd.facets)
            if (edge[].isSubsetOf(f[])) { hint = f.dup; break; }
        assert(hint.length == 4);

        foreach (slot; 0 .. SLIDE_SLOTS)
        {
            nProbed++;
            // Derive the frame separately, purely to measure how deep this
            // probe reaches: a derived frame means trySlideMove below really
            // applies the composite (and must therefore really roll it back).
            bool framed = false;
            {
                immutable c0 = slot / 6 == 0 ? edge[0] : edge[1];
                immutable c4 = slot / 6 == 0 ? edge[1] : edge[0];
                int[8] lb = 0;
                if (mfd.writeEdgeLinkCycle(edge[0], edge[1], hint, lb.ptr) == 3)
                {
                    int[3] lk = [lb[0], lb[1], lb[2]];
                    lk[].sort();
                    SlideFrame!int fr;
                    if (deriveSlideFrame(mfd, c0, c4, lk[picks[slot % 6][0]],
                                         lk[picks[slot % 6][1]], fr))
                        { nFramed++; framed = true; }
                }
            }

            real obj = mfd.objective(params);
            immutable objIn = obj;
            bool valid = false;
            real dS = real.nan;
            immutable committed = mfd.trySlideMove(obj, edge[0], edge[1],
                hint, slot, params, valid, null, null, null, null, null,
                SlideAccept.trialOnly, &dS);

            // trialOnly never commits, whatever the verdict
            assert(!committed);
            assert(obj == objIn, "trial pass leaked an objective change");
            assert(snapshot(mfd) == before, "slide trial did not roll back");
            assert(mfd.fVector == fvec0);
            // The incremental degree/ridge maps must survive the rollback
            // too, not just the facet set -- audit them wherever a composite
            // was actually applied and undone. (validateMaps is the cheap
            // per-probe check; the full topological audit runs once below,
            // since findProblems rebuilds a SimplicialComplex.)
            if (framed)
                assert(mfd.validateMaps is null);
            if (valid)
            {
                nClean++;
                assert(dS == dS, "a valid slide must report a finite dS");
            }
        }
    }
    assert(nProbed > 0, "no degree-3 chords to probe");
    assert(nFramed > 0, "no frame derived: the composite was never applied, "
                        ~ "so rollback exactness went untested");
    assert(mfd.findProblems.length == 0);
}

/******************************************************************************
Run one MCMC step using a unified proposal that naturally includes both
bistellar (Pachner) moves and 4-4 hinge moves.

Proposal: pick a random facet, pick a random sub-simplex (center). Check the
degree of the center. For edges (dim-1 center) with degree 4, propose a hinge
move instead of the (invalid) 3→2 bistellar move. All other cases follow the
standard Pachner move logic.

This avoids a separate hinge move probability parameter — hinge moves are
proposed at the natural rate determined by how many degree-4 edges exist.
*/
/// One contract/split proposal (see ContractSplitConfig). Returns true iff
/// a move was accepted. Proposal and Hastings bookkeeping:
///   contract x->y, merging (u,v) -> u (u < v):
///     q_fwd = 1/2 * deg(uv) / (6 f3(x))
///     q_rev = 1/2 * [deg3_y(u) / (4 f3(y))] * [1 / N_L(link_y(u))] * 1/2
///   split x->y at (w, gamma, side):
///     q_fwd = 1/2 * [deg3(w) / (4 f3(x))] * [1 / N_L(link(w))] * 1/2
///     q_rev = 1/2 * |gamma| / (6 f3(y))
/// where N_L counts simple cycles of length <= maxRing in the named link
/// (FK-catalog table when it matches, bounded DFS otherwise -- for the
/// merged link of a contraction it is always DFS on the PLANNED state).
/// The 1/2 direction factors cancel; the extra 1/2 is the side coin.
private bool tryContractSplit(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    ref Vertex[] unusedVertices,
    P params,
    ref ContractSplitConfig cs,
    VertexPotState!Vertex* potState,
    const(VertexPot)* pot,
    Deg3Set!Vertex* deg3Set,
    Deg3Set!Vertex* deg4Set)
{
    import link_cycles : CatalogClass, adjFromFaces, catalog, countCycles,
        kthCycle, matchCatalog, maxCycleLen, maxLinkVerts;
    import std.math : log;

    enum dim = 3;
    enum STAR_CAP = 96;

    immutable bool doContract = uniform(0, 2) == 0;
    auto facetR = mfd.randomFacetOfDim(dim);
    Vertex[4] facet;
    foreach (i; 0 .. 4)
        facet[i] = facetR[i];
    immutable long f3x = cast(long) mfd.fVector[dim];
    real baseObj = currentObjective
        - (potState !is null ? potState.total : 0.0L);

    // local-index link building shared by both directions
    Vertex[] linkVerts;
    int[Vertex] lidx;
    ubyte[3][] lfaces;
    bool linkOverflow = false;
    void addLinkFace(Vertex apex, Vertex[4] t)
    {
        ubyte[3] lf;
        int n = 0;
        foreach (x; t)
            if (x != apex)
            {
                auto p = x in lidx;
                int li;
                if (p is null)
                {
                    li = cast(int) linkVerts.length;
                    if (li >= maxLinkVerts)
                    {
                        linkOverflow = true;
                        return;
                    }
                    lidx[x] = li;
                    linkVerts ~= x;
                }
                else
                    li = *p;
                lf[n++] = cast(ubyte) li;
            }
        assert(n == 3);
        lf[].sort();
        lfaces ~= lf;
    }

    if (doContract)
    {
        static immutable int[2][6] pairIdx =
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
        immutable pi = pairIdx[uniform(0, 6)];
        Vertex u = facet[pi[0]];
        Vertex v = facet[pi[1]];
        if (u > v)
        {
            immutable t = u;
            u = v;
            v = t;
        }
        Vertex[2] uv = [u, v];
        immutable long rl = cast(long) mfd.degreeOrZero!1(uv[]);
        if (rl < 3 || rl > cs.maxRing)
        {
            cs.noValid++;
            return false;
        }

        Vertex[4][STAR_CAP] starBuf;
        immutable int nU = collectStar(mfd, u, facet, starBuf[], 0);
        immutable int nAll = collectStar(mfd, v, facet, starBuf[], nU);
        if (nAll >= STAR_CAP)
        {
            cs.noValid++;
            return false;
        }

        // classify: star(u), star(v) (collectStar dedups shared tets)
        Vertex[4][STAR_CAP] svBuf;
        int nSV = 0;
        Vertex[4][STAR_CAP] suBuf;
        int nSU = 0;
        foreach (i; 0 .. nAll)
        {
            bool hasU = false, hasV = false;
            foreach (x; starBuf[i])
            {
                hasU |= x == u;
                hasV |= x == v;
            }
            if (hasU)
                suBuf[nSU++] = starBuf[i];
            if (hasV)
                svBuf[nSV++] = starBuf[i];
        }

        // support = every vertex of star(u) + star(v); frozen rejection
        Vertex[] sup;
        foreach (i; 0 .. nAll)
            foreach (x; starBuf[i])
                if (!sup.canFind(x))
                    sup ~= x;
        if (mfd.anyFrozen(sup))
        {
            cs.noValid++;
            return false;
        }

        // link condition, locally: (a) common neighbours = ring vertices
        Vertex[] nbrsU, nbrsV, ringV;
        bool[Vertex[2]] ringEdges;
        foreach (i; 0 .. nSU)
            foreach (x; suBuf[i])
                if (x != u && !nbrsU.canFind(x))
                    nbrsU ~= x;
        foreach (i; 0 .. nSV)
            foreach (x; svBuf[i])
                if (x != v && !nbrsV.canFind(x))
                    nbrsV ~= x;
        foreach (i; 0 .. nSV)
        {
            bool hasU = false;
            foreach (x; svBuf[i])
                hasU |= x == u;
            if (!hasU)
                continue;
            Vertex[2] opp;
            int no = 0;
            foreach (x; svBuf[i])
                if (x != u && x != v)
                    opp[no++] = x;
            assert(no == 2);
            if (opp[0] > opp[1])
            {
                immutable t = opp[0];
                opp[0] = opp[1];
                opp[1] = t;
            }
            ringEdges[opp] = true;
            foreach (x; opp)
                if (!ringV.canFind(x))
                    ringV ~= x;
        }
        nbrsU.sort();
        nbrsV.sort();
        ringV.sort();
        auto common = setIntersection(nbrsU, nbrsV).array;
        if (common != ringV)
        {
            cs.noValid++;
            return false;
        }
        // (b) common link edges = ring edges exactly
        foreach (i; 0 .. common.length)
            foreach (j; i + 1 .. common.length)
            {
                immutable a = common[i], b = common[j];
                Vertex[3] tu = [u, a, b];
                tu[].sort();
                Vertex[3] tv = [v, a, b];
                tv[].sort();
                immutable bothTri = mfd.contains(tu[]) && mfd.contains(tv[]);
                Vertex[2] ab = [a, b];
                if (bothTri != ((ab in ringEdges) !is null))
                {
                    cs.noValid++;
                    return false;
                }
            }
        // (c) no common link triangles
        foreach (i; 0 .. common.length)
            foreach (j; i + 1 .. common.length)
                foreach (k; j + 1 .. common.length)
                {
                    Vertex[4] tu = [u, common[i], common[j], common[k]];
                    tu[].sort();
                    Vertex[4] tv = [v, common[i], common[j], common[k]];
                    tv[].sort();
                    if (mfd.contains(tu[]) && mfd.contains(tv[]))
                    {
                        cs.noValid++;
                        return false;
                    }
                }

        auto move = ContractMove!Vertex(u, v);
        mfd.planContractMove(move, svBuf[0 .. nSV]);

        // merged link of u on the PLANNED state: u-only old tets + added
        foreach (i; 0 .. nSU)
        {
            bool hasV = false;
            foreach (x; suBuf[i])
                hasV |= x == v;
            if (!hasV)
                addLinkFace(u, suBuf[i]);
        }
        foreach (ref t; move.addedFacets)
        {
            Vertex[4] tt = t;
            addLinkFace(u, tt);
        }
        if (linkOverflow)
        {
            cs.noValid++;
            return false;
        }
        uint[maxLinkVerts] adjY;
        adjFromFaces(lfaces, adjY);
        immutable long NY =
            countCycles(adjY[0 .. linkVerts.length], cs.maxRing)
            .total(cs.maxRing);
        assert(NY >= 1, "merged link lost its ring cycle");
        immutable long deg3yU = cast(long) lfaces.length;
        immutable long f3y = f3x + cast(long) move.addedFacets.length
            - cast(long) move.removedFacets.length;

        immutable real logQ =
            log(cast(real) deg3yU) - log(4.0L * f3y) - log(cast(real) NY)
            - log(2.0L) - (log(cast(real) rl) - log(6.0L * f3x));

        real dS = mfd.speculativeBlockDelta(move.removedFacets,
                                            move.addedFacets, baseObj, params);
        Vertex[1] removedV = [v];
        if (potState !is null)
            dS += mfd.potentialBlockDelta(move.removedFacets,
                move.addedFacets, [], removedV[], *potState, *pot, false);

        cs.contractTries++;
        immutable real logAlpha = -dS + logQ;
        if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
        {
            if (potState !is null)
                mfd.potentialBlockDelta(move.removedFacets, move.addedFacets,
                    [], removedV[], *potState, *pot, true);
            mfd.commitPlannedMove(move);
            unusedVertices ~= v;
            currentObjective += dS;
            cs.contractAccepts++;
            if (deg3Set !is null)
                reconcileDeg3(mfd, *deg3Set, sup);
            if (deg4Set !is null)
                reconcileDegSet(mfd, *deg4Set, 4, sup);
            return true;
        }
        return false;
    }
    else
    {
        immutable Vertex w = facet[uniform(0, 4)];
        Vertex[4][STAR_CAP] starBuf;
        immutable int nT = collectStar(mfd, w, facet, starBuf[], 0);
        if (nT >= STAR_CAP)
        {
            cs.noValid++;
            return false;
        }
        foreach (i; 0 .. nT)
            addLinkFace(w, starBuf[i]);
        if (linkOverflow)
        {
            cs.noValid++;
            return false;
        }
        immutable int nL = cast(int) linkVerts.length;

        Vertex[] sup = linkVerts.dup;
        sup ~= w;
        if (mfd.anyFrozen(sup))
        {
            cs.noValid++;
            return false;
        }

        uint[maxLinkVerts] adjX;
        adjFromFaces(lfaces, adjX);
        long NX;
        ubyte[maxCycleLen] cbuf;
        int clen = 0;
        ubyte[maxLinkVerts] perm;
        immutable ci = matchCatalog(nL, lfaces, perm);
        if (ci >= 0)
        {
            // FK-catalog fast path: O(1) count + transported cycle list
            const cat = &catalog(cast(CatalogClass) ci);
            NX = cat.counts.total(cs.maxRing);
            version (ExpensiveAsserts)
                assert(NX == countCycles(adjX[0 .. nL], cs.maxRing)
                    .total(cs.maxRing));
            immutable long k = uniform(0L, NX);
            auto cyc = cat.cycles[cast(size_t) k];
            clen = cast(int) cyc.length;
            foreach (i; 0 .. clen)
                cbuf[i] = perm[cyc[i]];
        }
        else
        {
            NX = countCycles(adjX[0 .. nL], cs.maxRing).total(cs.maxRing);
            if (NX == 0)
            {
                cs.noValid++;
                return false;
            }
            immutable long k = uniform(0L, NX);
            clen = kthCycle(adjX[0 .. nL], cs.maxRing, k, cbuf);
            assert(clen >= 3);
        }
        Vertex[] gamma;
        foreach (i; 0 .. clen)
            gamma ~= linkVerts[cbuf[i]];
        immutable bool freshSide0 = uniform(0, 2) == 0;

        if (unusedVertices.empty)
            unusedVertices ~= mfd.fVector[0].to!Vertex;
        immutable Vertex fresh = unusedVertices.back;
        Vertex[1] freshV = [fresh];
        if (mfd.contains(freshV[]))
        {
            // label-pool head collides with a live vertex (possible when the
            // pool was seeded from fVector[0] on a label-holey manifold);
            // graceful reject, mirroring hasValidMove's 1->4 path
            cs.noValid++;
            return false;
        }

        auto move = SplitMove!Vertex(w, fresh, gamma, freshSide0);
        mfd.planSplitMove(move, starBuf[0 .. nT]);
        immutable long f3y = f3x + cast(long) move.addedFacets.length
            - cast(long) move.removedFacets.length;

        immutable real logQ =
            log(cast(real) clen) - log(6.0L * f3y)
            - (log(cast(real) nT) - log(4.0L * f3x)
               - log(cast(real) NX) - log(2.0L));

        real dS = mfd.speculativeBlockDelta(move.removedFacets,
                                            move.addedFacets, baseObj, params);
        Vertex[1] createdV = [fresh];
        if (potState !is null)
            dS += mfd.potentialBlockDelta(move.removedFacets,
                move.addedFacets, createdV[], [], *potState, *pot, false);

        cs.splitTries++;
        immutable real logAlpha = -dS + logQ;
        if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
        {
            if (potState !is null)
                mfd.potentialBlockDelta(move.removedFacets, move.addedFacets,
                    createdV[], [], *potState, *pot, true);
            mfd.commitPlannedMove(move);
            unusedVertices.popBack;
            unusedVertices.assumeSafeAppend;
            currentObjective += dS;
            cs.splitAccepts++;
            sup ~= fresh;
            if (deg3Set !is null)
                reconcileDeg3(mfd, *deg3Set, sup);
            if (deg4Set !is null)
                reconcileDegSet(mfd, *deg4Set, 4, sup);
            return true;
        }
        return false;
    }
}

bool mcmcStep(Vertex, P)(
    ref Manifold!(3, Vertex) mfd,
    ref real currentObjective,
    ref Vertex[] unusedVertices,
    P params,
    real hingeMoveProb,  // ignored (kept for API compatibility)
    ref ulong hingeTries,
    ref ulong hingeAccepts,
    ref ulong[4] bistellarTries,
    ref ulong[4] bistellarAccepts,
    MoveCounters!Vertex* counters = null,
    GeometryLedger!Vertex* ledger = null,
    VertexPotState!Vertex* potState = null,
    const(VertexPot)* pot = null,
    CocycleState!Vertex* cocycle = null,
    SlideConfig* slide = null,
    NonlocalSlideConfig* nlSlide = null,
    Deg3Set!Vertex* deg3Set = null,
    WormConfig* worm = null,
    Deg3Set!Vertex* deg4Set = null,
    ContractSplitConfig* contractSplit = null)
{
    enum dim = 3;
    enum nVerts = dim + 1;
    enum maxMask = (1 << nVerts) - 1;
    alias BM = BistellarMove!(dim, Vertex);

    if (ledger !is null)
        ledger.clock++;               // one tick per attempted move

    // --- Knot-slide channel ---------------------------------------------
    // An INDEPENDENT move type, chosen with a STATE-INDEPENDENT probability,
    // so it never cannibalizes the unified proposal below. This structure is
    // load-bearing, not stylistic: selecting the slide *inside* the proposal
    // loop (i.e. "at a degree-3 edge, slide instead of doing the 3->2") would
    // throttle the 3->2 move by (1 - prob) at exactly the degree-3 edges a
    // slide competes for, while leaving untouched the 2->3 at a triangle
    // centre that CREATES such an edge. A 3->2 on a degree-3 edge is the
    // reverse of that 2->3, so suppressing one and not the other breaks
    // detailed balance and degree-3 defects accumulate without bound (n_ill
    // ran away 46 -> 232 at prob = 1 before this was restructured).
    //
    // Both branches satisfy detailed balance on their own -- the ordinary
    // proposal because it is untouched, the slide because its proposal is
    // symmetric (every degree-3 edge lies in exactly 3 facets, so a random
    // facet + random vertex pair is uniform over chords; n_3 and N_3 are
    // preserved on the clean class; k_f = k_r = 1; and the slot denominator
    // is the constant SLIDE_SLOTS) -- so mixing them with a state-independent
    // weight satisfies it too. A draw that lands on a non-chord, a malformed
    // slot or an unclean end state is simply a rejected step, which is why
    // the denominator can stay constant.
    if (slide !is null && slide.prob > 0 && uniform01 < slide.prob)
    {
        auto facet = mfd.randomFacetOfDim(dim);
        // one of the 6 vertex pairs of the facet, uniformly
        static immutable int[2][6] pairIdx =
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
        immutable pi = pairIdx[uniform(0, 6)];
        Vertex[2] chord = [facet[pi[0]], facet[pi[1]]];
        chord[].sort();
        if (mfd.degreeOrZero!1(chord[]) != 3)
            return false;             // not a chord: rejected step
        bool slideValid = false;
        immutable ok = trySlideMove(mfd, currentObjective,
            chord[0], chord[1], facet, uniform(0, SLIDE_SLOTS), params,
            slideValid, counters, ledger, potState, pot, cocycle,
            SlideAccept.metropolis, null, slide.cleanOnly, deg3Set);
        if (slideValid)
        {
            slide.tries++;
            if (ok) slide.accepts++;
        }
        if (ok && deg4Set !is null)
            rebuildDegSet(mfd, *deg4Set, 4);   // O(N); slides commit rarely
        return ok;
    }

    // --- Non-local slide channel (dim=3) -----------------------------------
    // Pick a degree-3 chord UNIFORMLY (1/n_3) from the live set, then
    // annihilate + re-create it up to maxStep tets down its BC chain. The
    // acceptance carries the 1/n_3 Hastings factor via n3Before = |deg3Set|;
    // direction is the slot's orientation so the step distribution is
    // symmetric about 0. deg3Set is maintained by reconcileDeg3 at every commit
    // site, so the draw and n_3 are O(1).
    // Not cocycle-safe (tryNonlocalSlide has no cocycle lockstep) and not
    // six-flip-safe (apply-then-decide: pre-move state unavailable at
    // commit) -- gated off under both, like the worm channel.
    if (nlSlide !is null && nlSlide.prob > 0 && deg3Set !is null
        && cocycle is null && !(ledger !is null && ledger.logSixFlips)
        && deg3Set.length > 0 && uniform01 < nlSlide.prob)
    {
        immutable long n3 = cast(long) deg3Set.length;
        immutable size_t di = uniform(0, deg3Set.length);
        Vertex[2] chord = deg3Set.edges[di];
        // slide hint = a facet containing the chord. Use the cached witness
        // (validated in O(1) via writeFaceApexes); only on a stale/unknown
        // witness fall back to the O(N) mfd.link, then cache the result.
        Vertex[4] hint = 0;
        {
            Vertex[4] w = deg3Set.wit[di];
            bool valid = false;
            if (w[0] >= 0 && Deg3Set!Vertex.edgeInFacet(chord, w))
            {
                Vertex c = -1, d = -1;
                foreach (v; w)
                    if (v != chord[0] && v != chord[1])
                    { if (c < 0) c = v; else d = v; }
                int[2] ap = 0;
                if (c >= 0 && d >= 0
                    && mfd.writeFaceApexes(chord[0], chord[1], c, ap.ptr) == 2)
                    valid = (ap[0] == d || ap[1] == d);
            }
            if (valid)
                hint = w;
            else
            {
                bool got = false;
                foreach (pr; mfd.link(chord[]))
                {
                    hint[0] = chord[0]; hint[1] = chord[1];
                    int i = 0;
                    foreach (v; pr) hint[2 + i++] = v;
                    got = true; break;
                }
                if (!got) return false;
                hint[].sort();
                deg3Set.wit[di] = hint;          // cache for next time
            }
        }
        immutable int slot = uniform(0, SLIDE_SLOTS);
        immutable int steps = uniform(1, nlSlide.maxStep + 1);
        bool nlValid = false;
        immutable ok = tryNonlocalSlide(mfd, currentObjective,
            chord[0], chord[1], hint[], slot, steps, params, nlValid,
            potState, pot, SlideAccept.metropolis, null, n3, null, null, null,
            deg3Set, ledger);
        if (nlValid)
        {
            nlSlide.tries++;
            if (ok) nlSlide.accepts++;
        }
        if (ok && deg4Set !is null)
            rebuildDegSet(mfd, *deg4Set, 4);   // O(N); commits are rare
        return ok;
    }

    // --- Deg-4 worm channel (dim=3) -----------------------------------------
    // Anchor uniform over the live deg-4 set (1/n_4 cancels: the move class
    // preserves the deg-4 count), then the catalysed 2-move transport step
    // with anchor-sum Hastings weights (see tryWormMove).  Not cocycle-safe:
    // gated off whenever a cocycle is attached.  Also gated off under
    // logSixFlips: six-flip records need PRE-move state, and the worm's
    // trial applies its moves before the accept decision (recordBistellar /
    // logEvent are post-move-safe and ARE emitted -- see tryWormMove).
    if (worm !is null && worm.prob > 0 && deg4Set !is null
        && cocycle is null && !(ledger !is null && ledger.logSixFlips)
        && deg4Set.length > 0 && uniform01 < worm.prob)
    {
        immutable size_t di = uniform(0, deg4Set.length);
        Vertex[2] e4 = deg4Set.edges[di];
        Vertex[4] hint = deg4Set.wit[di];
        bool wOk = hint[0] >= 0 && Deg3Set!Vertex.edgeInFacet(e4, hint);
        if (wOk)
        {
            Vertex[4] hs = hint; hs[].sort();
            wOk = mfd.contains(hs[]);
        }
        if (!wOk)
        {
            if (!findTetOfEdge(mfd, e4[0], e4[1], hint)) return false;
            hint[].sort();
            deg4Set.wit[di] = hint;
        }
        bool wormValid = false;
        immutable ok = tryWormMove(mfd, currentObjective, e4, hint[], params,
            wormValid, potState, pot, deg4Set, deg3Set,
            SlideAccept.metropolis, null, -1, null, null, ledger);
        if (wormValid)
        {
            worm.tries++;
            if (ok) worm.accepts++;
        }
        else
            worm.noCands++;
        return ok;
    }

    // --- Contract/split channel (dim=3) ---------------------------------
    // The only channel that changes f0 without a degree-4 vertex. Not
    // cocycle-safe (a vertex disappears; the lift cannot follow) and not
    // six-flip-safe -- gated off under both, like the worm channel.
    if (contractSplit !is null && contractSplit.prob > 0 && cocycle is null
        && !(ledger !is null && ledger.logSixFlips)
        && uniform01 < contractSplit.prob)
    {
        return tryContractSplit(mfd, currentObjective, unusedVertices,
            params, *contractSplit, potState, pot, deg3Set, deg4Set);
    }

    // Unified proposal loop
    while (true)
    {
        auto facet = mfd.randomFacetOfDim(dim);

        auto mask = uniform(1, maxMask + 1);
        Vertex[nVerts] centerBuf;
        int centerLen = 0;
        foreach (i; 0 .. nVerts)
        {
            if (mask & (1 << i))
                centerBuf[centerLen++] = facet[i];
        }
        auto center = centerBuf[0 .. centerLen];
        center.sort();

        auto centerDim = centerLen - 1;
        auto centerDeg = mfd.degree(center);

        // --- Edge of degree 4: propose hinge move ---
        if (centerDim == 1 && centerDeg == 4)
        {
            Vertex[2] edge = [center[0], center[1]];

            // Find a start vertex (any facet vertex not in the edge)
            Vertex startVertex;
            foreach (v; facet)
                if (v != edge[0] && v != edge[1]) { startVertex = v; break; }

            auto hm = mfd.hingeMove(edge, startVertex, uniform(0, 2));

            // Support = the 6 vertices whose stars the move touches.
            Vertex[6] hingeSupport;
            hingeSupport[0 .. 2] = hm.removedEdge[];
            hingeSupport[2 .. 6] = hm.linkCycle[];
            if (counters !is null)
                addSupport(counters.proposed, hingeSupport[]);

            // Frozen-region rejection: every facet this move adds/removes has
            // all its vertices in the support, so rejecting here preserves the
            // frozen set's closed star exactly (facets + hinge degrees).
            if (mfd.anyFrozen(hingeSupport[]))
                continue;

            if (!mfd.hasValidHingeMove(hm))
                continue;

            hingeTries++;
            if (counters !is null)
                addSupport(counters.valid, hingeSupport[]);
            // currentObjective includes the potential total (when enabled); the
            // base speculative delta needs the base-only objective as baseline.
            real baseObj = currentObjective
                - (potState !is null ? potState.total : 0.0L);
            real deltaObj = mfd.speculativeHingeDelta(hm, baseObj, params);
            if (potState !is null)
                deltaObj += mfd.potentialHingeDelta(hm, *potState, *pot, false);
            real logAlpha = -deltaObj;

            if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
            {
                if (potState !is null)
                    mfd.potentialHingeDelta(hm, *potState, *pot, true);
                if (ledger !is null && ledger.logSixFlips)
                    sixFlipsHinge(*ledger, mfd, hm);
                if (cocycle !is null)
                    cocycleHinge(*cocycle, hm);
                mfd.doHingeMove(hm);
                currentObjective += deltaObj;
                hingeAccepts++;
                if (deg3Set !is null)
                    reconcileDeg3(mfd, *deg3Set, hingeSupport[]);
                if (deg4Set !is null)
                    reconcileDegSet(mfd, *deg4Set, 4, hingeSupport[]);
                if (counters !is null)
                    addSupport(counters.acceptedHinge, hingeSupport[]);
                if (ledger !is null)
                {
                    if (ledger.trackRoles)
                        recordHinge(*ledger, hm.removedEdge, hm.addedEdge,
                                    hm.linkCycle);
                    if (ledger.logEvents)
                    {
                        // Rotate cycle so the diagonal is (labels[2], labels[4]).
                        Vertex[4] lc = hm.linkCycle;
                        immutable d0 = mkEdge(lc[0], lc[2]);
                        if (!(d0[0] == min(hm.addedEdge[0], hm.addedEdge[1])
                              && d0[1] == max(hm.addedEdge[0], hm.addedEdge[1])))
                            lc = [hm.linkCycle[1], hm.linkCycle[2],
                                  hm.linkCycle[3], hm.linkCycle[0]];
                        logEvent(*ledger, 4, hm.removedEdge[], lc[]);
                    }
                }
                return true;
            }
            return false;
        }

        // --- Standard bistellar move ---
        if (centerDeg + centerDim != dim + 1)
            continue;

        BM bm;
        if (centerDim == dim)
        {
            if (unusedVertices.empty)
                unusedVertices ~= mfd.fVector[0].to!Vertex;
            bm = BM(center, unusedVertices.back.only);
        }
        else
        {
            auto coCenter = mfd.coCenter(center, facet);
            bm = BM(center, coCenter[]);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        // Support = the 5 vertices of the bistellar ball (one glued 4-simplex).
        Vertex[nVerts + 1] ballBuf;
        Vertex[] ball;
        {
            size_t nb = 0;
            foreach (v; bm.center) ballBuf[nb++] = v;
            foreach (v; bm.coCenter) ballBuf[nb++] = v;
            ball = ballBuf[0 .. nb];
        }
        if (counters !is null)
            addSupport(counters.proposed, ball);

        // Frozen-region rejection (see hinge branch). For a 1->4 move the
        // coCenter vertex is new/unused, hence never frozen.
        if (mfd.anyFrozen(ball))
            continue;

        if (!mfd.hasValidMove(bm))
            continue;

        bistellarTries[bm.coCenter.length - 1]++;
        if (counters !is null)
            addSupport(counters.valid, ball);
        real baseObj = currentObjective
            - (potState !is null ? potState.total : 0.0L);
        real deltaObj = mfd.speculativeBistellarDelta(bm, baseObj, params);
        if (potState !is null)
            deltaObj += mfd.potentialBistellarDelta(bm, *potState, *pot, false);
        real logAlpha = -deltaObj;

        if (logAlpha >= 0 || uniform01 <= exp(logAlpha))
        {
            if (potState !is null)
                mfd.potentialBistellarDelta(bm, *potState, *pot, true);
            if (ledger !is null && ledger.logSixFlips)
                sixFlipsBistellar(*ledger, mfd, bm.center, bm.coCenter);
            if (cocycle !is null)
                cocycleBistellar(*cocycle, bm.center, bm.coCenter);
            mfd.doMove(bm);
            if (deg3Set !is null)
                reconcileDeg3(mfd, *deg3Set, ball);
            if (deg4Set !is null)
                reconcileDegSet(mfd, *deg4Set, 4, ball);
            if (counters !is null)
                addSupport(counters.acceptedBistellar, ball);
            if (ledger !is null)
            {
                if (ledger.trackRoles)
                    recordBistellar(*ledger, bm.center, bm.coCenter);
                if (ledger.logEvents)
                    logEvent(*ledger, cast(int) bm.coCenter.length - 1,
                             bm.center, bm.coCenter);
            }
            if (bm.coCenter.length == 1)
            {
                // Shrink, then tell the runtime we own the freed slot so the
                // next `~=` reuses the buffer instead of reallocating a fresh
                // int[] every accepted move (that churn was false-pinned by the
                // conservative GC, leaking ~0.14 MB/sweep). Mirrors ddg_capi.d.
                unusedVertices.popBack;
                unusedVertices.assumeSafeAppend;
            }
            if (bm.center.length == 1) unusedVertices ~= bm.center;
            currentObjective += deltaObj;
            bistellarAccepts[bm.coCenter.length - 1]++;
            return true;
        }
        return false;
    }
}

///
unittest
{
    // Integration test: run mixed MCMC for a number of steps, verify manifold integrity
    struct TestParams
    {
        int numFacetsTarget = 20;
        real hingeDegreeTarget = 4.5;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.2;
        real coDim3DegreeVarianceCoef = 0.1;
        real hingeDegreeTargetCoef = 0.15;
        real coDim3DegreeTargetCoef = 0.07;
        real coDim3DegreeTarget = 9.5;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();

    // Grow a bit first
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));

    auto currentObj = mfd.objective(params);
    int[] unusedVertices = [8];
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    MoveCounters!int mc;

    int accepted = 0;
    foreach (_; 0 .. 500)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts, &mc))
            accepted++;
    }

    assert(accepted > 0, "should have accepted some moves");
    assert(hingeTries > 0, "should have attempted some hinge moves");

    // Move-counter invariants: every event distributes total weight 1 over its
    // support, so ledger totals must reproduce the scalar tallies exactly.
    double total(double[int] aa)
    {
        double s = 0;
        foreach (v; aa.byValue) s += v;
        return s;
    }
    assert(isClose(total(mc.valid), cast(double)(hingeTries + bTries[].sum), 0.0, 1e-6),
        "valid ledger total != tries");
    assert(isClose(total(mc.acceptedBistellar), cast(double) bAccepts[].sum, 0.0, 1e-6),
        "acceptedBistellar ledger total != bistellar accepts");
    assert(isClose(total(mc.acceptedHinge), cast(double) hingeAccepts, 0.0, 1e-6),
        "acceptedHinge ledger total != hinge accepts");
    assert(total(mc.proposed) >= total(mc.valid) - 1e-9,
        "proposed must dominate valid");

    // Verify manifold integrity after many mixed moves
    assert(mfd.findProblems.length == 0,
        "manifold has problems after mixed MCMC: " ~ mfd.findProblems.to!string);

    // Verify objective tracking is consistent
    auto actualObj = mfd.objective(params);
    assert(isClose(currentObj, actualObj, 1e-4),
        "objective drift: tracked=%s actual=%s".format(currentObj, actualObj));
}

/// Contract/split channel integration: mixed MCMC with the channel enabled
/// (and the vertex potential on) must keep the manifold valid, keep the
/// tracked objective and potential counters drift-free, move f0 through
/// both channel directions, and roundtrip the vertex-label pool.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 30;
        real hingeDegreeTarget = 4.8;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.02;
        real hingeDegreeVarianceCoef = 0.05;
        real coDim3DegreeVarianceCoef = 0.02;
        real hingeDegreeTargetCoef = 0.0;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 0.0;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));

    VertexPot pot;
    pot.zlegCoef = 0.1;
    pot.impCoef = 0.2;
    VertexPotState!int potState;
    mfd.recomputeVertexPotState(potState, pot);

    auto currentObj = mfd.objective(params) + potState.total;
    int[] unusedVertices = [7];
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    ContractSplitConfig cs;
    cs.prob = 0.5;
    cs.maxRing = 6;

    int accepted = 0;
    foreach (step; 0 .. 3000)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts,
                null, null, &potState, &pot, null, null, null, null, null,
                null, &cs))
            accepted++;
    }

    assert(cs.contractTries > 0, "no contract proposals reached Metropolis");
    assert(cs.splitTries > 0, "no split proposals reached Metropolis");
    assert(cs.contractAccepts + cs.splitAccepts > 0,
        "contract/split channel never accepted");

    assert(mfd.findProblems.length == 0,
        "manifold has problems after contract/split MCMC: "
        ~ mfd.findProblems.to!string);

    // objective and potential-counter consistency (no drift)
    VertexPotState!int freshState;
    mfd.recomputeVertexPotState(freshState, pot);
    assert(isClose(potState.total, freshState.total, 1e-6, 1e-9),
        "potential drift: tracked=%s actual=%s"
        .format(potState.total, freshState.total));
    foreach (v, n6; freshState.n6)
        assert(potState.n6.get(v, 0) == n6, "n6 counter drift");
    foreach (v, m; freshState.mImp)
        assert(potState.mImp.get(v, 0) == m, "impurity counter drift");
    auto actualObj = mfd.objective(params) + potState.total;
    assert(isClose(currentObj, actualObj, 1e-4),
        "objective drift: tracked=%s actual=%s".format(currentObj, actualObj));
}

/// Contract/split Hastings quantities: the contract branch estimates the
/// reverse-split's inputs (merged-link cycle count, post-move vertex tet
/// degree, post-move f3) on the PLANNED state; verify them against ground
/// truth computed on the committed state, for every valid contraction of a
/// churned sphere.
unittest
{
    import link_cycles : adjFromFaces, countCycles, maxLinkVerts;
    alias BM = BistellarMove!3;

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([0,1,2], [5,6]));

    int nChecked = 0;
    foreach (ePair; mfd.simplices(1).map!array.array)
    {
        immutable int u = ePair[0], v = ePair[1];
        if (!mfd.hasValidContractMove(u, v))
            continue;

        // planned-state estimates, as tryContractSplit computes them
        auto move = ContractMove!int(u, v);
        mfd.planContractMove(move);
        int[] linkVerts;
        int[int] lidx;
        ubyte[3][] lfaces;
        void addLinkFace(int[4] t)
        {
            ubyte[3] lf;
            int n = 0;
            foreach (x; t)
                if (x != u)
                {
                    auto p = x in lidx;
                    int li;
                    if (p is null)
                    {
                        li = cast(int) linkVerts.length;
                        lidx[x] = li;
                        linkVerts ~= x;
                    }
                    else
                        li = *p;
                    lf[n++] = cast(ubyte) li;
                }
            lf[].sort();
            lfaces ~= lf;
        }
        foreach (f; mfd.star([u]))
        {
            int[4] t;
            int n = 0;
            foreach (x; f)
                t[n++] = x;
            bool hasV = false;
            foreach (x; t)
                hasV |= x == v;
            if (!hasV)
                addLinkFace(t);
        }
        foreach (ref t; move.addedFacets)
            addLinkFace(t);
        assert(linkVerts.length <= maxLinkVerts);
        uint[maxLinkVerts] adjP;
        adjFromFaces(lfaces, adjP);
        immutable plannedN =
            countCycles(adjP[0 .. linkVerts.length], 6).total(6);
        immutable plannedDeg = cast(long) lfaces.length;
        immutable plannedF3 = cast(long) mfd.fVector[3]
            + cast(long) move.addedFacets.length
            - cast(long) move.removedFacets.length;

        // ground truth on the committed state
        mfd.commitPlannedMove(move);
        int[] lv2;
        int[int] lidx2;
        ubyte[3][] lf2;
        foreach (f; mfd.star([u]))
        {
            ubyte[3] lf;
            int n = 0;
            foreach (x; f)
                if (x != u)
                {
                    auto p = x in lidx2;
                    int li;
                    if (p is null)
                    {
                        li = cast(int) lv2.length;
                        lidx2[x] = li;
                        lv2 ~= x;
                    }
                    else
                        li = *p;
                    lf[n++] = cast(ubyte) li;
                }
            lf[].sort();
            lf2 ~= lf;
        }
        uint[maxLinkVerts] adjA;
        adjFromFaces(lf2, adjA);
        immutable actualN = countCycles(adjA[0 .. lv2.length], 6).total(6);
        assert(plannedDeg == cast(long) lf2.length,
            "planned deg3_y mismatch: %s vs %s".format(plannedDeg, lf2.length));
        assert(plannedN == actualN,
            "planned merged-link cycle count mismatch: %s vs %s"
            .format(plannedN, actualN));
        assert(plannedF3 == cast(long) mfd.fVector[3], "planned f3 mismatch");
        mfd.undoContractMove(move);
        ++nChecked;
    }
    assert(nChecked > 0, "no valid contractions found to check");
}

/// Six-flip stream + disclination census: the flip records must exactly
/// account for the six-edge count change, and the census must satisfy its
/// internal identities on a churned manifold.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 64;
        real hingeDegreeTarget = 5.1;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;

    long countSixEdges()
    {
        long n = 0;
        foreach (e; mfd.simplices(1))
            if (mfd.degree(e) >= 6) ++n;
        return n;
    }

    GeometryLedger!int g;
    g.logSixFlips = true;
    g.sixBuf = new ubyte[1 << 20];

    immutable e6Before = countSixEdges();
    foreach (_; 0 .. 4000)
        mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
            hingeTries, hingeAccepts, bTries, bAccepts, null, &g);
    immutable e6After = countSixEdges();

    assert(!g.sixOverflow, "six-flip buffer overflowed in test");
    long net = 0;
    ulong lastClock = 0;
    for (size_t p = 0; p + sixRecordBytes <= g.sixUsed; p += sixRecordBytes)
    {
        import core.stdc.string : memcpy;
        ulong clk; int u, v, dir;
        memcpy(&clk, g.sixBuf.ptr + p, 8);
        memcpy(&u, g.sixBuf.ptr + p + 8, 4);
        memcpy(&v, g.sixBuf.ptr + p + 12, 4);
        memcpy(&dir, g.sixBuf.ptr + p + 16, 4);
        assert(dir == 1 || dir == -1, "bad flip direction");
        assert(u < v, "flip labels not sorted");
        assert(clk >= lastClock, "clock not monotone");
        lastClock = clk;
        net += dir;
    }
    assert(e6After == e6Before + net,
        "six-flip stream does not account for six-edge count change: "
        ~ "%s -> %s but net flips %s".format(e6Before, e6After, net));
    assert(e6After > 0, "test manifold never grew six-edges (weak test)");

    // Census identities on the churned manifold.
    auto c = mfd.disclinationCensus(0);
    c.nSixEdges.shouldEqual(e6After);
    long degTot = 0, vertTot = 0;
    foreach (k; 0 .. 8) { degTot += (k + 1) * c.netDegCensus[k]; vertTot += c.netDegCensus[k]; }
    vertTot.shouldEqual(c.nNetVerts);
    // netDegCensus is clamped at 8+, but a churned 64-facet manifold should
    // not produce network degree >= 8; if it does, skip the handshake check.
    if (c.netDegCensus[7] == 0)
        degTot.shouldEqual(2 * c.nSixEdges);
    (c.sumSegLen + c.sumLoopLen).shouldEqual(c.nSixEdges);
    c.cycleRank.shouldEqual(c.nSixEdges - c.nNetVerts + c.nComponents);
    assert(c.cycleRank >= 0);
    assert(c.giantSize >= c.secondSize);
    assert(c.giantDiameter <= c.giantSize);
    long compTot = 0;
    foreach (x; c.compSizeHist) compTot += x;
    compTot.shouldEqual(c.nComponents);

    // Valence census cross-check: rows n6 >= 1 must reproduce nNetVerts and
    // the six-edge handshake, with caps high enough to avoid clamping.
    enum CAP = 31;
    long[(CAP + 1) * (CAP + 1)] vc;
    mfd.valenceCensus(vc[], CAP, CAP);
    long f0Tot = 0, netTot = 0, handshake = 0;
    foreach (k; 0 .. CAP + 1)
        foreach (m; 0 .. CAP + 1)
        {
            immutable cnt = vc[k * (CAP + 1) + m];
            f0Tot += cnt;
            if (k >= 1) { netTot += cnt; handshake += k * cnt; }
        }
    f0Tot.shouldEqual(mfd.fVector[0]);
    netTot.shouldEqual(c.nNetVerts);
    handshake.shouldEqual(2 * c.nSixEdges);

    // Flattened layout must reproduce the struct (the C API path).
    long[disclinationCensusSlots] flat;
    flattenCensus(c, flat[]);
    flat[0].shouldEqual(c.nNetVerts);
    flat[1].shouldEqual(c.nSixEdges);
    flat[6].shouldEqual(c.cycleRank);
    flat[20].shouldEqual(c.eImpAny);
    flat[21].shouldEqual(0);   // reserved
    long flatDegTot = 0;
    foreach (k; 0 .. 8) flatDegTot += flat[24 + k];
    flatDegTot.shouldEqual(c.nNetVerts);
}

/// Frozen vertices: the sampler must preserve the frozen set's closed star
/// EXACTLY (facets and the degrees of every simplex meeting a frozen vertex)
/// while the rest of the manifold churns freely.
unittest
{
    import std.algorithm : canFind;

    struct TestParams
    {
        int numFacetsTarget = 48;
        real hingeDegreeTarget = 4.8;
        real numFacetsCoef = 0.05;
        real numHingesCoef = 0.0;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 12.0;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    // Grow so there is bulk to churn away from the frozen star.
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));
    mfd.doMove(BM([0,1,2,5], [8]));

    // Freeze vertex 7 and snapshot its closed star + incident edge degrees.
    mfd.setVertexFrozen(7, true);
    assert(mfd.vertexFrozen(7));
    assert(!mfd.vertexFrozen(3));
    assert(mfd.numFrozenVertices == 1);

    int[4][] frozenStar;
    foreach (f; mfd.facets)
        if (f.canFind(7))
        {
            int[4] ff;
            ff[] = f[0 .. 4];
            frozenStar ~= ff;
        }
    size_t[int[2]] edgeDegs;      // degrees of edges at the frozen vertex
    foreach (e; mfd.simplices(1))
        if (e.canFind(7))
        {
            int[2] ee;
            ee[] = e[0 .. 2];
            edgeDegs[ee] = mfd.degree(e);
        }
    immutable frozenVtxDeg = mfd.degree([7]);

    auto params = TestParams();
    auto currentObj = mfd.objective(params);
    int[] unusedVertices;
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;
    int accepted = 0;
    foreach (step; 0 .. 4000)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts))
            accepted++;
    }
    assert(accepted > 50, "sampler froze entirely; test not exercising moves");

    // The frozen closed star is exactly intact.
    foreach (ff; frozenStar)
        assert(mfd.contains(ff[]), "frozen-star facet destroyed");
    size_t starCount = 0;
    foreach (f; mfd.facets)
        if (f.canFind(7)) starCount++;
    starCount.shouldEqual(frozenStar.length);   // no facet added at 7 either
    mfd.degree([7]).shouldEqual(frozenVtxDeg);
    foreach (ee, d; edgeDegs)
        mfd.degree(ee[]).shouldEqual(d);        // hinge curvature at 7 pinned

    // Unfreezing restores full dynamics (flag actually clears).
    mfd.clearFrozenVertices();
    assert(!mfd.vertexFrozen(7));
    assert(mfd.numFrozenVertices == 0);
    assert(mfd.findProblems.length == 0);
}

/// Vertex 6-valence potential: incremental counter state and running total
/// stay exact through mixed MCMC (all bistellar types + hinge moves), and the
/// tracked objective (base + potential) matches a from-scratch recompute.
unittest
{
    struct TestParams
    {
        int numFacetsTarget = 24;
        real hingeDegreeTarget = 4.6;
        real numFacetsCoef = 0.1;
        real numHingesCoef = 0.05;
        real hingeDegreeVarianceCoef = 0.0;
        real coDim3DegreeVarianceCoef = 0.0;
        real hingeDegreeTargetCoef = 0.1;
        real coDim3DegreeTargetCoef = 0.0;
        real coDim3DegreeTarget = 0.0;
    }

    alias BM = BistellarMove!3;
    auto mfd = Manifold!3([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
    auto params = TestParams();
    mfd.doMove(BM([0,1,2,3], [5]));
    mfd.doMove(BM([0,1,2,4], [6]));
    mfd.doMove(BM([1,2,3,4], [7]));

    VertexPot pot;
    pot.zlegCoef = 0.3;
    pot.impCoef = 0.05;
    pot.tilt = [0.02L, 0, -0.04L, 0.01L, -0.03L];  // exercise every tilt slot
    assert(pot.enabled);

    // U spot checks: legal set flat at tilt, n6=1/5 pay 1, quadratic tail.
    assert(isClose(pot.U(0), 0.02L));
    assert(isClose(pot.U(1), 0.3L + 0.0L));
    assert(isClose(pot.U(2), -0.04L));
    assert(isClose(pot.U(4), -0.03L));
    assert(isClose(pot.U(5), 0.3L));
    assert(isClose(pot.U(7), 0.3L * 9));
    assert(isClose(pot.V(6), 0.05L * 36));

    VertexPotState!int st;
    mfd.recomputeVertexPotState(st, pot);

    // Counter invariant: sum_v n6 = 2 * #edges(deg >= 6), same for impurity.
    void checkInvariants()
    {
        long s6 = 0, sI = 0;
        foreach (v, k; st.n6) s6 += k;
        foreach (v, m; st.mImp) sI += m;
        long e6 = 0, eI = 0;
        foreach (e; mfd.simplices(1))
        {
            immutable d = cast(long) mfd.degree(e);
            if (d >= 6) e6++;
            if (d < 5 || d > 6) eI++;
        }
        assert(s6 == 2 * e6, "n6 counter sum mismatch");
        assert(sI == 2 * eI, "impurity counter sum mismatch");
    }
    checkInvariants();

    auto currentObj = mfd.objective(params) + st.total;
    int[] unusedVertices = [8];
    ulong hingeTries, hingeAccepts;
    ulong[4] bTries, bAccepts;

    int accepted = 0;
    foreach (step; 0 .. 600)
    {
        if (mfd.mcmcStep(currentObj, unusedVertices, params, 0.5,
                hingeTries, hingeAccepts, bTries, bAccepts, null, null,
                &st, &pot))
            accepted++;

        if (step % 100 == 99)
        {
            // Incremental state == from-scratch recompute
            VertexPotState!int fresh;
            mfd.recomputeVertexPotState(fresh, pot);
            assert(fresh.n6 == st.n6,
                "n6 drift at step %s".format(step));
            assert(fresh.mImp == st.mImp,
                "impurity drift at step %s".format(step));
            assert(isClose(fresh.total, st.total, 1e-9L, 1e-9L),
                "potential total drift: tracked=%s fresh=%s"
                .format(st.total, fresh.total));
            checkInvariants();

            // Tracked objective == base + potential, recomputed
            assert(isClose(currentObj, mfd.objective(params) + st.total, 1e-6L),
                "objective drift with potential at step %s".format(step));
        }
    }
    assert(accepted > 0, "should have accepted some moves");
    assert(hingeAccepts + bAccepts[].sum > 0);
    assert(mfd.findProblems.length == 0);
}


// ---------------------------------------------------------------------------
// f0 worm channel (scheme C, notes/bilocal-worm-design.md 3.2)
// ---------------------------------------------------------------------------
/*
Extended-ensemble vertex removal / insertion. Closed states T carry
weight e^-S; open states (T, head) carry zeta * e^{-S + U(star(head))}.
U (a frozen table over spoke-degree multisets + linear fallback) and
zeta shape only the auxiliary open sector: the closed-sector
conditional is e^-S exactly for ANY choice, so measurements gated on
closed states are unbiased. One episode = open -> local walk -> close;
a walk that hits the step cap is EXACTLY UNDONE (never committed).

Open-state kernel mixture (fixed weights; invalid draws auto-reject):
  ph   HEAD kernel: uniform over the enumerated class H of moves whose
       support CONTAINS the head. Only these can change a spoke degree
       (an edge (head,u) changes degree only when head is in the move
       support), so dU != 0 exactly on H.
         alpha = e^{-(dS - dU)} * |H| / |H'|
  pg   GLOBAL repair kernel: chooseRandomMove, auto-rejected when the
       head is in the support or the move is vertex-changing; plain
       Metropolis e^{-dS} (dU = 0 identically off the head; the same
       exact=false proposal treatment as run()).
  bcf  close-flag:  alpha = e^{-U} / (zeta f0) * aof/bcf
  bc4  close-41 (Z=4 only):
       alpha = e^{-dS41 - U} / zeta * aoi / (bc4 * f3' * pool')
Openings (one attempt per episode call, from closed):
  aof  open-flag:   alpha = zeta e^{U} f0 * bcf/aof
  aoi  open-insert: alpha = zeta e^{-dS14 + U(3,3,3,3)} f3 * pool
                            * bc4/aoi
The 1<->4 sector crossings maintain the caller's unusedVertices pool
exactly like the targeted-move capi path. Gates (enforced by the capi
wrapper): dim = 3, no cocycle, no geometry ledger (v2 adds CHANNEL_F0
brackets + label epochs per design doc 2.6).
*/

/// Frozen umbrella + kernel mixture for the f0 worm.
struct WormF0Params
{
    double zeta = 0.05;
    double aof = 0.5;        ///< open-flag share of opening attempts
    double ph = 0.45;        ///< head-kernel share of open steps
    double pg = 0.49;        ///< global-repair share
    double bcf = 0.01;       ///< close-flag share
    double bc4 = 0.05;       ///< close-41 share
    int maxstep = 100000;
    int lmax = 4096;
    double z0 = 12.0;
    double[6] ufb = 0.0;     ///< (n3,n4,n5,n6,n7plus,Zdef) fallback
    double ufbc = 0.0;
    double ucapHi = 35.0;    ///< U ceiling: soft confinement (a valid
    double ucapLo = -20.0;   ///< table choice; stops fallback runaway)
    double mu = 1.5;         ///< open-flag seed bias: p(v) ~ e^{-mu Z(v)}
    double zeta2 = 1.0;      ///< PAIR fugacity (two-ball episodes)
    bool zeta2Auto = false;  ///< derive zeta2 from the proposal density
    double pclTarget = 10.0; ///< auto p_close: mean episode length is
                             ///< maxstep/pclTarget, so abandonment is
                             ///< ~e^-pclTarget and the walk uses its
                             ///< budget instead of closing at once
    int chainK = 20;         ///< chord channel: chain windows searched
                             ///< around each defect for creation sites
    double bcp = 0.05;       ///< close-pair share of open steps
    double[ulong] utab;      ///< packed spoke multiset -> U (compiled)
    double[ulong] ctab;      ///< CHORD carrier: packed endpoint-spoke
                             ///< multiset -> U, replayed from measured
                             ///< catalysed paths (see build_chord_tube)
    double cfb = 0.0;        ///< chord off-tube flat value
    // --- AGGREGATION knobs (chord channel). Defaults reproduce the
    // certified behaviour bit-for-bit: regionMax = 0 is the original
    // region-clean test, aggBeta = 0 is the original uniform site draw.
    int regionMax = 0;       ///< other degree-3 edges tolerated in an
                             ///< EMPTY mark's region. The predicate stays
                             ///< symmetric ("whenever a mark is empty its
                             ///< region holds <= regionMax flickers"),
                             ///< applied identically at open and close --
                             ///< only the threshold moves. At 0 genuine
                             ///< aggregation is IMPOSSIBLE: the only
                             ///< tolerated neighbour is the adopted chord,
                             ///< i.e. the one leaving, so no partner ever
                             ///< remains (measured: 89.1% of sites have 0
                             ///< flickers near, 10.7% have exactly 1,
                             ///< 0.25% have >= 2 and can never pass).
    double aggBeta = 0.0;    ///< destination weight w = exp(aggBeta * n),
                             ///< n = distinct flickers whose SUPPORT
                             ///< shares a vertex with the site's (the
                             ///< exact condition under which bilocal
                             ///< factorization fails, so the only geometry
                             ///< in which the two interact). Lives purely
                             ///< in the proposal and cancels in the
                             ///< Hastings ratio, so any value samples the
                             ///< same equilibrium -- it selects the
                             ///< relaxation PATHWAY, not the fixed point.
    WF0Skel[ulong] skel;     ///< f-independent skeleton (see wf0Compile)
    double tubeF1 = 0;       ///< f-vector at skeleton build time
    double tubeF3 = 0;
}

/// f-independent tube skeleton entry: the corridor state's measured
/// cumulative dS at build time, plus its exact (f1, f3) offset from the
/// corridor start ('23' flip: +1,+1; spoke delete: -1,-1 per move).
struct WF0Skel
{
    double val = 0;          ///< cum dS measured at (tubeF1, tubeF3)
    double df1 = 0;
    double df3 = 0;
}

struct WormF0Result
{
    int opened;              ///< 0 rejected, 1 flag, 2 insert
    int head = -1;
    int steps;
    int closedHow;           ///< 0 none, 1 cf, 2 c4, 3 cap-undone
    double dS = 0.0;         ///< committed S change (undone => ~0)
    double umax = 0.0;
    long nH, accH, nG, accG;
    int zmin = 999;          ///< deepest link size reached
    long nZ4;                ///< steps spent at Z == 4
    long[4] df;              ///< per-episode net f change (census)
    // explicit zero-init: D default-initializes floating point to NaN,
    // and every comparison against NaN is false, so `if (d > dsArm[0])`
    // would never fire and the arm split silently reported NaN
    double[2] dsArm = [0.0, 0.0];  ///< strict: max accepted dS [head, global]
    long[2] nUpArm;          ///< strict: accepted uphill count per arm
}

/// One committed move (for exact cap-undo).
struct WF0Applied(Vertex)
{
    Vertex[4] cen;
    Vertex[4] coc;
    int cl, ccl;
}

/// Spoke multiset -> packed bucket counts (degrees 3..9+, 8 bits each).
private ulong wf0Key(scope const(size_t)[] degs) @nogc nothrow
{
    ulong k = 0;
    foreach (d; degs)
    {
        long b = cast(long) d - 3;
        if (b < 0) b = 0;
        if (b > 6) b = 6;
        k += 1UL << (8 * cast(int) b);
    }
    return k;
}

/// The global-pin part of the action as a function of (f1, f3): the only
/// action terms NONLINEAR in the f-vector (volume pin + edge pin). The
/// intensive variance terms are also nonlinear but are off (coef 0) in
/// worm runs; if on, tube flatness degrades slightly -- balance is
/// unaffected either way, since U is a free choice.
double wf0PinPart(P)(P params, real f1, real f3)
{
    immutable real dv = f3 - params.numFacetsTarget;
    immutable real dh = f1 - 6.0L * f3 / params.hingeDegreeTarget;
    return cast(double) (params.numFacetsCoef * dv * dv
        + params.numHingesCoef * dh * dh);
}

/// Compile the frozen scalar table for ONE episode from the stored
/// f-independent skeleton, at the episode-start f-vector (f1, f3):
///   U_k = cum_k + [g(f + d_k) - g(f)] - [g(f0 + d_k) - g(f0)]
/// (g = wf0PinPart, f0 = build-time f, d_k = the state's exact f offset).
/// The m^2 part of cum_k is f-independent (local, orbit-quantized), so
/// this reprices the corridor's pin content at current conditions -- one
/// skeleton serves the whole f0 range. U stays FROZEN within the episode
/// (the dU == 0 lemma for global repair moves survives verbatim); each
/// episode is balanced under its own compiled table and the closed
/// measure is U-independent, so per-episode recompiles are exactly
/// unbiased. Zero-delta skeletons (plain tables) compile to themselves.
void wf0Compile(P)(ref WormF0Params cfg, P params, long f1, long f3)
{
    if (cfg.skel.length == 0) return;
    immutable double gNow = wf0PinPart(params, f1, f3);
    immutable double gRef = wf0PinPart(params, cfg.tubeF1, cfg.tubeF3);
    foreach (k, sk; cfg.skel)
        cfg.utab[k] = sk.val
            + (wf0PinPart(params, f1 + sk.df1, f3 + sk.df3) - gNow)
            - (wf0PinPart(params, cfg.tubeF1 + sk.df1,
                          cfg.tubeF3 + sk.df3) - gRef);
}

private double wf0U(const ref WormF0Params cfg,
                    scope const(size_t)[] degs) nothrow
{
    immutable k = wf0Key(degs);
    if (auto p = k in cfg.utab)
    {
        double u = *p;
        if (u > cfg.ucapHi) u = cfg.ucapHi;
        if (u < cfg.ucapLo) u = cfg.ucapLo;
        return u;
    }
    double[5] n = 0.0;
    foreach (d; degs)
    {
        long b = cast(long) d - 3;
        if (b < 0) b = 0;
        if (b > 4) b = 4;
        n[cast(int) b] += 1.0;
    }
    double u = cfg.ufbc + cfg.ufb[5] * (cfg.z0 - cast(double) degs.length);
    foreach (i; 0 .. 5) u += cfg.ufb[i] * n[i];
    if (u > cfg.ucapHi) u = cfg.ucapHi;
    if (u < cfg.ucapLo) u = cfg.ucapLo;
    return u;
}

/// Head-class candidate (support contains the head).
private struct WF0Cand(Vertex)
{
    bool is23;
    Vertex[3] f;             ///< 2->3 center face | 3->2 link triangle
    Vertex[2] p;             ///< 2->3 apex pair   | 3->2 center edge
}

/// Star of the head from a known seed tet: fills tets, link vertices
/// and spoke degrees. Returns Z, or -1 on buffer overflow.
private int wf0Star(Vertex)(ref Manifold!(3, Vertex) mfd, Vertex v,
    Vertex[4] seed, Vertex[4][] tets, ref int nT,
    Vertex[] lk, size_t[] degs)
{
    nT = collectStar(mfd, v, seed, tets, 0);
    int z = 0;
    foreach (i; 0 .. nT)
        foreach (u; tets[i])
        {
            if (u == v) continue;
            bool seen = false;
            foreach (j; 0 .. z) if (lk[j] == u) { seen = true; break; }
            if (seen) continue;
            if (z >= lk.length) return -1;
            lk[z] = u;
            Vertex[2] e = v < u ? [v, u] : [u, v];
            degs[z] = mfd.degreeOrZero!1(e[]);
            z++;
        }
    return z;
}

/// Enumerate H at the head (optimistic validity: cheap static checks;
/// the apply step re-validates and a failed draw is a null step, the
/// same rule on both sides of any transition).
private int wf0EnumH(Vertex)(ref Manifold!(3, Vertex) mfd, Vertex v,
    scope Vertex[4][] tets, int nT, WF0Cand!Vertex[] outC)
{
    int n = 0;
    Vertex[3][160] seenF;
    int nSF = 0;
    Vertex[2][160] seenE;
    int nSE = 0;
    foreach (ti; 0 .. nT)
    {
        auto T = tets[ti];
        foreach (skip; 0 .. 4)
        {
            Vertex[3] face;
            int k = 0;
            foreach (i; 0 .. 4) if (i != skip) face[k++] = T[i];
            bool dup = false;
            foreach (i; 0 .. nSF) if (seenF[i] == face) { dup = true; break; }
            if (dup) continue;
            if (nSF < seenF.length) seenF[nSF++] = face;
            int[2] ap = 0;
            if (mfd.writeFaceApexes(face[0], face[1], face[2], ap.ptr) != 2)
                continue;
            Vertex[2] axis = ap[0] < ap[1]
                ? [cast(Vertex) ap[0], cast(Vertex) ap[1]]
                : [cast(Vertex) ap[1], cast(Vertex) ap[0]];
            bool hasV = false;
            foreach (x; face) if (x == v) hasV = true;
            foreach (x; axis) if (x == v) hasV = true;
            if (!hasV) continue;
            if (mfd.degreeOrZero!1(axis[]) != 0) continue;
            if (mfd.anyFrozen(face[])) continue;
            if (n < outC.length)
            {
                outC[n].is23 = true;
                outC[n].f = face;
                outC[n].p = axis;
                n++;
            }
        }
        foreach (a; 0 .. 4)
            foreach (b; a + 1 .. 4)
            {
                Vertex[2] e = T[a] < T[b] ? [T[a], T[b]] : [T[b], T[a]];
                bool dup = false;
                foreach (i; 0 .. nSE) if (seenE[i] == e) { dup = true; break; }
                if (dup) continue;
                if (nSE < seenE.length) seenE[nSE++] = e;
                if (mfd.degreeOrZero!1(e[]) != 3) continue;
                int[8] lkB = 0;
                if (mfd.writeEdgeLinkCycle(e[0], e[1], T[], lkB.ptr) != 3)
                    continue;
                Vertex[3] link = [cast(Vertex) lkB[0], cast(Vertex) lkB[1],
                                  cast(Vertex) lkB[2]];
                link[].sort();
                bool hasV = (e[0] == v || e[1] == v);
                foreach (x; link) if (x == v) hasV = true;
                if (!hasV) continue;
                if (n < outC.length)
                {
                    outC[n].is23 = false;
                    outC[n].f = link;
                    outC[n].p = e;
                    n++;
                }
            }
    }
    return n;
}

/// Debug: umbrella value at v's current star (capi probe).
double wormF0DebugU(Vertex)(ref Manifold!(3, Vertex) mfd,
    const ref WormF0Params cfg, Vertex v)
{
    Vertex[4] seed;
    bool ok = mfd.someFacetContaining(v, seed);
    if (!ok) return double.nan;
    Vertex[4][96] tets;
    int nT = 0;
    Vertex[48] lk;
    size_t[48] degs;
    immutable z = wf0Star(mfd, v, seed, tets[], nT, lk[], degs[]);
    if (z < 0) return double.nan;
    return wf0U(cfg, degs[0 .. z]);
}

/// One full f0-worm episode. Returns true if the episode CHANGED the
/// closed state (an f0-changing closure committed); the result struct
/// carries diagnostics either way. currentObjective is kept exact.
bool wormF0Episode(Vertex, P)(ref Manifold!(3, Vertex) mfd,
    ref real currentObjective, ref Vertex[] unusedVertices,
    P params, const ref WormF0Params cfg,
    VertexPotState!Vertex* potState, VertexPot* pot,
    scope WF0Applied!Vertex[] undoBuf, out WormF0Result res)
{
    import std.math : log, exp;
    import std.random : uniform, uniform01;
    alias BM = BistellarMove!(3, Vertex);

    real baseRun = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;
    int nApplied = 0;

    // -- shared apply/undo through the run() commit pipeline ------------
    real applyMove(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        auto bm = BM(cen, coc);
        immutable real dBase =
            mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        return dBase + dPot;
    }

    void record(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        assert(nApplied < undoBuf.length, "wormF0 undo buffer overflow");
        undoBuf[nApplied].cl = cast(int) cen.length;
        undoBuf[nApplied].ccl = cast(int) coc.length;
        undoBuf[nApplied].cen[0 .. cen.length] = cen[];
        undoBuf[nApplied].coc[0 .. coc.length] = coc[];
        nApplied++;
    }

    void unwindAll()
    {
        foreach_reverse (k; 0 .. nApplied)
        {
            auto cen = undoBuf[k].coc[0 .. undoBuf[k].ccl];
            auto coc = undoBuf[k].cen[0 .. undoBuf[k].cl];
            if (undoBuf[k].cl == 4 && undoBuf[k].ccl == 1)
            {
                // undoing an opening 1->4: give the label back
                unusedVertices ~= undoBuf[k].coc[0];
            }
            applyMove(cen, coc);
        }
        nApplied = 0;
    }

    // -- head star / H-class cache --------------------------------------
    Vertex[4][96] starT;
    int nStar = 0;
    Vertex[48] lk;
    size_t[48] degs;
    int z = 0;
    double uCur = 0.0;
    static WF0Cand!Vertex[512] hBuf;
    int nH = 0;
    Vertex head;
    Vertex[4] headSeed;
    // vertices whose state the H cache depends on ({v} u lk u apexes):
    // a move whose support misses this set cannot invalidate the cache
    Vertex[160] hDep;
    int nDep = 0;

    bool refreshHead()
    {
        z = wf0Star(mfd, head, headSeed, starT[], nStar, lk[], degs[]);
        if (z < 0) return false;
        uCur = wf0U(cfg, degs[0 .. z]);
        nH = wf0EnumH(mfd, head, starT[], nStar, hBuf[]);
        nDep = 0;
        void dep(Vertex x)
        {
            foreach (i; 0 .. nDep) if (hDep[i] == x) return;
            if (nDep < hDep.length) hDep[nDep++] = x;
        }
        dep(head);
        foreach (i; 0 .. z) dep(lk[i]);
        foreach (i; 0 .. nH)
        {
            foreach (x; hBuf[i].f) dep(x);
            foreach (x; hBuf[i].p) dep(x);
        }
        return true;
    }

    // exact-restore snapshot of the head cache (rejected H proposals
    // restore the state bit-for-bit, so the cache is restored too)
    Vertex[4][96] snapT;
    int snapNStar, snapZ, snapNH, snapNDep;
    Vertex[48] snapLk;
    size_t[48] snapDegs;
    double snapU;
    static WF0Cand!Vertex[512] snapH;
    Vertex[160] snapDep;
    Vertex[4] snapSeed;

    void saveCache()
    {
        snapT[0 .. nStar] = starT[0 .. nStar];
        snapNStar = nStar;
        snapLk[0 .. z] = lk[0 .. z];
        snapDegs[0 .. z] = degs[0 .. z];
        snapZ = z;
        snapU = uCur;
        snapH[0 .. nH] = hBuf[0 .. nH];
        snapNH = nH;
        snapDep[0 .. nDep] = hDep[0 .. nDep];
        snapNDep = nDep;
        snapSeed = headSeed;
    }

    void restoreCache()
    {
        starT[0 .. snapNStar] = snapT[0 .. snapNStar];
        nStar = snapNStar;
        lk[0 .. snapZ] = snapLk[0 .. snapZ];
        degs[0 .. snapZ] = snapDegs[0 .. snapZ];
        z = snapZ;
        uCur = snapU;
        hBuf[0 .. snapNH] = snapH[0 .. snapNH];
        nH = snapNH;
        hDep[0 .. snapNDep] = snapDep[0 .. snapNDep];
        nDep = snapNDep;
        headSeed = snapSeed;
    }

    // -- opening attempt -------------------------------------------------
    immutable long f0 = cast(long) mfd.fVector[0];
    immutable long f3 = cast(long) mfd.fVector[3];
    immutable double aoi = 1.0 - cfg.aof;

    if (uniform01 < cfg.aof)
    {
        // open-flag: seed biased toward low-Z vertices, p(v) = w/W with
        // w = e^{-mu Z(v)}. Z is O(1) from the tet-degree (link-sphere
        // Euler: Z = 2 + D/2). The W/w Hastings factor replaces f0; the
        // close-flag reverse uses the SAME W (T unchanged by flag moves).
        res.opened = 0;
        double W = 0.0;
        foreach (sv; mfd.simplices(0))
        {
            immutable zz = 2.0
                + 0.5 * cast(double) mfd.degreeOrZero!0(sv);
            W += exp(-cfg.mu * zz);
        }
        immutable double rPick = uniform01 * W;
        double acc = 0.0;
        Vertex v = Vertex.init;
        double wv = 0.0;
        bool got = false;
        foreach (sv; mfd.simplices(0))
        {
            immutable zz = 2.0
                + 0.5 * cast(double) mfd.degreeOrZero!0(sv);
            immutable w = exp(-cfg.mu * zz);
            acc += w;
            if (acc >= rPick) { v = sv.front; wv = w; got = true; break; }
        }
        if (!got) return false;
        // seed tet: any facet containing v (O(N), once per episode)
        immutable bool seedOk = mfd.someFacetContaining(v, headSeed);
        if (!seedOk) return false;
        head = v;
        if (!refreshHead()) return false;
        immutable double la = log(cfg.zeta) + uCur + log(W / wv)
            + log(cfg.bcf / cfg.aof);
        if (!(la >= 0 || uniform01 <= exp(la))) return false;
        res.opened = 1;
    }
    else
    {
        // open-insert: uniform tet, uniform pool label, 1->4
        res.opened = 0;
        auto facet = mfd.randomFacetOfDim(3);
        Vertex[4] tet;
        int i = 0;
        foreach (x; facet) tet[i++] = x;
        tet[].sort();
        // label ~ uniform over [0, lmax); an occupied draw is a null
        // attempt (the same rule prices the reverse side, so the
        // proposal factor is 1/lmax on both legs -- no pool counting)
        auto vn = cast(Vertex) uniform(0, cfg.lmax);
        auto bm14 = BM(tet[], vn.only);
        if (!mfd.hasValidMove(bm14)) return false;
        immutable real dB = mfd.speculativeBistellarDelta(bm14, baseRun,
                                                          params);
        real dP = 0.0L;
        if (potState !is null)
            dP = mfd.potentialBistellarDelta(bm14, *potState, *pot, false);
        size_t[4] freshDegs = [3, 3, 3, 3];
        immutable double u14 = wf0U(cfg, freshDegs[]);
        immutable double la = log(cfg.zeta) - cast(double)(dB + dP) + u14
            + log(cast(double) f3) + log(cast(double) cfg.lmax)
            + log(cfg.bc4 / aoi);
        if (!(la >= 0 || uniform01 <= exp(la))) return false;
        Vertex[1] vnA = [vn];
        applyMove(tet[], vnA[]);
        record(tet[], vnA[]);
        // consume the label from the tracked pool if present
        foreach (j; 0 .. unusedVertices.length)
            if (unusedVertices[j] == vn)
            {
                unusedVertices[j] = unusedVertices[$ - 1];
                unusedVertices = unusedVertices[0 .. $ - 1];
                unusedVertices.assumeSafeAppend;
                break;
            }
        head = vn;
        headSeed = [tet[0], tet[1], tet[2], vn];
        headSeed[].sort();
        if (!refreshHead()) { unwindAll(); return false; }
        res.opened = 2;
    }
    res.head = cast(int) head;
    res.umax = uCur;

    // -- open-sector walk ------------------------------------------------
    immutable double pcum1 = cfg.ph;
    immutable double pcum2 = cfg.ph + cfg.pg;
    immutable double pcum3 = cfg.ph + cfg.pg + cfg.bcf;

    while (true)
    {
        if (res.steps >= cfg.maxstep)
        {
            unwindAll();
            res.closedHow = 3;
            res.dS = cast(double) deltaTotal;   // must be ~0
            currentObjective += deltaTotal;
            return false;
        }
        res.steps++;
        if (z < res.zmin) res.zmin = z;
        if (z == 4) res.nZ4++;
        immutable double r = uniform01;
        if (r < pcum1)
        {
            // HEAD kernel
            res.nH++;
            if (nH == 0) continue;
            auto c = hBuf[uniform(0, nH)];
            BM bm = c.is23 ? BM(c.f[], c.p[]) : BM(c.p[], c.f[]);
            if (!mfd.hasValidMove(bm)) continue;
            immutable int nH0 = nH;
            immutable double u0 = uCur;
            saveCache();
            immutable real d = applyMove(c.is23 ? c.f[] : c.p[],
                                         c.is23 ? c.p[] : c.f[]);
            // head seed survives iff it still exists; recompute from
            // the move: pick any post-move tet containing the head
            Vertex[4] newSeed = headSeed;
            {
                // post tets of the move that contain the head
                Vertex[5] sup;
                int nsup = 0;
                foreach (x; c.f) sup[nsup++] = x;
                foreach (x; c.p) sup[nsup++] = x;
                bool inSup = false;
                foreach (x; sup[0 .. nsup]) if (x == head) inSup = true;
                if (inSup)
                {
                    if (c.is23)
                    {
                        // new tets: (f_i, f_j, p0, p1)
                        outer23: foreach (a; 0 .. 3) foreach (b; a + 1 .. 3)
                        {
                            Vertex[4] t = [c.f[a], c.f[b], c.p[0], c.p[1]];
                            bool has = false;
                            foreach (x; t) if (x == head) has = true;
                            if (has) { t[].sort(); newSeed = t; break outer23; }
                        }
                    }
                    else
                    {
                        // new tets: (link tri, e0) and (link tri, e1)
                        foreach (bb; 0 .. 2)
                        {
                            Vertex[4] t = [c.f[0], c.f[1], c.f[2], c.p[bb]];
                            bool has = false;
                            foreach (x; t) if (x == head) has = true;
                            if (has) { t[].sort(); newSeed = t; break; }
                        }
                    }
                }
                headSeed = newSeed;
            }
            if (!refreshHead())
            {
                // pathological: reject and restore
                applyMove(bm.coCenter, bm.center);
                restoreCache();
                continue;
            }
            immutable double la = -cast(double) d + (uCur - u0)
                + log(cast(double) nH0) - log(cast(double)(nH == 0 ? 1 : nH));
            if (nH != 0 && (la >= 0 || uniform01 <= exp(la)))
            {
                res.accH++;
                record(c.is23 ? c.f[] : c.p[], c.is23 ? c.p[] : c.f[]);
            }
            else
            {
                applyMove(bm.coCenter, bm.center);
                restoreCache();   // state restored exactly => cache too
            }
            if (uCur > res.umax) res.umax = uCur;
        }
        else if (r < pcum2)
        {
            // GLOBAL repair kernel (plain Metropolis; dU = 0 off-head)
            res.nG++;
            Vertex fresh = unusedVertices.length > 0
                ? unusedVertices[$ - 1]
                : cast(Vertex) mfd.fVector[0];
            auto bm = mfd.chooseRandomMove(fresh, params);
            if (bm.center.length == 1 || bm.coCenter.length == 1)
                continue;                    // no 1<->4 inside episodes
            bool touches = false;
            foreach (x; bm.center) if (x == head) touches = true;
            foreach (x; bm.coCenter) if (x == head) touches = true;
            if (touches) continue;           // head moves belong to H
            immutable real dB = mfd.speculativeBistellarDelta(bm, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm, *potState, *pot,
                                                 false);
            immutable double la = -cast(double)(dB + dP);
            if (la >= 0 || uniform01 <= exp(la))
            {
                res.accG++;
                applyMove(bm.center, bm.coCenter);
                record(bm.center, bm.coCenter);
                // refresh only if the move can see the H cache's deps
                bool hit = false;
                foreach (x; bm.center)
                    foreach (i; 0 .. nDep) if (hDep[i] == x) hit = true;
                foreach (x; bm.coCenter)
                    foreach (i; 0 .. nDep) if (hDep[i] == x) hit = true;
                if (hit) refreshHead();
            }
        }
        else if (r < pcum3)
        {
            // close-flag: reverse of the biased open-flag draw
            double W = 0.0;
            foreach (sv; mfd.simplices(0))
            {
                immutable zz = 2.0
                    + 0.5 * cast(double) mfd.degreeOrZero!0(sv);
                W += exp(-cfg.mu * zz);
            }
            immutable double zh = 2.0
                + 0.5 * cast(double) mfd.degreeOrZero!0(head.only);
            immutable double wv = exp(-cfg.mu * zh);
            immutable double la = -uCur - log(cfg.zeta)
                - log(W / wv) + log(cfg.aof / cfg.bcf);
            if (la >= 0 || uniform01 <= exp(la))
            {
                res.closedHow = 1;
                res.dS = cast(double) deltaTotal;
                currentObjective += deltaTotal;
                return nApplied > 0;
            }
        }
        else
        {
            // close-41
            if (z != 4) continue;
            Vertex[4] nb = [lk[0], lk[1], lk[2], lk[3]];
            nb[].sort();
            auto bm41 = BM(head.only, nb[]);
            if (!mfd.hasValidMove(bm41)) continue;
            immutable real dB = mfd.speculativeBistellarDelta(bm41, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm41, *potState, *pot,
                                                 false);
            immutable long f3after = cast(long) mfd.fVector[3] - 3;
            immutable double la = -cast(double)(dB + dP) - uCur
                - log(cfg.zeta) + log(aoi / cfg.bc4)
                - log(cast(double) f3after)
                - log(cast(double) cfg.lmax);
            if (la >= 0 || uniform01 <= exp(la))
            {
                Vertex[1] hA = [head];
                applyMove(hA[], nb[]);
                unusedVertices ~= head;
                res.closedHow = 2;
                res.dS = cast(double) deltaTotal;
                currentObjective += deltaTotal;
                return true;
            }
        }
    }
    assert(0);
}

/*
BILOCAL (two-ball) episodes.
=============================

The single-head episode above opens ONE ball, walks it under the frozen
umbrella, and closes. A bilocal move is the same structure with TWO
balls, and it exists because a pair whose NET f-change vanishes pays no
global pin at all: the pin cost is g(f + df_A + df_B) - g(f), identically
zero when df_B = -df_A, at ANY separation. Measured: the single-head 1->4
pays +22.75 in pin alone at gap +10, the conserving pair pays 0.0e+00.

Roles. Head 0 is opened by a 1->4 (a fresh vertex) and closes by being
DROPPED (it stays); head 1 is opened by FLAGGING an existing vertex and
closes by 4->1 (it is deleted). Net over the episode: one vertex created
at 0, one destroyed at 1, f0 unchanged -- vertex transport. Running the
family with the roles swapped gives the reverse transition, so the family
is its own inverse and the two closures pair crosswise, exactly as
open-flag/close-flag and open-insert/close-41 do for one head.

Factorization. Everything factorizes when the two balls' vertex supports
are DISJOINT, and only then: measured residual 1e-13 at graph distance
>= 1 versus up to +29.4 when the stars share a vertex. Hence
  U(h0,h1) = U(h0) + U(h1)   -- the SAME single-head tube serves both,
  alpha_open  = alpha_openinsert(0) * alpha_openflag(1)  with one zeta2,
  alpha_close = alpha_close41(1)    * alpha_closeflag(0) with one zeta2.
Disjointness is enforced as a STATE FUNCTION (reject any move whose
support meets the other head's closed star), so it never breaks balance;
no minimum separation is needed, and none is imposed -- transport cost is
range-independent (measured flat out to 60 chain steps).

The dU == 0 lemma survives with two exclusions: a global repair move
whose support contains neither head leaves both umbrella values fixed, so
that kernel stays plain Metropolis.
*/

/// Per-head cache for a two-ball episode (the single-head episode keeps
/// the same state in locals; here it is indexed so two can coexist).
private struct WF0Ball(Vertex)
{
    Vertex head;
    Vertex[4] seed;
    int nT;
    int z;
    double u = 0.0;
    int nH;
    int nDep;
    bool isInsert;           ///< opened by 1->4 (role 0) vs flagged
}

/// One bilocal (two-ball) episode. Returns true if the closed state
/// changed. currentObjective is kept exact; a capped walk is unwound
/// exactly, like the single-head episode.
bool wormPairEpisode(Vertex, P)(ref Manifold!(3, Vertex) mfd,
    ref real currentObjective, ref Vertex[] unusedVertices,
    P params, const ref WormF0Params cfg,
    VertexPotState!Vertex* potState, VertexPot* pot,
    scope WF0Applied!Vertex[] undoBuf, out WormF0Result res)
{
    import std.math : log, exp;
    import std.random : uniform, uniform01;
    alias BM = BistellarMove!(3, Vertex);

    // big per-ball buffers live in TLS, not on the stack
    static Vertex[4][96][2] bTets;
    static Vertex[48][2] bLk;
    static size_t[48][2] bDegs;
    static WF0Cand!Vertex[512][2] bH;
    static Vertex[160][2] bDep;
    WF0Ball!Vertex[2] ball;

    // ROLE-SIGNED UMBRELLA. The tube is a one-way ramp (U +3.85 at an
    // ordinary z=11 star rising to +20.60 at the collapsed (3,3,3,3)),
    // which is right for a single collapsing ball but wrong for a PAIR:
    // transport needs one star to collapse and the other to GROW, and a
    // single ramp drives both to z=4. Measured consequence with one
    // ramp: 58045 visits to z=4, 100% of opened episodes hitting the
    // step cap, and ZERO closes -- the close subtracted u_del + u_keep =
    // 20.60 + 20.60 while the open had only added 20.60 + 0.85, a ~19.8
    // deficit it can never recover. Giving the CREATED ball -U and the
    // ADOPTED ball +U makes the two balls' umbrellas simply SWAP over
    // the episode (A goes 11->4 as B goes 4->11), so the U content is
    // conserved open-to-close and both legs sit at O(1).
    double[2] usgn = [1.0, 1.0];

    real baseRun = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;
    int nApplied = 0;

    real applyMove(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        auto bm = BM(cen, coc);
        immutable real dBase =
            mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        return dBase + dPot;
    }

    void record(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        assert(nApplied < undoBuf.length, "wormPair undo buffer overflow");
        undoBuf[nApplied].cl = cast(int) cen.length;
        undoBuf[nApplied].ccl = cast(int) coc.length;
        undoBuf[nApplied].cen[0 .. cen.length] = cen[];
        undoBuf[nApplied].coc[0 .. coc.length] = coc[];
        nApplied++;
    }

    void unwindAll()
    {
        foreach_reverse (k; 0 .. nApplied)
        {
            auto cen = undoBuf[k].coc[0 .. undoBuf[k].ccl];
            auto coc = undoBuf[k].cen[0 .. undoBuf[k].cl];
            if (undoBuf[k].cl == 4 && undoBuf[k].ccl == 1)
                unusedVertices ~= undoBuf[k].coc[0];
            applyMove(cen, coc);
        }
        nApplied = 0;
    }

    bool refresh(int i)
    {
        ball[i].z = wf0Star(mfd, ball[i].head, ball[i].seed,
                            bTets[i][], ball[i].nT, bLk[i][], bDegs[i][]);
        if (ball[i].z < 0) return false;
        // signed by role, so the walk's +(u - uOld) drive pushes the
        // adopted ball DOWN the ramp and the created ball UP it
        ball[i].u = usgn[i] * wf0U(cfg, bDegs[i][0 .. ball[i].z]);
        ball[i].nH = wf0EnumH(mfd, ball[i].head, bTets[i][], ball[i].nT,
                              bH[i][]);
        ball[i].nDep = 0;
        void dep(Vertex x)
        {
            foreach (k; 0 .. ball[i].nDep) if (bDep[i][k] == x) return;
            if (ball[i].nDep < bDep[i].length)
                bDep[i][ball[i].nDep++] = x;
        }
        dep(ball[i].head);
        foreach (k; 0 .. ball[i].z) dep(bLk[i][k]);
        foreach (k; 0 .. ball[i].nH)
        {
            foreach (x; bH[i][k].f) dep(x);
            foreach (x; bH[i][k].p) dep(x);
        }
        return true;
    }

    /// closed star of ball i (head + its link) -- the disjointness set
    bool inStar(int i, Vertex x)
    {
        if (x == ball[i].head) return true;
        foreach (k; 0 .. ball[i].z) if (bLk[i][k] == x) return true;
        return false;
    }

    /// a move's support may not meet the OTHER ball's closed star
    bool supportClear(int self, scope const(Vertex)[] a,
                      scope const(Vertex)[] b)
    {
        immutable o = 1 - self;
        foreach (x; a) if (inStar(o, x)) return false;
        foreach (x; b) if (inStar(o, x)) return false;
        return true;
    }

    // -- biased-seed normalizer (shared by open-flag and close-flag) -----
    double seedWeightTotal()
    {
        double W = 0.0;
        foreach (sv; mfd.simplices(0))
            W += exp(-cfg.mu
                     * (2.0 + 0.5 * cast(double) mfd.degreeOrZero!0(sv)));
        return W;
    }

    // -- OPEN the pair ---------------------------------------------------
    // p_close is the ACTUAL probability of the close branch in the step
    // mixture below; it must be the same number in both acceptance
    // ratios or open/close balance is off by their ratio. (cfg.bcp is
    // the requested share; the config validates ph + pg + bcp == 1, so
    // this equals it -- computed here so the two can never drift.)
    // p_close IS cfg.bcp, with ph/pg rescaled to make room for it. It
    // used to be 1 - ph - pg, which made bcp inert (the in-code claim
    // that "the config validates ph + pg + bcp == 1" was never true --
    // the f0 config validates ph + pg + bcf + bc4 <= 1 and knows
    // nothing about bcp). The value matters: at 1 - ph - pg = 0.10 the
    // mean episode is 10 steps, so the close always fires while the
    // created ball is still the bare z=4 spike it opened as, and EVERY
    // closure is a roundtrip. Transport needs the corridor to collapse
    // the adopted ball first, which takes hundreds of steps.
    immutable double pcl = cfg.bcp;
    if (!(pcl > 0.0 && pcl < 1.0)) return false;
    if (!(cfg.ph + cfg.pg > 0.0)) return false;
    immutable double pmix = (1.0 - pcl) / (cfg.ph + cfg.pg);
    immutable double phE = cfg.ph * pmix;
    immutable double pgE = cfg.pg * pmix;
    // FREE MODES: which ball creates (1->4) and which adopts (flag) is
    // drawn here, and which ball is DELETED is drawn again at the close.
    // Both are uniform over the two balls, so the 1/2's cancel between
    // the forward close and the reverse open -- the family is then
    // manifestly self-inverse instead of self-inverse-by-relabelling.
    // Deleting the ADOPTED ball is transport (a vertex moves); deleting
    // the CREATED one is a round trip that commits only corridor work.
    immutable int cre = uniform(0, 2);
    immutable int adp = 1 - cre;
    usgn[cre] = -1.0;
    usgn[adp] = 1.0;
    immutable long f3 = cast(long) mfd.fVector[3];
    res.opened = 0;

    // ball 0: 1->4 at a uniform tet with a uniform label
    auto facet = mfd.randomFacetOfDim(3);
    Vertex[4] tet;
    {
        int i = 0;
        foreach (x; facet) tet[i++] = x;
    }
    tet[].sort();
    auto vn = cast(Vertex) uniform(0, cfg.lmax);
    auto bm14 = BM(tet[], vn.only);
    if (!mfd.hasValidMove(bm14)) return false;

    // ball 1: biased seed over existing vertices, p(v) = w/W
    immutable double W = seedWeightTotal();
    immutable double rPick = uniform01 * W;
    double acc = 0.0, wv = 0.0;
    Vertex vf = Vertex.init;
    bool got = false;
    foreach (sv; mfd.simplices(0))
    {
        immutable double w = exp(-cfg.mu
            * (2.0 + 0.5 * cast(double) mfd.degreeOrZero!0(sv)));
        acc += w;
        if (acc >= rPick) { vf = sv.front; wv = w; got = true; break; }
    }
    if (!got) return false;
    // the flagged vertex must not sit inside the insert's own tet
    foreach (x; tet) if (x == vf) return false;

    immutable real dB14 = mfd.speculativeBistellarDelta(bm14, baseRun,
                                                        params);
    real dP14 = 0.0L;
    if (potState !is null)
        dP14 = mfd.potentialBistellarDelta(bm14, *potState, *pot, false);
    size_t[4] freshDegs = [3, 3, 3, 3];
    // the fresh vertex belongs to the CREATED ball, so it carries that
    // role's sign -- it opens at the BOTTOM of its own ramp and the
    // walk is rewarded for growing its star back out
    immutable double u0 = -wf0U(cfg, freshDegs[]);

    // ball 1's umbrella BEFORE any move (flagging changes nothing)
    {
        immutable bool seedOk =
            mfd.someFacetContaining(vf, ball[adp].seed);
        if (!seedOk) return false;
        ball[adp].head = vf;
        ball[adp].isInsert = false;
        if (!refresh(adp)) return false;
    }
    immutable double u1 = ball[adp].u;

    // alpha_open = alpha_openinsert(0) * alpha_openflag(1), one zeta2
    immutable double laOpen = log(cfg.zeta2) - cast(double)(dB14 + dP14)
        + u0 + u1 + log(cast(double) f3) + log(cast(double) cfg.lmax)
        + log(W / wv) + log(pcl);
    if (!(laOpen >= 0 || uniform01 <= exp(laOpen))) return false;

    Vertex[1] vnA = [vn];
    applyMove(tet[], vnA[]);
    record(tet[], vnA[]);
    foreach (j; 0 .. unusedVertices.length)
        if (unusedVertices[j] == vn)
        {
            unusedVertices[j] = unusedVertices[$ - 1];
            unusedVertices = unusedVertices[0 .. $ - 1];
            unusedVertices.assumeSafeAppend;
            break;
        }
    ball[cre].head = vn;
    ball[cre].isInsert = true;
    ball[cre].seed = [tet[0], tet[1], tet[2], vn];
    ball[cre].seed[].sort();
    if (!refresh(cre) || !refresh(adp)) { unwindAll(); return false; }
    // the two balls must be support-disjoint from the start
    foreach (k; 0 .. ball[cre].z)
        if (inStar(adp, bLk[cre][k])) { unwindAll(); return false; }
    if (inStar(adp, ball[cre].head)) { unwindAll(); return false; }
    res.opened = 2;
    res.head = cast(int) ball[adp].head;
    res.umax = ball[cre].u + ball[adp].u;
    double laBest = -1e300;   // best close log-alpha (diagnostic)

    // -- two-ball walk ---------------------------------------------------
    immutable double pcum0 = 0.5 * phE;
    immutable double pcum1 = phE;
    immutable double pcum2 = phE + pgE;

    while (true)
    {
        if (res.steps >= cfg.maxstep)
        {
            unwindAll();
            res.closedHow = 3;
            res.dS = cast(double) deltaTotal;
            currentObjective += deltaTotal;
            return false;
        }
        res.steps++;
        // track the ADOPTED ball: transport closes by deleting IT, so
        // whether it reaches z=4 is the question. (This used to track
        // ball[1], an arbitrary index -- with the roles drawn per
        // episode that was the created ball half the time, which sits
        // at z=4 from birth and made the counter look healthy while
        // transport was never actually reachable.)
        if (ball[adp].z < res.zmin) res.zmin = ball[adp].z;
        if (ball[adp].z == 4) res.nZ4++;
        immutable double r = uniform01;

        if (r < pcum1)
        {
            // HEAD kernel on ball i
            immutable int i = (r < pcum0) ? 0 : 1;
            res.nH++;
            if (ball[i].nH == 0) continue;
            auto c = bH[i][uniform(0, ball[i].nH)];
            if (!supportClear(i, c.f[], c.p[])) continue;
            BM bm = c.is23 ? BM(c.f[], c.p[]) : BM(c.p[], c.f[]);
            if (!mfd.hasValidMove(bm)) continue;
            immutable int nH0 = ball[i].nH;
            immutable double uOld = ball[i].u;
            auto seedOld = ball[i].seed;
            immutable real d = applyMove(c.is23 ? c.f[] : c.p[],
                                         c.is23 ? c.p[] : c.f[]);
            // re-seed from a post-move tet containing this head
            {
                bool inSup = false;
                foreach (x; c.f) if (x == ball[i].head) inSup = true;
                foreach (x; c.p) if (x == ball[i].head) inSup = true;
                if (inSup)
                {
                    if (c.is23)
                    {
                        outer23: foreach (a; 0 .. 3) foreach (b2; a + 1 .. 3)
                        {
                            Vertex[4] t = [c.f[a], c.f[b2], c.p[0], c.p[1]];
                            foreach (x; t) if (x == ball[i].head)
                            { t[].sort(); ball[i].seed = t; break outer23; }
                        }
                    }
                    else
                    {
                        foreach (bb; 0 .. 2)
                        {
                            Vertex[4] t = [c.f[0], c.f[1], c.f[2], c.p[bb]];
                            bool has = false;
                            foreach (x; t) if (x == ball[i].head) has = true;
                            if (has) { t[].sort(); ball[i].seed = t; break; }
                        }
                    }
                }
            }
            if (!refresh(i))
            {
                applyMove(bm.coCenter, bm.center);
                ball[i].seed = seedOld;
                refresh(i);
                continue;
            }
            immutable double la = -cast(double) d + (ball[i].u - uOld)
                + log(cast(double) nH0)
                - log(cast(double)(ball[i].nH == 0 ? 1 : ball[i].nH));
            if (ball[i].nH != 0 && (la >= 0 || uniform01 <= exp(la)))
            {
                res.accH++;
                record(c.is23 ? c.f[] : c.p[], c.is23 ? c.p[] : c.f[]);
            }
            else
            {
                applyMove(bm.coCenter, bm.center);
                ball[i].seed = seedOld;
                refresh(i);
            }
            // (umax is repurposed as the best close log-alpha)
        }
        else if (r < pcum2)
        {
            // GLOBAL repair: dU = 0 iff the support misses BOTH heads
            res.nG++;
            Vertex fresh = unusedVertices.length > 0
                ? unusedVertices[$ - 1]
                : cast(Vertex) mfd.fVector[0];
            auto bm = mfd.chooseRandomMove(fresh, params);
            if (bm.center.length == 1 || bm.coCenter.length == 1) continue;
            bool touches = false;
            foreach (x; bm.center)
                if (x == ball[0].head || x == ball[1].head) touches = true;
            foreach (x; bm.coCenter)
                if (x == ball[0].head || x == ball[1].head) touches = true;
            if (touches) continue;
            immutable real dB = mfd.speculativeBistellarDelta(bm, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm, *potState, *pot,
                                                 false);
            immutable double la = -cast(double)(dB + dP);
            if (la >= 0 || uniform01 <= exp(la))
            {
                res.accG++;
                applyMove(bm.center, bm.coCenter);
                record(bm.center, bm.coCenter);
                foreach (i; 0 .. 2)
                {
                    bool hit = false;
                    foreach (x; bm.center)
                        foreach (k; 0 .. ball[i].nDep)
                            if (bDep[i][k] == x) hit = true;
                    foreach (x; bm.coCenter)
                        foreach (k; 0 .. ball[i].nDep)
                            if (bDep[i][k] == x) hit = true;
                    if (hit) refresh(i);
                }
            }
        }
        else
        {
            // CLOSE the pair: 4->1 at one ball (drawn), drop the other.
            // alpha = alpha_close41(del) * alpha_closeflag(keep), one
            // zeta2. Either ball may be the deleted one: deleting the
            // adopted ball transports a vertex, deleting the created
            // ball is a round trip that keeps only the corridor work.
            immutable int del = uniform(0, 2);
            immutable int keep = 1 - del;
            if (ball[del].z != 4) continue;
            Vertex[4] nb = [bLk[del][0], bLk[del][1],
                            bLk[del][2], bLk[del][3]];
            nb[].sort();
            Vertex[1] h1 = [ball[del].head];
            if (!supportClear(del, h1[], nb[])) continue;
            auto bm41 = BM(h1[], nb[]);
            if (!mfd.hasValidMove(bm41)) continue;
            immutable real dB = mfd.speculativeBistellarDelta(bm41, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm41, *potState, *pot,
                                                 false);
            // the reverse open would re-insert at ball 1's site and flag
            // ball 0's vertex: f3 and the seed normalizer are POST-close
            immutable long f3after = cast(long) mfd.fVector[3] - 3;
            immutable double Wc = seedWeightTotal();
            immutable double wc = exp(-cfg.mu * (2.0 + 0.5
                * cast(double) mfd.degreeOrZero!0(ball[keep].head.only)));
            // ROLE SWAP. Detailed balance pairs this close with the
            // REVERSE episode's open, and in that episode the roles are
            // exchanged: the ball deleted here is the one re-CREATED
            // there, and the ball kept here is the one ADOPTED there.
            // So the U content to cancel is (-1)*Uraw(del) +
            // (+1)*Uraw(keep), and la subtracts it. Pricing under the
            // forward roles instead would put the close at +16.75 where
            // the reverse open sits at -16.75 -- balance broken by twice
            // the umbrella. Uraw = usgn * ball.u recovers the unsigned
            // tube value (usgn^2 == 1).
            immutable double uRawDel = usgn[del] * ball[del].u;
            immutable double uRawKeep = usgn[keep] * ball[keep].u;
            immutable double la = -cast(double)(dB + dP)
                + uRawDel - uRawKeep - log(cfg.zeta2)
                - log(cast(double) f3after) - log(cast(double) cfg.lmax)
                - log(Wc / wc) - log(pcl);
            // DIAGNOSTIC: the best close attempt, and its ingredients.
            // Kept in its OWN accumulator: res.umax is already seeded at
            // the open with ball[cre].u + ball[adp].u, so folding close
            // attempts into it reports max(la, -19.75) and buries every
            // attempt worse than the open's U sum.
            if (la > laBest)
            {
                laBest = la;
                res.df[0] = cast(long)(cast(double)(dB + dP) * 100.0);
                res.df[1] = cast(long)(uRawKeep * 100.0);
                res.df[2] = cast(long)(log(Wc / wc) * 100.0);
                res.df[3] = ball[keep].z;
                res.umax = laBest;
            }
            if (la >= 0 || uniform01 <= exp(la))
            {
                applyMove(h1[], nb[]);
                unusedVertices ~= ball[del].head;
                res.closedHow = (del == adp) ? 4 : 6;   // 4 = transport
                res.dS = cast(double) deltaTotal;
                currentObjective += deltaTotal;
                return true;
            }
        }
    }
    assert(0);
}

/*
CHORD (2<->3) bilocal channel -- the flicker carrier.
=====================================================

A single 2->3 on pristine crystal IS the elementary (3,4,4) knot
(f = (5,10,9,3) always, ~95% of thermal births at lam=0.40), so the
"create" mode here is literally the flicker. The channel gives two
physical readings from one kernel, told apart only by which closure
fires:

  create -> annihilate at the SAME ball, with work in between
      the flicker appears, enables moves that were not available, and
      vanishes: CATALYSIS (it does not appear in the net reaction);
  create -> drop at ball 0 + adopt -> annihilate at ball 1
      the flicker relocates: TRANSPORT of a disclination quantum.

Net f-change of the transport pairing is (0,+1,+2,+1) + (0,-1,-2,-1) = 0
exactly, so it is pin-free at ANY separation, like the vertex pairing.

Unlike the vertex carrier there is NO umbrella: a chord's closure
condition (degree exactly 3) is what a 2->3 creates, so it holds the
instant the flicker is born -- nothing to flatten. Hence U == 0 here,
dU == 0 for EVERY move, and the walk is plain Metropolis. The only
moves excluded are those that would annihilate a head chord outside the
close (i.e. centered on it), which would destroy the head.

Proposal densities. A face is drawn as (uniform tet, uniform face index):
each face lies in exactly two tets, so it is hit with probability
2/(4 f3) = 1/f2 -- uniform over faces, using f2 = 2 f3 exactly. The
adopted chord is drawn uniformly from the degree-3 edges (count n3). The
reverse of the close is an open that CREATES at ball 1's site and ADOPTS
ball 0's chord, so the close carries the post-move f2' and n3'.
*/

/// Link vertices of an edge (its degree, and the link cycle's vertices).
/// Returns the count, or -1 on overflow.
private int wf0EdgeLink(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex a, Vertex b, Vertex[] outv)
{
    Vertex[2] e = a < b ? [a, b] : [b, a];
    int n = 0;
    foreach (f; mfd.star(e[]))
        foreach (x; f)
        {
            if (x == a || x == b) continue;
            bool seen = false;
            foreach (i; 0 .. n) if (outv[i] == x) { seen = true; break; }
            if (seen) continue;
            if (n >= outv.length) return -1;
            outv[n++] = x;
        }
    return n;
}

/// Count degree-3 edges; if pick >= 0, also return the pick-th one.
private long wf0Deg3Scan(Vertex)(ref Manifold!(3, Vertex) mfd,
    long pick, Vertex* outA, Vertex* outB)
{
    long n = 0;
    foreach (ed; mfd.simplices(1))
    {
        Vertex[2] e;
        int i = 0;
        foreach (x; ed) { if (i < 2) e[i] = x; i++; }
        if (i != 2) continue;
        if (mfd.degreeOrZero!1(e[]) != 3) continue;
        if (pick >= 0 && n == pick && outA !is null)
        { *outA = e[0]; *outB = e[1]; }
        n++;
    }
    return n;
}

/// Enumerate valid 2->3 sites lying within `kmax` chain windows of SOME
/// degree-3 chord, and return the count; with pick >= 0 also return the
/// pick-th site. This is the chord channel's creation-site proposal.
///
/// A flicker born in pristine crystal has nothing to act on -- the
/// crystal is a local minimum, so it simply relocates (measured: 42 of
/// 44 commits were pure transport). Born next to a defect it can react.
/// The creation-site distribution therefore SELECTS WHICH PHYSICS the
/// channel samples: pristine -> transport, defect-adjacent -> catalysis.
///
/// The set is deliberately keyed on "near ANY defect", not on the chord
/// this episode adopted: that makes it a pure STATE FUNCTION, so the
/// count is recomputable identically after the move and the Hastings
/// ratio |C|/|C'| is exact. Walking outward from the adopted chord would
/// make the set proposal-dependent, and its reverse would depend on the
/// chain still connecting the two sites after both ends were modified.
///
/// Sites reachable from several chords are enumerated once per route;
/// that multiplicity is a proposal weight, identical on both legs, so it
/// cancels. Cost is n3 * 3 * kmax windows, each O(1).
private long wf0ChainSites(Vertex)(ref Manifold!(3, Vertex) mfd,
    int kmax, long pick, Vertex[3]* outFace, Vertex[2]* outAx)
{
    long n = 0;
    // Two phases. Phase 1 sweeps facets and COLLECTS the start windows
    // (a tet qualifies iff one of its six edges has degree 3); phase 2
    // walks them. They must not interleave: writeFaceApexes performs a
    // ridge-link lookup, and doing that inside a simplices(3) traversal
    // corrupts the traversal. mfd.star() is avoided entirely -- it
    // allocates, and this runs twice per episode.
    static Vertex[4][4096] starts;   // zero-init, never = void
    int nStart = 0;
    foreach (f; mfd.simplices(3))
    {
        Vertex[4] w;
        {
            int k = 0;
            foreach (x; f) { if (k < 4) w[k] = x; k++; }
            if (k != 4) continue;
        }
        w[].sort();                           // canonical start window
        bool near = false;
        foreach (i2; 0 .. 4)
            foreach (j2; i2 + 1 .. 4)
            {
                Vertex[2] e = [w[i2], w[j2]];
                if (mfd.degreeOrZero!1(e[]) == 3) near = true;
            }
        if (!near) continue;
        if (nStart >= starts.length) break;
        starts[nStart++] = w;
    }
    foreach (si; 0 .. nStart)
    {
        Vertex[4] w = starts[si];
        {
            foreach (step; 0 .. kmax)
            {
                int[2] ap = 0;
                if (mfd.writeFaceApexes(w[1], w[2], w[3], ap.ptr) != 2)
                    break;
                Vertex[2] axis = ap[0] < ap[1]
                    ? [cast(Vertex) ap[0], cast(Vertex) ap[1]]
                    : [cast(Vertex) ap[1], cast(Vertex) ap[0]];
                // CANONICALIZE: the window is in chain-walk order, but
                // moves require sorted center/coCenter -- an unsorted
                // face corrupts the ridge bookkeeping downstream.
                Vertex[3] face = [w[1], w[2], w[3]];
                face[].sort();
                if (mfd.degreeOrZero!1(axis[]) == 0
                    && !mfd.anyFrozen(face[]))
                {
                    if (pick >= 0 && n == pick && outFace !is null)
                    { *outFace = face; *outAx = axis; }
                    n++;
                }
                // advance: drop w[0], adopt the apex opposite it
                Vertex nxt = (ap[0] == w[0]) ? cast(Vertex) ap[1]
                                             : cast(Vertex) ap[0];
                if (ap[0] != w[0] && ap[1] != w[0]) break;
                w = [w[1], w[2], w[3], nxt];
            }
        }
    }
    return n;
}

/// The chord's local signature: the degrees of every edge at either
/// endpoint. This is the chord analogue of the vertex head's spoke
/// multiset, and it is the right key because U must change exactly when
/// the wide head class fires: an edge (x,u) changes degree only if a tet
/// on it is touched, which needs BOTH x and u in the support, so U
/// changes iff an endpoint is in the support. That restores the
/// dU == 0 lemma for every move outside the head class.
private int wf0ChordDegs(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] chord, scope Vertex[4][] tets, int nT, size_t[] outD,
    out int n3reg)
{
    // The region's degree-3 edge count -- i.e. IS A HELPER PRESENT. The
    // degree multiset alone cannot tell a bare flicker from a catalysed
    // end-state (measured: 18 of 40 fresh flickers aliased onto tube
    // entries, inheriting catalysed prices as low as -5.96), and that
    // aliasing is what forced one U value to serve two incompatible
    // roles. This is exactly the variable the catalysed path turns on.
    n3reg = 0;
    {
        Vertex[2][512] seenE;
        int nSE = 0;
        foreach (ti; 0 .. nT)
            foreach (a2; 0 .. 4)
                foreach (b2; a2 + 1 .. 4)
                {
                    Vertex[2] e = tets[ti][a2] < tets[ti][b2]
                        ? [tets[ti][a2], tets[ti][b2]]
                        : [tets[ti][b2], tets[ti][a2]];
                    bool dup = false;
                    foreach (i; 0 .. nSE)
                        if (seenE[i] == e) { dup = true; break; }
                    if (dup) continue;
                    if (nSE < seenE.length) seenE[nSE++] = e;
                    if (mfd.degreeOrZero!1(e[]) == 3) n3reg++;
                }
    }
    Vertex[128] nb;
    int nn = 0;
    foreach (ti; 0 .. nT)
        foreach (x; tets[ti])
        {
            if (x == chord[0] || x == chord[1]) continue;
            bool seen = false;
            foreach (i; 0 .. nn) if (nb[i] == x) { seen = true; break; }
            if (seen) continue;
            if (nn >= nb.length) return -1;
            nb[nn++] = x;
        }
    int n = 0;
    foreach (i; 0 .. nn)
        foreach (c; chord)
        {
            Vertex[2] e = c < nb[i] ? [c, nb[i]] : [nb[i], c];
            immutable d = mfd.degreeOrZero!1(e[]);
            if (d == 0) continue;
            if (n >= outD.length) return -1;
            outD[n++] = d;
        }
    return n;
}

/// Umbrella for the chord carrier: the replayed catalysed-path table.
/// wf0Key uses buckets 0..6 (bits 0..55); the top byte is free, so the
/// region's degree-3 count rides there and cannot collide with it.
private ulong wf0ChordKey(scope const(size_t)[] degs, int n3reg) @nogc nothrow
{
    return wf0Key(degs)
        + (cast(ulong)(n3reg > 255 ? 255 : n3reg) << 56);
}

private double wf0ChordU(const ref WormF0Params cfg,
                         scope const(size_t)[] degs, int n3reg) nothrow
{
    if (cfg.ctab.length == 0) return 0.0;
    immutable k = wf0ChordKey(degs, n3reg);
    if (auto p = k in cfg.ctab)
    {
        double u = *p;
        if (u > cfg.ucapHi) u = cfg.ucapHi;
        if (u < cfg.ucapLo) u = cfg.ucapLo;
        return u;
    }
    return cfg.cfb;
}

/// Collect the union of the two endpoint stars -- NOT just the tets
/// containing the whole chord. A move that touches only ONE endpoint
/// cannot alter the chord's degree (see wf0ChordEnumH), and those are
/// exactly the moves that do work while the flicker stays closable, so
/// they must be in the region. Uses collectStar (O(deg^2),
/// allocation-free) seeded from a tet on the chord, which dedups across
/// the two calls because it appends through a running count.
private int wf0ChordStar(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] chord, Vertex[4][] outT)
{
    Vertex[4] seed;
    bool got = false;
    foreach (f; mfd.star(chord[]))
    {
        int k = 0;
        foreach (x; f) { if (k < 4) seed[k] = x; k++; }
        if (k == 4) { got = true; }
        break;
    }
    if (!got) return -1;
    int n = collectStar(mfd, chord[0], seed, outT, 0);
    if (n < 0) return -1;
    n = collectStar(mfd, chord[1], seed, outT, n);
    if (n < 0) return -1;
    foreach (i; 0 .. n) outT[i][].sort();
    return n;
}

/// Enumerate moves whose support contains BOTH chord endpoints -- the
/// chord analogue of wf0EnumH, and by the same lemma: only such a move
/// can change the chord's own local configuration, so every other move
/// leaves this head invariant.
///
/// This is what lets a flicker CATALYSE anything. The global repair
/// kernel proposes uniformly over the whole manifold, so it lands in the
/// flicker's ~5-vertex support about once in 300 draws; targeting the
/// creation site fixed WHERE the flicker is born but not where the walk
/// looks. Proposals drawn from this set are, by construction, exactly
/// the moves the flicker can affect.
///
/// The 3->2 centered on the chord itself is EXCLUDED: that is the
/// episode's close, not a head move.
private int wf0ChordEnumH(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] chord, scope Vertex[4][] tets, int nT,
    WF0Cand!Vertex[] outC)
{
    int n = 0;
    Vertex[3][1024] seenF;
    int nSF = 0;
    Vertex[2][1024] seenE;
    int nSE = 0;
    // A move's support must MEET the chord -- one endpoint is enough,
    // and one endpoint is what we want. The chord's degree is the number
    // of tets containing BOTH endpoints, and every tet a move touches
    // has all its vertices in the support; so a move omitting an
    // endpoint cannot change the degree, and the flicker stays closable
    // while work happens beside it. Demanding BOTH endpoints (the
    // vertex-head analogue) guarantees the degree changes and so
    // guarantees the closure condition is destroyed -- measured 0/125
    // vs 2186/2286 moves leaving the chord closable.
    bool meetsChord(scope const(Vertex)[] x, scope const(Vertex)[] y)
    {
        foreach (q; x) if (q == chord[0] || q == chord[1]) return true;
        foreach (q; y) if (q == chord[0] || q == chord[1]) return true;
        return false;
    }
    foreach (ti; 0 .. nT)
    {
        auto T = tets[ti];
        foreach (skip; 0 .. 4)
        {
            Vertex[3] face;                  // T is sorted => face is too
            int k = 0;
            foreach (i; 0 .. 4) if (i != skip) face[k++] = T[i];
            bool dup = false;
            foreach (i; 0 .. nSF) if (seenF[i] == face) { dup = true; break; }
            if (dup) continue;
            if (nSF < seenF.length) seenF[nSF++] = face;
            int[2] ap = 0;
            if (mfd.writeFaceApexes(face[0], face[1], face[2], ap.ptr) != 2)
                continue;
            Vertex[2] axis = ap[0] < ap[1]
                ? [cast(Vertex) ap[0], cast(Vertex) ap[1]]
                : [cast(Vertex) ap[1], cast(Vertex) ap[0]];
            if (!meetsChord(face[], axis[])) continue;
            if (mfd.degreeOrZero!1(axis[]) != 0) continue;
            if (mfd.anyFrozen(face[])) continue;
            if (n < outC.length)
            { outC[n].is23 = true; outC[n].f = face; outC[n].p = axis; n++; }
        }
        foreach (a2; 0 .. 4)
            foreach (b2; a2 + 1 .. 4)
            {
                Vertex[2] e = T[a2] < T[b2] ? [T[a2], T[b2]]
                                            : [T[b2], T[a2]];
                if (e == chord) continue;        // that is the CLOSE
                bool dup = false;
                foreach (i; 0 .. nSE) if (seenE[i] == e) { dup = true; break; }
                if (dup) continue;
                if (nSE < seenE.length) seenE[nSE++] = e;
                if (mfd.degreeOrZero!1(e[]) != 3) continue;
                int[8] lkB = 0;
                if (mfd.writeEdgeLinkCycle(e[0], e[1], T[], lkB.ptr) != 3)
                    continue;
                Vertex[3] link = [cast(Vertex) lkB[0], cast(Vertex) lkB[1],
                                  cast(Vertex) lkB[2]];
                link[].sort();
                if (!meetsChord(link[], e[])) continue;
                if (n < outC.length)
                { outC[n].is23 = false; outC[n].f = link; outC[n].p = e; n++; }
            }
    }
    return n;
}

/// Isolation wrapper for the site enumerator (see ddg_sampler_chain_sites).
long wf0ChainSitesProbe(Vertex)(ref Manifold!(3, Vertex) mfd, int kmax)
{
    return wf0ChainSites!Vertex(mfd, kmax, -1, null, null);
}

/// AGGREGATION-WEIGHTED creation-site enumeration for the chord channel.
///
/// Same candidate set as wf0ChainSites (valid 2->3 sites within `kmax`
/// chain windows of some degree-3 chord), but each site i carries
/// w_i = exp(beta * n_i) with n_i the number of DISTINCT flickers whose
/// support shares a vertex with site i's support -- exactly the
/// condition under which bilocal factorization fails (measured 1e-13
/// when disjoint), hence the only geometry in which the two interact.
/// The adopted chord (exA, exB) is excluded from n_i: it is about to be
/// annihilated, so counting it would reward the flicker for staying put.
///
/// Returns the total weight W. With `pick` in [0, W) the selected site
/// is written to outFace/outAx. With (tgA, tgB) >= 0 the weight of the
/// site whose AXIS equals that pair is written to *outTgW (0 if absent)
/// -- that is what the close needs for the reverse proposal density.
///
/// beta == 0 short-circuits to the uniform path, so the default is
/// bit-for-bit the certified behaviour and W == the plain site count.
/// One enumerated creation site, for the COLLECT mode below.
struct WF0Site(Vertex)
{
    Vertex[3] face;
    Vertex[2] axis;
    double wt = 0.0;
}

/// ditto (COLLECT mode). Pass `outAll`/`outNAll` to have every accepted
/// site appended in enumeration order along with its weight. The caller
/// can then make the weighted draw itself over the cached array instead
/// of paying a SECOND full enumeration for it -- phase 1 is O(N), and
/// the open used to run it twice purely because the draw needs the total
/// weight before it can be scaled. Selecting the first i with
/// pick < sum(wt[0..i+1]) over the cache is exactly the in-loop rule, so
/// this is bit-identical (same distribution, same single uniform01).
/// Collection stops silently if `outAll` fills; `*outNAll` reports how
/// many were stored and the caller must fall back if it is short of the
/// enumerated count.
private double wf0ChainSitesW(Vertex)(ref Manifold!(3, Vertex) mfd,
    int kmax, double beta, Vertex exA, Vertex exB,
    double pick, Vertex[3]* outFace, Vertex[2]* outAx,
    Vertex tgA, Vertex tgB, double* outTgW,
    WF0Site!Vertex[] outAll = null, int* outNAll = null)
{
    import std.math : exp;
    // ---- flicker supports, as a per-vertex bitmask (up to 64 tracked;
    // beyond that the weight is a well-defined function of the first 64
    // in enumeration order, which is deterministic, so balance holds).
    static ulong[8192] fmask;
    static Vertex[512] touched;
    int nTouch = 0;
    int nF = 0;
    if (beta != 0.0)
    {
        foreach (ed; mfd.simplices(1))
        {
            Vertex[2] e;
            int i = 0;
            foreach (x; ed) { if (i < 2) e[i] = x; i++; }
            if (i != 2) continue;
            if (mfd.degreeOrZero!1(e[]) != 3) continue;
            if ((e[0] == exA && e[1] == exB) || (e[0] == exB && e[1] == exA))
                continue;                       // the adopted chord
            if (nF >= 64) break;
            Vertex[8] lk;
            immutable nl = wf0EdgeLink(mfd, e[0], e[1], lk[]);
            if (nl != 3) continue;
            immutable ulong bit = 1UL << nF;
            Vertex[5] sup = [e[0], e[1], lk[0], lk[1], lk[2]];
            foreach (x; sup)
            {
                if (x < 0 || x >= fmask.length) continue;
                if (fmask[x] == 0 && nTouch < touched.length)
                    touched[nTouch++] = x;
                fmask[x] |= bit;
            }
            nF++;
        }
    }

    double W = 0.0;
    bool tookPick = false;
    if (outTgW !is null) *outTgW = 0.0;
    if (outNAll !is null) *outNAll = 0;

    static Vertex[4][4096] starts;
    int nStart = 0;
    foreach (f; mfd.simplices(3))
    {
        Vertex[4] w;
        {
            int k = 0;
            foreach (x; f) { if (k < 4) w[k] = x; k++; }
            if (k != 4) continue;
        }
        w[].sort();
        bool near = false;
        foreach (i2; 0 .. 4)
            foreach (j2; i2 + 1 .. 4)
            {
                Vertex[2] e = [w[i2], w[j2]];
                if (mfd.degreeOrZero!1(e[]) == 3) near = true;
            }
        if (!near) continue;
        if (nStart >= starts.length) break;
        starts[nStart++] = w;
    }
    foreach (si; 0 .. nStart)
    {
        Vertex[4] w = starts[si];
        foreach (step; 0 .. kmax)
        {
            int[2] ap = 0;
            if (mfd.writeFaceApexes(w[1], w[2], w[3], ap.ptr) != 2) break;
            Vertex[2] axis = ap[0] < ap[1]
                ? [cast(Vertex) ap[0], cast(Vertex) ap[1]]
                : [cast(Vertex) ap[1], cast(Vertex) ap[0]];
            Vertex[3] face = [w[1], w[2], w[3]];
            face[].sort();
            if (mfd.degreeOrZero!1(axis[]) == 0 && !mfd.anyFrozen(face[]))
            {
                double wt = 1.0;
                if (beta != 0.0)
                {
                    ulong acc = 0;
                    foreach (x; face)
                        if (x >= 0 && x < fmask.length) acc |= fmask[x];
                    foreach (x; axis)
                        if (x >= 0 && x < fmask.length) acc |= fmask[x];
                    int n = 0;
                    for (ulong t = acc; t; t &= t - 1) n++;
                    wt = exp(beta * cast(double) n);
                }
                if (tgA >= 0 && axis[0] == tgA && axis[1] == tgB
                    && outTgW !is null)
                    *outTgW = wt;
                if (!tookPick && pick >= 0 && pick < W + wt
                    && outFace !is null)
                {
                    *outFace = face; *outAx = axis; tookPick = true;
                    if (outTgW !is null && tgA < 0) *outTgW = wt;
                }
                if (outAll !is null && outNAll !is null
                    && *outNAll < outAll.length)
                {
                    outAll[*outNAll].face = face;
                    outAll[*outNAll].axis = axis;
                    outAll[*outNAll].wt = wt;
                    ++*outNAll;
                }
                W += wt;
            }
            Vertex nxt = (ap[0] == w[0]) ? cast(Vertex) ap[1]
                                         : cast(Vertex) ap[0];
            if (ap[0] != w[0] && ap[1] != w[0]) break;
            w = [w[1], w[2], w[3], nxt];
        }
    }
    foreach (i; 0 .. nTouch) fmask[touched[i]] = 0;
    return W;
}

/// One CHORD bilocal episode (see the note above). Returns true if the
/// closed state changed. currentObjective stays exact; a capped walk is
/// unwound exactly.
bool wormChordPairEpisode(Vertex, P)(ref Manifold!(3, Vertex) mfd,
    ref real currentObjective, P params, const ref WormF0Params cfg,
    VertexPotState!Vertex* potState, VertexPot* pot,
    scope WF0Applied!Vertex[] undoBuf, out WormF0Result res)
{
    import std.math : log, exp;
    import std.random : uniform, uniform01;
    alias BM = BistellarMove!(3, Vertex);

    real baseRun = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;
    int nApplied = 0;

    real applyMove(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        auto bm = BM(cen, coc);
        immutable real dBase =
            mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        return dBase + dPot;
    }

    void record(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        assert(nApplied < undoBuf.length, "chord undo buffer overflow");
        undoBuf[nApplied].cl = cast(int) cen.length;
        undoBuf[nApplied].ccl = cast(int) coc.length;
        undoBuf[nApplied].cen[0 .. cen.length] = cen[];
        undoBuf[nApplied].coc[0 .. coc.length] = coc[];
        nApplied++;
    }

    void unwindAll()
    {
        foreach_reverse (k; 0 .. nApplied)
            applyMove(undoBuf[k].coc[0 .. undoBuf[k].ccl],
                      undoBuf[k].cen[0 .. undoBuf[k].cl]);
        nApplied = 0;
    }

    // mixture: cfg.ph funds the two chord head kernels (split evenly),
    // cfg.pg the global repair, the remainder the close
    immutable double pcl = 1.0 - cfg.ph - cfg.pg;
    if (!(pcl > 0.0 && pcl < 1.0)) return false;
    immutable double pch0 = 0.5 * cfg.ph;
    immutable double pch1 = cfg.ph;
    immutable double pgTop = cfg.ph + cfg.pg;

    // per-head caches: tets containing the chord, and the candidate set
    static Vertex[4][256][2] hTets;
    static WF0Cand!Vertex[1024][2] hCand;
    Vertex[2][2] chords;
    int[2] hNT = 0;
    int[2] hNH = 0;
    double[2] hU = 0.0;
    static size_t[256][2] hDegs;

    int refreshChord(int i)
    {
        hNT[i] = wf0ChordStar(mfd, chords[i], hTets[i][]);
        if (hNT[i] < 0) { hNH[i] = 0; return -1; }
        hNH[i] = wf0ChordEnumH(mfd, chords[i], hTets[i][], hNT[i],
                               hCand[i][]);
        int n3reg;
        immutable nd = wf0ChordDegs(mfd, chords[i], hTets[i][], hNT[i],
                                    hDegs[i][], n3reg);
        hU[i] = (nd < 0) ? cfg.cfb
                         : wf0ChordU(cfg, hDegs[i][0 .. nd], n3reg);
        return hNH[i];
    }

    /// a head move may not reach into the other head's tets
    bool chordClear(int self, scope const(Vertex)[] a,
                    scope const(Vertex)[] b)
    {
        immutable o = 1 - self;
        foreach (ti; 0 .. hNT[o])
            foreach (x; hTets[o][ti])
            {
                foreach (q; a) if (q == x) return false;
                foreach (q; b) if (q == x) return false;
            }
        return true;
    }

    // -- OPEN --------------------------------------------------------------
    // ball 0: CREATE a chord (2->3) at a site drawn uniformly from the
    // chain-targeted candidate set -- valid 2->3 sites within chainK
    // windows of some defect, so the flicker is born where it has
    // something to act on (see wf0ChainSites).
    immutable long nSite = wf0ChainSites!Vertex(mfd, cfg.chainK, -1,
                                                null, null);
    if (nSite <= 0) return false;
    Vertex[3] face = -1;
    Vertex[2] axis = -1;
    wf0ChainSites!Vertex(mfd, cfg.chainK, uniform(0, nSite),
                         &face, &axis);
    // the pick MUST have been written: the count pass and the pick pass
    // enumerate the same state in the same order. Guard anyway -- a
    // silent miss would hand BM() a degenerate simplex.
    if (face[0] < 0 || axis[0] < 0) return false;

    // ball 1: ADOPT a degree-3 chord, uniform over the n3 of them
    immutable long n3 = wf0Deg3Scan(mfd, -1, null, null);
    if (n3 <= 0) return false;
    Vertex bA, bB;
    wf0Deg3Scan(mfd, uniform(0, n3), &bA, &bB);
    Vertex[2] chordB = bA < bB ? [bA, bB] : [bB, bA];
    Vertex[16] lkB;
    immutable int nlkB = wf0EdgeLink(mfd, chordB[0], chordB[1], lkB[]);
    if (nlkB != 3) return false;

    // supports must be disjoint (the factorization condition)
    {
        Vertex[5] supA = [face[0], face[1], face[2], axis[0], axis[1]];
        Vertex[5] supB = [chordB[0], chordB[1], lkB[0], lkB[1], lkB[2]];
        foreach (x; supA) foreach (y; supB) if (x == y) return false;
    }

    auto bm23 = BM(face[], axis[]);
    if (!mfd.hasValidMove(bm23)) return false;
    immutable real dB23 = mfd.speculativeBistellarDelta(bm23, baseRun,
                                                        params);
    real dP23 = 0.0L;
    if (potState !is null)
        dP23 = mfd.potentialBistellarDelta(bm23, *potState, *pot, false);
    // U of the pair AFTER the open; the created chord's own value is
    // only known once it exists, so the ratio is assembled below.
    immutable double laOpenBase = cfg.zeta2 - cast(double)(dB23 + dP23)
        + log(cast(double) nSite) + log(cast(double) n3) + log(pcl);
    // The created chord's umbrella only exists once the chord does, so
    // apply first, measure both, then decide and roll back on reject.
    applyMove(face[], axis[]);
    Vertex[2] chordA = axis;
    chords[0] = chordA;
    chords[1] = chordB;
    if (refreshChord(0) < 0 || refreshChord(1) < 0)
    { applyMove(axis[], face[]); return false; }
    immutable double laOpen = laOpenBase + hU[0] + hU[1];
    if (!(laOpen >= 0 || uniform01 <= exp(laOpen)))
    { applyMove(axis[], face[]); return false; }
    record(face[], axis[]);
    res.opened = 3;                      // chord-pair open
    res.head = cast(int) chordB[0];

    // -- walk: plain Metropolis (U == 0, so dU == 0 for every move) -------
    while (true)
    {
        if (res.steps >= cfg.maxstep)
        {
            unwindAll();
            res.closedHow = 3;
            res.dS = cast(double) deltaTotal;
            currentObjective += deltaTotal;
            return false;
        }
        res.steps++;
        immutable double rr = uniform01;
        if (rr < pch1)
        {
            // CHORD HEAD kernel: propose only moves the flicker can
            // affect. U == 0 here, so the ratio is pure Metropolis with
            // the per-class count correction |H|/|H'|.
            immutable int i = (rr < pch0) ? 0 : 1;
            res.nH++;
            if (hNH[i] == 0) continue;
            auto c = hCand[i][uniform(0, hNH[i])];
            if (!chordClear(i, c.f[], c.p[])) continue;
            auto bmh = c.is23 ? BM(c.f[], c.p[]) : BM(c.p[], c.f[]);
            if (!mfd.hasValidMove(bmh)) continue;
            immutable int nH0 = hNH[i];
            immutable double uPre = hU[0] + hU[1];
            immutable real dh = applyMove(c.is23 ? c.f[] : c.p[],
                                          c.is23 ? c.p[] : c.f[]);
            if (refreshChord(0) < 0 || refreshChord(1) < 0)
            {
                applyMove(bmh.coCenter, bmh.center);
                refreshChord(0); refreshChord(1);
                continue;
            }
            immutable double lah = -cast(double) dh
                + (hU[0] + hU[1] - uPre)
                + log(cast(double) nH0)
                - log(cast(double)(hNH[i] == 0 ? 1 : hNH[i]));
            if (hNH[i] != 0 && (lah >= 0 || uniform01 <= exp(lah)))
            {
                res.accH++;
                record(c.is23 ? c.f[] : c.p[], c.is23 ? c.p[] : c.f[]);
            }
            else
            {
                applyMove(bmh.coCenter, bmh.center);
                refreshChord(0); refreshChord(1);
            }
            continue;
        }
        if (rr < pgTop)
        {
            res.nG++;
            Vertex fresh = cast(Vertex) mfd.fVector[0];
            auto bm = mfd.chooseRandomMove(fresh, params);
            if (bm.center.length == 1 || bm.coCenter.length == 1) continue;
            // never annihilate a head chord outside the close
            // with a chord umbrella, U changes iff an endpoint is in
            // the support, so those moves belong to the head kernel and
            // the dU == 0 lemma covers everything left here
            {
                bool touches = false;
                foreach (x; bm.center)
                    foreach (cc; chords)
                        if (x == cc[0] || x == cc[1]) touches = true;
                foreach (x; bm.coCenter)
                    foreach (cc; chords)
                        if (x == cc[0] || x == cc[1]) touches = true;
                if (touches) continue;
            }
            immutable real dB = mfd.speculativeBistellarDelta(bm, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm, *potState, *pot,
                                                 false);
            immutable double la = -cast(double)(dB + dP);
            if (la >= 0 || uniform01 <= exp(la))
            {
                res.accG++;
                applyMove(bm.center, bm.coCenter);
                record(bm.center, bm.coCenter);
                refreshChord(0);
                refreshChord(1);
            }
            continue;
        }

        // -- CLOSE: annihilate ball 1's chord (3->2), drop ball 0's -------
        // both heads must be back at the state a 2->3 creates: degree 3
        if (mfd.degreeOrZero!1(chordA[]) != 3) continue;
        if (mfd.degreeOrZero!1(chordB[]) != 3) continue;
        Vertex[16] lk;
        immutable int nlk = wf0EdgeLink(mfd, chordB[0], chordB[1], lk[]);
        if (nlk != 3) continue;
        Vertex[3] tri = [lk[0], lk[1], lk[2]];
        tri[].sort();
        auto bm32 = BM(chordB[], tri[]);
        if (!mfd.hasValidMove(bm32)) continue;
        immutable real dB = mfd.speculativeBistellarDelta(bm32, baseRun,
                                                          params);
        real dP = 0.0L;
        if (potState !is null)
            dP = mfd.potentialBistellarDelta(bm32, *potState, *pot, false);
        // The reverse open draws ball 1's site from the POST-close
        // candidate set and adopts ball 0's chord from the post-close
        // degree-3 set, so both counts are taken after the move. The
        // site set has no cheap local delta (a 3->2 can move defects in
        // and out of chain range), so it is recomputed on the applied
        // state and the move is rolled back if the draw rejects.
        immutable long n3after = wf0Deg3Scan(mfd, -1, null, null)
            + wf0Deg3Delta(mfd, chordB, tri);
        if (n3after <= 0) continue;
        // measure the post-close candidate count on the applied state,
        // then restore so the decision is taken in a clean state and the
        // move is re-applied only on acceptance
        applyMove(chordB[], tri[]);
        immutable long nSiteAfter = wf0ChainSites!Vertex(mfd, cfg.chainK,
                                                         -1, null, null);
        applyMove(tri[], chordB[]);
        if (nSiteAfter <= 0) continue;
        immutable double la = -cast(double)(dB + dP) - cfg.zeta2
            - hU[0] - hU[1]
            - log(cast(double) nSiteAfter) - log(cast(double) n3after)
            - log(pcl);
        if (la >= 0 || uniform01 <= exp(la))
        {
            applyMove(chordB[], tri[]);
            res.closedHow = 5;           // chord pair closed
            res.dS = cast(double) deltaTotal;
            currentObjective += deltaTotal;
            return true;
        }
    }
    assert(0);
}

/// Change in the degree-3 edge count caused by the 3->2 on (chord, tri).
/// The three tets round the chord become two, so: the chord itself
/// leaves the set; each LINK-triangle edge gains a tet (it was in one of
/// the three, it is in both survivors), so one at degree 3 leaves; each
/// of the six SPOKE edges (chord endpoint, link vertex) loses a tet, so
/// one at degree 4 arrives at 3. Edge degrees are >= 3 in a closed
/// manifold, so nothing can enter the set from below.
private long wf0Deg3Delta(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] chord, Vertex[3] tri)
{
    long d = -1;                          // the chord goes away
    foreach (i; 0 .. 3)
        foreach (j; i + 1 .. 3)
        {
            Vertex[2] e = tri[i] < tri[j] ? [tri[i], tri[j]]
                                          : [tri[j], tri[i]];
            if (mfd.degreeOrZero!1(e[]) == 3) d--;   // 3 -> 4, leaves
        }
    foreach (c; chord)
        foreach (t; tri)
        {
            Vertex[2] e = c < t ? [c, t] : [t, c];
            if (mfd.degreeOrZero!1(e[]) == 4) d++;   // 4 -> 3, joins
        }
    return d;
}

/*
STRICT-CLOSURE chord channel.
=============================

The transport channel above has a WEAK closure condition -- degree 3,
which a 2->3 satisfies the instant the flicker is born -- so it commits
relocations and almost nothing else. Trying to fix that with the
umbrella cannot work: the open-sector weight is zeta e^{-S+U}, so a high
U both DRAWS the walk to a state and, through alpha_close ~ e^{-dS-U},
BLOCKS closing from it. Same parameter, same sign, both effects.
Measured: large positive tube entries sent abandonment to 101/200, and
clamping them merely traded that back for lost opens. There is no U that
makes a state unattractive AND hard to close from.

So the discrimination moves out of U and into the CLOSURE CRITERION,
where it belongs, expressed as an absolute local configuration (a state
function) rather than a comparison with the start (history, which would
break balance).

Both marks are pure flags -- no move at open or close -- so the whole
f-change comes from the walk, and the sector boundary costs nothing.
The open-sector state is the UNORDERED pair of marked vertex-pairs, so
the roles need not be tracked and the reverse is automatic:

  open   pick e_full  (an existing degree-3 chord, uniform over n3)
         pick e_empty (an ABSENT pair from the chain-targeted sites)
         gate: e_empty's region carries no degree-3 edge
  close  fires when one mark is degree 3 and the other is ABSENT, and
         the absent one's region carries no degree-3 edge

At open the close condition already holds, so a null closure is
available; that is harmless (it commits nothing). A PRODUCTIVE episode
swaps them -- the old chord is gone, a new one stands at the target, and
the gate says the old neighbourhood came out clean, i.e. the reaction
happened rather than a bare relocation. The same gate sits on both legs,
which is what keeps the close and its reverse open in balance.
*/

/// Does the region around a vertex pair carry no degree-3 edge? The
/// region is the union of the two vertices' stars, so this is defined
/// whether or not the pair is currently an edge.
private bool wf0RegionClean(Vertex)(ref Manifold!(3, Vertex) mfd,
    Vertex[2] mark, Vertex[2] ignore, int maxAllowed = 0)
{
    // count DISTINCT degree-3 edges in the union of the mark's two
    // vertex stars, excluding `ignore`, and pass iff at most maxAllowed.
    // maxAllowed = 0 is the original all-or-nothing test.
    Vertex[2][64] found;
    int nf = 0;
    foreach (v; mark)
    {
        Vertex[4] seed;
        if (!mfd.someFacetContaining(v, seed)) return false;
        static Vertex[4][256] st;
        immutable n = collectStar(mfd, v, seed, st[], 0);
        if (n < 0) return false;
        foreach (i; 0 .. n)
            foreach (a2; 0 .. 4)
                foreach (b2; a2 + 1 .. 4)
                {
                    Vertex[2] e = st[i][a2] < st[i][b2]
                        ? [st[i][a2], st[i][b2]] : [st[i][b2], st[i][a2]];
                    if (e == ignore) continue;   // the other mark
                    if (mfd.degreeOrZero!1(e[]) != 3) continue;
                    bool seen = false;
                    foreach (j; 0 .. nf)
                        if (found[j] == e) { seen = true; break; }
                    if (seen) continue;
                    if (nf >= found.length) return false;   // way over
                    found[nf++] = e;
                    if (nf > maxAllowed) return false;
                }
    }
    return true;
}

/// One STRICT-CLOSURE chord episode (see the note above). Returns true
/// if the closed state changed.
bool wormChordStrictEpisode(Vertex, P)(ref Manifold!(3, Vertex) mfd,
    ref real currentObjective, P params, const ref WormF0Params cfg,
    VertexPotState!Vertex* potState, VertexPot* pot,
    scope WF0Applied!Vertex[] undoBuf, out WormF0Result res)
{
    import std.math : log, exp;
    import std.random : uniform, uniform01;
    alias BM = BistellarMove!(3, Vertex);

    real baseRun = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;
    int nApplied = 0;

    real applyMove(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        auto bm = BM(cen, coc);
        immutable real dBase =
            mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        return dBase + dPot;
    }

    void record(scope const(Vertex)[] cen, scope const(Vertex)[] coc)
    {
        assert(nApplied < undoBuf.length, "strict undo buffer overflow");
        undoBuf[nApplied].cl = cast(int) cen.length;
        undoBuf[nApplied].ccl = cast(int) coc.length;
        undoBuf[nApplied].cen[0 .. cen.length] = cen[];
        undoBuf[nApplied].coc[0 .. coc.length] = coc[];
        nApplied++;
    }

    void unwindAll()
    {
        foreach_reverse (k; 0 .. nApplied)
            applyMove(undoBuf[k].coc[0 .. undoBuf[k].ccl],
                      undoBuf[k].cen[0 .. undoBuf[k].cl]);
        nApplied = 0;
    }

    // CATALYSIS AUDIT. Every walk move here is priced alone (head is
    // -dS + Hastings counts, global is plain -dS, and U == 0), so no
    // move can spend action released by another. The observable
    // consequence is that accepted uphill moves must follow the bare
    // Metropolis tail e^-X -- so the largest accepted dS over N
    // accepted uphill moves should sit near log N, NOT beyond it.
    // res.umax = largest accepted single-move dS; res.nZ4 = number of
    // accepted UPHILL moves; res.zmin = 100 * max running excursion
    // max_k (S_k - S_0), the barrier the episode actually crossed.
    // The GLOBAL kernel is the built-in control: it is plain thermal
    // Metropolis over the whole manifold and is explicitly forbidden
    // from touching either mark, while the HEAD kernel works right at
    // the flicker. If catalysis were real the head arm would accept
    // uphill moves the global arm cannot. Split the statistics.
    res.umax = 0.0;
    res.nZ4 = 0;
    res.zmin = 0;
    res.dsArm[] = 0.0;
    res.nUpArm[] = 0;
    void noteAccept(double d, bool isHead)
    {
        if (d > res.umax) res.umax = d;
        if (d > 0.0) res.nZ4++;
        if (isHead)
        {
            if (d > res.dsArm[0]) res.dsArm[0] = d;
            if (d > 0.0) res.nUpArm[0]++;
        }
        else
        {
            if (d > res.dsArm[1]) res.dsArm[1] = d;
            if (d > 0.0) res.nUpArm[1]++;
        }
        immutable double exc = cast(double) deltaTotal;
        if (exc * 100.0 > cast(double) res.zmin)
            res.zmin = cast(int)(exc * 100.0);
    }

    double pcl = 1.0 - cfg.ph - cfg.pg;
    if (!(pcl > 0.0 && pcl < 1.0)) return false;
    double phE = cfg.ph, pgE = cfg.pg;
    if (cfg.zeta2Auto && cfg.maxstep > 0 && cfg.pclTarget > 0)
    {
        // AUTO p_close. The close is a geometric trial, so the mean
        // episode length is 1/p_close; tying it to the step budget as
        // maxstep/pclTarget makes abandonment ~e^-pclTarget and stops
        // the walk from closing before it has done anything. p_close
        // is NOT a free knob once zeta2 is auto: it already sits inside
        // zeta2* = -(log n3 + log nSite + log p_close), so the two must
        // be derived together or they fight. Frozen for the episode,
        // like zeta2.
        double want = cfg.pclTarget / cast(double) cfg.maxstep;
        if (want > pcl) want = pcl;
        if (want < 1e-12) want = 1e-12;
        immutable double scale = (1.0 - want) / (cfg.ph + cfg.pg);
        phE = cfg.ph * scale;
        pgE = cfg.pg * scale;
        pcl = want;
    }

    // -- OPEN: two pure flags, no move --------------------------------
    immutable long n3 = wf0Deg3Scan(mfd, -1, null, null);
    if (n3 <= 0) return false;
    Vertex fa, fb;
    wf0Deg3Scan(mfd, uniform(0, n3), &fa, &fb);
    Vertex[2] mFull = fa < fb ? [fa, fb] : [fb, fa];

    // AGGREGATION-WEIGHTED destination draw. The adopted chord is
    // excluded from the neighbour count (it is leaving). At aggBeta = 0
    // the weights are all 1, W == the old nSite, and the draw is the
    // original uniform one.
    // ONE enumeration, cached. The draw needs the total weight before it
    // can be scaled, which is why this used to run the (O(N)) enumeration
    // twice; collecting the sites instead and drawing over the cache is
    // the same rule -- first i with pick < sum(wt[0 .. i+1]) -- so it is
    // bit-identical, and it removes one of the episode's three sweeps.
    static WF0Site!Vertex[16384] siteBuf;   // zero-init, never = void
    int nSite2 = 0;
    immutable double Wopen = wf0ChainSitesW!Vertex(mfd, cfg.chainK,
        cfg.aggBeta, mFull[0], mFull[1], -1.0, null, null, -1, -1, null,
        siteBuf[], &nSite2);
    if (!(Wopen > 0)) return false;
    Vertex[3] face = -1;
    Vertex[2] mEmpty = -1;
    double wPick = 0.0;
    {
        immutable double pk = uniform01 * Wopen;
        double acc = 0.0;
        foreach (i; 0 .. nSite2)
        {
            if (pk < acc + siteBuf[i].wt)
            {
                face = siteBuf[i].face;
                mEmpty = siteBuf[i].axis;
                wPick = siteBuf[i].wt;
                break;
            }
            acc += siteBuf[i].wt;
        }
    }
    if (mEmpty[0] < 0 || !(wPick > 0)) return false;
    // the two marks must be disjoint, and the target region clean
    foreach (x; mEmpty) foreach (y; mFull) if (x == y) return false;
    if (!wf0RegionClean(mfd, mEmpty, mFull, cfg.regionMax))
        return false;

    // AUTO zeta2. Both marks are pure flags, so the open carries no dS
    // term at all and the balanced fugacity is simply minus the log
    // proposal density -- no probe needed, unlike the vertex carrier
    // where dS14 has to be measured. Frozen at the OPEN and reused at
    // the close, so it is a constant of the episode (state-dependent
    // within an episode would break balance); recomputing it per
    // episode is exactly unbiased by the same argument that licenses
    // retubing. At zeta2 = 0 the open saturated at +9.1 and the close
    // died at e^-9.1: nothing ever closed.
    immutable double dens = log(cast(double) n3)
        + log(Wopen / wPick) + log(pcl);
    immutable double z2 = cfg.zeta2Auto ? -dens : cfg.zeta2;
    immutable double laOpen = z2 + dens;
    if (!(laOpen >= 0 || uniform01 <= exp(laOpen))) return false;
    res.opened = 4;                       // strict chord open
    res.head = cast(int) mFull[0];
    immutable long[4] fOpen = [cast(long) mfd.fVector[0],
        cast(long) mfd.fVector[1], cast(long) mfd.fVector[2],
        cast(long) mfd.fVector[3]];

    Vertex[2][2] marks = [mFull, mEmpty];

    // -- walk ----------------------------------------------------------
    static Vertex[4][256][2] mTets;
    static WF0Cand!Vertex[1024][2] mCand;
    int[2] mNT = 0;
    int[2] mNH = 0;

    int refreshMark(int i)
    {
        Vertex[4] seed;
        if (!mfd.someFacetContaining(marks[i][0], seed))
        { mNT[i] = 0; mNH[i] = 0; return -1; }
        int n = collectStar(mfd, marks[i][0], seed, mTets[i][], 0);
        if (n < 0) { mNH[i] = 0; return -1; }
        n = collectStar(mfd, marks[i][1], seed, mTets[i][], n);
        if (n < 0) { mNH[i] = 0; return -1; }
        foreach (k; 0 .. n) mTets[i][k][].sort();
        mNT[i] = n;
        mNH[i] = wf0ChordEnumH(mfd, marks[i], mTets[i][], n, mCand[i][]);
        return mNH[i];
    }
    if (refreshMark(0) < 0 || refreshMark(1) < 0) return false;

    immutable double pm0 = 0.5 * phE;
    immutable double pm1 = phE;
    immutable double pgTop = phE + pgE;

    while (true)
    {
        if (res.steps >= cfg.maxstep)
        {
            unwindAll();
            res.closedHow = 3;
            res.dS = cast(double) deltaTotal;
            currentObjective += deltaTotal;
            return false;
        }
        res.steps++;
        immutable double rr = uniform01;
        if (rr < pm1)
        {
            immutable int i = (rr < pm0) ? 0 : 1;
            res.nH++;
            if (mNH[i] == 0) continue;
            auto c = mCand[i][uniform(0, mNH[i])];
            auto bmh = c.is23 ? BM(c.f[], c.p[]) : BM(c.p[], c.f[]);
            if (!mfd.hasValidMove(bmh)) continue;
            immutable int nH0 = mNH[i];
            immutable real dh = applyMove(c.is23 ? c.f[] : c.p[],
                                          c.is23 ? c.p[] : c.f[]);
            if (refreshMark(0) < 0 || refreshMark(1) < 0)
            {
                applyMove(bmh.coCenter, bmh.center);
                refreshMark(0); refreshMark(1);
                continue;
            }
            immutable double lah = -cast(double) dh
                + log(cast(double) nH0)
                - log(cast(double)(mNH[i] == 0 ? 1 : mNH[i]));
            if (mNH[i] != 0 && (lah >= 0 || uniform01 <= exp(lah)))
            {
                res.accH++;
                noteAccept(cast(double) dh, true);
                record(c.is23 ? c.f[] : c.p[], c.is23 ? c.p[] : c.f[]);
            }
            else
            {
                applyMove(bmh.coCenter, bmh.center);
                refreshMark(0); refreshMark(1);
            }
            continue;
        }
        if (rr < pgTop)
        {
            res.nG++;
            Vertex fresh = cast(Vertex) mfd.fVector[0];
            auto bm = mfd.chooseRandomMove(fresh, params);
            if (bm.center.length == 1 || bm.coCenter.length == 1) continue;
            bool touches = false;
            foreach (x; bm.center)
                foreach (mk; marks) if (x == mk[0] || x == mk[1])
                    touches = true;
            foreach (x; bm.coCenter)
                foreach (mk; marks) if (x == mk[0] || x == mk[1])
                    touches = true;
            if (touches) continue;
            immutable real dB = mfd.speculativeBistellarDelta(bm, baseRun,
                                                              params);
            real dP = 0.0L;
            if (potState !is null)
                dP = mfd.potentialBistellarDelta(bm, *potState, *pot,
                                                 false);
            immutable double la = -cast(double)(dB + dP);
            if (la >= 0 || uniform01 <= exp(la))
            {
                res.accG++;
                noteAccept(cast(double)(dB + dP), false);
                applyMove(bm.center, bm.coCenter);
                record(bm.center, bm.coCenter);
                refreshMark(0); refreshMark(1);
            }
            continue;
        }

        // -- CLOSE: one mark full, the other empty, empty region clean --
        immutable long d0 = cast(long) mfd.degreeOrZero!1(marks[0][]);
        immutable long d1 = cast(long) mfd.degreeOrZero!1(marks[1][]);
        int full = -1, empty = -1;
        if (d0 == 3 && d1 == 0) { full = 0; empty = 1; }
        else if (d1 == 3 && d0 == 0) { full = 1; empty = 0; }
        else continue;
        if (!wf0RegionClean(mfd, marks[empty], marks[full],
                            cfg.regionMax)) continue;
        // no move at the close; the counts are the CURRENT ones
        immutable long n3c = wf0Deg3Scan(mfd, -1, null, null);
        // the REVERSE episode adopts marks[full] and must draw
        // marks[empty] as its destination, so its density needs that
        // site's weight in the post-move state. wRev == 0 means the
        // reverse cannot propose this move at all -- decline.
        double wRev = 0.0;
        immutable double Wc = wf0ChainSitesW!Vertex(mfd, cfg.chainK,
            cfg.aggBeta, marks[full][0], marks[full][1], -1.0, null, null,
            marks[empty][0], marks[empty][1], &wRev);
        if (n3c <= 0 || !(Wc > 0) || !(wRev > 0)) continue;
        immutable double la = -z2
            - log(cast(double) n3c) - log(Wc / wRev)
            - log(pcl);
        if (la >= 0 || uniform01 <= exp(la))
        {
            foreach (k; 0 .. 4)
                res.df[k] = cast(long) mfd.fVector[k] - fOpen[k];
            res.closedHow = 7;            // strict close
            res.dS = cast(double) deltaTotal;
            currentObjective += deltaTotal;
            return nApplied > 0;
        }
    }
    assert(0);
}

// ---------------------------------------------------------------------------
// Move selection
// ---------------------------------------------------------------------------

BistellarMove!(dim, Vertex) chooseRandomMove(int dim, Vertex, P)(
    ref Manifold!(dim, Vertex) manifold, Vertex newVertex, P parameters)
{
    alias BM = BistellarMove!(dim, Vertex);
    enum nVerts = dim + 1;
    enum maxMask = (1 << nVerts) - 1; // 2^(dim+1) - 1

    while(true)
    {
        auto facet = manifold.randomFacetOfDim(dim);

        // Pick a random non-empty subset via bitmask (avoids materializing all subsets)
        auto mask = uniform(1, maxMask + 1);
        Vertex[nVerts] centerBuf;
        int centerLen = 0;
        foreach (i; 0 .. nVerts)
        {
            if (mask & (1 << i))
                centerBuf[centerLen++] = facet[i];
        }
        auto center = centerBuf[0 .. centerLen];
        center.sort();

        auto centerDim = centerLen - 1;
        auto centerDeg = manifold.degree(center);

        if (centerDeg + centerDim != dim + 1)
            continue;

        BM bm;
        if (centerDim == dim)
        {
            bm = BM(center, newVertex.only);
        }
        else
        {
            auto coCenter = manifold.coCenter(center, facet);
            bm = BM(center, coCenter[]);
        }

        if (uniform01 > 2.0 / centerDeg)
            continue;

        if (!manifold.hasValidMove(bm))
            continue;

        return bm;
    }
}
///
@safe unittest
{
    // Smoke test: chooseRandomMove should return without hanging
    auto rp3 = Manifold!3([[1, 2, 3, 7], [1, 2, 3, 11], [1, 2, 6, 9], [1,
            2, 6, 11], [1, 2, 7, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 7,
            10], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 8, 9], [1, 4, 8, 10], [1,
            5, 6, 8], [1, 5, 6, 11], [1, 5, 8, 10], [1, 6, 8, 9], [2, 3, 4, 8],
            [2, 3, 4, 11], [2, 3, 7, 8], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4,
            8, 10], [2, 5, 7, 8], [2, 5, 7, 9], [2, 5, 8, 10], [2, 5, 9, 10],
            [2, 6, 9, 10], [3, 4, 5, 9], [3, 4, 5, 11], [3, 4, 8, 9], [3, 5, 9,
            10], [3, 6, 7, 8], [3, 6, 7, 10], [3, 6, 8, 9], [3, 6, 9, 10], [4,
            5, 6, 7], [4, 5, 6, 11], [4, 5, 7, 9], [4, 6, 7, 10], [5, 6, 7, 8]]);
}

// ---------------------------------------------------------------------------
// Unused vertex management
// ---------------------------------------------------------------------------

Vertex[] getUnusedVertices(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd, Vertex[] initialVertices)
{
    Vertex[] unusedVertices;
    // all gaps in list of vertices should be unused vertices
    if (initialVertices.front != 0)
    {
        unusedVertices ~= initialVertices.front.iota.array;
    }
    foreach (i; 0 .. initialVertices.length - 1)
    {
        if (initialVertices[i] + 1 != initialVertices[i + 1])
        {
            unusedVertices ~= iota(initialVertices[i] + 1, initialVertices[i + 1]).array;
        }
    }
    assert(unusedVertices.all!(v => !mfd.contains(v.only)));
    return unusedVertices;
}

/// Convenience overload: compute initial vertices from the manifold.
Vertex[] getUnusedVertices(int dim, Vertex)(const ref Manifold!(dim, Vertex) mfd)
{
    auto verts = mfd.simplices(0).joiner.array.dup.sort.array;
    if (verts.length == 0) return [];
    return getUnusedVertices(mfd, verts);
}

// ---------------------------------------------------------------------------
// FPKMC slide-graph scan (notes/FPKMC_DESIGN.md, M1/M2)
// ---------------------------------------------------------------------------

/// Decode a slide slot at chord (a,b): the sorted 3-vertex link, the
/// (c2,c3) pick, and the derived frame. Shared by trySlideMove-style
/// callers and the graph scan. Returns false if the chord is not a
/// degree-3 edge or the frame does not derive.
bool slideDecode(Vertex)(const ref Manifold!(3, Vertex) mfd,
    Vertex a, Vertex b, const(Vertex)[] hintTet, int slot,
    out SlideFrame!Vertex f, out Vertex[3] link)
{
    int[8] linkBuf = 0;
    auto nl = mfd.writeEdgeLinkCycle(a, b, hintTet, linkBuf.ptr);
    if (nl != 3) return false;
    link = [cast(Vertex) linkBuf[0], cast(Vertex) linkBuf[1],
            cast(Vertex) linkBuf[2]];
    link[].sort();
    immutable int orient = slot / 6;
    immutable int pick = slot % 6;
    immutable Vertex c0 = orient == 0 ? a : b;
    immutable Vertex c4 = orient == 0 ? b : a;
    static immutable int[2][6] picks =
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]];
    return deriveSlideFrame(mfd, c0, c4, link[picks[pick][0]],
                            link[picks[pick][1]], f);
}

/// Apply the four Pachner moves of the slide at (a,b)/slot and KEEP them,
/// exporting the per-move records for exact rollback with slideRollback.
/// Uses the same speculative-delta path as trySlideMove; potState (if any)
/// advances with the moves. Returns false (nothing applied) if the slot
/// does not form a legal slide. dS receives the total action change.
bool slideApplyKeep(Vertex, P)(
    ref Manifold!(3, Vertex) mfd, real currentObjective,
    Vertex a, Vertex b, const(Vertex)[] hintTet, int slot, P params,
    VertexPotState!Vertex* potState, const(VertexPot)* pot,
    out SlideRec!Vertex[4] recs, out real dS)
{
    alias BM = BistellarMove!(3, Vertex);
    SlideFrame!Vertex f;
    Vertex[3] link;
    if (!slideDecode(mfd, a, b, hintTet, slot, f, link)) return false;

    real baseRun = currentObjective
        - (potState !is null ? potState.total : 0.0L);
    real deltaTotal = 0.0L;
    int nApplied = 0;

    void rollback()
    {
        foreach_reverse (k; 0 .. nApplied)
        {
            auto inv = BM(recs[k].coCenter[0 .. recs[k].coCenterLen],
                          recs[k].center[0 .. recs[k].centerLen]);
            if (potState !is null)
                mfd.potentialBistellarDelta(inv, *potState, *pot, true);
            mfd.doMove(inv);
        }
    }

    bool step(scope const(Vertex)[] center, scope const(Vertex)[] coCenter)
    {
        auto bm = BM(center, coCenter);
        if (!mfd.hasValidMove(bm)) { rollback(); return false; }
        real dBase = mfd.speculativeBistellarDelta(bm, baseRun, params);
        real dPot = 0.0L;
        if (potState !is null)
            dPot = mfd.potentialBistellarDelta(bm, *potState, *pot, true);
        mfd.doMove(bm);
        baseRun += dBase;
        deltaTotal += dBase + dPot;
        recs[nApplied].centerLen = cast(int) center.length;
        recs[nApplied].coCenterLen = cast(int) coCenter.length;
        recs[nApplied].center[0 .. center.length] = center[];
        recs[nApplied].coCenter[0 .. coCenter.length] = coCenter[];
        nApplied++;
        return true;
    }

    Vertex[9] supBuf = 0;
    int nsup = 0;
    void addSup(Vertex v)
    {
        foreach (i; 0 .. nsup) if (supBuf[i] == v) return;
        supBuf[nsup++] = v;
    }
    addSup(f.c0); addSup(f.c2); addSup(f.c3); addSup(f.c4);
    addSup(f.c5); addSup(f.c6); addSup(f.c7); addSup(f.c8);
    foreach (v; link) addSup(v);
    if (mfd.anyFrozen(supBuf[0 .. nsup])) return false;

    Vertex[2] m1c = [f.c0, f.c4]; m1c[].sort();
    if (!step(m1c[], link[])) return false;
    Vertex[3] m2c = [f.c3, f.c4, f.c5]; m2c[].sort();
    Vertex[2] m2cc = [f.c2, f.c6]; m2cc[].sort();
    if (!step(m2c[], m2cc[])) return false;
    Vertex[3] m3c = [f.c5, f.c6, f.c7]; m3c[].sort();
    Vertex[2] m3cc = [f.c4, f.c8]; m3cc[].sort();
    if (!step(m3c[], m3cc[])) return false;
    Vertex[2] m4c = [f.c2, f.c6]; m4c[].sort();
    {
        if (mfd.degreeOrZero!1(m4c[]) != 3) { rollback(); return false; }
        Vertex[4] hint = 0;
        {
            int[2] ap = 0;
            if (mfd.writeFaceApexes(f.c2, f.c6, f.c3, ap.ptr) != 2)
            { rollback(); return false; }
            hint = [f.c2, f.c6, f.c3, cast(Vertex) ap[0]];
            hint[].sort();
        }
        int[8] lb4 = 0;
        auto n4 = mfd.writeEdgeLinkCycle(m4c[0], m4c[1], hint[], lb4.ptr);
        if (n4 != 3) { rollback(); return false; }
        Vertex[3] m4cc = [cast(Vertex) lb4[0], cast(Vertex) lb4[1],
                          cast(Vertex) lb4[2]];
        m4cc[].sort();
        if (!step(m4c[], m4cc[])) return false;
    }
    {
        Vertex[2] arrival = [f.c4, f.c8]; arrival[].sort();
        if (mfd.degreeOrZero!1(arrival[]) != 3) { rollback(); return false; }
    }
    dS = deltaTotal;
    return true;
}

/// Exact rollback of a kept slide (inverse moves in reverse order,
/// potState advanced in lockstep).
void slideRollback(Vertex)(ref Manifold!(3, Vertex) mfd,
    ref SlideRec!Vertex[4] recs,
    VertexPotState!Vertex* potState, const(VertexPot)* pot)
{
    alias BM = BistellarMove!(3, Vertex);
    foreach_reverse (k; 0 .. 4)
    {
        auto inv = BM(recs[k].coCenter[0 .. recs[k].coCenterLen],
                      recs[k].center[0 .. recs[k].centerLen]);
        assert(mfd.hasValidMove(inv), "slideRollback: inverse invalid");
        if (potState !is null)
            mfd.potentialBistellarDelta(inv, *potState, *pot, true);
        mfd.doMove(inv);
    }
}
