# 4-4 Hinge Moves: Implementation and Mixing Analysis

Date: 2026-03-15

## What are 4-4 moves?

A 4-4 (hinge) move operates on a degree-4 edge in a 3-manifold triangulation.
The edge is shared by 4 tetrahedra whose link forms a quadrilateral cycle of 4
vertices. The move removes the edge and replaces it with the opposite diagonal
of the quadrilateral, producing a different set of 4 tetrahedra filling the
same region.

Key properties:
- **Not a Pachner move**: involves 6 vertices (not 5), so it lies outside the
  Pachner move framework.
- **Preserves the f-vector**: the number of vertices, edges, triangles, and
  tetrahedra are all unchanged.
- **Symmetric proposal**: forward and reverse proposal probabilities are equal
  (same number of containing facets, same edge count per facet, same f-vector),
  so no Hastings correction is needed beyond the objective delta.
- **O(1) construction**: the link cycle is walked via 3 ridge-link lookups.
- **Very high acceptance rate**: ~90-99% depending on size and energy function,
  because the moves are small local rearrangements.

## Implementation

Files modified:
- `source/manifold_moves.d` -- `HingeMove` struct (oldFacets, newFacets, inverse)
- `source/manifold.d` -- `hingeMove()` construction, `doHingeMove()`/`undoHingeMove()`
- `source/sampler.d` -- `speculativeHingeDelta()`, `tryProposeHingeMove()`, `mcmcStep()`
- `bench_sampler.d` -- comprehensive relaxation benchmark

Also fixed: `chooseRandomMove` was taking the manifold by value (copying all
hash maps on every call). Changed to `ref` -- this was causing a spurious 3x
"speedup" from hinge moves that was actually just avoiding the copy overhead.

## Mixing experiments

### Setup

Seeds equilibrated at sizes 500--10000 with physically motivated parameters:
- Target mean edge degree = 2pi/arccos(1/3) ~ 5.104 (flat Regge geometry)
- numFacetsCoef = 0.1, numHingesCoef = 0.1, coDim3DegreeVarianceCoef = 0.2
- 4-4 moves used during seed creation (30% attempt rate)

Six perturbation scenarios tested at each size:
1. Target mean edge degree +0.1 / -0.1
2. Target volume +5% / -5%
3. Vertex degree variance penalty 2x / 0.5x

Each scenario: 20 sweeps of relaxation, compared bistellar-only vs
bistellar + 30% hinge moves. Metric: lag-1 autocorrelation of the objective
trajectory (lower = faster decorrelation = better mixing).

### Results

**Wall-clock time**: Essentially identical between bistellar-only and
bistellar+hinge at all sizes (speedup ~ 1.0x). Hinge moves are cheap but
displace bistellar moves from the mix.

**Where hinge moves help (lower autocorrelation):**

| Size | Scenario           | Bistellar AC | +Hinge AC | Improvement |
|------|--------------------|-------------|-----------|-------------|
| 1002 | vtx_var_beta 2x    | 0.918       | 0.660     | strong      |
| 2004 | vtx_var_beta 2x    | 0.853       | 0.492     | strong      |
| 1002 | edge_deg -0.1      | 0.635       | 0.446     | moderate    |
| 10010| volume +5%         | 0.398       | 0.059     | strong      |

**Where hinge moves are neutral:**
- Most edge degree and volume perturbations at 5K--10K tets: both modes
  achieve low autocorrelation (<0.4) and similar final objectives.
- vtx_var_beta 2x at 5K: both stuck at AC ~ 0.83.

**Where hinge moves hurt (higher autocorrelation):**

| Size  | Scenario           | Bistellar AC | +Hinge AC | Effect      |
|-------|--------------------|-------------|-----------|-------------|
| 506   | vtx_var_beta 0.5x  | 0.320       | 0.491     | worse       |
| 2004  | vtx_var_beta 0.5x  | 0.312       | 0.498     | worse       |
| 5008  | vtx_var_beta 0.5x  | 0.143       | 0.454     | worse       |
| 10010 | edge_deg -0.1      | 0.028       | 0.300     | worse       |

### Interpretation

**Hinge moves assist descent toward regularity but resist ascent away from it.**

When the energy function demands *reducing* vertex degree variance (2x beta),
hinge moves help by locally homogenizing the degree distribution -- swapping
edges redistributes degree among the 6 involved vertices. The effect is
strongest at moderate sizes (1K-2K) where the per-move impact is significant
relative to the system size.

When the energy function *relaxes* the variance constraint (0.5x beta), hinge
moves are counterproductive. The system needs to explore states with higher
degree variance, which requires large-scale topological rearrangements
(adding/removing vertices via 1-4/4-1 moves). Hinge moves consume 30% of the
move budget with cheap local swaps that don't drive the system toward the new,
more disordered equilibrium.

**Volume perturbations** show an interesting size-dependent effect: at 10K tets,
hinge moves dramatically improve decorrelation for +5% volume growth (AC
0.40 -> 0.06) but are neutral for -5% shrinkage. This may be because growing
the manifold creates many new degree-4 edges (from 1-4 stellar subdivisions),
and hinge moves efficiently smooth out the resulting degree inhomogeneities.

**Edge degree perturbations** show mixed results. The target mean edge degree
is controlled by the global curvature penalty (numHingesCoef), which acts on
the total edge count rather than individual edge degrees. Hinge moves don't
change the edge count, so they can only affect this indirectly through degree
redistribution.

## Observations about 3-manifold triangulation structure

1. **Vertex degree variance grows with system size.** Equilibrated seeds show
   deg_var ~ 200 at 500 tets, ~450 at 1K, ~580 at 2K, ~975 at 5K, ~1320 at
   10K. This is roughly linear in system size, suggesting a constant
   per-vertex contribution to degree variance -- i.e., the degree distribution
   has a fixed shape that doesn't narrow with size.

2. **Mean edge degree converges to ~5.1.** All equilibrated seeds have mean
   edge degree between 5.11 and 5.18, close to the flat Regge target of
   2pi/arccos(1/3) ~ 5.104. The sampler successfully maintains this geometric
   constraint.

3. **Degree-4 edges are abundant.** Hinge move proposal acceptance rates are
   88-99%, meaning most random edges have degree 4 and most diagonals are
   valid. This implies degree-4 edges are a substantial fraction of all edges
   in equilibrated 3-manifold triangulations near the flat geometry.

4. **The energy landscape has different characteristic scales for different
   observables.** Volume (f-vector) relaxes quickly (high acceptance of 1-4 and
   4-1 moves), while vertex degree variance relaxes slowly (requires many
   coordinated local rearrangements). Edge degree is intermediate. This
   separation of timescales is a fundamental challenge for MCMC sampling of
   triangulations.

## Recommendations

- **Use hinge moves when the energy function penalizes degree variance.**
  A 30% attempt rate works well. The benefit is strongest at moderate sizes
  (1K-5K tets).
- **Don't use hinge moves when relaxing toward higher disorder.** They slow
  convergence by consuming move budget without driving the needed topological
  changes.
- **Adaptive mixing**: an ideal sampler would adjust the hinge move rate based
  on whether the current degree variance is above or below its target
  equilibrium value.
- **The `chooseRandomMove` by-value bug was a significant performance issue.**
  Any function taking a `Manifold` should use `ref` to avoid copying hash maps.

## Open questions

- How does the degree distribution shape depend on the target mean edge degree?
  Is it universal (up to scaling) or does it change qualitatively?
- What is the autocorrelation time for vertex degree variance as a function of
  system size? The data suggests it grows, but the functional form is unclear.
- Do higher hinge move rates (50%, 70%) help more for the vtx_var_beta 2x
  scenario, or is there a diminishing return?
- Are there other non-Pachner local moves (e.g., 2-3-2 composite moves) that
  could improve mixing for different observables?
