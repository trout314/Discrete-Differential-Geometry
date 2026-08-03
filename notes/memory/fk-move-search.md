---
name: fk-move-search
description: Machinery + key insight for the exhaustive FK->FK block-move search (constraint-surface moves)
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-22T02:45:48.896Z
---

Goal: find LOCAL moves FK-legal -> FK-legal (analog of Pachner moves on the
constraint surface). Motivated by: single Pachner moves provably CANNOT preserve
FK-legality (every one creates/needs an illegal edge -- see [[defect-kinetics-run5h]]),
so FK-legal states are isolated under the sampler's moves -> the glassiness/no-
transport we measured. A local FK->FK move algebra would give constraint-surface
dynamics (the discrete physical diffeomorphisms / the transport behind the
gravitational-wave question) and an ergodic sampler for the constrained liquid
needed for the TT/radiative test [[curvature-length-scale]].

MOVE = a pair of FK-legal fillings of the same FRAMED BOUNDARY. Framed boundary =
boundary triangulation + per-boundary-edge FRAMING (its degree WITHIN the piece).
Same framing => regluing preserves every boundary degree => FK-legal with NO
illegal intermediate (bypasses the single-move obstruction by acting at block
scale = your "cut a union of FK neighborhoods, glue a different one").

Core machinery: scripts/defect_dynamics/fk_moves.py (validated, 3 self-tests
pass). analyze() (interior/boundary, framing, FK check), canonical() (colored-
graph individualization-refinement), framed_boundary_canon (boundary TYPE),
ball_canon (abstract iso), filling_canon (REL FIXED boundary), find_moves.

KEY INSIGHT (from the octahedron self-test): move detection is an ENUMERATION
problem, NOT a bucketing problem. Fillings related by a boundary symmetry are
abstractly ISOMORPHIC balls yet are DISTINCT moves in a fixed context (the
exterior pins the boundary labels) -- e.g. the octahedron's 3 diagonal fillings.
So: (1) a single perfect crystal shows NO moves (every copy of a hole is filled
identically); (2) you find moves only by ENUMERATING alternative fillings of a
fixed framed boundary, or by mining a symmetry-broken source (the FK glass, a
phase boundary, defected states). The 2-3 Pachner move is NOT a framed-boundary
match -- it shifts the framing -- which is WHY FK moves must be composite/block.

NEXT PIECE (not built yet): the filling ENUMERATOR -- given a fixed framed
boundary, enumerate all FK-legal fillings with <= k interior vertices; >=2 => a
move. This is where the "most matches ruled out cheaply" pruning lives (partial
fillings die fast on the FK edge/vertex constraints). Open decision: enumerate
abstract framed boundaries (build from the 4 FK coordination polyhedra) vs mine
framed boundaries from the glass/multi-phase library then enumerate their fills.

Glass producer ready: scripts/defect_dynamics/glass_quench.py (melt fix: pin BOTH
f3 and mean edge degree via num_hinges_coef=2, hinge_degree_target_coef=0; V held
at 10176, melts to disclination fluid edv~7.4). BUT the quench ARRESTS: legaledge
plateaus ~0.52, legalvert stays 0.000 -- no FK glass forms (same glass wall the
VDV route hit; forming FK from a melt via local moves is the SAME hardness as the
frozen defects). So no disordered-FK mining source. Set aside.

ALPHABET EXTENSION (done): fk_moves.interior_ok(a, alphabet) -- 'fk' (strict,
species vac) or 'deg4' (FK + one fully-interior deg-4 edge = the minimal
disclination unit; species vac/deg4). glass_mine takes alphabet; moves classify
vac<->deg4=CREATION, deg4<->deg4=TRANSPORT. deg4 mine on run5h snapshot: 0 (the
isolated deg-4 edge is too rare/clustered statically).

EVENT REPLAY (built + VERIFIED): scripts/defect_dynamics/event_replay.py.
apply_event() reconstructs the triangulation exactly from the accepted-move
stream (event_changes per move type derived from move_geometry.replay_role_counts);
verified == sampler.facets() across 550k moves, all 5 types. Reusable.
scripts/defect_dynamics/deg4_moves.py: replay a live constrained R chain, catch
isolated-deg-4 birth/hop/death, capture minimal 2-vertex block before/after.
FINDING: isolated deg-4 edges DO occur transiently (up to 13 at once) -- refines
"size-2 never exists" (that was STABLE complexes). BUT can't capture a CLEAN
single-deg-4 block: dense regime -> defects clump, block contaminated; perfect
crystal -> near-frozen (60 moves/2500 sweeps), too few. Catch-22.

ALPHABET now generalized: fk_moves.ALPHABETS = {fk:{}, deg3:{3:1}, deg4:{4:1},
e34:{3:1,4:1}}. interior_ok(a, alphabet) caps illegal-edge multiset; deg-3 = the
2<->3 quantum, deg-4 = the 4<->4 quantum, e34 = both (the (3,4,4)-knot's
constituents) so block moves compose with the sampler's native flicker/rattle.

WHY REALIZED-CAPTURE IS BLOCKED IN R (3 independent confirmations): (1) dense
regime -> defects cluster, block contaminated; (2) dilute-from-perfect -> near
frozen (60 moves/2500 sweeps); (3) R has ZERO all-deg-6 triangles (in fact no
triangle has even 2 deg-6 edges -- the disclination lines are too spread), so a
2->3 always makes a deg-3+deg-4 CLUSTER, never a clean isolated deg-3. So the
minimal defects are never cleanly isolated by R's dynamics.

CONCLUSION / NEXT (firm): the minimal move is real but NOT cleanly realized by
the dynamics. The ONLY route to it is the filling ENUMERATOR (still deferred):
construct one FK-legal framed boundary, enumerate fillings with the e34 alphabet;
the vacuum filling and a one-defect filling of the SAME framed boundary ARE the
minimal creation move, guaranteed clean + source-independent. Build this next.
