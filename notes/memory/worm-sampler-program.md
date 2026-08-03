---
name: worm-sampler-program
description: "Worm-sampler program (user-driven, 2026-07-24): move library worm_moves.py + stage-1 creation catalog done — (3,4,4) knot = unique cheapest excitation (+17.1 lam-units); next: orbit refinement, two-star stage, exhaustive non-crystal search"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T01:55:07.347Z
---

Goal: worm sampler for the legal manifold (create -> walk -> annihilate
excursions), with a move alphabet whose symmetric structure (inverse closure +
space-group orbits) supports a later momentum/Bloch analysis (dispersion of
the worm head from the orbit-resolved hop table).

Alphabet = V-preserving bistellar moves = complete worm grammar:
2-3 = pair creation (always makes a deg-3 edge; face type fixes species),
4-4 = translate deg-4 illegality (cannot fire in a legal region),
3-2 = annihilate deg-3. Inverse-closed. Hinge moves deliberately excluded
for now (ergodicity question deferred).

Code (commit 0b1f3a2): scripts/defect_dynamics/worm_moves.py = REUSABLE
crystal-independent library (bare facet arrays; user explicitly wants a
future EXHAUSTIVE non-crystal search to reuse it — keep all move arithmetic
there). worm_catalog.py = r-crystal driver, canonical S3xS2 signatures,
exact verification (apply + global recount) per class.

Move arithmetic (exact, verified): 2-3 on face abc apexes d,e: new edge de
deg 3, face edges -1, six side edges +1 (valid iff de absent). 3-2 inverse
(valid iff face abc absent). 4-4 on deg-4 ef, link cycle abcd, diagonal ac:
ef removed, ac created deg 4, equator +1, poles-to-b,d -1 (valid iff ac
absent; two diagonal choices).

STAGE-1 RESULTS (r reference, N3=24462): only face types (5,5,5)/(5,5,6)
exist — no two six-edges share a face in r. All 17982 (5,5,6) sites create
the (3,4,4) knot (11664 clean, rest +7s); (5,5,5) creates (3,4,4,4)+-7s at
+25..+51. KNOT = UNIQUE CHEAPEST EXCITATION: dS_shape +17.1..+18.3 (lam=1
units; x lam for acceptance — matches observed birth rates at lam=0.4).
Create+step totals exactly degenerate at +20.608 across 6 creation classes
= first path-independence evidence (groupoid structure). Some 4-4 steps
downhill (-8..-9 from dense creations). A 4-4 with the deg-3 head on its
equator converts it 3->4 (head species changes mid-walk — sampler must
track). Catalog: 47 signature classes (fine env splits ±0.6 = Wyckoff-orbit
structure, not yet registry-resolved).

STAGE 2 (2026-07-24 eve, commits 76ac1a0/e545899): targeted moves exposed in
C API (do/has_bistellar_move, do/has_hinge_move — FULL explicit validation,
asserts compiled out in release; crossval_moves.py = regression test, python
worm_moves == D core exact on 290 moves). State queries: illegal_edges (from
D incremental degree map), edge_link_cycle (O(deg) ridge walk, 11us, hint
tet + slow fallback). worm_cycles.py v2 = apply/undo DFS on ONE manifold
(bistellar undo = center/cocenter swap; 4-4 undo = diagonal's interleaved
cycle), Overlay net-diff via exact deltas, facet-diff hashing, --memo
pruning. ~2400x speedup (walks-2: 35min -> 0.9s); profiler built in.

STAGE-2 FINDINGS: minimal nontrivial legal->legal worm cycle needs >= 7
MOVES — at walks<=4 (503k sequences, 1M moves, 15.7k legal closures) EVERY
closure is the identity. The (3,4,4) worm is rigid: no local rearrangement
of the legal manifold in <= 6 moves. Path-independence structure rich:
44880 multi-path states at walks-4 (commuting distant 4-4s, complementary-
diagonal relations, create_A;walks = create_B conjugations). walks-5 run
launched (memo). WARNING for worm sampler: targeted moves bypass sampler
cocycle tracking (docstring warns).

RIGIDITY BOUNDS — **PARTLY RETRACTED 2026-08-02, see below.**
(2026-07-24 night, as claimed then): single-worm grammar exhausted
through 8 moves (walks-6, 11.2M seq memo) — ALL identity. Interleaved grammar
(mid-path 2-3/3-2, budget 6) exhausted through 6 moves (1.37M) — ALL identity;
7-move interleaved run in background.

>  **!! THOSE RUNS WERE NOT EXHAUSTIVE.** `worm_cycles.py` seeded one DFS per
>  `canon_sig` bucket. That signature is only an Aut invariant, so its buckets
>  are UNIONS of orbits: on R it gives 47 buckets where there are 102 face
>  orbits, and the available FOLLOW-UP MOVE SET differs between orbits inside
>  28 of the 34 fused buckets. So 55 starting orbits were never explored, and
>  their move trees are demonstrably not copies of the explored ones. "No
>  nontrivial legal->legal cycle within 8 moves" was therefore never
>  established. Fixed 2026-08-02 (seeds per exact orbit now); see
>  [[crystal-symmetry-group]].
>
>  **RE-VERIFIED so far, with all 102 seeds (2026-08-02, R m2):**
>  through **6 moves** (walks<=4, 1,080,498 sequences, 34162 legal endpoints)
>  — 0 nontrivial cycles. Note 1.08M vs the old 503k at the same depth, the
>  expected ~2.17x = 102/47.
>
>  **STILL OPEN:** 7 moves (walks-5, ~14M seq, ~30 min) and 8 moves (walks-6,
>  ~180M seq, ~6 h) — costs measured 2026-08-02, growth ~13x per walk level.
>  Until those run, quote the bound as SIX moves, not eight.
>
>  Everything downstream that rests on the 8-move bound inherits the caveat:
>  the "contractible worms are all trivial, only torus-wrapping worms act"
>  hypothesis, and the motivation for the whole BC-helix worm program. Matveev-Piergallini guarantees cycles
exist at fixed V, so this = LOCAL RIGIDITY of the r legal manifold, dual to
enumerate_fillings/glass_mine finding no small two-filling regions.
HYPOTHESIS (spin-ice pattern): contractible worms are ALL trivial; nontrivial
transformations require torus-WRAPPING worms (change web winding W).

BC-HELIX WORM (user idea, from trout314/quantum-random-walks; commit b839fff):
sliding-window BC chain (drop oldest vertex, exit opposite face, adopt apex —
face_apexes O(1) query) is deterministic+invertible => pure cycles, and they
WRAP: r-crystal orbits len 2286 winding (0,0,10) (~229 tets/wrap), len 2439
winding (-43,-29,-8). Chain-aligned worms are CHAIN-INTERNAL (face between
consecutive chain tets has both apexes on chain) => exact chain-relative
codes => motif = path segment between equal codes at different offsets.
MOTIF FOUND (worm_helix.py, depth-10 DFS, 120s): period 4 chain steps, 6
moves/period (3x 2-3 + 3x 3-2), steady 10-edge worm code. NEXT: propagate
motif around one wrap (~229 steps ~ 344 moves; replay via short per-period
DFS re-derivation, robust to env variation), close, measure net transform +
dW. If nontrivial: first legal cycle = sector-changing operator (W-program
milestone). See [[six-web-gauge]], [[fk-move-search]], [[five-illegal-knot]].

initObjective-NaN RESOLVED (2026-07-24): every ddg_sampler_set_* poisons
currentObjective to NaN (lazy invalidation) and only run() healed it;
ddg_sampler_current_objective leaked the sentinel — so any sampler built by
the python ctor with fixed-target coefs (ALL production runs) read NaN
before its first run(). Fix in ddg_capi.d: single recomputeObjective(state)
helper (dim dispatch + n6-pot total) replaces 3 duplicated recompute blocks;
current_objective now lazily recomputes. Validated: tracked==fresh to 1e-13
after 20k moves; D-SIDE dS CROSS-CHECK NOW CLOSED — 72/72 slide directions
on lam40x_snap15000, |dS_between - sampler objective diff| <= 2e-14 (params:
edge coef ESTAR/6 @ target ESTAR, n6 pot (0.6,1.0), lam=1 units). Snapshot
lives in OLD session scratchpad /private/tmp/claude-501/...1ed7f5d0.../
scratchpad/mgas/ — at risk of tmp cleanup; copy out if still needed.

SLIDE IS A NATIVE D MOVE TYPE (2026-07-24). Option B (user-approved): the
move class is CLEAN (species-preserving) slides ONLY. Gate verified first —
exhaustive inverse-closure on 4 states (crystal+1knot, lam40x_snap15000,
ab1_snap14000, b2e3_snap20000): 58/58 clean transitions inverse-closed,
k_f = k_r = 1 on EVERY one, 0 local-vs-global cleanliness disagreements. So
alpha = min(1, exp(-dS)) — plain Metropolis, NO Hastings correction.

source/sampler.d: SLIDE_SLOTS=12, SlideConfig{prob,tries,accepts},
SlideAccept{metropolis,trialOnly,force}, deriveSlideFrame, trySlideMove.
Two-pass: TRIAL applies the 4 Pachner moves for real (dS accumulated through
the SAME speculativeBistellarDelta path as every other move), tests
cleanliness on the O(1) changed-edge support (frame + chord link, <=9 verts
=> <=36 pairs), then rolls back exactly; COMMIT replays with full
instrumentation in lockstep (pot/cocycle/sixFlips/ledger), so an accepted
slide is indistinguishable from 4 ordinary accepted bistellar moves.
mcmcStep branch: deg-3 edge center + uniform01 < slide.prob (a deg-3 edge
sits in exactly 3 facets => the facet proposal is UNIFORM over chords, no
global list needed); invalid/unclean slot => continue (redraw), same as any
invalid proposal. Also checks arrival chord (c4,c8) has degree 3 ("landed").

capi/python: ddg_sampler_set_slide_prob / _slide_stats / _slide_at(a,b,slot,
mode,out_dS) -> s.set_slide_prob(p), s.slide_stats(), s.slide_at(a,b,slot,
commit=) [scripted/crossval only, bypasses MH]. ddg.SLIDE_SLOTS exported.

Regressions (worm_slide.py, the oracle is now ONLY a test reference):
  --dcross SNAP    slot-for-slot D vs oracle: verdict, dS, end facet set
  --dsampler SNAP  live sampler: maps validate, N3 preserved, cocycle
                   attached, tracked objective vs from-scratch
Results: lam40x 72 slots / 27 clean, b2e3 36 slots / 15 clean — 0 verdict
mismatches, 0 end-state mismatches, |dS_D - lam*dS_oracle| <= 1.2e-14.
Live: 200 sweeps @ slide_prob=1.0 -> 620 tries / 545 accepts (88%), cost
+1.4..3.4% wall clock (vs ~4 s/proposal for the old python oracle path),
objective drift 2.7e-12. mobile_gas.py: nslide=K REPLACED by slideprob=P.

CAVEAT: `meson test` still cannot LINK the unittest binary on macOS (ldc
rpath quirk, see [[build-state]]) — sampler.d unittests compile but do not
run, so the Python crossval above IS the regression for this move type.

DEFECT TRACKING REBUILT (2026-07-25, commits dce7332 + 716eb79).

`scripts/defect_dynamics/defect_state.py` is now THE defect census: one
incremental implementation over the accepted-move stream, 0.39 ms/sample vs
715 ms for the from-scratch rebuild every script was doing (1846x). Seeded
O(N) once, then 59 us/event. `audit()` compares the whole state against the
old path AND checks n6 == Z-12 on all-legal vertices (link is a triangulated
2-sphere: 2Z-4 triangles, degree sum 6Z-12).

BROADENED DEFECT = incident to an illegal edge (deg not 5,6) OR non-FK
coordination (n6 not in {0,2,3,4}). The second clause closes a real blind
spot: the FK potential penalises n6=5 but such a vertex lands in NO complex
and mobile_gas's legalvert scores it LEGAL. n6=1 (line terminating in bulk,
topologically forbidden) is never observed among all-legal vertices; n6>=5 is
rare (0/2/7 vs 46-63 illegal-edge vertices) and ALWAYS attaches to an existing
complex -- component count unchanged 8->8, no percolation. Labels keep the old
edge-degree signature as species: ((3,4,4), ()) bare, ((3,4,4), (5,)) with a
5-fold node; `nodes` counts only all-legal vertices.

DECIDED AGAINST D-side defect tracking. Prototyped a vertex-anchor map in
Manifold (star() is an O(N) filter, no incidence map exists) then REVERTED it:
the only consumer needing fine resolution already replays the event stream,
and the 8 producers pay just ~3% for the full rebuild at TS=150. Everything is
reconstructable from the move record. Would only revisit if the SAMPLER itself
needs defect structure (e.g. defect-targeted move classes) -- replay is
inherently after the fact.

knot_worldlines.py rebuilt on it with PER-MOVE identity over whole complexes
(a Pachner move touches <=5 vertices, so unmeeting complexes are untouched;
the rest relabel by majority inheritance). This KILLED an artifact: the old
version tracked a knot by its chord (a vertex pair) whose lift position is
constant by construction, so "MSD=0 with slides off" was a tautology. Real
numbers, lam=0.40, 300 sweeps, cells^2: dt=1 0.0030 off vs 0.0082 on (2.77x),
dt=3 0.0149 vs 0.0289, dt=8 0.0412 vs 0.0515. Thermal transport is REAL;
slides ~2-3x the short-time MSD while being only ~10% of steps.
CAVEATS: large-dt is survivorship-biased (non-monotonic past dt~8); centroid
MSD conflates translation with shape fluctuation, so trust the RATIO.

CORRECTION (same day, commit 7274167): pass1_kinematics is NOT retired. I
claimed its overlap linkage (>=2 shared or Jaccard>=0.3 across 150-sweep
frames) invents deaths where defects moved, one-sidedly suppressing the
slides-on arm. `linkage_bias.py` tests both linkages on identical data:
false-death rate is 0.0% (slides off) and 1.1% (slides on). The hypothesis is
REFUTED; pass1's kinematics stand and its "blinking" (82% of deaths followed
by nearby rebirth) is corroborated by exact identity.

Why both agree is the real finding: DEFECT COMPLEXES ARE TRANSIENT at
lam=0.40. Only ~20% survive a 150-sweep frame (19-22 continuations vs 73-89
deaths); the exact tracker sees 38-72 merges and 41-72 splits per ~10^4 moves.
The population is a reconnecting network, not a gas of persistent particles.
This BOUNDS the MSD result above: tracks live ~1 sweep, so large-dt points come
only from rare survivors. dt=1-3 (n~10^4) stand; do NOT read the curve as a
diffusion coefficient without conditioning on complex size/lifetime.

STILL TO DO: point
the 8 producers at defect_state (emit BOTH old imp>0 counts and broadened ones
so existing time series stay comparable). Leave the static structural
analyzers (pass2_structure, survivor_anatomy, carrier_gr, sl_verdict, web_*)
alone -- full scans are right for structure-at-a-state.

DEFECT SPECIES CLASSIFICATION (2026-07-25, commits a9a7267..fc5f4a1).

`species_report.py --group {sig,induced,decorated,full}` over an ensemble,
with the CONVENTIONS.md state-report block. Four keys, increasing resolution:
  sig        illegal-edge signature (historical, e.g. "(3,4,4)")
  induced    isomorphism class of the induced subcomplex (bare shape)
  decorated  + ambient edge degrees on edges
  full       + n6 on vertices  <- pins Z, hence curvature, at every vertex

D side: Manifold.inducedSubcomplex / closedStar (commit f64baea), both taking
a transient vertexFacetIndex so one build serves many extractions; pinned to a
scan-everything oracle in simplicial_complex.d by unittest. Python mirrors in
defect_state (induced/star/induced_facets/decorations), brute-force verified.

KEY RESULTS at lam=0.40 (79 snapshots, 15 chains, 708 complexes):
* (3,4,4) is the most common species, 33.8%, and RIGID: it always induces
  f=(5,10,9,3), 1-skeleton = K5, canonical key == three tets of the boundary
  of the 4-simplex glued in a row -- a stacked 3-ball, chi=1. Identical in all
  instances. The closed star is also a ball (f=(35,143,190,81), chi=1).
* The induced subcomplex is DENSE, not curve-like. The curve is the
  ILLEGAL-EDGE subgraph: for (3,4,4) that is 3 edges, degseq (2,1,1,1,1) = a
  2-path plus a disjoint edge (2 components).
* Edge decoration PURIFIES the knot class (242 mixed -> 239 pure (3,4,4)); the
  shape and signature labellings agree once shape carries the physics.
* THE KNOT IS A CHARGE LADDER. Adding n6 to the key splits it into rungs
  indexed by TOTAL COORDINATION sumZ over its 5 vertices:
     sumZ 67/68/69/70/71/72 -> Q_c -0.0074/-0.5586/-1.1099/-1.6612/-2.2125/
     -2.7638 rad, spacing exactly pi-3*theta = -0.5512856. sumZ=70 dominant.
  Three of the five vertices always carry n6=(2,3,3); only two vary.
* USE RAW CHARGES (no background subtraction). Q_c = sum q_R is reference-free
  so identical complexes give bit-identical values: sd collapses to 4.4e-16
  (machine epsilon). Subtracting the state mean qbar ties the measurement to
  the whole configuration and manufactures sd ~ 1.8e-3 -- I first mistook that
  artifact for an irreducible residual. Report the relative value too, but raw
  is primary. CLOSED FORM (exact on all 60 classes):
     Q_c = (sum Z)*(pi - 3*theta) + 6*n_verts*theta
  so charge is fixed by ONE integer, the total coordination.
  => CONVENTIONS.md's "knot Q_c ~ -1.19 +- 0.18" is averaging over rungs, not
  one value; the +-0.18 is narrower than the 0.551 rung spacing.

USEFUL IDENTITIES (all verified on 8000 vertices, legal and illegal alike;
they follow from the vertex link being a triangulated 2-sphere):
    sum_{e at v} deg(e) = 6Z - 12
    n6(v) = Z + 5*imp(v) - sum(d_illegal at v) - 12
    q_R(v) = Z*(pi - 3*theta) + 6*theta      <- depends on Z ALONE
So q_R is redundant given Z; n6 and Z are interchangeable given the edge
decoration; and every illegal edge lies INSIDE an induced defect subcomplex
(both endpoints have imp>0), so edge colours already fix imp and sum d_ill.
