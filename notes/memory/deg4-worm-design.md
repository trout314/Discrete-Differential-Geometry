---
name: deg4-worm-design
description: "deg-4 worm move program — tip motif, design principles, search status (retract NOT found at depth ≤6 on first targets; sweep + deep probe running)"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-30T14:54:23.470Z
---

Goal: a slide/worm move for DEG-4 edges (analog of the validated deg-3 knot
slide, [[nonlocal-slide-move]]) that gives small dS at typical tips of the
linear defect complexes.

**Established structure (2026-07-29, λ=0.35 R m4 gas snapshots):**
- Deg-4 disclination lines are short open filaments (1–5 edges); complexes are
  topologically thin (induced E/V ≈ 2 vs bulk 6.7).
- UNIVERSAL tip motif: every line-end vertex has star {1×deg-4, 10×deg-5,
  rest deg-6}, charge exactly 12 — 100% of tips in both dilute and
  over-defected states (figs data/rxn_lam035_m4/deg4_line_ends.png,
  complex_shape_topo.png). No deg-3/deg-7 baggage anywhere on the lines.
- 4→4 hinge flip baseline (hinge_dS_by_position.png): isolated dS~1.9
  (⟨acc⟩0.34), tips 6.4 (reach downhill), interior 4.3 (never downhill),
  branch 7.7. All flips raise n_ill (+1..+3) — it roughens, doesn't transport.

**Design principles (agreed with user):** small dS ⟺ endpoint congruence
(translate the tip motif; dS depends only on endpoints); search between KNOWN
endpoints; allow dirty intermediates; primitive pair = EXTEND/RETRACT
(inverses, ±1 line length); frame/slots from the universal motif via
chain_walk (constant slot count like SLIDE_SLOTS). 4→4 = (2→3)∘(3→2), so a
2→3/3→2 search engine covers hinge moves at depth 2.

**Search verdict (2026-07-29, POSITIVE-CONTROLLED):**
- Positive control PASSES: planted 5-edge cluster (worm_deg4 CREATE on pristine
  m2) → search finds TWO 3-move full retracts (dS −10) in 9s. Search logic
  (goal/focus/caps/dedup) verified to find real solutions.
- Thermal states: NO retract anywhere. Lone[0] depth ≤7 exhausted +
  depth-8 partial (40M nodes, ILLX=5): zero. Sweep 73+/97 lone+tip targets at
  depth ≤6: ALL zero. Local annihilation of thermal deg-4s DOES NOT EXIST at
  small depth/radius — they are structurally different from planted clusters
  (non-locally compensated / paired remnants).

**Move-record provenance (deg4_provenance.py, 1000 sw on the gas, all audits
pass incl. hinge diagonal=(c0,c2) convention; data/rxn_lam035_m4/deg4_provenance.jsonl):**
- 4→4 HINGE FLIPS ARE 73% OF ALL ACCEPTED MOVES (18/sw); 2→3=3→2≈3.3k/1000sw;
  volume moves zero. Median lone-deg-4 residence = 1 sweep (p90 7).
- Wild deg-4 "deaths" are RELOCATIONS/EXCHANGES, never annihilation: example
  hinge event kills 2 deg-4s and creates 2 elsewhere in one move. Fates of
  lone deg-4s: 55% hinge-relocated, 31% healed by neighbor's flip, 11%
  healed by 3→2 (flicker absorption), 2% demoted to 3.
- Resolution of the paradox: deg-4 edges are simultaneously locally
  indestructible (search), fast-churning edge-identity-wise (records), and
  caged as coarse objects ([[defect-travel]]) — the hinge walk is recurrent.

**HINGE FLIPS CANNOT TRANSPORT (user caught this):** the orbit of iterated
4→4 flips is exactly the 3 diagonals of one 6-vertex octahedral cavity —
proven by the flip algebra. Sequences of flips never leave the cavity. Wild
churn = flips rattling in cavities while rare 2→3/3→2 deform the cavities.

**Full retract sweep FINAL: 0/97** (all lone + line tips, depth ≤6; deep
probe depth ≤7 + partial 8: zero). No local annihilation, universal.

**LONG-LIVED STRUCTURES (longlived_structures.py; tracks × canonical_key):**
lifetime ranking (med sw): 6-edge deg-4 bundle f=(10,15,5,0) 58.5 > 8-edge
42 > 10-edge 31 > DECORATED lone deg-4+n6-node f=(3,3,1,0) 23.9 (2.4k trk)
> shared-vertex deg-4 pair (4,4,1,0) 18.8 (max 1530) > 4-edge bundles ~14 >
bare pair (3,2,0,0) 9.5 > n6-monomer (1,0,0,0) 4.7 (26k trk) > flicker
(5,10,9,3) 2.8 (72k trk, 60%). KEY: long-lived complexes = BUNDLES of
disjoint 1-3-edge deg-4 fragments stitched by LEGAL deg-5/6 edges + n6
nodes; induced complexes sparse (top: 0 tets); NO deg-3 anywhere in the
long-lived sector (why deg-3 slide can't equilibrate it). tracks.jsonl only
has finalized tracks (death/merge) — lifetimes censored, maxima are lower
bounds.

**FINAL VERDICT (2026-07-29, committed 1761202):** CLEAN DEG-4 TRANSLATION
DOES NOT EXIST. Fixed-goal translate search (other illegal edges must
survive): 0/12 bundle fragment tips at depth ≤6, excess ≤4. Wild record:
0 clean escapes + only 4 clean in-cavity hinge relocations in 24,558
accepted moves — nature never does it either. Deg-4 charge moves ONLY via
dressed multi-edge reactions.

**DISCOVERED REACTION (fusion_verify.py, D-core verified, dS=−0.73
DOWNHILL, f-vector neutral, 2 moves):** 3-path + 2-path → closed deg-4
TRIANGLE + lone deg-3 exhaust. First closed disclination loop seen; loop
closure is FAVORABLE (tips with their 10×deg-5 caps are the expensive
feature — continuum "lines close or end on defects" made discrete). Rings
have no tips ⇒ likely extra-stable species; check census for rings.
GOTCHA: illegal-set goals compare edge SETS — degree changes among
retained illegal edges (e.g. 4→3 demotion) pass silently; strict goals
must compare edge→degree maps.

**PAIR CENSUS (deg4_pair_census.py, commit ca717d5): 2,024,911 reactions,
106 tips.** Vocabulary: mostly expensive double-2→3 junk (dS 20–45);
fusion-assisted transport (dS −5..−0.7, consumes partner, self-limiting);
and THE WORM PRIMITIVE: deg-3-CATALYZED content-neutral transport
(tip(4,5) | 0>3 3>0 [4>5 5>4 5>4 variants]) — n=367, median dS 0.000,
100% escaped landings ~0.25 cells, needs a deg-3 dipole in range
(fired at 8/8 such tips), ~50 slot variants/tip. D-CORE VERIFIED: one
2→3 + one 3→2, f-vector + illegal composition exactly preserved; deg-4
hops out of its cavity while the catalyst is consumed & re-emitted.

**THE WORM WORKS (worm_walk.py, commit 20d290e):** 8-step ballistic walk
on lam35r_snap15000 from (3282,3318): 35–68 candidates at EVERY head
(self-composable), net 2.758 cells at alignment +0.61..+1.00 (vs caged
thermal max ~1.26 cells/lifetime), per-step dS ∈ {−0.42,0,+0.42} (sum
+0.21), composition {4:121,3:6} preserved (catalyst travels with head),
all 16 moves undone → EXACT facet restoration (invertible). First real
deg-4 transport. NOTE: walker is chart-greedy (maximize harmonic-chart
displacement dot direction), NOT chain-constrained — whether steps track
a BC spiral is a separate measurement; landings are often vertex-disjoint
from the head (~0.3–0.44 cells jumps).

**ORACLE COMPLETE (worm_deg4_slide.py, commits 1690c4f + 816e09b) — full
worm_slide.py validation parity:**
- Move class: one 2→3 + one 3→2; anchor deg-4 healed/removed; k∈{1,2}
  healed/born pairs (k=1 subclass is closed but IMMOBILE — ping-pongs
  between two pivot edges; transport = the k=2 fragment-as-unit family);
  ≤1 catalyst (deg-3) relocation; symmetric octa escape; strict illegal
  edge→degree MAP preservation elsewhere.
- Hastings: k=2 breaks anchor uniqueness ⇒ anchor-sum
  q=(1/n4)Σ_a k_f(a)/n_f(a), n4 cancels (class preserves it).
- Closure test 8 anchors: 69/69 transitions inverse-closed, dS exactly
  antisymmetric, k_f/k_r ∈ 1..4. 4/8 anchors zero candidates =
  catalyst gating (correct rejection branch).
- MH integration vs live sampler: tracked objective == oracle dS to
  3.6e-14 on every accept; cocycle attached; composite kernel drift
  5.7e-14; 31/40 catalyst-gated, 5/9 live-proposal acceptance.
**D-CORE PORT DONE (2026-07-30):** sampler.d: WormConfig + WormCand +
wormEnumerate (star-BFS regions, class checks, canonical net-key) +
tryWormMove (anchor-sum Hastings by enumeration at all gone4/new4 anchors,
speculative+potential lockstep dS, exact rollback); mcmcStep channel
(anchor 1/n4 from a second Deg3Set holding deg-4 edges, reconciled at every
commit site; slide/nonlocal commits rebuild it O(N)); NOT cocycle-safe
(gated off when cocycle attached). capi set_worm_prob/worm_stats/worm_at;
Python set_worm_prob/worm_stats/worm_enum/worm_at. VALIDATION: end-state
crossval vs oracle PASS (all oracle end states reachable, dS exact; D is a
strict superset — no radius-2 patch cutoff); live channel 67–74% acceptance
among live proposals (~5.5 proposals/sweep at prob 1e-4), objective drift
0.0 vs rebuild, maps validate, census audit clean with worms on. GOTCHAS:
collectStar's visited set must be per-call (shared-array BFS blocking bug);
the "landing" label is arbitrary for k=2 double-escape candidates (oracle
dict order vs D AA order) — crossval must compare END STATES, not landings.

**RUNNING: worms-on ensemble, ATTEMPT 3** data/rxn_lam035_m4_worm w0–w7
(8×5h, worm 1e-4, slide 0, seeds 201–208, 4 dilute + 4 over starts) vs
baseline slide-off arm (data/rxn_lam035_m4 c4–c7): stationarity slopes /
melting drift, k2/k1, compound survival. reaction_census.py has
--worm-prob. Attempt 1 KILLED (GC leak, see [[embedded-gc-rt-init]];
partial in _killed_partial). Attempt 2 ran 5h but ALL NUMBERS GARBAGE:
tryWormMove didn't mirror its 2 moves into the event log, census tracker
diverged (AUDIT FAILED, fictitious nill≈5000 melt) — kept in
_audit_broken. FIX: tryWormMove ledger param (post-hoc recordBistellar +
logEvent, safe post-move; sixFlips NOT — worm gated off under
logSixFlips); verified worm 1e-3 run: 1726 accepts, audit_failures=[].

**TEMPLATE CENSUS (2026-07-31, lam35r_snap15000, all 121 anchors): NO
small template library.** 454 candidates at 24 live anchors (97 catalyst-
gated) -> 137 distinct move classes (exact decorated-iso certificates,
patch+move-pair canonicalized jointly). FLAT distribution: top class 3.5%,
top-20 only 34% coverage. Link-pair domain nearly all singletons: 115
distinct decorated star(a)∪star(b) classes among 121 anchors; 190/454
candidates' second move reaches OUTSIDE the link-pair patch. User's
discretization point confirmed: dS quantization (λζ washboard) coexists
with huge geometric diversity. VERDICT: enumeration IS the production
algorithm; template distillation dead for the thermal gas (possible only
as clean-tip sub-library for the lift). Script:
scratchpad/template_census.py, results template_census_results.json.

**WORMS-ON ENSEMBLE RESULT (2026-07-30, attempt 3, audits clean):** worm
1e-4 does NOT cure the λ=0.35 m4 melting drift. Complete chains: w2
(dilute, 12.9k sw) late-window slope +18.4±2.0 /1000sw, w5 (over, 18.6k
sw) +4.9±1.2 — same range as baseline c4–c7 (+1.6 to +28.9). No R̂_q
pass anywhere (baseline 4-chain R̂=1.77 ESS 17; worm 2-chain R̂=2.98).
Channel active & healthy: 72–74% acceptance, ~1.5 tries/sweep, ~74% of
firings catalyst-gated (no candidates). Chemistry (merge/split ~0.98)
unchanged. CAVEAT: 6/8 worm chains killed at ~10–15k sw for user
experiments (no .json ts; events/tracks jsonl valid) → arm comparison
power limited. Data: data/rxn_lam035_m4_worm (attempt-3);
_killed_partial and _audit_broken dirs are unusable siblings.

**NEXT: BILOCAL WORM (design written 2026-07-30,
notes/bilocal-worm-design.md):** joint-balance two-region proposer —
unbalanced half-moves at A and B (the pair census's "junk" sector is the
alphabet), admitted on the JOINT ledger, paired by pluggable kernels (P1
chain-walk with bidirectional-agreement admission, P2 uniform control, P3
defect-to-vacuum later). Key corrected physics: Σ(6−d)=12 per vertex is
link-Gauss-Bonnet, KINEMATIC — no charge obstruction to unbalanced
patches; local-worm balance was only the empirical small-depth search
negative. dS additive when patches disjoint (endpoint congruence ⇒
distance-independent cost). Phase 0 measurements first: chain-reach
census + half-move type census. q_R = curvature, not "charge"
(terminology note for CONVENTIONS.md pending user OK).

**RE-VERIFIED POST-BUGFIX (2026-08-03, lam35r_snap15000, independent of
wormEnumerate's own filters — everything re-derived from raw facet sets):**
- LOCAL: a committed candidate changes exactly 10 facets (= 2->3 plus 3->2)
  and 20 edge degrees, max chart distance among them 0.431 cells. Zero
  action at a distance.
- ESCAPES THE HINGE CAVITY: 255/255 candidates satisfy the D core's own
  condition (some landing outside the anchor octahedron) AND the symmetric
  end-state condition, checked independently. Strict readings: 84.3% have
  EVERY born edge outside; 90.6% have the born PAIR outside the union of the
  healed pair's octahedra (the fragment-as-unit reading — the right one for
  k=2, which is 253/255 of candidates). The other 9.4% are in-cavity
  rearrangements and are exactly the short-displacement tail (median 0.121
  vs 0.216 cells) — filter on displacement if transport is the goal.
  Iterated 4->4 flips are confined to one cavity, so outside = unreachable
  by any number of hinge flips.
- ROLLBACK EXACT: trial every candidate with commit=False, facet set
  unchanged, 0 corruptions.
- Sampler construction does not perturb labels (facet set and edge-degree
  map identical to a plain Manifold.load).
- **ITERATES BALLISTICALLY:** greedy directional march, 20 steps, dS = 0 at
  every step, illegal multiset + f-vector preserved at every step:
  net **7.27 cells** with path 7.68, **straightness 0.94-0.95** (2 of 4 tips
  tested). Box is 4x4x4 cells, so it wrapped the torus ~twice; caged thermal
  max is 1.26 cells/lifetime. Step size median 0.39 cells. This reproduces
  and extends worm_walk's net 2.758 over 8 steps (~0.39*8*0.95 = 2.9).
- The other 2 of 4 tips PING-PONG: net exactly 0.000 over 20 steps,
  straightness 0.00, anchor alternating between two edges — the documented
  k=1 immobile subclass. A directional driver must forbid returning to the
  previous anchor.

**MEASUREMENT GOTCHA (cost two wrong readings this session):** the midpoint
of an edge straddling the periodic boundary is NOT (X[a]+X[b])/2 — that lands
in the middle of the box and fabricates ~2-cell "jumps" (half of a 4x4x4 box
= the maximum minimum-image distance, so it looks like a huge but plausible
displacement). Correct form is mid = X[a] + 0.5*minimg(X[b]-X[a]), and net
displacement must accumulate minimum-image STEP vectors in the covering
space rather than wrapping the endpoint. worm_walk's ref-based `mid()` was
always right; only the scratch checks were wrong. Chart itself is clean —
all 68167 real edges are 0.10-0.46 cells, none over 0.5.

**WORM AVAILABILITY MEASURED (2026-08-03, see [[lifted-ecmc-transport]]):**
all worm candidates are dS = 0 EXACTLY under a pure edge action (forced:
illegal-degree-map preservation + f-vector neutrality pins n5 and n6, so the
whole degree multiset is preserved); INTERIOR deg-4 edges are 0% live
(0/18 pooled) so the worm is strictly a TIP move; live fraction 20% (lam.35)
/ 45% (lam.40), gated by graph distance to the nearest deg-3 catalyst
(dist 1 -> 73-100% live, dist >=8 -> 0%).

STILL AHEAD: event-chain lift (momentum sector); cocycle support if needed.

Related: [[reaction-census-lam035-m4-results]] (lone deg-4 monomers are the
dominant blinking species), [[five-illegal-knot]].
