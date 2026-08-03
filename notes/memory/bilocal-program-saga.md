---
name: bilocal-program-saga
description: "Master narrative of the bilocal move program — two carriers built, what was tried and refuted, the structural impossibility lessons, measured constants, bug taxonomy, and the full script inventory"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-08-02T14:53:47.073Z
---

The synthesis memory for the bilocal-move program (through 2026-08-02).
Detail lives in the linked memories; this is the arc, the dead ends, and
where things are.

## The goal

MCMC moves acting at TWO separated locations at once, so the sampler can
(a) transport a defect non-locally instead of diffusing it, and (b) reach
rearrangements a local walk is barrier-blocked from. Design doc:
`notes/bilocal-worm-design.md`.

## The two carriers

**CHORD carrier** — `wormChordStrictEpisode` (sampler.d). Relocates a
degree-3 edge (a "flicker"). Both marks are pure flags; the walk is
2<->3 only; the close fires when one mark is a degree-3 chord, the other
absent, and the vacated region is clean.
STATUS: **works, balanced, fast, and does not help.** Two-sided
certificate passed twice (2026-07-31, re-verified 2026-08-02). Made
15.7x faster on 2026-08-02. Then MEASURED not to accelerate relaxation
— [[percolation-ab-test]], [[strict-chord-channel]].

**VERTEX (pair) carrier** — `wormPairEpisode`. One 1->4 at an opened
location + one 4->1 at another, so it conserves f0 exactly.
STATUS: **opens and closes on ROUNDTRIPS, transport never fires.**
94/200 closes at mu=0, median close log alpha -0.26, after five bug
fixes. But the adopted ball's zmin never moves off its starting value,
even at a 60,000-step cap: U is a function of the spoke MULTISET while a
collapse is a PATH, the single-orbit tube covers 1 of 1533 vertices, and
the linear fallback has per-step residuals to 7.52.
— [[pair-carrier-calibration]]

**Both carriers conserve f0 EXACTLY**, so neither can move the f0
read-out alone. The pure-f0 composite would be E = vertex - 3*chord,
giving Delta f = (+1, 0). Not built.

## Structural results — things that CANNOT work (all proved, not guessed)

1. **The umbrella is self-cancelling.** Open weight ~ zeta e^{-S+U}, close
   ~ e^{-dS-U}: a high U draws the walk INTO a state and BLOCKS closing
   from it, same parameter, same sign. No U makes a state attractive to
   reach and easy to close from. Measured: big tube entries drove
   abandonment to 101/200.

2. **Selectivity cannot live in the closure criterion.** A criterion is a
   predicate on the open state, and the reverse of the close IS the open,
   so both gates must be the SAME predicate — hence the close condition
   already holds at the open and a null closure is always available.
   Selectivity belongs in the MEASUREMENT: post-selecting on the
   per-episode `df` census is exactly unbiased.

3. **The invariance theorem.** log alpha_open + log alpha_close = -Delta,
   so P(open)*P(close) <= e^-Delta where Delta = bracket_fwd -
   bracket_rev. **No zeta can fix a bracket mismatch**; only reducing
   Delta helps. This kills tuning-based rescues generally.

4. **Sequential Metropolis cannot transfer action between moves.**
   Annihilating a flicker releases its dS at its own step; the next
   uphill move still pays e^-X in full. So no per-move-accepted channel
   can lower a barrier. The only mechanism that would is BUNDLING —
   pricing the annihilation and the barrier in ONE acceptance, e^-(X-Q)
   instead of e^-Q * e^-X. Never built. — [[flicker-catalysis]]

5. **Concentration entropy, four times over.** Scheme B's guidance /
   acceptance tension (e^-70), the anti-umbrella impossibility, the
   "battery" idea's n3/N_s cancellation, and clustering entropy are all
   the same lesson: *any scheme that concentrates a diffuse resource pays
   its concentration entropy.* — [[pair-carrier-calibration]]

6. **Tier 1 vs Tier 2** (design doc 2.5): self-mirror alphabets work
   (6/6, knot slide is the existence proof); deep/dressed halves do not
   (0/24) — "no frame fixes that".

## Measured physics (constants worth not re-deriving)

- **Flicker annihilation quantum = 5.930 EXACTLY** (sd 0.000 over all 17
  in quench_down5q_wOFF), on a (5,5,6) face in a pure d5/d6 region.
  Independently equals the design doc's self-mirror value +-5.93.
- **NOT universal**: fresh-site creation 6.191, recreation at a vacated
  site 5.930; the 0.261 gap is exactly one step of the quadratic pins
  (volume 0.20 + edge 0.07).
- **NOT an ideal gas**: the quadratic volume pin couples flickers
  all-to-all through f3; the (m+1)-th costs 0.1(2m+1) more.
- **Catalysis is REAL but oracle-dependent**: one flicker adjacent to a
  collapse target drops the barrier 29.40 -> 17.87, each further flicker
  exactly -5.93, four flatten it to 0.00. Free-energy barrier 20.60 ->
  ~9.82 after clustering entropy. BUT every measurement used the PLANNER
  as an oracle to place the flicker; delivery by the channel is unsolved,
  and BC chains show no enrichment (0.52x, chance).
- **Bilocal factorization is exact iff supports are VERTEX-DISJOINT**
  (1e-13), not a function of distance. — [[bilocal-factorization]]
- **A lone flicker on pristine crystal** = 5 illegal edges
  {3:1, 4:3, 7:1}, which form illegal-edge-graph components {1, 3, 1}
  (the deg-3 apex chord shares NO vertex with its three deg-4 face
  edges). Pristine C15 m4 is exactly {5: 9216, 6: 1024}.

## Bug taxonomy — the expensive ones

**Silent-kill class** (channel reports a clean zero):
- `set_worm_pair(zeta2=NaN)` -> auto-stored 0 -> log(0) = -inf -> 0 opens
  in 2000 episodes. Now a loud capi error.
- `cfg.bcp` was inert; p_close was 1-ph-pg = 0.10, so every closure was a
  roundtrip.
- D default-initialises floating point to **NaN**: `double[2] dsArm;`
  made `if (d > res.dsArm[0])` never fire. Same hazard class as the
  standing "never `= void`" rule — zero-initialise explicitly.

**Wrong-model class**:
- One-way umbrella ramp applied to BOTH balls (close subtracted
  20.60+20.60 where the open added 20.60+0.85) -> 100% cap-undone. Fixed
  by the role-signed umbrella: created ball -U, adopted +U, and the close
  must price both under the REVERSE episode's roles.
- `calib_zeta2` used the neighbour count (13.4); the engine uses
  `_vertexDegrees` = tets in the star (22.8). Threw log(W/wv) off by ~17.
- Diagnostics aliased `res.umax`, which the pair open seeds with
  ball[cre].u + ball[adp].u = -19.75, flooring every readout.

**Measurement-artifact class** (I believed a result that wasn't there):
- `w(f) = [face has an edge of degree <= 4]` "refuted" — I computed the
  signature while the flicker was still placed, so I read POST-placement
  degrees (a 2->3 drops all three face edges by one).
- "Ideal lattice gas, 3% closure" was luck: (5,5,6) is a superset with
  sd 3.259, not the exact class.
- The disjointness prediction was vacuous: every candidate placement
  touched the target vertex.
- First chain-symmetry test created a second flicker instead of
  RELOCATING one.
- "One ring out" aggregation baseline (48%) was an artifact of a ~50-vertex
  neighbourhood; the support-overlap criterion gives 17.8%.
- The catalysis claim itself was published, then RETRACTED, then
  re-established on a different (barrier, not endpoint) measurement.

**Performance class** (see [[percolation-ab-test]] for full numbers):
- `Manifold.star()` is `facets.filter(...)` and `facets()` GC-dups AND
  SORTS all N facets before filter runs => ~1.4 ms per call at N=8704.
  Five sampler call sites used it just to grab ONE seed tet and `break`.
  Fixed by `_vertexWitness` (13x on the episode).
- `wf0ChainSites` phase 2 scales with **defect density** (start windows
  ~ n3), not just N. Low-density smokes underestimate badly — cost a 3x
  runtime miss and a 3x ep_frac miscalibration in one run.

## Script inventory

**D core** (`source/sampler.d`): `wormChordStrictEpisode`,
`wormChordPairEpisode`, `wormF0Episode`, `wf0ChainSites`,
`wf0ChainSitesW` (+ WF0Site collect mode), `wf0StartWindows`,
`wf0RegionClean`, `wf0ChordDegs/Key/U`, `Deg3Set`, `collectStar`.
`source/manifold.d`: `someFacetContaining` / `_vertexWitness`.

**Committed drivers** (`scripts/defect_dynamics/`):
- `f0_worm.py` — the main worm driver (both carriers, all knobs)
- `f0_channel.py`, `harvest_f0.py`, `harvest_f0_verdict.py`
- `twosided_chord.py` — **the balance certificate**; env START OUT NEP SEED NPREP
- `link_planner.py` — (link, collar) 2D planner, `optimize='z'|'cost'`
- `catalysis_search.py`, `catalysis_wide.py`, `catalysis_audit.py`
- `flicker_quantum.py`, `flicker_site.py`, `flicker_spectrum.py`,
  `flicker_subclass.py`, `flicker_catalysis_barrier.py`
- `chord_tube.py`, `dressed_generators.py`, `tip_retract_search.py`
- `worm_deg4_slide.py` (the Python `Live` index), `worm_slide.py`,
  `deg4_pair_census.py`, `deg4_provenance.py`

**Scratchpad, uncommitted but worth keeping** (session
d562feea-0887-419b-9cd9-9899d180c7c2):
- `perc_ab.py` — the A/B driver with honest work accounting
- `perc_pilot.py` — density ladder for percolation
- `deg4_observable.py`, `deg4_pilot.py`, `deg4_density_scan{,2}.py`
- `ep_cost.py` — episode cost vs plain sampler, with N-scaling
- `chainsites_validate.py` — enumeration-count equality across builds
- `planner_timing{,2}.py`, `planner_profile.py`, `planner_hoist.py`
- `agg_knobs_test.py`, `aggregation_baseline.py`
- `bilocal_factorize.py`, `bilocal_range.py`, `flicker_ladder.py`,
  `corridor_support.py`, `wf_vs_planner.py`

## Where it stands

A **correct, certified, fast channel that does not do its job.** The
chord carrier is balanced and cheap; the vertex carrier cannot transport.
The catalysis physics is real but only reachable with an oracle the
channel does not have.

Open threads:
- BUNDLING is the one untried mechanism the structure permits (item 4).
- The pure-f0 composite E = vertex - 3*chord is unbuilt.
- The vertex carrier needs its adopted-ball seed biased toward
  flicker-adjacent vertices, or it will never see a cheap corridor.
- `notes/STATE_AUDIT.md` — duplicated-state audit, inventory only;
  `Deg3Set` is the highest-risk item and would break balance SILENTLY.

METHOD LESSON that cost the most time: measure before optimising, and
test end-to-end behaviour before trusting a component-level check. Two
hot-spot diagnoses without measurement produced one reverted regression
and one near-miss; a commit went in on an enumeration-count check that
structurally could not detect what it changed.
