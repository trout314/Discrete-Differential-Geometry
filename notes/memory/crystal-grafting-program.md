---
name: crystal-grafting-program
description: "Grafting lumps of one TCP crystal into another via exact decorated-boundary-surface matching: L1/L2/L3 match hierarchy (L3 = zero-cost seam), graft_signature.py library validated, C36<->C14 slab control"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1fc2ced4-2271-4a60-ba17-bc68ad6d30ad
  modified: 2026-08-06T18:23:24.920Z
---

**Goal (Aaron, 2026-08-06):** graft a piece of one TCP crystal into another by
finding surfaces with isomorphic decorated graphs. Reframed: search over
LUMPS (tet sets) and match their decorated boundary surfaces.

**Match hierarchy (each exactly checkable, proofs in graft_signature.py docstring):**
- L1: simplicial iso of boundary 2-complexes → glued object is automatically a
  closed 3-manifold (disc ∪ disc = S² links).
- L2: + per-edge (total degree, lump-interior tet count d_in) preserved → ALL
  edges of the graft stay in {5,6}.
- L3: + per-vertex (cn, n6, n_interior_edges, n_interior_6edges) preserved →
  disclination web reconnects (no Z-class change, no broken lines; "lines
  can't end" is a theorem in the pure-{5,6} sector, see fk_skeleton.py).
  **An L3 graft has exactly zero local energy cost under any degree-only
  action — the seam is invisible to the sampler.**
- Vertex Z classes alone (the original idea) are NOT sufficient; edge splits +
  web continuity are the real conditions.

**Physics motivation:** the pin gap (e* − e_native) is extensive with opposite
signs across the library ([[crystal-library-gas-campaign]]); a mixed-sign L3
composite at the right volume fraction has ZERO net pin gap — no pure crystal
except R can do that.

**Infrastructure:**
- `scripts/graft_signature.py` — CrystalContext, lump_boundary (decorations +
  surface-manifoldness checks), EXACT canonical certificates by
  combinatorial-map BFS over all start darts (no WL hashing — exact, and the
  canonical vertex orders double as φ candidates), match_phis, graft surgery
  (fresh interior ids + compaction), validate_facets (link sum rule assert,
  {5,6} census, broken-line count, χ, orientability). `--selftest` = Z12
  star-swap on in-memory c15 m3: PASSES (12 decorated link symmetries,
  nontrivial φ, Z census preserved).
- `scripts/graft_c36_control.py` — positive control: C36 = hc Laves polytype
  interleaving C14-like/C15-like stacking. Enumerates basal slabs (cuts
  between atomic z-planes; c36 m4 has 16 planes/cell, c14 m4 has 8),
  certificates, within-crystal groups (same boundary + different interior =
  stacking-swap candidates), cross-crystal join, then performs + validates the
  best c14→c36 graft, saving to data/grafts/.
- Orientation gotcha: gluing along both slab boundary tori needs φ
  orientation-compatible on BOTH simultaneously or the result is
  non-orientable — handled by iterating match_phis candidates until
  validation passes.

**CONTROL RESULT (2026-08-06, m=4, kmin 3 kmax 12):** c36: 156 slabs → 6
distinct L3 boundary classes; c14: 78 slabs → the SAME 6 classes (identical
certs) — c36 and c14 share their entire basal boundary-class set, and
L1=L2=L3=78 cross pairs (for Laves basal slabs the surface complex already
determines the decorations). Every class contains slabs with different
interiors (stacking-swap candidates). Performed a THICKNESS-CHANGING graft:
12-plane c36 slab (1632 tets) out, 3-plane c14 slab (544 tets) in, first φ
validated: all {5,6}, 0 broken lines, χ=0, orientable, Z census pure
Z12:Z16 = 896:448 = 2:1 → a new 7616-tet mixed-stacking Laves polytope-box,
saved data/grafts/T3_C36m4_graftC14_0.8905-1.6405_f37616.mfd (+ summary JSON
c36_control_m4.json). ~1.5 min/crystal census, single core.

**C15 CROSS RESULT (graft_c15_cross.py):** epitaxy arithmetic: the (111) wrap
of a cubic m-box is a 2m x 2m hex torus of the Laves net, so **c15 m=2 is the
epitaxial partner of c14/c36 m=4**. c15 m2 (111): 12 slabs → 4 L3 classes.
Cross joins: c15×c14 L3=216 (L1=344), c15×c36 L3=432 (L1=688) — first case
where **L1 > L3**: some surfaces match topologically but NOT with decorations
(would be defective if glued; c36↔c14 basal had L1=L3). Grafted a cubic C15
(111) slab into BOTH hex hosts, first φ each: all {5,6}, 0 broken lines,
χ=0, orientable, Z12:Z16=2:1. Both grafts CHANGE total V by −192
(768→576 c14 host; 1536→1344 c36 host) with a zero-defect seam —
**vertex-count-changing grafts do NOT force non-FK seam vertices** (Aaron's
suspicion refuted within the Laves family; thickness-changing slab swaps
change V freely at L3). Saved: data/grafts/T3_C14m4_graftC15m2_f33264.mfd,
T3_C36m4_graftC15m2_f37616.mfd, c15_cross.json.

**BALL-GRAFT SEARCH (Aaron's ask: dV != 0 ball grafts, e.g. 3-vertex star
union → 2-interior filling).** Reuses the FK-move-search formalism
([[fk-move-search]]): ball graft = pair of FK-legal fillings of one FRAMED
boundary (= L2 decoration on S²) with different V_int. New
`scripts/graft_ball_search.py`: best-first refilling enumeration =
enumerate_fillings.py moves + the missing 1→4/4→1 channel (vertex-preserving
flips can NEVER change V_int) + exact framing bookkeeping (the old filler's
2→3 can silently change boundary-edge framing — my acceptance requires
frame_dev==0; NOTE this means enumerate_fillings.py itself may collect
framing-broken "fillings"). Lazy priorities: children scored by exact delta
bookkeeping, analyze only on pop. Acceptance = L2-exact; full L3/Z check by
validate_facets on the performed graft.
**c15 m3 results (max_bad 12, grow 12):** first flip out of a pristine
crystal filling costs ~8 (block-scale echo of "no single Pachner move
preserves FK"). Vertex stars (Z12, Z16) and edge cavities (edge1/2):
search EXHAUSTED → provably UNIQUE FK filling in radius. Face cavities
(V_int=3, 47-62 tets): cap-limited (100k pops); no dV hit yet, but best
V_int=2 filling reached total violation 3 (face2). Deep runs (400k
bad-first; 200k vfirst) launched 2026-08-06 evening.

**EDGE-CONTRACTION MOVE (Aaron's donor-free scheme, 2026-08-06 evening).**
"Delete two adjacent stars, cone the boundary from one new vertex" = EDGE
CONTRACTION (u,v)→w; validity = link condition = lump_boundary's pinch check.
`scripts/contract_relax.py`: exact local delta bookkeeping (verified vs full
recompute), apply, and A/B relaxation driver. **Measured on c15 m3:**
contracting a Z12-Z12 deg-5 edge costs remarkably little — CN(w)=17 with ALL
17 edges legal (m(w)=0, a pure-{5,6} "Z17", n6=5), just 3 ring edges go
deg-4; Δ(f0,f1,f3) = (−1,−6,−5), Σm² +10. Deg-6 ring edges drop to LEGAL 5.
**A/B verdicts:** (1) at c_imp≤0.5 the gas isn't f0-stuck (sampler sheds
inserted vertices in <100 sweeps; at 0.2 equilibrium f0 < crystal). (2) at
c_imp=3.0 f0 is FROZEN: prep-at-inflated-f3-target then quench → sampler
stuck at f0=650; +contraction channel reaches exactly 648 instantly (guard),
BUT its debris can't heal at 3.0 (S_B > S_A). (3) UNGUARDED contraction
ratchets BELOW crystal f0 (volume pin pays for it) and strands defect-rich —
the reverse (vertex-split) channel is needed for equilibrium; guard
(f0 > f0_ref) is the annealing stopgap. D-port design: split = choose vertex
+ splitting cycle in link; detailed balance needs per-vertex cycle counts
(bounded-length cycle enumeration in links, or restrict both channels to
deg(uv) ≤ L).

**D-SIDE CONTRACT/SPLIT CHANNEL (built 2026-08-06 late, all meson tests
pass):**
- `source/link_cycles.d`: bounded simple-cycle DFS on adjacency bitmasks
  (allocation-free countCycles + rank-select kthCycle for uniform sampling
  without enumeration), FK-catalog (Z12/14/15/16 canonical polyhedra
  embedded; counts cross-validated vs independent Python enumeration: Z12
  5-cycles 72, 6-cycles 240, total<=8 2702 etc.), cycle lists length-sorted
  (prefix property → O(1) draws under a cap), exact face-based dart-
  propagation matchCatalog with vertex perm for cycle transport.
- `manifold_moves.d/manifold.d`: ContractMove/SplitMove; hasValidContractMove
  = full 3-level link condition (vertex/edge/triangle — vertex sets alone are
  NOT sufficient); plan/commit split (planContractMove/planSplitMove fill
  removed/addedFacets with NO mutation → speculative deltas; commitPlannedMove
  applies; undo exact). Split side convention: side 0 = component containing
  lexicographically smallest link face; freshSide0 flag. Tests: ∂Δ⁴ rejects
  all contractions (level c); contraction inverts 1→4 exactly; split w/
  triangle γ = the 1→4; exhaustive split→contract roundtrips.
- `sampler.d`: generic block-move deltas (speculativeBlockDelta /
  potentialBlockDelta driven by planned facet lists — reusable for ANY block
  move); tryContractSplit with exact Hastings (contract fwd deg(uv)/(6f3),
  reverse needs cycle count of MERGED link by DFS on planned state; split fwd
  deg3(w)/(4f3)·1/N_L·1/2side, catalog fast path + transported cycles);
  channel gate in mcmcStep (new trailing ContractSplitConfig* param), gated
  off under cocycle/six-flips like the worm. maxRing caps BOTH directions
  (must be capped together or DB breaks). Label pool: contract pushes v,
  split pops (mirrors 1↔4).
- capi `ddg_sampler_set_contract_split(prob, max_ring)` +
  `..._contract_split_stats`; Python `set_contract_split` /
  `contract_split_stats`. Integration unittest: 3000 mixed steps with
  potential on — zero objective/counter drift, both directions accept.
- Smoke on c15 m3 @ cimp 0.5: 30 sweeps, contract 4/2678, split 1/2690
  accepts, f0 648→645, manifold valid.
- `scripts/validate_contract_split.py`: VALIDATION-DESIGN SAGA (three traps,
  all diagnosed):
  (1) The naive A/B equilibrium test (channel off/on) is INVALID: the
  sampler's own 1↔4 bistellar pair carries an uncorrected O(1) proposal
  asymmetry (factor ~2 — THAT's why run_exact and importance_weight exist),
  so baseline A does not sample exp(−S); mixing in a CORRECT channel shifts
  B anyway. A/B "FAIL" ≠ channel bug.
  (2) Pair test with a weak invariant key (degree profiles) AGGREGATES
  isomorphism classes: many distinct γ's → key-equal targets, measured
  Σ P(x→y′) vs single reverse mimics a "missing 1/N" violation with
  aggregated forward rate EXCEEDING the true single-transition proposal
  bound (the 25× tell). Forensics variant-table (predict absolute rates per
  candidate mis-formula) exposed this.
  (3) EXACT test: labeled pair test on split-forward pairs — X labels
  0..n−1 compact, split creates fresh = n deterministically, reverse
  contraction of (w,n) keeps w<n and restores X BIT-FOR-BIT, so
  sorted-facet-bytes keys give single-transition rates. 150k-trial run:
  ratio tests pass (z −0.25 clean, +3.9 borderline); absolute rates exceed
  single-path theory by ~4× → (4) PROPOSAL MULTIPLICITY, the final piece:
  a labeled transition can be realized by SEVERAL (w,γ,side) paths (1→4-
  type splits by 4, one per tet vertex), each in BIJECTION with a reverse
  contract-edge path — and the per-path acceptance pairs exactly the
  corresponding reverse path, so each sub-kernel satisfies DB individually
  (composite-kernel MH, saturation included) ⇒ THE CHANNEL IS EXACT AS
  IMPLEMENTED. The multiplicity only makes single-path rate predictions
  underestimates; ratios are unaffected. Forensics counts labeled path
  multiplicity. **FINAL VERDICT (4 pairs × 150k trials, power-aware pair
  selection = smallest |dS| among trafficked probe targets): ALL PASS,
  z = −0.05, +1.06, −0.13, −0.19; multiplicity = 4 confirmed in each.
  Channel validated to exact detailed balance. All meson tests + graft
  selftest green.**
  Also fixed en route: split fresh-label collision on label-holey manifolds
  (pool seeded from fVector[0]) — graceful reject added (was a release-mode
  corruption risk); pair-test perf: cache the start Manifold, sampler
  copies it (~1.3 ms/trial).

**FIRST UNGUARDED QUENCH (native channel, c15 m3, cimp 3.0, committed
cb5d139):** prep at f3_ref+150/cimp 0.2 → quench. Channel OFF: frozen at
f0=653 (5 above crystal) for 1600 sweeps, S≈3950 static. Channel ON
(prob 0.1, NO guard): 45 accepted moves in the first cycle, overshoot to
f0=628, then detailed balance brings it back UP 628→634 with two-way
traffic (74 contract/51 split accepts) — no ratchet. S=2481 (37% below
control), n_ill 337 vs 368. **500-cycle run
(10k sweeps, prob 0.02): the quenched glass at cimp 3.0 PREFERS
f0 ≈ 635-636, THIRTEEN BELOW crystal 648** — f3 rides the volume target
exactly; S = 1979 vs frozen control's 3310 (40% lower), n_ill ~297 vs 335;
plateau flat for 200+ cycles with sparse two-way traffic (58 contract/41
split accepts over 727k proposals ≈ 1.4e-4 acceptance). PROVISIONAL
(single chain; f0 dynamics itself may be glassy-slow at this coupling) —
certification needs multi-chain R̂ on f0 per conventions. COST NOTE: thermal
links are non-catalog so every split proposal pays the ~0.5 ms DFS; at
prob 0.1 the channel dominated runtime (321 s vs 18 s) — prob 0.02 is
still acceptance-limited and ~5× cheaper; catalog fast path only pays in
near-crystal states.

**c_imp SCAN (native A/B quenches, c15 m3, 6k sweeps, committed 53235ef):**
two-tier stuckness. Tier 1 (surplus shedding above crystal): A relaxes
fully at c_imp 1.6, strands +2 at 2.0/2.4, +4 at 3.0 — wall opens below
~2. Tier 2 (crystal → vertex-poor): B's equilibrium f0 ≈ 635-638 at EVERY
coupling 1.6-3.0, unreachable by A at any c_imp (needs contraction of
pristine regions). **KEY IDENTITY: f1 = f0 + f3 (Euler, closed 3-mfd), so
at pinned f3 the flat pin IS an f0 pin with optimum f0* = f3(6/e* − 1) ≈
643.7 — the crystal's 648 exceeds it by exactly the pin-gap debt (4.3):
the channel pays the forced-defect debt by DELETING VERTICES instead of
holding defect structure.** Equilibrium lands ~8 BELOW f0* (636 vs 644):
defect-gas/entropy correction beyond the naive pin optimum. S_B vs S_A:
1230/1733 (1.6), 1606/1969 (2.0), 1497/2828 (2.4), 1979/3310 (3.0).
Single chains, uncertified.

**Next steps:** cross-library graft compatibility matrix beyond the Laves
family (a15, sigma, z, mu, p, delta, r — expect the vertex-density/web
obstruction to appear as NON-EXISTENCE of shared L3 surfaces, degrading to
L2/L1 with quantifiable seam debt); thermal stability of a grafted seam under
the minimal 3-term action; mixed-sign pin-gap composite.
