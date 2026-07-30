# Bilocal worm: joint-balance two-region move proposer — design

*Status: design (2026-07-30). Prerequisite reading: the deg-4 worm move
(sampler.d `wormEnumerate`/`tryWormMove`), the pair census
(`deg4_pair_census.py`), the nonlocal slide, `notes/CONVENTIONS.md`.*

## 1. Goal and physical picture

A Metropolis–Hastings channel whose elementary event is a **pair of local
move composites applied at two separated regions A and B**, admitted by a
**joint** postcondition ledger, so that a completed move is the transfer of
a definite balanced unit between *here* and *there* — primarily along a
Boerdijk–Coxeter chain, so the event carries a direction. This is the
long-range generalization of the deg-4 worm: the worm is the degenerate
case where A and B are one octahedron apart and the balance is enforced
inside a single patch.

Why the split beats composing local hops:

- **Direction is intrinsic** (chain + orientation), not reconstructed
  post-hoc from chart displacement.
- **Cost is distance-independent**: with the intermediate chain legal and
  untouched, dS = dS_A + dS_B exactly (support disjointness); endpoint
  congruence means a 20-cell transfer costs the same as a 1-cell hop.
- **Barrier tunneling**: creation-at-B + annihilation-at-A is evaluated as
  one MH event at the *net* dS; the expensive intermediate (cf.
  open–walk–close) is never occupied.
- **ECMC-ready**: in the lifted dynamics the chain is the persistent
  direction and the compensating far end is the natural event handoff.

### 1.1 No conservation-law obstruction (corrected understanding)

The per-vertex identity Σ_{e∋v}(6−d(e)) = 12 is Gauss–Bonnet on the vertex
link (Σ over link vertices of (6−deg_lk) = 6χ(S²) = 12, with deg_lk(u) =
d(uv)). It is **kinematic** — satisfied by every configuration — and
forbids nothing. Deg-5 edges are legal deficit carriers, so illegal-sector
structures condense from and dissolve into the 5/6 sea freely. The reason
the *local* worm class demands per-patch balance is empirical, not
fundamental: exhaustive search found no small-depth composite that removes
a thermal deg-4 (or retracts a tip) leaving all other illegal edges
intact. The bilocal class **lifts precisely that restriction**: each end
contributes an *unbalanced* half-move; only the joint ledger balances.
The pair census's "junk" sector (unbalanced creations at +20…+45 and
downhill fusions) is exactly the half-move alphabet — junk locally,
currency bilocally.

Terminology (to be added to CONVENTIONS.md): q_R(v) is **curvature**
(never "charge"); "the 12" is the link-GB **deficit-sum constraint**;
n₆ = Z − 12 is its corollary on legal vertices (web valence); "charge" is
reserved for dynamical, empirically conserved quantities (Q_c, rung Q).

## 2. Move class (v1)

A **half-move** at region X = a composite of 1 or 2 elementary bistellar
moves (2→3 / 3→2, any mix — note: *not* restricted to one-of-each as in
the local worm) with support inside a radius-R_BALL patch around a seed,
recorded as:

- net facet diff (canonical, cancellation-aware) — for end-state keys;
- **net illegal ledger** Δ_X: the map {edge → (d_before, d_after)}
  restricted to transitions where either side is illegal;
- **type** τ_X: Δ_X with positions erased — the signed multiset of
  illegal transitions, e.g. τ = {4→0 ×1, 0→3 ×1} ("absorbs a deg-4,
  emits a deg-3");
- dS_X (exact, from the speculative + potential lockstep), support
  vertex set.

A **bilocal candidate** = (half-move at A, half-move at B, pairing
witness) admitted iff:

1. **Joint balance**: across A ∪ B, #gone-4 = #new-4 ≤ K (K = 2 in v1),
   #gone-3 = #new-3 ≤ 1 *per end*; every other edge transition in either
   patch is legal→legal. (v1 enforces exact species balance → the move is
   content-neutral transport and n_ill, n₄, n₃ are preserved, so the
   site-selection normalizations cancel in Hastings, as 1/n₄ did for the
   worm. An "approximate balance" relaxation — transport + small reaction
   in one event — is a v2 dial; it forfeits the cancellation and must
   carry the count factors explicitly.)
2. **Genuine bilocality**: patch supports vertex-disjoint (guarantees
   exact dS additivity; verified anyway by lockstep on the applied
   composite). Short-range pairs where supports collide are the local
   worm's business — the channel rejects them.
3. **Pairing witness valid**: the pairing kernel (§3) connects A and B in
   both directions in the pre- and post-state (reverse reachability).

The half-move alphabet at each end is found by **exhaustive
generate-and-test**, exactly as in the worm: gather elementary moves
touching the seed (stage 1) and the disturbed region (stage 2), apply for
real, read the net ledger, roll back. No pattern matching; postconditions
only. Depth ≤ 2 per end in v1 (the census showed depth-2 already contains
creation, absorption, and fragment-shift types); depth ≤ 3 is a
measured-cost extension.

### 2.1 Matching without the n_A × n_B blowup

Half-moves are indexed by type τ. Joint balance is a constraint on
(τ_A, τ_B) alone, so admissible pairs are found by hashing B's half-moves
by the **complement type** needed to balance A's: O(n_A + n_B) with a
per-(τ, τ̄) verification of conditions 2–3. Typical n per end is ~10²
(census), so this is cheap even before caps.

## 3. Pairing kernels (modular — answers "ends of BC spiral, or any two regions?")

The region-pairing rule is a pluggable kernel. MH correctness requires
only that the kernel's selection probability is computable in both
directions; the chain is not needed for correctness. Build three:

- **P1 — chain-paired (primary).** Seed A = an illegal edge drawn
  uniformly (n_ill preserved by the class → factor cancels); a chain slot
  (bounded slot count, as SLIDE_SLOTS) picks the BC helix and orientation
  through A's dock; `chain_walk` (already in D, ~1 ms/site) follows the
  helix deterministically until the stopping rule fires → B. v1 stopping
  rule: *first tet whose edges include an illegal edge* (defect-to-defect
  transport). The geometry guarantees termination: every helix leaving a
  complex runs through legal crystal until it meets another complex (or
  wraps the torus — walk capped at L_max, wrap = no candidate).
  **Reversibility**: the walk's interior is untouched by the move (patch
  disjointness + walk-interior-legality are admission conditions), so the
  reverse walk from B's post-state dock along the reversed slot must
  reach A; this bidirectional consistency is *checked at admission time*
  (walk both ways, require agreement), which makes the kernel's reverse
  probability well-defined by construction rather than by hope.
- **P2 — uniform pair (control + validation).** A, B = uniform unordered
  pair of illegal edges with disjoint patches. Trivially symmetric.
  Physics use: comparing P1 vs P2 acceptance and dS spectra *measures
  whether the chain matters* — the rung-ladder result (slide dS = c·ΔQ,
  free ⟺ same rung, per-rung crystal-spanning webs) predicts P1 finds
  exactly-compensating pairs at a much higher rate on same-rung
  connections. If P1 ≈ P2, the chain is irrelevant to statics and its
  role is purely the lift's direction label — either outcome is a result.
- **P3 — defect-to-vacuum (later).** B = a pristine chain site at a
  walk distance drawn from a symmetric computable distribution
  (geometric in steps). Pairs an absorption at A with a *creation into
  empty crystal* at B — pair-creation/annihilation transport, the
  generalization of the nonlocal slide beyond the (3,4,4) species.
  Deferred until P1/P2 physics is understood.

## 4. Hastings

Forward weight of a candidate M with end state s′:

$$
q_f(M) \;=\; \sum_{\text{paths} \to s'} P(\text{seed}) \, P(\text{slot}) \, \frac{k(\text{pair})}{N_{\text{pairs}}(A,B)}
$$

summed over all (seed, slot, ordering) generation paths that produce the
same end state — the worm's anchor-sum, promoted to site pairs. All terms
are enumerable at proposal time: seeds are the illegal edges of the
half-moves' ledgers (bounded), slots are bounded, and N_pairs comes from
the matching table. Reverse weight q_r identically in s′ (the class is
inverse-closed by construction: each half-move's inverse is a half-move
with negated ledger, and the pairing witness survives because the chain
interior is untouched). Acceptance:

$$
\alpha \;=\; \min\!\left(1,\; e^{-\Delta S}\,\frac{q_r}{q_f}\right),
\qquad \Delta S = \Delta S_A + \Delta S_B .
$$

v1 proposes uniformly over admissible pairs and lets Metropolis find the
low-dS ones. **Low-dS targeting** without breaking balance is a v1.5
option: propose pair i with probability ∝ e^{−ΔS_i/2} — legal because
both directions enumerate the full pair set anyway, so the weight enters
q_f and q_r exactly. Adopt only if uniform acceptance is poor.

## 5. Implementation phases

**Phase 0 — measurements (before any move code; each is hours, not days):**

- *Chain-reach census* (`chain_reach_census.py`): for every illegal edge
  in the reference snapshots, walk all slots both ways; record distance
  to the first illegal edge, rung class Q of the traversed chain, whether
  the far dock admits stage-1 moves, catalyst presence at both ends.
  Deliverables: distance distribution, same-rung fraction, dockable
  fraction. Decides P1's expected candidate density and typical range.
- *Half-move type census* (`halfmove_census.py`): run the per-end
  enumerator (no joint check) at tips/monomers of the reference
  snapshots; histogram types τ with dS spectra. Deliverable: the type
  algebra table — which complementary (τ, τ̄) pairs exist, at what
  combined dS. (Reuses `deg4_pair_census.py` machinery; much of this is
  reprocessing.)

Go/no-go: Phase 1 proceeds if complementary type pairs with
ΔS_A + ΔS_B ≲ 1 exist at ≥ O(10⁻²) per site pair.

**Phase 1 — Python oracle (`worm_bilocal.py`), validation parity with the
worm oracle:**

1. Live-based half-move enumerator (depth ≤ 2, net ledgers, types,
   supports).
2. Pairing kernels P1/P2; bidirectional walk check.
3. Joint class check; exact dS lockstep; apply/undo with facet-set
   restoration tests.
4. **Closure test**: inverse-closure over a snapshot's admissible
   candidates (dS antisymmetric, k_f/k_r consistent) — the worm's 69/69
   pattern.
5. **MH integration**: composite kernel (thermal + bilocal) against the
   live sampler objective to ≲1e-13; drift audit.
6. Labs: (a) pristine two-knot teleport — plant two (3,4,4) knots on one
   chain, verify the bilocal move reproduces the nonlocal slide's
   physics as a special case; (b) thermal snapshot: acceptance, dS
   spectra, accepted-distance histogram; (c) P1 vs P2 head-to-head (the
   does-the-chain-matter measurement).

**Phase 2 — D port (only after oracle physics looks worthwhile):**

- `bilocalEnumerate` (shares the worm's gather machinery; typed
  half-moves; type-hash matching), `tryBilocalMove` (joint anchor-sum
  Hastings), `BilocalConfig` channel in mcmcStep.
- Lessons already paid for, applied from day one: allocation-free hot
  path (fixed scratch, no per-call AAs — the GC/rt_init episode);
  event-log mirroring of every committed elementary move in applied
  order (the audit-divergence episode); gate off under cocycle and
  logSixFlips; reconcile deg3/deg4 sets over both patches at commit.
- Crossval vs oracle: end-state superset + dS multiset (the worm's
  standing pattern, `worm_crossval2.py` template). Rollback-exactness
  unittest.

**Phase 3 — physics:**

- Ensemble arm: thermal + worm + bilocal vs the current worms-on arm
  (same starts/seeds protocol, 8 chains, over-dispersed starts); compare
  equilibration via the split-R̂_q/ESS gate.
- Transfer statistics: accepted-move distance histogram, direction
  persistence, rung dependence of transfer rate between complex pairs —
  the first direct measurement of the chain-mediated defect–defect
  coupling (momentum-sector phenomenology).

## 6. Risks and mitigations

- **Reverse-walk subtleties** (kernel P1): mitigated by making
  bidirectional walk agreement an *admission condition* checked on both
  pre- and post-states, and by keeping the walk interior strictly outside
  both patches. Unit-tested with adversarial cases (walk grazing a patch,
  wrap-around, chain through a third complex).
- **Combinatorial growth** of depth-2-per-end alphabets: type-indexed
  matching keeps it linear; caps with logged truncation (no silent
  coverage loss); Phase-0 census may justify restricting to the observed
  useful types.
- **Low acceptance** (compensation rarer than hoped): Phase 0 measures
  this before code is written; v1.5 dS-weighted proposal is the fallback;
  same-rung restriction of P1 is a second fallback (the free-slide
  network says same-rung is where exact compensation lives).
- **n_ill drift under approximate balance** (v2 only): keep exact balance
  until the count-factor bookkeeping is explicitly built and tested.
- **Short-pair degeneracy**: pairs whose patches touch are rejected
  (condition 2); the local worm channel owns that regime, so the two
  channels partition cleanly by range.

## 7. Relation to existing machinery (reuse map)

| Need | Existing piece |
|---|---|
| chain following | `chain_walk` (D, capi + `fpkmc.py`), ~1 ms/site |
| per-end move gathering | `wormEnumerate` stage-1/2 gather + rollback |
| exact dS | speculative + potential lockstep (worm/slide pattern) |
| net-ledger keys | worm's canonical net facet diff + `canonical_key` |
| half-move vocabulary priors | `deg4_pair_census.py` results |
| special-case ground truth | nonlocal slide (= bilocal move for (3,4,4)) |
| convergence gate | `convergence.py` quantized_split_rhat / ESS |
