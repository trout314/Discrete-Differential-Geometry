# Bilocal moves: general joint-proposal design (v2)

*Status: design v2 (2026-07-30). v1 (same file, git history) restricted the
class to content-neutral deg-4 transport; v2 makes the GENERAL class
primary after two experimental anchors landed (§1.2). Prerequisites: the
deg-4 worm (sampler.d), the pair census (`deg4_pair_census.py`), the
channel-bracket event annotations (EVT_CHANNEL_*, commit 1467b24),
`notes/CONVENTIONS.md`.*

## 1. Idea and status

A Metropolis–Hastings channel whose elementary event is a **pair of local
move composites applied at two separated regions A and B**, priced and
admitted **jointly**. Because the patches are disjoint, dS = dS_A + dS_B
exactly; because each end enumerates its local moves exhaustively, the
proposal can be dS-targeted with *exact* Hastings bookkeeping. A completed
move transfers something definite between *here* and *there* — with the
pairing kernel (§4) supplying direction when it follows a BC chain.

### 1.1 v2 shift: conservation is emergent, not imposed

v1 defined the class by a conservation ledger (balanced deg-4 exchange,
≤1 catalyst, all else preserved). v2 drops that: the admission criterion
is purely structural (§2), and the proposal's dS weighting supplies the
selectivity that the conservation rules used to. The old classes survive
as **sectors** of the general class — the deg-4 worm is the f-neutral
deg-4-balanced sector, the nonlocal slide is the (3,4,4) sector,
defect↔vacuum pair creation is another — but which sectors carry
acceptance mass is now *measured* (Phase 0 census, §6), not designed.
The pair census already showed the half-move vocabulary is richer than
any hand-written class; the "junk" sector (unbalanced creations and
absorptions) is exactly the currency that pairs across distance.

### 1.2 Experimental anchors (both 2026-07-30)

- **Transport is orthogonal to the global sector — measured.** Quench
  experiment (m²-only C15 gas, ±5-quanta target quench, worms on/off,
  same seeds): ~2,700 accepted worm hops produced ZERO change in the
  pin-gap or density relaxation, both quench directions. f-neutral
  transport cannot touch a global-sector perturbation; a channel that
  should is precisely the f-changing pair sector of this design (§3).
- **Composite-move instrumentation exists.** The channel-bracket event
  records (BEGIN/EDGE/END with identity sets and per-move dS) let the
  defect tracker treat any composite as one atomic transaction — worm
  validation: transport events == accept counter exactly, 44% of
  apparent births were relabeling artifacts, median transport |dS| = 0.
  Any new bilocal channel MUST emit these brackets from day one.

## 2. The general class

A **half-move** at region X: a composite of 1..D elementary bistellar
moves (D = 2 in v1 implementation) supported in a patch around a seed,
recorded with its net illegal ledger Δ_X, its net f-vector change Δf_X,
its exact dS_X (speculative + potential lockstep), its support vertex
set, and its position-erased **type** τ_X.

A **bilocal candidate** = (half at A, half at B, pairing witness),
admitted iff:

1. both halves are individually valid move sequences;
2. **patch disjointness**: the two support vertex sets (including every
   simplex whose degree/counter changes) are disjoint — this is what
   makes dS_A + dS_B exact and the two enumerations independent. Keep
   this even in the most general version: it is the definition of
   *bi*-local (pairs whose patches collide belong to single-patch
   channels);
3. the pairing witness validates in both directions (§4).

No balance conditions. Species-changing, f-changing, and asymmetric
pairs are all admissible; Metropolis prices them.

### 2.1 dS-targeted proposal (core in v2, was v1.5 option)

Each end's enumeration yields every half-move WITH its dS. Propose pair
(i, j) with probability ∝ exp(−(dS_i + dS_j)/2) (or enumerate only
pairs below a logged dS cutoff — no silent caps). Exactness: both the
forward and reverse states enumerate the same way, so the weights enter
q_f and q_r exactly. Junk pairs are not rejected; they are barely
proposed. Without conservation rules this weighting is the class's
selectivity mechanism, not an optimization.

### 2.2 Pricing factorizes along the action's structure

ΔS_glob (volume pin, f₁ pin — all f-vector terms, including the affine
part of per-edge terms via Σd = 6f₃) is a function of the pair's joint
type (Δf_A + Δf_B) alone: one table lookup, position-blind. The local
nonlinear sector (washboard) is the per-end dS from enumeration. So
matching is type-indexed (hash B's halves by complement/companion type,
O(n_A + n_B)) and only the local scores need per-instance computation.

### 2.3 Hastings

q_f(M) = Σ over generation paths of P(sites) · w(pair)/W(A,B), with
w the dS weight and W the normalization over admissible pairs at
(A, B); reverse identically in the end state. With unbalanced halves,
the site-selection normalizations (1/n_ill, 1/n₄ …) NO LONGER cancel —
carry them explicitly (they are known integers on both sides). Each
half's inverse is a half at the same site, so inverse-closure is
structural; reverse enumeration runs at the post-state ends (the worm's
anchor-sum pattern, generalized).

## 3. The f-changing sector, and where it can actually be tested

Pairs with joint Δf ≠ 0 collect the global-sector price/bonus at the
pair level while placing each ledger where it is locally cheapest —
manufacturing the coincidences that sequential dynamics must wait for.

**Testbed requirement (learned the hard way):** before crediting any
mover with accelerating a global mode, MEASURE that the mode is
kinetically limited. The m²-gas quench showed its f₃ sector relaxes in
~10 sweeps (nothing to accelerate) — but the same experiment exposed a
genuinely frozen mode:

**The f₀ benchmark.** Under a down-quench of the f₁ pin, the state can
satisfy both global pins only by shedding vertices (f₁ = f₀ + f₃ at
χ = 0; e.g. f₀ 1536 → ~1522 in the C15 m4 down5q quench, which would
take BOTH pins to ~0 instead of the observed compromise at gap ≈ +10,
cost ≈ +100). The fast sector equilibrates exactly to the
constrained-at-fixed-f₀ optimum (predicted f₃ = 8722, gap +10.5;
observed 8722–8726, +10) — but f₀ is quasi-conserved because 4→1 needs
a coordination-4 vertex, which this landscape essentially never forms:
the volume-move channel measures ZERO accepts across all campaigns. A
composite/bilocal sequence that performs an effective vertex removal
(drive a vertex's star toward Z = 4, then collapse; possibly paired
with a compensating half elsewhere) has a crisp, quantitative success
criterion: close the +10 gap, f₀ 1536 → 1522, in a state where we KNOW
the target configurations exist. This is the first concrete target for
the f-sector.

## 4. Pairing kernels (unchanged from v1, one addition)

- **P1 chain-paired** (primary): seed + slot → deterministic
  `chain_walk` to a stopping site; bidirectional walk agreement is an
  admission condition. With unbalanced halves the stopping rule must
  tolerate ends whose defect content the move itself changes — the
  subtlest reversibility point in v2; needs a dedicated closure test.
- **P2 uniform pair** (control): measures whether the chain matters
  (rung-congruence prediction: it does, for exact compensation).
- **P3 defect↔vacuum**: subsumed by the general class (a creation half
  pairs with anything); keep as a kernel only if vacuum sites need
  their own selection distribution.

## 5. Implementation phases (v2)

**Phase 0 — measurements (unchanged, now richer):** chain-reach census;
half-move type census WITH full dS spectra and Δf per type at both
ends. Output IS the proposal weight table. Go/no-go: complementary
low-|dS| pairs at ≥ O(10⁻²) per site pair.

**Phase 1 — Python oracle:** general enumerator (net ledgers, Δf, dS);
P1/P2 kernels; dS-weighted joint proposal; explicit count factors;
closure test incl. content-changed ends; MH lockstep ≤ 1e-13; labs:
(a) two-knot teleport reproduces the nonlocal slide; (b) worm sector
emerges as the dominant f-neutral acceptance mass (consistency with
v1); (c) P1 vs P2; (d) the f₀ benchmark (§3).

**Phase 2 — D port:** allocation-free from day one (GC episode);
channel brackets from day one (EVT_CHANNEL_BILOCAL: add channel id 4);
gated off under cocycle and logSixFlips like the other composite
channels; crossval vs oracle (end-state superset + dS multiset).

**Phase 3 — physics:** transfer statistics (distance/direction/rung
dependence); f₀-benchmark quench races; λ-gas arm (creation-dominated
drift is the other known kinetically-limited mode).

## 6. Risks

- **Acceptance dilution** without conservation rules: mitigated by the
  dS-weighted proposal; Phase 0 measures the achievable enriched
  acceptance before any move code is written.
- **Count-factor bookkeeping** (no cancellations): explicit integers,
  but every one omitted is a silent detailed-balance bug — closure
  tests must cover species-changing pairs specifically.
- **Reverse-walk validity when ends change content** (§4).
- **Proposal cost**: two enumerations + reverse pass, ~worm-scale
  (tens of ms); dS cutoff bounds the pair table.

## 7. Reuse map

| Need | Existing piece |
|---|---|
| chain following | `chain_walk` (D, capi + `fpkmc.py`) |
| per-end gather + rollback | `wormEnumerate` stages |
| exact dS | speculative + potential lockstep |
| atomic tracking + transport events | EVT_CHANNEL_* brackets (1467b24) |
| half-move vocabulary priors | `deg4_pair_census.py` |
| ground truth special cases | nonlocal slide; deg-4 worm |
| convergence gate | `convergence.py` R̂_q / ESS |
| f₀ benchmark states | quench_* recorder streams (scratchpad) |
