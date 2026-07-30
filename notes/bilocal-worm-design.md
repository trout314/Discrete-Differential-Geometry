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
- **Dressed generators exist and can be favorable.** Targeted edge
  removal: 32/32 across d₀ = 3..6 on the quenched gas (2–4 s each,
  crystallographically quantized costs; deg-3 removal DOWNHILL at
  −5.9). First effective VERTEX removal: Z = 10 collapsed to 4 in 12
  moves (+23.4 = the whole barrier), 4→1 refund −35.0, NET dS = −11.6 —
  a composite channel would accept it outright. The f₀-frozen direction
  is crossable in one dressed proposal. Machinery:
  `dressed_generators.py` (deterministic engine, commit a8a1162).

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

### 2.4 Reversibility ledger for DEEP (dressed) generators

Shallow composites (the worm) use enumerate-then-draw: the constructor
factor is 1/n(x) and multiplicity is the anchor-sum. At depth ~13,
enumeration is impossible and a randomized search has an intractable
path-sum. The exact scheme for deep generators counts FOUR factors and
adds one check:

1. **Menu factor** p_g — and the inverse generator p_ḡ must be in the
   menu (a removal channel without its insertion partner is
   irreversible by construction).
2. **Seed factor** 1/|C_g(x)| forward, 1/|C_ḡ(x′)| reverse. These do
   NOT cancel for f-changing moves (the move changes the counts). Any
   state-computable seed set is valid — failed seeds are rejected
   proposals; C's definition is a rate knob, never a correctness knob.
3. **Constructor factor** — probability the search emits THIS composite
   given (x, seed). The determinism contract (canonical orderings and
   tie-breaks; never dict iteration order, which is
   insertion-history-dependent and breaks balance invisibly) collapses
   this to an indicator: 1 for the unique output, else 0. Bit-identical
   rerun tests are the certificate.
4. **Multiplicity** k_f = #{seeds in C_g(x) whose output is exactly this
   transition}, k_r likewise at x′. Vertex removal: k = 1 structurally
   (different seeds delete different labels ⇒ different end states).
   Edge removal: collateral collisions possible; defense is bounded —
   only seeds inside the realized support can collide, so run the
   constructor for each (tens) and count.

**The reverse-check** (support symmetry): forward and reverse
constructors are independent searches, so nothing guarantees the
inverse constructor at x′ reproduces x. Run it as part of the FORWARD
proposal; if it does not return the exact inverse composite, reject
unconditionally. This prunes the kernel to a symmetric support (one
extra construction per proposal), resolves the insertion-goal ambiguity
(balance only credits mutual-inverse pairs — the check IS the usable
class), and creates healthy design pressure toward mirror-image
constructors. Acceptance:

    alpha = min(1, exp(-dS) * [p_ḡ k_r / |C_ḡ(x′)|]
                            / [p_g  k_f / |C_g(x)|])   given check passes.

NOT counted: search effort (deterministic constructor ⇒ irrelevant);
other generators reaching the same x′ (each generator pair is its own
channel; fixed mixtures of balanced channels are balanced);
intermediate states (scaffolding — their dS cancels telescopically in
the lockstep total). Label bookkeeping caveat: vertex labels are
recycled by the capi; the offline 1↔4 path must assign labels
deterministically too, since labels are part of the chain's state.

**MEASURED (Phase-1 lab, 2026-07-30,
scratchpad/phase1_reversibility_lab.py):** naive independent mirror
goals give reverse-check pass rates of **6/6 for the self-mirror
sector** (d₀ = 3 edges: removal = inverse of creation structurally;
dS exactly antisymmetric ±5.93) and **0/24 for everything deeper**
(d₀ ≥ 4 edges: the mirror "creates e" the cheap way — one 2→3 at
degree 3 in the mangled environment — never the original-degree
rebuild; vertex tier: legalize gets stuck rebuilding Z = 4 → 10).
Conclusion: deterministic + reverse-check survives as the CHEAP FAST
PATH for structurally self-mirror sectors only. The general deep
sector needs **scheme B — configurational-bias (Rosenbluth)
constructors**: build the composite stepwise, sampling the (symmetric)
alphabet with computable stagewise weights; the reverse probability of
the recorded inverse path is computed by RETRACING it under the
inverse goal's bias. Alphabet symmetry (every 2→3 has its 3→2)
guarantees nonzero reverse support, eliminating the exact-match
brittleness; the Hastings ratio carries the two Rosenbluth products.
Forward-direction bonuses from the same lab: the canonical engine cut
the vertex-removal barrier to +16.40 (net −22.82), and costs are
CRYSTAL-ORBIT-QUANTIZED (two same-orbit vertices: bit-identical
trajectories) — generator weight tables in ordered states are tiny.

## 3. The f-changing sector, and where it can actually be tested

Pairs with joint Δf ≠ 0 collect the global-sector price/bonus at the
pair level while placing each ledger where it is locally cheapest —
manufacturing the coincidences that sequential dynamics must wait for.

### 3.0 Menu completeness (do not curate by ensemble)

Elementary Pachner moves have HARD availability preconditions (3→2
needs a deg-3 edge; 4→1 needs a Z=4 vertex), so ensembles can be
**proposal-starved** in some f-directions — categorically different
from Metropolis suppression and invisible in acceptance statistics.
Our usual conditions starve exactly the Δf₀ directions ("zero volume
moves ever" is a property of the gating, not the physics), which made
**f₀ a hidden ensemble parameter of every run to date** — all sampling
has been implicitly microcanonical in vertex count.

Design rule: the move menu must contain an *unconditionally available*
dressed generator for each simple Δf consistent with χ = 0, regardless
of what fires under any particular ensemble. On T³ (f₁ = f₀ + f₃,
f₂ = 2f₃) the lattice is ℤ² in (Δf₀, Δf₃); the thermal bath supplies
(0, ±1) abundantly, so ONLY Δf₀ ≠ 0 needs engineering — one dressed
generator per sign, composition with thermal moves reaches the rest.

Taxonomy: (1) **dressed generators** — single-region deep composites,
the unconditional f-movers (and the building blocks generally);
(2) **bilocal pairings** — complementary halves give the f-NEUTRAL
transport sector; the pairing's role for f-movers is cost-sharing
(finance an uphill half with a downhill half elsewhere, e.g. deg-3
removals at −5.9 in quenched states); (3) elementary moves — the fast
bath, gating accepted since (1) covers what it starves.

NOTE (corrected 2026-07-30, user): enabling f₀ generators does NOT
change the ensemble. The full Pachner alphabet is ergodic across f₀
(1→4 is always proposable at finite dS, and the pins already legislate
f₀ = f₁ − f₃ softly), so π has always been the action's distribution
with its f₀ marginal included — the sequential chain just cannot reach
it (spontaneous realization of the 13-move collapse path costs
~e^(−barrier) ≈ e^(−23) per attempt; the dressed proposal converts
that rare fluctuation into one O(1)-acceptance move). Generators are
pure equilibration technology: mixing changes, the measure does not.
The operational caveat that survives: previously certified states were
certified as stationary CONDITIONAL on their kinetically frozen f₀
slice — correct metastable-state measurements that may relax to a
different f₀ (and different observable values) once the menu is
completed.

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
| dressed constructors + goals | `dressed_generators.py` (a8a1162) |
| edge surgery + cost survey | `edge_removal.py` (de7cc4b) |
| convergence gate | `convergence.py` R̂_q / ESS |
| f₀ benchmark states | quench_* recorder streams (scratchpad) |
