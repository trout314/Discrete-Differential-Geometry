# First-passage kinetic Monte Carlo for the defect gas — design

Status: DESIGN (2026-07-25). Nothing here is implemented yet.
Companion measurements: no-halo theorem (defect_halo), washboard + S-matrix
(knot_collider / crossing_collider / knot_smatrix_sweep), contact algebra
(compound_lab, harvest_collider), caging (defect_travel).

## 0. Goals

Three samplers, one infrastructure, in order of construction:

  HB    segment heat-bath — exact equilibrium accelerator for the FULL
        model (thermal + flicker + strings), no time variable. Unlocks:
        true equilibrium g(r) of complexes, association/dissociation
        equilibrium of the compound family, species distribution.
  FP    first-passage KMC — exact kinetics of the slide channel in a
        frozen background (v1), thermally dressed (v2). Unlocks: encounter
        rates, association/dissociation kinetics, reaction-network
        trajectories ("worldlines interacting at vertices").
  EC    lifted event-chain — non-reversible ballistic variant; the
        momentum-sector candidate. Deferred; only its infrastructure
        requirements are recorded here so nothing built now blocks it.

Performance stance (user requirement): any component with even a small
chance of needing performance goes on the D side behind the C API. Python
keeps orchestration, linear algebra on tiny matrices, and analysis.

The same D-side additions also fix the S-matrix sweep's cost problem
(currently ~8 s/collision in Python); the sweep is a consumer of this
infrastructure, not a separate investment.

## 1. Physical foundations (what licenses what)

F1  NO-HALO (measured, exact): a defect complex alters NOTHING outside its
    own vertices — every legal vertex beyond tet-adjacency of a defect
    vertex has exactly its reference-crystal class and charge. Hence:
    (a) the action of a multi-defect configuration is additive until
        complexes are tet-adjacent;
    (b) a knot's slide dynamics inside a clear chain segment is an
        autonomous 1D process with site energies given by the PRISTINE
        washboard;
    (c) two configurations identical within a region behave identically
        under any move whose frame lies in that region.

F2  WASHBOARD (measured): the single-knot site energy S1(j) along a BC
    chain takes values in a small set of translation classes (17.1–51.0 in
    λ=1 units on R). It is cacheable by the translation class of the
    site's local ball.

F3  SLIDE ALPHABET IS INVERSE-CLOSED (verified in worm_slide): every clean
    slide j→j+4 has a clean inverse j+4→j. Consequence: the number of
    clean forward slots out of site j equals the number of clean backward
    slots out of site j+4 (call both n(j, j+4)). This is what makes the
    1D chain satisfy detailed balance with the Boltzmann weights e^{−λS1}
    (see §3.1 invariant I1 — VERIFY, do not assume).

F4  CONTACT ALGEBRA (measured): interactions are contact-only —
    vertex-disjoint touch is free (V = 0), vertex-sharing carries a
    quantized penalty (+4.4 … +16 observed so far), all products unfuse by
    one slide. The S-matrix table predicts dock outcomes but the samplers
    below never rely on it for correctness — contact is always resolved by
    explicit moves.

## 2. State space and channels

State: the triangulation T (fixed topology T³), action

    S = c_N (N3 − N*)² + λ [ (e*/6) Σ_e (deg e − e*)² + zleg Σ_v U(n6)
                             + cimp Σ_v m² ]

Channels of the mixed sampler:

  C-th   ordinary thermal Pachner channel (existing mcmcStep).
  C-sl   per-step Metropolis slide channel (existing; k_f = k_r = 1).
  C-hb   NEW segment heat-bath jump (HB mode).
  C-fp   NEW first-passage flight (FP mode; replaces C-sl for tracked
         knots).

HB correctness note: a heat-bath jump is a conditional resampling of one
knot's position given everything else, using the CURRENT segment; it is a
valid MH move on the full state space no matter how the background evolves
between proposals. Therefore C-th + C-sl + C-hb is exact for the full
model with no further conditions. This is the low-risk first deliverable.

FP correctness is subtler (frozen background exact; thermal dressing
approximate) — see §5.

## 3. The 1D chain: exact formulas

### 3.1 Rates

Site j (spacing 4 chain steps). Per attempted sampler move, the
probability of proposing THIS knot's forward slide is

    p(j→j+4) = p_slide · (3 / N3) · (1/6) · (n(j, j+4) / 12)

(3 facets contain a deg-3 chord; 6 vertex pairs per facet; 12 slots, of
which n are the clean ones toward j+4 — read the exact constants OFF
sampler.d at implementation time and treat the formula above as a
hypothesis to validate, not a spec). Acceptance is Metropolis on
ΔS = λ(S1(j+4) − S1(j)) plus, near contact, interaction terms — which is
why segments end before contact. Chain rates:

    q(j→j±4) = ν · n(j, j±4) · min(1, e^{−λ(S1(j±4) − S1(j))})

with ν the constant prefactor.

INVARIANT I1 (verify computationally over whole orbits before anything
else is built): n(j, j+4) as counted forward equals n as counted backward
from j+4. If I1 holds, detailed balance w.r.t. π(j) ∝ e^{−λ S1(j)} holds
on the chain, and everything below is exact.

### 3.2 Splitting probabilities (which end of the segment)

Segment sites j = 0..n (0 and n are the docking sites; interior clear).
Standard birth–death exit problem. With Metropolis rates and I1, the
ratio q(i→i−1)/q(i→i+1) telescopes to pure Boltzmann factors. Define

    φ_k = exp( +λ S̄_k ),   S̄_k = the "edge potential": for the edge
    (k, k+1), S̄_k = max(S1(k), S1(k+1))  (Metropolis: the uphill side)

then the exit-right probability from start site s is

    P_R(s) = [ Σ_{k=0}^{s−1} φ_k n_k^{−1} ] / [ Σ_{k=0}^{n−1} φ_k n_k^{−1} ]

with n_k = n(k, k+1). (Derive once more carefully in the implementation
notebook — the φ form depends on the Metropolis min(1,·) convention; the
generic correct statement is P_R(s) = Σ_{k<s} r_k / Σ_{k<n} r_k with
r_k = Π_{i≤k} q(i→i−1)/q(i→i+1), and the φ expression is its
simplification. TEST against direct simulation regardless.)

### 3.3 First-passage time

The generator Q is tridiagonal (≤ ~200 sites/segment). Exact options:

  mean/variance   closed-form sums (standard birth–death FPT formulas) —
                  enough for time bookkeeping at leading order;
  full FPT law    eigendecomposition of Q restricted to the interior
                  (numpy, ~200×200, once per segment CLASS, cached), or
                  uniformization for direct exact sampling.
  interruption    if an external clock (flicker, §5) fires at t < T_exit,
                  the knot's position is drawn from the
                  conditioned-on-no-exit propagator exp(Q_int t) row —
                  same eigendecomposition, no new machinery.

All linear algebra stays in Python/numpy: matrices are tiny and per-class
cached; this is NOT a performance path.

### 3.4 Time calibration

FP time is in attempted-move units (the same clock as sweeps: 1 sweep =
N3 attempts), so FP trajectories and brute-force trajectories are
directly comparable. All rates inherit ν from §3.1; ν errors cancel in
equilibrium quantities and scale all kinetics uniformly — still validate
ν explicitly (V2 test below) so absolute rates are trustworthy.

## 4. Segments

For a knot with chord at chain position j0 on chain c:

  blocked(k) = ANY of:
    (a) the knot-at-k 5-vertex complex would be tet-adjacent to a vertex
        of another defect complex (contact);
    (b) a slide frame into or out of k is template-invalid in the current
        background;
    (c) k is beyond a chain wrap into already-scanned territory (cyclic
        orbits: cap segment length at L_orbit).

  segment(j0) = maximal clear run around j0; docking sites = first
  blocked site in each direction.

Flicker policy (both modes): segments are computed against the долгоживущие
complexes ONLY (age or size threshold, as in species_interactions);
flicker is handled by the §5 clock in FP mode and needs no handling in HB
mode (HB is exact regardless of what the segment is, as long as the SAME
segment rule is applied when proposing and when reversing — the segment
depends only on the background, not on the knot's position in it, so the
proposal is symmetric by construction; a knot INSIDE another complex's
contact range has no segment and simply gets no HB proposal).

## 5. FP mode: background hierarchy (be honest about exactness)

  v1 FROZEN     background fixed (no thermal channel). FP is EXACT for
                the slide-channel kinetics in a quenched background. This
                already answers: encounter rates on the chain network,
                association/dissociation of specific pairs, washboard-
                limited transport statistics. This is the deliverable
                that replaces the abandoned 40-day sweep economics.
  v2 DRESSED    thermal channel represented by: (i) a Poisson clock for
                knot–flicker encounters with rate per site calibrated
                from measured absorption statistics (species_interactions:
                120 events / 2500 sweeps / known geometry), (ii) measured
                effective slide acceptance in thermal background
                (worm_slide thermal test) replacing the clean washboard
                acceptances. Documented APPROXIMATION; label all outputs
                PROVISIONAL-DRESSED.
  v3 FULL       true concurrent thermal field — protective-domain FPKMC.
                Not planned; revisit only if v2 vs brute-force disagrees
                where it matters.

## 6. D-core additions (C API)

Read-only surveys (no sampler-state mutation; safe, high value):

  ddg_manifold_chain_walk(h, window[4], n, out_verts)
      sliding-window BC walk, n steps (replaces Python bc_orbit /
      walk_stretch inner loops; also used by the sweep).
  ddg_manifold_site_survey(h, chain*, len, estar, zleg, cimp, out)
      per chain site k: (dS_create(k), n_clean(k, fwd), n_clean(k, bwd),
      dS_slide_fwd(k), dS_slide_bwd(k), template_valid flags).
      One call per segment; this is the washboard + I1 data in one shot.
      Implementation reuses the sampler's speculative-delta machinery.
  ddg_manifold_segment_scan(h, chord[2], defect_verts*, nd, max_extent,
      out_blocked, out_reason)
      walks both directions from the chord, applying §4's blocked()
      predicate with the D-side degree maps and tet adjacency.

Targeted mutations (needed for HB jumps and contact resolution while
keeping SAMPLER bookkeeping consistent — the current targeted-move path
is Manifold-level only):

  ddg_sampler_do_bistellar(h, center*, nc, cocenter*, nco)
      apply a specified bistellar move through the sampler so its
      incremental action/degree bookkeeping stays valid (the speculative
      delta code path already computes everything needed).
  (slide execution already exists: ddg_sampler_slide_at.)

Everything else — segment linear algebra, event queue, caches, drivers,
analysis — stays in Python.

Sweep synergy: chain_walk + site_survey + (existing) slide_at reduce the
S-matrix sweep's per-collision cost from ~8 s (Python edeg_dict/census
rebuilds) to milliseconds-scale D calls plus one incremental census; the
sweep should be rewritten on this API before any relaunch.

## 7. Caching and translation classes

  site class key   relkey (sorted rounded relative positions) of the ball
                   of radius 1.0 cells around the site — the provably
                   sufficient radius (no-halo + frame reach).
  washboard cache  site class → (S1, n_fwd, n_bwd). Populated lazily via
                   site_survey; expected only O(10²–10³) classes per
                   crystal (cf. 680 A-classes, 47 flicker species).
  segment class    tuple of its sites' class keys → eigendecomposition,
                   splitting tables, FPT samplers.
  invalidation     none needed in frozen mode; in mixed/HB mode the
                   caches key on PRISTINE classes (F1 guarantees validity
                   whenever the segment is clear — that is what "clear"
                   means).

## 8. Contact resolution

When a flight docks (or an HB jump lands next to its segment boundary and
the thermal/slide channels take over):

  - switch to EXPLICIT dynamics within a local window (existing channels;
    ordinary Metropolis) until the knot either re-enters a clear segment
    (new flight) or merges (worldline bookkeeping via defect_state).
  - the S-matrix table is used for PREDICTION and VALIDATION, never for
    dynamics: dock outcomes in simulation must reproduce the deterministic
    collider table — a free continuous integrity test. Log disagreements
    loudly; they mean the background near the dock was not clear.

## 9. Validation matrix (all against brute force, small boxes)

  V1  I1 slot symmetry: survey whole orbits on R m2/m3; assert
      n(j,j+4)_fwd == n(j,j+4)_bwd at every edge.
  V2  rate constant ν: brute-force slide-channel simulation of one knot
      on a clear orbit (m2), fit hop rate per attempt vs §3.1 formula.
  V3  HB equilibrium: single knot, m2, HB+thermal vs long thermal+slide:
      site-occupancy histograms must agree (χ² on translation classes);
      two-knot: pair-separation histogram (the first true equilibrium
      g(r) measurement).
  V4  FP kinetics (frozen): exit-site fractions and FPT distributions vs
      brute-force slide-only runs from identical starts (KS tests).
  V5  contact round-trip: dock → explicit resolution products vs the
      collider S-matrix table.
  V6  bookkeeping: sampler-integrated targeted moves vs dS_between on
      every audit interval (reuse the Ledger audit pattern).

## 10. Milestones

  M0  this document reviewed; formulas §3.2/3.3 re-derived and checked
      symbolically (SymPy) before coding.
  M1  D: chain_walk + site_survey (+ tests co-located, D unittest).
      Python: washboard/site-class cache; I1 verification (V1); ν (V2).
  M2  D: segment_scan + sampler_do_bistellar. Python: HB channel driver;
      V3. First physics: equilibrium two-knot g(r), association constant
      of the compound family.
  M3  FP frozen driver (event loop, FPT sampling, contact resolution);
      V4 + V5. First kinetics: encounter-rate and dissociation-rate
      tables; compare against detailed-balance predictions from HB.
  M4  FP dressed (flicker clock, dressed acceptances); labeled
      approximation.
  M5  sweep rewrite on the M1 API (the exhaustive S-matrix resumes with
      new economics).
  M6  EC lifted variant (momentum sector) — separate design note when
      reached.

## 11. Risks / open questions

  R1  I1 could fail at exotic sites (slot asymmetry) → chain detailed
      balance needs slot-count corrections in π; formulas generalize but
      MUST be derived before caching anything (this is why V1 is first).
  R2  proposal constants in §3.1 depend on sampler.d internals; wrong ν
      silently rescales all kinetics — hence V2 before any physics claim.
  R3  self-interaction / chain revisits (RESOLVED 2026-07-25, second
      revision after discussion). Three layered facts:
      (i)   No physical self-interaction: no-halo means a knot's energy
            never depends on its history; a lone complex cannot contact
            itself. Segment scans need no self-checks for single knots.
      (ii)  No double cover: the BC walk is a reversible bijection on
            directed windows (unique continuation both ways, ANY
            triangulation — user's argument), so chains never merge or
            branch; and MEASURED on pristine R m4, the face -> apex-pair
            map is INJECTIVE (115968 faces, 115968 distinct pairs, 0
            duplicates), so no chain visits the same physical chord
            twice. The 1D line model keyed by sorted chords is exact on
            this crystal. The duplicate-apex-pair configuration IS
            combinatorially possible in general (two 2-simplices whose
            edge-links are the same vertex pair) — just absent here.
            GUARD (required, M1): the infrastructure runs the injectivity
            scan on the REFERENCE crystal at startup and hard-errors on
            any duplicate. By no-halo, every clear region of every state
            is exactly the reference crystal, so reference-level
            injectivity propagates to all segments of all configurations
            — one scan covers the samplers permanently; defect-adjacent
            zones are blocked sites and need no guarantee. ROBUSTNESS
            FALLBACK (documented, not built): if a crystal ever fails the
            scan, key states by (chord, link-face) instead of bare chord
            and use the §3 sparse-graph generalization — no redesign, a
            different key.
      (iii) The one open finite question — the SLOT CENSUS (in V1): the
            slide proposal re-derives its frame from the 12-slot menu
            using only (chord + unordered link); the two true chain
            continuations always validate, and the question is whether
            any site class validates a THIRD clean frame, which would be
            a legal chain-to-chain hop through a crossing (NOT a chain
            merge — chains stay distinct; it is a property of the move
            set, not the walk). Slot validity is class-deterministic
            (no-halo), so a census over the few hundred site classes
            decides it exactly: all classes = 2 -> per-crystal proof the
            slide move is chain-confined and §3 stays tridiagonal; any
            class > 2 -> use the (already specified) sparse-graph
            generalization and report the hop channel as physics (it
            would help R4 ergodicity).
      Genuine self-contact of EXTENDED objects (a string curling onto
      itself = loop closure) is real physics handled by the contact
      algebra as a species change — a measurement target for the string
      program, not a sampler pathology.
  R4  HB ergodicity: jumps move knots along their current chain only;
      chain-changing remains with the thermal channel. Mixing claims must
      cite both channels; measure chain-change rate (cheap worldline
      observable) before claiming equilibrated pair correlations.
  R5  dressed-mode honesty: v2 is an approximation with measured inputs;
      keep it clearly labeled and validated against v1 limits.
