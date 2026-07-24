# Notation, terminology, and reporting conventions

Standards for communication in research sessions (human ↔ assistant) and for
figures/reports produced by scripts. Rules here are referenced from CLAUDE.md;
edit freely — this file is the single source of truth.

## 1. Standard notation

| symbol | meaning | units / value |
|---|---|---|
| N_k | N_k(T) = # of k-dimensional simplices in triangulation T | count |
| N, E, V | N = N₃(T) tetrahedra, E = N₁(T) edges, V = N₀(T) vertices | counts |
| N\* | target number of tetrahedra | count |
| E\* | target number of edges; set by ē\* and N\* via E\* = 6N\*/ē\* | count |
| v̄, ē | mean vertex-degree, mean edge-degree | — |
| v̄\*, ē\* | target mean vertex-degree, target mean edge-degree | — |
| ē_f, v̄_f | flat edge/vertex degree (zero angle-defect): ē_f = 2π/θ ≈ 5.1044, v̄_f ≈ 22.80 | — |
| c_N | coefficient of the objective term (N − N\*)² | — |
| c_E | coefficient of the objective term (E − E\*)² | — |
| c_VV | coefficient of the vertex-degree variance term | — |
| c_EV | coefficient of the edge-degree variance term | — |
| c_ET | coefficient of the mean quadratic edge-degree term (with floor) | — |
| c_VT | coefficient of the mean quadratic vertex-degree term (with floor) | — |
| S[T] | sampler objective function (≡ ADM Hamiltonian ≡ Euclidean action; see registers below) | — |
| m | number of crystal domains along each axis (box = m³ cells; "m4") | — |
| λ | coefficient scaling all S[T] terms except the volume pin | — |
| λ_EDQ | as λ, but for the R²-only ablation (n6 off). *Never numerically comparable to λ.* | — |
| Δē\* | a change in the target mean edge degree | — |
| ⟨X⟩ | ensemble average of X (e.g. ⟨ē⟩) | — |
| θ | regular-tetrahedron dihedral angle, arccos(1/3) | rad |
| δ(e) | angle defect at an edge, δ(e) = 2π − θ·deg(e) | rad |
| q_R(v) | curvature "charge" at vertex v, q_R(v) = ½ Σ_{e∋v} δ(e) | rad |
| ρ_R(v) | curvature charge density, ρ_R(v) = q_R(v)/deg(v) | rad |
| δq(v) | curvature deviation from mean, q_R(v) − q̄_R (state the reference — usually cell-mean) | rad |
| w(e) | link (disclination) charge of edge e, w(e) = 6 − deg(e); flat at deg 6 | — |
| n₆(v) | # of degree-6 edges incident to v | — |
| imp(v) | # incident illegal edges (deg ∉ {5,6}) at v | — |
| n_ill | **census**: # vertices with imp > 0 | count |
| N_cx, s | # complexes; complex size (vertices) | counts |
| Q_c | total δq of a complex (knot ≈ −1.19 ± 0.18) | rad |
| f_FK | FK-legal vertex fraction (`legalvert`) | % |
| f_e | legal-edge fraction (`legaledge`) | % |
| f_reg | in-registry (interior-crystalline) vertex fraction | % |
| N_gr | # crystalline grains (≥ 10 vertices) | count |
| f_G1 | fraction of vertices in the largest grain | % |
| τ | lag / elapsed time | sweeps |
| Δr²(τ) | defect MSD | cells² |
| a_v | vertex spacing (r phase: 159 vertices/cell) | 0.185 cell |
| J(τ) | illegal-set Jaccard overlap at lag τ | — |
| R̂_q | quantized split R-hat (`quantized_split_rhat`; never hand-rolled) | — |
| ESS | effective sample size | — |
| ω | integer cocycle (the T³ **chart**) | raw units, 10⁶/cell |
| φ | oriented six-flux (cycle-space projection) | — |
| D_v | six-flux divergence at v (monopole charge) | — |
| W | six-web Z³ homology winding | quanta |
| ΔP | charge-polarization increment (per move / per window) | rad·cell |
| J_Q | monopole current, dP/dt | rad·cell/sweep |

**Two curvature charges — do not conflate.** *q* (from δ(e)) is the **ambient**
hinge curvature, flat at ē_f = 2π/θ ≈ 5.104: degree-5 edges are near-neutral
(δ ≈ +0.13), degree-6 edges strongly charged (δ ≈ −1.10). *w* = 6 − deg is the
**link** curvature (2D Gauss–Bonnet on the vertex's S² link), flat at degree 6,
so degree-6 edges are the *w*-zero species and Σ_{e∋v} w = 6χ(S²) = 12 at every
vertex (legal or not). The six-web / n6 / FK-legality story is a *w* statement;
physical curvature and its transport are *q* statements. State which when you
say "charge."

## 2. Standard terminology

| term | definition | acceptable synonyms | avoid |
|---|---|---|---|
| full action | edge pin + n6 legality terms (+ volume pin, always on) | — | "the objective" |
| EDQ-only | full action with n6 off (the R² ablation) | R²-only | bare "EDQ" as an ensemble name |
| complex | connected set of imp>0 vertices | illegal complex, defect complex | "defect" for a single vertex |
| knot | complex with (3,4,4)-type edge signature — quantized full-action species | — | — |
| worm | open chain of deg-4 edges — EDQ-only species | 4-segment | — |
| monopole | edge-illegal vertex = divergence site of the six-flux | — | "line endpoint" (impossible at legal vertices — fullerene theorem) |
| halo | FK-legal but out-of-registry vertices dressing a defect | — | — |
| host | the crystalline grain containing the defect gas | — | — |
| census | the n_ill time series | population | — |
| blinker | track living ≤ 1 frame (150 sw) | — | — |
| walker | long-lived track with net displacement ≳ a_v (centroid-measured) | — | walker claims from visitation counts |
| rattler | long-lived but caged (MSD plateau ≲ a_v) | — | — |
| turnover | population replacement (Jaccard decay) | churn | "transport" unqualified |
| tracer transport | displacement of tracked complexes (MSD) | self-diffusion | conflation with collective |
| collective transport | conserved-quantity flow (ΔP, W-drift) | coherent channel | — |
| densifying | census rising while the host stays intact | — | "melting" |
| melt | registry collapse (f_reg ≪ 1) | — | — |
| crumple | ⟨e⟩ runaway — geometric collapse | crumple-melt | — |
| stationary | late-window slope consistent with 0 (block-bootstrap σ) | leveled out | claims from full-series OLS |
| certified | passed the two-sided R̂_q gate | — | "equilibrated" for uncertified data |
| provisional | any ensemble statement not yet certified | uncertified | — |
| chart | the cocycle ω | cocycle | — |
| detached | cocycle no longer tracking the manifold (see memory: detachment bug) | — | — |

## 3. Standard state-report block

Every report on a crystal-derived state includes:

- **Provenance**: start state + segment lineage (e.g. "m4 r-crystal → lam35 →
  snap11000 → lam35c"), seeds, # chains.
- **Action**: full action / EDQ-only, λ or λ_EDQ, b (and zleg/cimp if
  nonstandard).
- **Legality**: f_FK (and f_e when relevant).
- **Crystallinity**: N_gr, f_G1, f_reg, host phase.
- **Curvature**: e\* vs ⟨e⟩.
- **Census**: n_ill (and N_cx when species matter).
- **Certification status**: certified (gate numbers) or provisional.
- Additional observables as relevant, if not disproportionately expensive.

## 4. Figure rules

- Stamp the **action and couplings** (λ or λ_EDQ, b), the **box/host**
  ("m4 r-crystal start"), and #chains/seeds on every figure (title, subtitle,
  or caption).
- Mark **provisional data** explicitly whenever certification matters for the
  claim the figure supports.
- Distances in **cells** with vertex-spacing equivalents (a_v = 0.185 cell)
  where displacement is the point; time axes in sweeps (or sweeps/1000,
  labeled).
- Print the **full file path** of every generated figure.

## 5. Claim standards

- **Stationarity**: late-window slope ± block-bootstrap σ, window stated.
  Full-series OLS alone is not evidence (transients dominate it).
- **Transport**: name the channel — tracer (MSD) / population (Jaccard) /
  collective (ΔP, W). Visitation counts are upper bounds, never rates.
- **Ensemble averages**: quote R̂_q, window, and ESS when the average carries
  the claim; label provisional data.
- **Species**: complex censuses report size spectrum + illegal-edge signature
  + Q_c (cell-mean reference).

## 6. Run-launch checklist

- SNAP ≡ 0 (mod TS) — otherwise snapshots silently skip.
- Instrumentation stated at launch: cents (on by default), event log, flip log.
- Resume only from **verified-consistent** snapshot pairs (edge-set check;
  regenerate the chart via `regen_cocycle.py` when broken).
- `caffeinate -i` on every launch.
- Chain naming: `<box><phase>-<action>-<λ>-<segment><seed>` preferred for new
  campaigns (e.g. `m4r-fa0.35-s2c`); legacy names stay as-is.
- Certification campaigns: both sides (below + above inits) planned from the
  start.
