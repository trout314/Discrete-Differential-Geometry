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
| δq(v) | curvature deviation from a reference, q_R(v) − q̄_R (state the reference). **Not** for complex charges — see Q_c | rad |
| w(e) | link (disclination) charge of edge e, w(e) = 6 − deg(e); flat at deg 6 | — |
| n₆(v) | # of degree-6 edges incident to v | — |
| imp(v) | # incident illegal edges (deg ∉ {5,6}) at v | — |
| Z(v) | coordination = # edges at v = (#tets at v + 4)/2 | count |
| ΣZ | total coordination of a complex, Σ_{v∈c} Z(v) — **the integer that fixes Q_c** | count |
| n_ill | **census**: # vertices with imp > 0 | count |
| N_cx, s | # complexes; complex size (vertices) | counts |
| Q_c | **raw** curvature charge of a complex, Σ_{v∈c} q_R(v) — no background subtraction. See §1.1 | rad |
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

### 1.1 Curvature charge is fixed by coordination — use raw sums

Because the link of a vertex is a triangulated 2-sphere on Z vertices (so it
has 2Z − 4 triangles and its degrees sum to 6Z − 12), and a link-vertex's
degree is the ambient degree of the corresponding edge,

$$\sum_{e \ni v} \deg(e) = 6Z - 12
\qquad\Longrightarrow\qquad
q_R(v) = Z(\pi - 3\theta) + 6\theta$$

so **q_R depends on the coordination alone** — for every vertex, legal or not,
since only Euler is used. Summing over a complex of n vertices:

$$Q_c = \sum_{v\in c} q_R(v) = \left(\textstyle\sum Z\right)(\pi - 3\theta) + 6n\theta$$

A complex's curvature charge is therefore set by **one integer**, its total
coordination ΣZ, and is quantized in steps of

$$\pi - 3\theta = -0.5512856\ \text{rad per unit of } \textstyle\sum Z.$$

Two related identities, useful because they say what is *not* independent
(all three verified exactly on 8000 vertices):

- n₆(v) = Z + 5·imp(v) − Σ(illegal degrees at v) − 12, so n₆ and Z are
  interchangeable once edge degrees are known;
- every illegal edge lies **inside** a defect complex (both endpoints have
  imp > 0), so edge degrees already determine imp(v) and that degree sum.

**Report Q_c raw.** Do not subtract a background mean. Raw Q_c depends only on
the complex's own vertices, so identical complexes give identical values; with
the state mean q̄_R subtracted, every measurement is tied to the whole
configuration and acquires a spurious spread. Measured over 60 species classes
at λ = 0.40: sd(Q_c) is 4.4 × 10⁻¹⁶ (machine epsilon) raw, against 1.8 × 10⁻³
background-subtracted. Report the relative value in addition if it is useful,
never instead.

**The (3,4,4) knot is a charge ladder**, not a single value. Its rungs, indexed
by ΣZ over its five vertices (ΣZ = 70, i.e. mean Z = 14, is the most common):

| ΣZ | 67 | 68 | 69 | 70 | 71 | 72 |
|---|---|---|---|---|---|---|
| Q_c (rad) | −0.0074 | −0.5586 | −1.1099 | −1.6612 | −2.2125 | −2.7638 |

An earlier figure of Q_c(knot) ≈ −1.19 ± 0.18 in this file was an average over
rungs, not a measurement of one: the quoted ±0.18 is narrower than the 0.551
rung spacing. Quote the rung (or ΣZ) rather than a mean over the ladder.

## 2. Standard terminology

| term | definition | acceptable synonyms | avoid |
|---|---|---|---|
| full action | edge pin + n6 legality terms (+ volume pin, always on) | — | "the objective" |
| EDQ-only | full action with n6 off (the R² ablation) | R²-only | bare "EDQ" as an ensemble name |
| defect vertex | imp > 0 **or** non-FK coordination (n₆ ∉ {0,2,3,4}) | — | equating it with imp > 0 alone |
| complex | connected set of defect vertices | illegal complex, defect complex | "defect" for a single vertex |
| induced subcomplex | every simplex whose vertices all lie in a complex — what the defect **is**. Dense, not curve-like | — | confusing it with the closed star |
| closed star | every facet incident to a complex — the region it **occupies**; its boundary is the flux surface | — | confusing it with the induced subcomplex |
| knot | complex whose induced subcomplex is 3 tets of ∂Δ⁴ (f = (5,10,9,3), 1-skeleton K₅, χ = 1); edge signature (3,4,4) | — | treating its Q_c as single-valued — it is a ladder, §1.1 |
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
- **Curvature**: e\* vs ⟨e⟩ (and Q_c **raw**, with ΣZ, when complex charges are
  reported — see §1.1).
- **Census**: n_ill (and N_cx when species matter). When species are the point,
  say which key: edge signature, induced subcomplex, or that subcomplex
  decorated by edge degree ± n₆ (`species_report.py --group`).
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
- **Species**: complex censuses report size spectrum + the grouping key used
  (illegal-edge signature, induced subcomplex, or decorated subcomplex) + Q_c
  **raw** with ΣZ. Never a mean over ladder rungs — see §1.1.

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
