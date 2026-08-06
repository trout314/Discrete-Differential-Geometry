---
name: crystal-library-gas-campaign
description: "Library-wide defect-gas survey (2026-08-06): all 10 TCP crystals under the minimal 3-term action (volume pin + FLAT mean-edge-degree pin + c_imp*sum m^2); the pin gap is an EXTENSIVE forced-defect debt whose SIGN splits the library in two"
metadata: 
  node_type: memory
  type: project
  originSessionId: 93e02d27-494a-461f-a866-e57819478b52
  modified: 2026-08-06T13:01:43.061Z
---

**Question (Aaron, 2026-08-06):** for EVERY crystal in the library, can the
minimal action hold a stable dilute defect gas at the FLAT edge-degree target,
and what do the resulting defect species look like (HTML catalog)?

**Action (nothing else on):**
S = 0.1(f3−f3_ref)² + 1.0(f1 − 6f3/e*)² + c_imp·Σ_v m(v)²,
e* = 2π/arccos(1/3) = 5.1042993. Wire-up: `SamplerParams(num_facets_coef=0.1,
hinge_degree_target=e*, num_hinges_coef=1.0, …_variance_coef=0,
hinge_degree_target_coef=0, codim3_degree_target_coef=0)` +
`set_n6_potential(0.0, c_imp, tilt=[0]*5)`. `num_hinges_coef` IS the mean-edge
pin: globalCurvPenalty = (f1 − 6f3/e*)² (source/sampler.d:37). Variance
couplings default ON — must be explicitly zeroed ([[vdv-hdv-conflicts-with-tcp]]).

**KEY STRUCTURAL FACT (derived, then measured).** The two pins are not jointly
satisfiable by the pristine crystal. e_native = 6 − 12/CN, so the pin gap
gap = f1 − 6f3/e* is EXTENSIVE, and the crystal must pay it in defects:
~|gap|/0.1755 forced 2→3-equivalent moves (Δgap: 2→3 = −0.1755, 1→4 = +0.4734).
Per 1000 vertices the forced debt is
  r 0.6 | mu 6.4 | z,p,delta 15.6 | c15,c14,c36 32.0 | sigma 34.7 | a15 51.3
and the SIGN of (e* − e_native) splits the library:
  e_nat < e* (pays in 2→3): c15, c14, c36, mu, r
  e_nat > e* (pays in the reverse): a15, sigma, z, p, delta
**R is the flat crystal** — e_nat = 5.1042254 vs e* = 5.1042993, gap 0.12 f1 at
m2, i.e. <1 forced move in 1272 vertices. Its gas is purely thermal.

**Infrastructure built (this session):**
- `scripts/defect_dynamics/crystal_gas_scan.py` — one (crystal, c_imp) cell;
  `LIBRARY` dict fixes the reference cell per crystal, all f0 ∈ 1000–2500 for
  comparability. Uses `Recorder` (maximalist recording); adds `.obs.jsonl` with
  vertex-based complex sizes, max_m, Σm², turnover. `--cocycle` attaches the
  harmonic cocycle AT SWEEP 0 (post-hoc regen FAILS on these gases — monodromy).
  GOTCHA: `recording._components` returns component sizes in EDGES; every defect
  census in the tree counts VERTICES — hence local `vertex_components()`.
- `crystal_gas_sweep.py` — process-pool driver over the crystal × c_imp grid.
- `crystal_gas_report.py` — verdict table; tests DISPERSED / DILUTE / STABLE
  separately, stationarity via late-window slope ± moving-block bootstrap σ,
  with `turnover` alongside to catch glassy arrest.

**References built:** a15 m5, c14 m5, c36 m4, z m6, mu m4, p m3, delta m3,
sigma m4 added to `data/tcp_reference/` (all VALIDATED: all edges {5,6},
Z-census exact, χ=0, orientable).

**Benchmark:** ~10 ms/sweep at f3=3672; c15 m4 at c_imp=0.5 equilibrates within
50 sweeps to n_ill≈104 in 41 complexes, top1=4, max_m=3, turnover 0.95/50sw —
a dispersed mobile gas. The f1 pin is NOT driven to zero (sits ~+5): paying it
costs m² energy, so equilibrium is a compromise between the two terms.

Related: [[m2-only-gas]], [[edq-only-melting]], [[defect-viewer-tools]],
[[maximalist-recording]], [[reporting-conventions]].
