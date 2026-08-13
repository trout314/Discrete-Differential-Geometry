---
name: cocycle-coordinate-quality
description: "MEASURED: harmonic-gauge cocycle coords are a NON-FOLDING embedding (sum|V|/det B = 1 to 1e-6) with tet min/max edge ratio 0.83 median, mode exactly 2/sqrt6; CONVENTIONS' '~25% registry spread' was a lattice-fractional artifact (true Cartesian spread 7%)"
metadata: 
  node_type: memory
  type: project
  originSessionId: a2047098-7e89-4d8e-b1c5-3231dfa0bf45
  modified: 2026-08-13T16:27:36.732Z
---

2026-08-13. Quantified how well-behaved the cocycle coordinates are
(scratch scripts in `$CLAUDE_JOB_DIR/tmp/{aspect,run_pristine,run_states,run_defects,run_fold}.py`;
figure `data/figs/cocycle_tet_aspect.png`). Metric used: per-tet
ℓ_min/ℓ_max, with edge vectors read straight off the cochain (harmonic
gauge ⇒ pos(v)−pos(u) = ω_h(u→v), so no min-imaging/basis is involved).

**Gauge subtlety that dominates every number.** `harmonic_gauge` solves the
Laplacian per column, so it is GL(3)-equivariant on the Z³ value space: the
Euclidean structure is a FREE MODULUS. Canonical fix = whiten by G^{-1/2},
G = ⟨ω ω^T⟩_edges. Validated: on pristine non-cubic crystals whitening
reproduces the true Cartesian embedding (frac @ L) to ~1%, without knowing L.
**Always whiten before quoting any length/angle in cocycle coordinates.**

**Results (whitened):**
- Pristine crystals: ratio median 0.82–0.87, min 0.73–0.80, cv(edge) 6.5–7.6%.
  Mode is EXACT: 2/√6 = √(2/3) = 0.81650 = ratio of the two TCP bond lengths
  (1/(2√2) and √3/4). Harmonic gauge ≈ a no-op on a pristine crystal
  (|ω_h − ω|/⟨|ω|⟩ = 0 for cubic a15/c15/z, 2–4% for sigma/r/c14).
- Sampled states (crystal_gas defect gases, mgas m4, run5h): median 0.82–0.85,
  1st pct 0.60–0.77, 0.01–0.04% below 1/2, essentially none below 1/3.
- **All of the tail is defects**: clean tets median 0.831 (min 0.43) vs tets
  touching a deg≠5,6 edge median 0.706 (min 0.219). Two disjoint populations.
- **No folding**: Σ|V_tet| / det(B) = 1.000000–1.000184 across states
  (GL(3)-invariant) — the harmonic embedding is a genuine bijection onto T³,
  not merely a lift. Degenerate tets are ~1e-5 of the population.

**APPLIED (same day).**
- `cocycle.whitening_transform(gram)` added; `torus_positions(..., whiten=True)`
  is now the DEFAULT. `frac` is untouched by it (whitening acts on the value
  space), only `basis` moves, and it is unimodular so det(B) is preserved.
  Tests `test_whitening_*` in tests/test_hyperuniformity.py; 239 pass.
- notes/CONVENTIONS.md §6 table rewritten + new rule 6.

**WHY IT MATTERS — the winding lattice in cochain units is CUBIC for every
host** (coordinates are fractional), so the old `basis` reported a cubic cell
for R/σ/μ/Z/C14/C36/P/δ and pushed that into every |k| bin. Whitening recovers
the true aspect without being told L: R 1.774 → 1.780, C14 1.633 → 1.633,
σ 1.932 → 1.911, μ 5.379 → 5.525; cubic hosts unchanged.
**Impact measured** (crude single-snapshot log-log fit of S/S_null over
k ≤ 2.5, curvature_charge): r_gas exponent 1.579 → 2.241, mgas ab1 (R m4)
−0.160 → −0.043, c15_gas 0.828 → 0.824 (cubic, null). **So any pre-2026-08-13
S(k) result on an R host may shift — re-check [[statics-hu-verdict]] and
[[defect-density-hu]]; A15/C15 results stand.**

See [[intrinsic-geometry]] (§6 policy) and [[cocycle-vertex-lift]] (how the
lift is maintained).
