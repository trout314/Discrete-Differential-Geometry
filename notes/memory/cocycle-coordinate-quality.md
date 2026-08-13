---
name: cocycle-coordinate-quality
description: "MEASURED: harmonic-gauge cocycle coords are a NON-FOLDING embedding (sum|V|/det B = 1 to 1e-6) with tet min/max edge ratio 0.83 median, mode exactly 2/sqrt6; CONVENTIONS' '~25% registry spread' was a lattice-fractional artifact (true Cartesian spread 7%)"
metadata: 
  node_type: memory
  type: project
  originSessionId: a2047098-7e89-4d8e-b1c5-3231dfa0bf45
  modified: 2026-08-13T17:27:47.288Z
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
**EXACT SCOPE OF THE FIX (verified):** S_obs and S_null are BIT-IDENTICAL per
mode whitened vs not (they are built from `frac` and integer n, both untouched);
only the |k| LABEL of each mode moves (max shift 2.4 units, rank corr 0.90).
So a derived number moves only insofar as it aggregates over modes by |k|.
Measured on R m4: shell-binned low-k ratios essentially unchanged (0.0017 vs
0.0016); a narrow-window per-mode log-log fit swings a lot (r_gas 1.579 →
2.241) purely because window membership changes. Don't quote the latter kind
of fit without saying which modes it used.

**R-hosted re-run done same day** — see [[defect-density-hu]]: the Poisson
verdict SURVIVES (lam=0.40 rigid_charge 0.999±0.026 → 0.980±0.018). The metric
fix is a few-percent effect on those estimators; the big shifts that showed up
in the re-run were an unrelated script drift.

**SECOND, PRE-EXISTING FRAME BUG — FOUND AND FIXED (same day).**
`lattice_basis` returns a Euclidean-reduced basis that is NOT diagonal (R m4
gives [[4,4,0],[0,4,0],[0,0,4]]), so `frac * np.diag(basis)` is not a Cartesian
position and per-axis min-imaging is not the min-image rule. Two distinct
errors, both independent of whitening:
 (a) diag-only distances on R m4 are wrong by >10% for **68%** of pairs and
     >30% for 18%, range 0.62–1.62x — that SMEARS a correlation function, it
     does not merely shift it;
 (b) even with the full basis, `d -= round(d)` gives the primary image, not
     the nearest: wrong for **26%** of random pairs, up to 2.4x too far.
Both closed by `cocycle.min_image(dfrac, basis)` (minimizes over neighbouring
images; validated against a ±3 brute force). Consumers keep FRACTIONAL
positions and let the metric enter only there.

**Audit outcome, script by script:**
- `sl_verdict.py` — **NO BUG** (I was wrong to list it). Its `P` cancels
  identically: the centroid is frac[0] + mean(wrap(frac − frac[0])) and a
  per-axis scale commutes with the per-axis wrap; S_N/S_Q use fractional coords
  against integer n with |n| shells. Verified bit-identical whitened vs not.
  Note in the file's docstring so nobody "fixes" it.
- `pass2_structure.py` — fixed. Conclusions unchanged: Rg ~ N^nu 0.49 → 0.50,
  blinker templating still null (KS p 0.52 → 0.35). Only the Rg CONSTANTS move,
  down 10–13% (0.31 → 0.27 cells for the largest bin) — exactly the
  bi-Lipschitz expectation (exponents survive, constants don't).
- `carrier_gr.py` — fixed, and **the conclusion changes**. See
  [[defect-density-hu]] §10.

See [[intrinsic-geometry]] (§6 policy) and [[cocycle-vertex-lift]] (how the
lift is maintained).
