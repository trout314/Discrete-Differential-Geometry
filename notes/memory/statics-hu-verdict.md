---
name: statics-hu-verdict
description: "BIG RESULT — lam=0.40 melt's TOTAL curvature field is crystal-grade class-I HU (exp 1.98 vs pristine 2.01) even though complexes cluster; screening is essentially complete; pair/orientation correlators inconclusive at n=94"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T22:22:29.154Z
---

Measured 2026-07-25. Tools: `scripts/defect_dynamics/defect_statics.py`
(complex-level pair stats + ball variance by flicker attribution) and
`curvature_balls.py` (full vertex-field ball variance with pristine + shuffled
controls). Data: 11 lam40*/lam40x snapshots (2 chains), R m4, box = 4 cells;
positions = cocycle lift (crystal chart, real Euclidean balls). PROVISIONAL.
JSONs: scratchpad `statics_lam40.json`, `curvballs_lam40.json`; figure
`data/figs/defect_statics_lam40.png`.

**The headline (point 2 of the ADM discussion): the TOTAL curvature field of
the melt is exactly as hyperuniform as the perfect crystal.** sigma^2_Q(R) of
dq(v) = q_R(v) − qbar over ALL vertices: melt tracks pristine to within a few
percent at every R (3.63 vs 3.53 at R=0.3 ... 177 vs 182 at R=2.0); scaling
exponents melt 1.98, pristine 2.01 (class-I HU surface law), shuffled-charge
ceiling 2.55 and 4–16x higher. The defect gas does NOT degrade crystal-grade
hyperuniformity — equilibrium screening by the legal matrix compensates the
sources almost completely. Organized defect MOTION is not needed for HU
curvature; this rehabilitates the HU goal without a worm sampler (contra the
worry in the ADM discussion; see [[defect-travel]] for why motion couldn't do
it anyway).

**Meanwhile the complex-level SOURCE process is NOT hyperuniform**: g(r) shows
a shell at ~1–1.4 cells (g ~ 1.3–1.6) and depletion at ~1.75 (0.74±0.11), and
the complex-point-charge ball variance is super-Poisson (ratio 1.14 at
R>=1.2). Sources cluster; the response field cancels them. Classic Coulomb
phenomenology.

**Inconclusive at n=94 complexes**: connected C_QQ(r) and the chord-director
nematic P2(r) (null = +0.118 from crystal-axis anisotropy, cross-snapshot
pairs) both fluctuate around 0 without coherent trend. Need ~10x more
snapshots for the directional-interaction question. F/P cross-covariance in
balls ~ matched-Poisson (no vacuum-screens-strings signal resolvable).

**Bugs caught (methodology notes)**:
- Raw <dQ_i dQ_j> is positive at ALL r because <dQ> = −1.72 != 0 — must use
  the CONNECTED correlator.
- P2 null is NOT 0: chords lie on crystallographic axes (discrete axis set)
  — null from cross-snapshot pairs.
- g(r) hole below r~0.4 is partly definitional (adjacent complexes merge).
- reference_frac_positions returns PERIOD-M (cell-unit) coords — do NOT
  multiply by mcell (aliasing gave pristine variance above the coherent
  Cauchy-Schwarz bound, which is how it was caught).

Caveats: box only 4 cells (R <= 2), R=2.0 ball is over half the box; melt
positions are lift/reference-chart (combinatorial curvature in crystal chart —
right chart for this question); single lambda; uncertified.

Next steps flagged: more snapshots for C_QQ/P2; bigger box (m6?) to extend R;
two-knot controlled interaction experiment (frozen vertices + directed slides)
as the decisive directional test. See [[flicker-background]],
[[lifetime-vs-charge]], [[six-web-gauge]].
