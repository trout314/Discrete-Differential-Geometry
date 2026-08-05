---
name: defect-density-hu
description: "RETRACTION + RESULT: the 'crystal-grade HU' headline is a normalization artifact (the old estimator calls RANDOMISED defects hyperuniform too); shuffling WHOLE DEFECTS rigidly gives 0.999+-0.026 = exactly Poisson; complexes carry residual NEGATIVE charge Q=-0.97, 77% negative"
metadata:
  type: project
---

Measured 2026-08-05 after the user flagged that "our conclusions about
curvature fluctuations may be misleading because the crystal itself is
hyperuniform." Confirmed, with a mechanism, and replaced with the right
measurement. Tools: `scripts/defect_dynamics/defect_density_hu.py` (+
`_figure.py`). Data: `data/defect_hu/*.json`, figure
`data/figs/defect_density_hu.png`. R m4 crystal on T^3, V=10176, box 4 cells.

    ensemble                  snaps n_cplx perm-null vtx-reloc vertex centroid RIGID-Q RIGID-n
    lam=0.40 (8 chains)         55    9.3   0.0015    0.523     7.31   1.07   0.999   1.003
    lam=0.35 (3 chains)         18   21.3   0.0038    0.546     8.20   0.93   0.958   0.865
    CONTROL defects randomised  55   60.6   0.0215    0.988     1.03   0.95   0.971   0.947
    sem                                     .0001-.0013 .03-.07  .03-.4  .03-.04 .026-.053

## 1. THE RETRACTION -- the old estimator is normalization-limited

[[statics-hu-verdict]]'s "melt is as hyperuniform as the perfect crystal"
and `sk_torus --field curvature_charge` low-k ratio ~2e-3 do NOT constrain
the defect arrangement. Two independent reasons:

* **Normalization.** S(k) is divided by S2 = sum_v dq_v^2, and the periodic
  crystal carries 97.7-99.2% of it (defect share 0.76% at lam=0.40, 2.27% at
  lam=0.35) while contributing nothing sub-Bragg. **The measured ratio tracks
  the defect share almost exactly** -- 0.76% -> 0.0015 and 2.27% -> 0.0038,
  a 3.0x rise for a 3.0x share. It is a normalization readout, not physics.
* **The null.** The charge-permutation null scrambles the crystal too, so the
  comparison is melt-vs-amorphous and any crystal-like state passes.

**Decisive control**: relocate every defect uniformly at random and re-run
the OLD estimator -- it still returns **0.0213 +/- 0.0011**, still "strongly
hyperuniform", on a configuration with no arrangement at all.

## 2. THE RESULT -- defect complexes are POISSON

The right null (user's, 2026-08-05): shuffle WHOLE DEFECTS rigidly, not the
individual vertex charges. With F_i(k) = sum_{v in i} w_v exp(ik.x_v) the
complex's own amplitude and independent uniform torus translations t_i,
every cross term carries E[exp(ik.(t_i - t_j))] = 0 at nonzero commensurate
k, so the null is ANALYTIC:

    S_obs/S_null = |sum_i F_i(k)|^2 / sum_i |F_i(k)|^2

Each complex's form factor AND net charge divide out; only the relative
PHASES -- the arrangement -- survive. As k->0, F_i -> Q_i, so it becomes
the charge-weighted centroid S(k). Complexes are translated, NOT reoriented
(the crystal locks defect orientation; the P2 null is +0.118, not 0), so
this is positional arrangement alone. Validated by an MC self-test that
applies actual random translations: 0.997 (must be 1.000).

**lam=0.40: 0.999 +/- 0.026 charge-weighted, 1.003 +/- 0.030
count-weighted.** Flat across every shell. The defect arrangement is
EXACTLY POISSON -- not hyperuniform, not clustered. lam=0.35 gives 0.958
+/- 0.053 and 0.865 +/- 0.039 against a control baseline of 0.971/0.947,
so at most a ~2 sigma hint of mild suppression, consistent with hard-core
exclusion between complexes; not a claim.

charge- and count-weighting agree at lam=0.40 (0.999 vs 1.003), so there is
no charge-position coupling either -- big-|Q| complexes are not placed
differently from small ones.

This also CORRECTS [[statics-hu-verdict]]'s "sources cluster, super-Poisson
1.14": with 6x the data and a form-factor-free estimator there is no
clustering. The old g(r) shell at 1-1.4 cells is short-range packing/merge
structure, not a long-wavelength effect.

## 2b. COMPLEXES CARRY RESIDUAL NEGATIVE CURVATURE (user's prediction, confirmed)

Net charge per complex Q_i = sum_{v in i} dq_v:

    lam=0.40   Q = -0.965 +/- 1.347   77% negative   |Q|/sum|dq| = 0.308
    lam=0.35   Q = -0.895 +/- 1.706   70% negative   |Q|/sum|dq| = 0.301
    CONTROL    Q = -0.148 +/- 0.945   52% negative   |Q|/sum|dq| = 0.988

So ~2/3 of a complex's per-vertex charge magnitude cancels internally, but a
definite one-signed NEGATIVE residual survives -- the user's reasoning that
defects persist in order to satisfy the mean-edge-degree pin, and so must
carry residual negative curvature. The control decomposes it: scattering the
SAME charges into singletons gives mean -0.148 per defect vertex and 52%
negative, and -0.148 x <m>=6.9 = -1.02 ~ the measured -0.965. The complex
residual is the per-vertex bias accumulated over ~7 vertices; the 77%-negative
consistency at complex level is that weak per-vertex bias aggregating.

Consequence: the defect gas is a population of LIKE-CHARGED objects whose
total is (softly) pinned by the action. That is exactly the setting in which
hyperuniformity could have appeared -- and it does not. The pin fixes the
total charge without inducing any spatial correlation.

## 3. The per-VERTEX relocation number, and what it actually measures

Defect SOURCE charge vs a PER-VERTEX relocation null: **0.523 +- 0.030**
(lam=0.40), 0.546 +- 0.070 (lam=0.35), control 0.988 +- 0.034.

This is NOT an arrangement measurement -- scattering individual vertex
charges destroys the intra-complex +/- cancellation as well as the
positions, so it measures how NEUTRAL a complex is (cf. sec 2b: neutrality
fraction 0.31, and 0.31 ~ 0.52 in power-like units). Reported for the
record; the rigid shuffle (sec 2) is the arrangement estimator. The
mechanism is INTRA-COMPLEX charge cancellation: the same estimator on the
indicator field (all values = 1) gives 7.31, on the signed charge 0.523 --
a factor 14 apart on identical positions. Each complex is a partially
neutral multipole; the null scatters its +/- vertices to independent sites
and destroys the cancellation. **Local neutrality alone suppresses low-k
charge power regardless of arrangement** -- exactly the reason the Level-2
"S(k->0) -> 0 => Hamiltonian constraint" test of
`notes/ensemble-as-wavefunction-program.md` is passed trivially by any
locally-neutral medium and is NOT evidence of a constraint. Only partial
(0.5, not ->0) because knots carry net charge (Q ~ -5, [[lifetime-vs-charge]]).

## 4. THE FORM-FACTOR TRAP

lam=0.40 defects come in complexes of mean size 6.9 (9.3 complexes per
snapshot). A compact-cluster process has S(k->0)/S_Poisson = <m^2>/<m> = 7.5
even with perfectly Poisson centres, so the raw defect-VERTEX ratio 7.31 is
essentially ALL form factor. Consistency identity centroid x form_factor vs
vertex holds to ~20% (8.88 vs 7.31), the gap being the cluster form factor
P(k) < 1 at finite k -- right sign, right size.

## 5. THE BIASED NULL THE CONTROL CAUGHT (methodology)

First attempt at the charge null was a TRANSPOSITION: swap defect charges
with those of random crystal sites, keeping the crystal. It returned
**0.61 +- 0.06 on randomised defects where it must be 1.0** -- a 40% "low-k
suppression" that was pure construction, because a transposition perturbs
BOTH ends (3 disturbed sites in the null vs 2 in the observed). Fix:
relocate a ZERO-BACKGROUND field so observed and null differ only in the
positions of an identical value multiset. Same shape as the centred-smoother
bug in [[ecmc-blob-ab]] and the aliasing in [[statics-hu-verdict]]:
**always run the estimator on a state where you know the answer.**

## 6. WHAT LIMITS A k->0 CLAIM

Only ~9 complexes per snapshot in a 4-cell box => mean spacing 1.9 cells =>
only the |k|=1 shell (9 modes) is genuinely long-wavelength. A real
hydrodynamic-limit statement needs a bigger box (m6 = 3.4x volume, ~32
complexes) or higher density. lam=0.35 has 2.3x the density but is
NON-STATIONARY (slowly densifying, [[mobile-gas-liquid]]) -- cross-check
only. Datasets with cocycle lifts at all: `data/mgas` (126 pairs),
`data/run5h` (32); nothing else has them.

## 7. What is still OPEN

The TOTAL curvature field (sources + screening halo on legal vertices) is
NOT measured here -- separating the halo from the crystal needs a background
model. So "is the melt's total curvature field hyperuniform?" remains
genuinely open; what is settled is that the published number never answered
it. Note also [[six-web-gauge]] is a COMPETING explanation for any curvature
correlation (L_vec(k) = S_mono(k)/k^2 exactly).

Related: [[statics-hu-verdict]], [[curvature-length-scale]],
[[no-halo-verdict]], [[ecmc-blob-ab]].
