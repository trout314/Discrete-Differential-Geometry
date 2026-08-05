---
name: defect-density-hu
description: "RETRACTION + RESULT: the 'crystal-grade HU' headline is a normalization artifact (the old estimator calls RANDOMISED defects hyperuniform too); the defect complex arrangement is exactly POISSON, and the source-charge suppression 0.515 is intra-complex neutrality, not order"
metadata:
  type: project
---

Measured 2026-08-05 after the user flagged that "our conclusions about
curvature fluctuations may be misleading because the crystal itself is
hyperuniform." Confirmed, with a mechanism, and replaced with the right
measurement. Tools: `scripts/defect_dynamics/defect_density_hu.py` (+
`_figure.py`). Data: `data/defect_hu/*.json`, figure
`data/figs/defect_density_hu.png`. R m4 crystal on T^3, V=10176, box 4 cells.

    ensemble                   snaps  n_cplx   perm-null  reloc-null  vertex  centroid
    lam=0.40 (8 chains)          55     9.3    0.0015     0.515       7.31    1.07+-0.04
    lam=0.35 (3 chains)          18    21.3    0.0038     0.548       8.20    0.93+-0.04
    CONTROL defects randomised   55    60.5    0.0213     0.997       1.02    0.92+-0.03

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

Form-factor-free complex-centroid point process (one point per connected
defect complex, circular torus mean, vs the random-crystal-site null), flat
at every accessible shell:

    |k|          1.0    2.0    3.0    4.0    5.0
    lam=0.40   1.055  1.022  1.026  0.994  0.980   (+-0.05 .. +-0.01)
    lam=0.35   0.971  0.906  0.985  0.948  0.983

No hyperuniformity (that would fall toward 0), no clustering. Taking the
control's 0.92 as the estimator's residual scale, the honest statement is
**S_centroid(k->0) = 1.0 +- 0.1 = Poisson**. This CORRECTS
[[statics-hu-verdict]]'s "sources cluster, super-Poisson 1.14" -- with 6x
the data and a form-factor-free estimator there is no clustering; the old
g(r) shell at 1-1.4 cells is short-range packing/merge structure, not a
long-wavelength effect.

## 3. The one real signal, and why it is NOT order

Defect SOURCE charge (anomalous dq on impurity vertices, zero elsewhere) vs
the matched relocation null: **0.515 +- 0.029** (lam=0.40) and 0.548 +-
0.069 (lam=0.35), control 0.997 +- 0.039. A genuine 12-sigma factor-2
suppression that reproduces across densities.

But the positions are Poisson (sec 2), so it cannot be spatial order. The
mechanism is INTRA-COMPLEX charge cancellation: the same estimator on the
indicator field (all values = 1) gives 7.31, on the signed charge 0.515 --
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
