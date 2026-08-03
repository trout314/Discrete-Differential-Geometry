---
name: curvature-length-scale
description: "Curvature fluctuations decay as 1/R^2 (class-I HU) in constrained states — measured with real Euclidean balls, not the graph proxy"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T19:53:06.359Z
---

Measured how curvature-charge (q_R = ½Σδ_e) fluctuations decay with length
scale (2026-07-21). Prediction (Torquato-Stillinger, d=3): if S(k)~k² at small
k (Gauss-law / Hamiltonian constraint), the charge variance in a window of size R
grows as the SURFACE R² (not volume R³), i.e. RMS of the *mean* curvature over a
region ~ R⁻². Poisson/unconstrained: Var~R³, RMS~R^{-3/2}.

CONFIRMED for the sampler's constrained states (run5h above/below-native R
crystal, pooled snapshots): Var_Q ~ R^1.93 ≈ R², RMS ~ R^{-2.03}. Charge-shuffle
null (positions kept, charges permuted): Var_Q ~ R^2.8 ≈ R³, RMS ~ R^{-1.6} ≈
R^{-3/2}. Perfect crystal r_m4: R^1.74 / R^{-2.13} (surface-like with lattice
commensurability wobble). The obs (~-2.0) vs null (~-1.6) gap is the
hyperuniformity signature.

KEY GOTCHA: the graph-BFS window proxy (graph_hyperuniformity.window_variances)
CANNOT see this — graph balls have rough/fractal boundaries whose surface is an
O(1) fraction of the ball, swamping the R² surface law (everything came out
β≈0.87, RMS≈-1.7, crystal ≈ its own shuffle). MUST use real Euclidean balls on
the metric torus via cocycle.torus_positions. Also: a MELT's cocycle embedding
COLLAPSES (needs crystalline connectivity to spread in 3D) -> real-space balls
meaningless (shuffle slope came out 1.12, unphysical); establish melt non-HU via
S(k) instead (r_m3_melt low-k ratio 1.95, clustered).

DIRECT RECIPROCAL alpha FIT (S(k)~k^alpha) IS NOT RESOLVABLE at this box size
(2026-07-21 follow-up). r_m4 has only ~4 commensurate k-shells below the first
CHARGE Bragg peak (~k=5 in box-mode units); the smallest accessible mode (k=1)
already sits at a low-k local max (~k=1.4) with a pre-Bragg dip -- mid-k
structure, NOT the k->0 tail. A naive sub-Bragg fit gives a spurious NEGATIVE
slope (-0.76+/-0.07) from that shape. Do NOT report it as the exponent. Perfect
crystal: S(k)=1e-30 sub-Bragg = STEALTHY hyperuniform (undefined power law).
Pinning alpha=2 directly needs a LARGER defected crystal (fresh sampler run,
more sub-Bragg modes).

What IS solid (real space, decisive control): mean-subtracted charge (point
pattern removed) still gives RMS-density slope -2.01+/-0.01, number-variance
exponent p=1.97+/-0.03 (=surface law R^2), vs its own charge-SHUFFLE null at
-1.59+/-0.03 (Poisson). So the CHARGE ARRANGEMENT is genuinely class-I
hyperuniform (alpha>1), not just the crystalline vertex positions. Surface law
saturates for any alpha>1 so it confirms alpha>1 / consistent with the k^2
Hamiltonian-constraint expectation but cannot pin alpha=2.

Scripts: scratchpad/{curv_scale_real, sk_exponent, reconcile, final_fig}.py.
Figures: curv_charge_surfacelaw.png (headline), sk_exponent.png (why direct fit
fails). See [[hyperuniformity-refactor]], [[five-illegal-knot]].
