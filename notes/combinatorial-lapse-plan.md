# The combinatorial lapse: foliation freedom in the sampler

*2026-07-15. Status: proposal (adopted after discussion with Aaron). Companion documents:
`ensemble-as-wavefunction-program.md` (the |Ψ|² program this extends),
`sampler-dynamics-matter-coupling-plan.md` (quench/dynamics tooling),
`physics-review-2026-07.md` §IV (background and current spectroscopy results).*

---

## 0. The question this answers

The RK construction gives, from the sampler's stationary measure π(T) ∝ e^{−S(T)}, a
Hamiltonian H = 𝟙 − D^{1/2}WD^{−1/2} with ground state ψ₀ = √π. Two standing discomforts:

1. **Non-uniqueness.** Many Hamiltonians share this ground state — the excited spectrum
   (everything dynamical) depends on the move set and kernel, not on S alone. Which H, if any,
   is "the physical one"?
2. **Bespoke time.** Monte Carlo sweep time is an externally imposed clock. What relation does
   it bear to time in general relativity, where "time" is a gauge choice (the lapse N and shift
   Nⁱ of the ADM decomposition t^μ = N n^μ + N^μ)?

**Central claim: these are the same puzzle.** The kernel freedom at fixed π *is* a foliation
freedom, and the lapse is the dictionary entry that makes this precise.

## 1. Refinements to the setup (for the record)

- Ground-state expectation matching ⟨ψ₀|Ô|ψ₀⟩ = Σ_T π(T)O(T) holds for **diagonal**
  observables (functions of the configuration). ψ₀ = √π > 0 is a *stoquastic* ground state;
  off-diagonal operators are determined too, but matching is only automatic on the diagonal
  algebra.
- The non-uniqueness class is larger than "kernels with stationary π": any positive
  semidefinite H′ with H′√π = 0 qualifies, including non-stochastic ones. The invariant
  content of the construction is exactly the measure.
- The chain is already the **Euclidean** side: it furnishes e^{−Hτ} for one member of the
  family, with τ = sweep time. Whether any Wick rotation e^{−iHt} of any member is "physical
  dynamics" is the open question — this plan is about organizing that question correctly.

## 2. The identification: kernel freedom = lapse freedom

Metropolis with local (Pachner) moves satisfies detailed balance **move by move**, so the RK
Hamiltonian is a sum of local positive terms each annihilating the ground state:

  H = Σ_m h_m,  h_m ⪰ 0,  h_m √π = 0  for every elementary move m.

Modulate the *attempt rates* by a positive weight N(x) on the triangulation (attempt moves in
region x at rate ∝ N(x)). Detailed balance w.r.t. π is preserved term by term, and

  H[N] = Σ_x N(x) h_x.

Compare ADM: H_GR[N] = ∫ d³x N(x) 𝓗_⊥(x). Conclusions:

- Every member of {H[N] : N > 0} has ground state √π. **The Hamiltonian non-uniqueness (at
  least the physically implementable part of it) is parameterized by a lapse.**
- Standard uniform-sweep MCMC is the **unit-lapse gauge** N ≡ 1 (the synchronous / Gaussian
  slicing of the sampler). Our bespoke time is bespoke exactly the way a coordinate time is.
- There is a shift analog: composing updates with relabeling/transport flows tangent to the
  triangulation (moving the "attention frame") also fixes π — the momentum-constraint
  multipliers.
- **Gauge-invariant content = the measure** (equal-time correlations: roundness, structure
  factors, degree statistics). **Lapse-dependent content = all relaxation dynamics** (every τ,
  the GEVP spectra, the gap taxonomy — all computed in the N ≡ 1 gauge and transforming under
  lapse changes).
- Local transformation law = the proper-time relation dτ_eff = N(x)dt: an observable localized
  near x relaxes with τ_x → τ_x / N(x). ("Gravitational redshift" for the sampler.)

**Lapse collapse.** In GR, singularity-avoiding slicings (maximal, 1+log) work by driving
N → 0 where geometry degenerates, so proper time stops accumulating there. In the sampler, the
measured local acceptance rate collapses in stiff/glassy regions — the chain's proper time
stalls exactly where the landscape is singular. The glass wall (g ≳ 0.01) is, in this language,
a region where **no positive lapse keeps evolution ergodic**.

## 3. Tier (a): lapse as gauge data — the consistency probe

*Cheap; requires only weighted move selection in the sampler.*

Add static vertex weights N_v > 0; select move sites with probability ∝ (local mean of N over
the simplices the move touches). No physics changes — that is the point; it is a falsifiable
test of the entire time-as-gauge picture.

Predictions to verify:
1. **π-invariance**: certification gate passes identically under strongly nonuniform N
   (e.g. hemispherical profile with N ratio 5:1).
2. **Redshift law**: local autocorrelation times scale as τ_v ∝ 1/N_v; cross-region
   correlators show the mismatch dictated by dτ = N dt.
3. **Foliation independence**: evolve half the manifold "fast" then the other half (and
   vice versa); the equal-time ensemble reached is identical — discrete many-fingered
   evolution consistency.

Failure of (1) would mean the move-selection change broke detailed balance (bug). Failure of
(2)/(3) with (1) intact would be genuinely interesting — a breakdown of locality of the
h_x decomposition.

## 4. Tier (b): lapse as measured geometry — the history is a 4D spacetime

A Pachner move in dimension 3 is the boundary trace of **gluing one 4-simplex** onto the
triangulation (2–3 move = attaching a 4-simplex across two adjacent tets; 1–4/4–1 and 3–2
likewise). The Monte Carlo history is therefore a 4-dimensional simplicial **cobordism**
between initial and final 3-triangulations — the chain *builds* a discrete spacetime, one
4-simplex per accepted move. (Precedent: spin foams as histories of evolving spin networks;
Markopoulou–Smolin causal spin-network evolution [MS97].)

In ADM, √(−g) = N√h: the lapse is **4-volume deposited per unit coordinate time per unit
3-volume**. Combinatorial transcription — a *measured* lapse field requiring no new machinery:

  N(v) ≡ (accepted moves involving vertex v per sweep) / (local 3-volume share D_v/4),

with D_v the vertex degree (Σ_v D_v = 4f₃, so D_v/4 is v's tet-share). Under uniform attempt
rates this is the **local acceptance density**. Its spatial structure is *output*, not gauge
input: the sampler exhibits a dynamically determined foliation.

Measurements (see §7 for the concrete plan):
- N(v) vs local curvature data (vertex degree/deficit, star VDV, incident edge degrees) across
  phases: extended (g = 2–4e−3), hub (β = 0 flat-pinned), glass-adjacent (8e−3).
- Spatial correlations of δN on the 1-skeleton: correlation length ξ_N vs g — approaching the
  glass this is exactly **dynamical heterogeneity** (χ₄-style physics) in the geometry
  ensemble.
- Split kinematic vs energetic suppression: a move can fail because *invalid* (hard
  combinatorial constraint) or *Metropolis-rejected* (energetic). Two distinct lapse-collapse
  mechanisms; record attempted-valid and accepted separately.

## 5. Tier (c): lapse as dynamical data — imposing a constraint

The ADM lesson: the lapse matters because the action is **linear in N**, which is what makes
𝓗_⊥ a constraint rather than a dynamical field. Transcription: extend configurations to
(T, N) with weight e^{−∫N·𝓒(T)} for a candidate scalar-constraint density 𝓒; integrating over
the lapse (over the correct contour — Teitelboim's proper-time representation; the contour
subtleties are the causal-vs-acausal propagator issue [T82, HH91]) *imposes* 𝓒(T) = 0.

> **Adding a combinatorial lapse dynamically = choosing and imposing a discrete Hamiltonian
> constraint. The lapse is scaffolding; the physics question is what 𝓒(T) is** — the discrete
> analog of 𝓗_⊥ = ⁽³⁾R + K² − K_ij K^ij.

Two measured facts already constrain this program:

1. **The glass is an obstruction to sharp constraints.** Our quadratic penalties are soft
   constraints; the constraint reading is their β → ∞ limit, and we have measured that limit
   to fail ergodically (the glass wall): *the constraint surface is disconnected under Pachner
   moves*. Ergodicity breaking as an obstruction to constraint quantization on discrete
   superspace — a physical statement, not a numerical nuisance.
2. **Gravity-likeness requires criticality.** Gravity's signature is not having *a*
   constraint but the **hypersurface-deformation algebra**: [H[N], H[M]] closing on the
   momentum constraint with metric-dependent structure functions, and an enormous constraint
   surface (extensively many zero modes = many physical states). Our extended phase is gapped
   with a *unique* ground state — the constraint surface is a point, and the on-shell algebra
   closes vacuously. A necessary condition for the gravitational reading is therefore **gap
   closing / extensive near-zero modes somewhere in (g, d̄, k)**. The S(k→0) hyperuniformity
   test and a gap search along the d̄ axis are the hunt. Note also that discreteness itself
   resists: vertex-translation symmetry (the classical residue of lapse freedom) is generically
   broken by discretization and restored only at special "perfect" actions [D08, BD09] — we
   should *measure the violation*, not assume it away.

Long-range target: a numerical **deformation-algebra test** — whether the commutator of two
locally-modulated evolutions, ([H[N], H[M]] applied via short weighted-sweep protocols), acts
like a transport/shift generator with the right N∂M − M∂N structure. Formulable via response
functions; hard; the payoff is a direct measurement of how close the ensemble's time structure
is to GR's.

## 6. Falsifiers and decision points

| test | confirms | falsifies |
|---|---|---|
| (a1) π-invariance under nonuniform N | implementation + DB | (bug, not physics) |
| (a2) τ_v ∝ 1/N_v redshift | locality of h_x decomposition; time-as-gauge | local RK picture |
| (b) N(v) structure vs curvature | dynamically-determined foliation; heterogeneity → glass | featureless N(v) would weaken the lapse-collapse story |
| (c) gap closing on some locus | possibility of constraint-like continuum | gapped everywhere ⇒ RK-fluid reading only |
| deformation algebra ~ Dirac | gravitational time structure | generic non-closure ⇒ sampler time ≠ GR time beyond analogy |

## 7. Concrete first steps (ordered)

1. **Measured lapse from existing machinery (tier b, zero D-core risk):** instrument
   `sampler.d` with per-vertex accepted/attempted-valid counters (increment for vertices in
   each move's support; O(f₀) memory), expose via C API like `degree_histogram`; script
   `scripts/local_lapse.py` to map N(v) on certified seeds across the three phases, correlate
   with degree/deficit, and extract ξ_N(g).
2. **Weighted move selection (tier a):** static N_v vertex weights, weighted site selection;
   run predictions (a1)–(a3) on small certified families.
3. **Write up + decide on 𝓒 candidates (tier c):** survey discrete 𝓗_⊥ forms (Regge-style
   deficit combinations vs our penalty densities); design the lapse-integration toy on a small
   N where exact enumeration is possible.

## 8. References

- **[ADM62]** R. Arnowitt, S. Deser, C. Misner, *The dynamics of general relativity*, in
  Gravitation: An Introduction to Current Research (1962); reprinted GRG 40 (2008).
- **[T82]** C. Teitelboim, *Quantum mechanics of the gravitational field*, Phys. Rev. D 25
  (1982) — proper-time/lapse-integration representation of the propagator.
- **[HH91]** J. Halliwell, J. Hartle, *Wave functions constructed from an invariant sum over
  histories satisfy constraints*, Phys. Rev. D 43 (1991) — lapse contours, causal vs acausal.
- **[MS97]** F. Markopoulou, L. Smolin, *Causal evolution of spin networks*, Nucl. Phys. B 508
  (1997) — local graph moves stacked into discrete spacetime; closest realization of tier (b).
- **[D08]** B. Dittrich, *Diffeomorphism symmetry in quantum gravity models*,
  arXiv:0810.3594 — broken vertex-translation (lapse) symmetry in discrete gravity.
- **[BD09]** B. Bahr, B. Dittrich, *Improved and perfect actions in discrete gravity*,
  Phys. Rev. D 80 (2009).
- **[PW86]** T. Piran, R. Williams, *Three-plus-one formulation of Regge calculus*,
  Phys. Rev. D 33 (1986) — lapse/shift as vertex-displacement freedom, simplicially.
- **[H09]** P. Hořava, *Quantum gravity at a Lifshitz point*, Phys. Rev. D 79 (2009) — a z = 2,
  RK-adjacent gravity; the projectable-vs-nonprojectable lapse issue as a cautionary tale.
- **[PW81]** G. Parisi, Y.-S. Wu, *Perturbation theory without gauge fixing*, Sci. Sin. 24
  (1981) — stochastic quantization: the rival reading of MCMC time (fictitious 5th time).
- **[PaW83]** D. Page, W. Wootters, *Evolution without evolution*, Phys. Rev. D 27 (1983) —
  relational time from a timeless state; the matter-clock alternative to a lapse.
- Background on RK structure: D. Rokhsar, S. Kivelson, PRL 61 (1988); C. Henley, J. Phys.:
  Condens. Matter 16 (2004); C. Castelnovo et al., Ann. Phys. 318 (2005). ADM/lapse pedagogy:
  R. Wald, *General Relativity*, Ch. 10; E. Gourgoulhon, *3+1 Formalism in General Relativity*
  (Springer 2012). Dynamical heterogeneity (for §4's ξ_N): L. Berthier et al. (eds.),
  *Dynamical Heterogeneities in Glasses, Colloids, and Granular Media* (OUP 2011).
