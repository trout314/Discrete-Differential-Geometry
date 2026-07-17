# The Equilibrium Ensemble as |Ψ[h]|²: the Sampler as a Wavefunction of the Universe

*Program document, drafted 2026-07-11. Companion to
`notes/sampler-dynamics-matter-coupling-plan.md` (Stage-0 dynamics program) and
`notes/emergent-nonlocality-research-plans.md`. This document develops framing
(a) from the ADM discussion: drop the claim that sampler time is physical time,
and interpret the equilibrium ensemble itself as the squared wavefunctional of a
quantum state of 3-geometry.*

---

## 1. The proposal in one paragraph

Canonical quantum gravity on a closed spatial manifold Σ = S³ describes states
by wavefunctionals Ψ[h] on **superspace** — the space of 3-geometries
Riem(Σ)/Diff(Σ) — and |Ψ[h]|² is a probability measure over 3-geometries. Our
sampler produces exactly such an object in discrete form: a normalizable
probability measure P[T] ∝ e^(−S[T]) over equilateral triangulations T of S³
modulo relabeling (= discrete superspace). We propose to take the identification
**P = |Ψ|², i.e. Ψ[T] = e^(−S[T]/2)**, seriously as a *definition* of a
candidate quantum state of 3-geometry, and to ask, quantitatively: *for which
couplings (edge-degree target, β/N, HDV) does this state look like a physical
state of gravity — a semiclassical universe with the right fluctuation
statistics?* The sampler is then not a model of dynamics at all; it is a
**preparation device** for the state and a **spectroscope** for an associated
Hamiltonian (§3). Everything physical lives in equilibrium correlations, which
is where our certified-equilibrium machinery is strongest.

## 2. The kinematic dictionary

| Canonical gravity | Our system |
|---|---|
| spatial slice Σ = S³ | fixed S³ topology (enforced by move set + init) |
| 3-metric h_ij mod diffeos | triangulation T mod isomorphism |
| superspace Riem/Diff | set of combinatorial triangulations of S³ |
| wavefunctional Ψ[h] | Ψ[T] = e^(−S[T]/2) |
| \|Ψ[h]\|² measure | equilibrium measure P[T] ∝ e^(−S[T]) |
| curvature invariants of h | degree fields: edge deficit δ(e) = 2π − deg(e)·arccos(1/3); vertex deficits |
| volume of Σ | N (pinned by facet term — a Lagrange multiplier / lapse-like chemical potential) |
| choice of state (e.g. Hartle–Hawking) | choice of couplings (edge target, β/N, HDV/N) |

Notes on the dictionary:

- **Measure conventions.** Our chain samples *labeled* triangulations uniformly
  within isomorphism class, which weights unlabeled T by N_labelings ∝ 1/|Aut(T)|
  — the standard DT measure. Aut(T) is trivial for generic large T, so this is a
  detail, but it should be stated in anything written up.
- **Fixed edge lengths** mean we access the *combinatorial* (conformal-like /
  scalar) sector of geometry only. Spin-2 (TT) metric fluctuations are not
  directly represented; see §6 (Limitations) and §5 Part V for the Regge-length
  extension that would add them.
- **The facet pin is a lapse-like multiplier:** we work canonically at fixed
  volume rather than with a cosmological-constant ensemble; converting between
  the two is a Legendre transform (as in fixed-N vs fixed-Λ DT ensembles).
- **The Hartle–Hawking analogy is structural:** Ψ_HH[h] = ∫ Dg e^(−S_E[g]) over
  4-geometries filling in h. Our Ψ[T] = e^(−S[T]/2) with *local* S is what such
  a path integral would produce if the bulk integral localized to a local
  boundary functional — we do not claim this; we treat Ψ as a variational ansatz
  whose physicality is to be *tested*, not derived.

## 3. The exact statement: the sampler is imaginary-time QM on discrete superspace

This is the technical heart, and it is exact, not analogy.

The Metropolis chain has transition matrix W(T→T′) (Pachner-move proposals ×
acceptance min(1, e^(−ΔS))), satisfying **detailed balance**:
P(T) W(T→T′) = P(T′) W(T′→T). Define D = diag(P) and symmetrize:

  **H := 𝟙 − D^(1/2) W D^(−1/2).**

Then H is symmetric and positive semi-definite, and

  **H ψ₀ = 0 with ψ₀(T) = √P(T) = e^(−S[T]/2) — the ground state of H is exactly our Ψ.**

The chain evolved for t sweeps is e^(−tH) up to the discrete-time embedding: the
sampler literally performs **imaginary-time Schrödinger evolution on discrete
superspace**, with ground state Ψ and excitation energies E_i equal to the
inverse autocorrelation times of the chain. Every τ our convergence gate already
measures is a **1/E_gap of this Hamiltonian**. Equilibrium correlators of
observables are ground-state expectation values ⟨ψ₀|Ô₁Ô₂|ψ₀⟩; time-displaced
correlators C_O(t) = Σ_i |⟨ψ_i|Ô|ψ₀⟩|² e^(−E_i t) are spectral decompositions.

Three consequences:

1. **This is the Rokhsar–Kivelson (RK) construction.** Any classical equilibrium
   measure is the ground state of an RK-type Hamiltonian built from its
   reversible dynamics. Our state is an RK/Jastrow-type wavefunction on
   triangulations — "log-local" (log Ψ = −S/2 is a sum of local terms). Known RK
   phenomenology gives sharp expectations: *equal-time* correlations of the
   quantum state equal classical critical correlations, and RK dynamics
   generically has **dynamical exponent z ≈ 2** (diffusive). The Stage-0
   measurement of z is therefore also a test of the RK character of the state.
2. **Scheme dependence is confined to excited states.** Different move sets /
   proposal distributions give different H with the *same* ground state. So the
   state Ψ is well-defined by S alone; the excitation spectrum is
   scheme-dependent except for universal content (exponents, multiplet
   structure). Claims should be organized accordingly: ground-state properties
   are physics of the state; spectral properties are physics of the (chosen)
   stochastic quantization.
3. **Wheeler–DeWitt is aspirational, not automatic.** H ≥ 0 with HΨ = 0 has
   *the form* of a constraint equation ("the state is annihilated by a positive
   Hamiltonian"), which is suggestive of WdW, but H is built from our move set,
   not derived from the Einstein–Hilbert action. The scientifically honest
   question is not "is H the WdW operator" but: **"does the ground state Ψ have
   the fluctuation structure that solutions of the gravitational constraints
   are known to have?"** That is testable (§4).

## 4. What "looks like a state of gravity" means — three testable levels

**Level 1 — Semiclassical background.** The typical geometry drawn from |Ψ|²
should be an extended, homogeneous, isotropic round S³ with self-averaging
fluctuations (δρ/ρ → 0 as N → ∞). *Status: strong evidence in the low-VDV,
edge-pinned corner (Δ_S3 → 0.02 at N=1e5, d_H → 3); crumpled/hub phases fail
it.* The couplings are the state-selection knobs.

**Level 2 — Constraint-like fluctuation statistics.** In canonical gravity the
Hamiltonian and momentum constraints eliminate the scalar gravitational mode:
long-wavelength *scalar* fluctuations of a physical state are suppressed
relative to a generic local field theory (only the two TT modes propagate).
Since our accessible sector is exactly the scalar one, this is the *right*
first quantitative test, not a compromise:

> Measure the curvature structure factor S(k) — variance of the deficit-angle
> field projected on long-wavelength modes. A generic massive-scalar ensemble
> gives S(k→0) → const (finite compressibility). A constraint-satisfying
> gravitational state should show **anomalous suppression S(k→0) → 0**.

The closest precedent is CDT's matching of spatial-volume fluctuations to a
minisuperspace wavefunction (the "de Sitter from CDT" program) — ours is the
within-slice, fixed-volume analog.

**Level 3 — Emergent isometry multiplets.** If the state describes small
fluctuations on a round S³, both equal-time correlations and the excitation
spectrum of H should organize into **SO(4) representation multiplets**: the
eigenmodes of the fluctuation two-point kernel should come in degenerate
families with the dimensions of S³ harmonic spaces, and the low-lying
relaxation modes of the chain should inherit the same degeneracy pattern.
Finding harmonic multiplets in the *relaxation spectrum of an MCMC chain* would
be a striking and, to our knowledge, novel demonstration that the stochastic
dynamics has organized itself around an emergent symmetric space.

## 5. Investigation plan mapped to existing code

**Assets already in hand:**
- `seeds/`: ~260 certified families × 32 independent draws each = i.i.d. samples
  of |Ψ|² across a 3-parameter state family (N up to 5.6e4 and growing).
- Chain time series: per-sample observables (vdv, edge_deg, hdv, num_facets) and
  degree histograms (`--record-histograms`) retained in `data/*/staging/*.csv`
  from every produce run — thousands of equilibrium trajectories.
- Geometry pipeline: `distance_distribution.py`, `roundness_analysis.py`,
  `scale_curvature.py`; `Manifold` API (facets, simplices, degree_histogram,
  save_edge_graph/save_dual_graph); checkpointed chains (`--checkpoint-every`)
  that can dump manifold snapshots along a trajectory.
- Convergence machinery that already computes τ (= 1/E_gap) per observable.

### Part I — Formalize the dictionary (paper math + sanity checks; ~days)
1. Write the H = 𝟙 − D^(1/2) W D^(−1/2) construction for our exact move set
   (Pachner + hinge moves, proposal probabilities as implemented in `sampler.d`);
   note where proposal asymmetries enter the acceptance ratio.
2. Empirical detailed-balance/self-adjointness check: time-reversal symmetry of
   cross-correlators, C_AB(t) = C_BA(t), on existing chain CSVs. Any violation
   = sampling bug, so this doubles as QA.

### Part II — State atlas: where is the semiclassical corner? (existing scripts; ~week)
Systematize Level 1 over the library: for each family, (Δ_S3, d_H, spectral
dimension d_s, correlation length of the deficit field, self-averaging exponent
of δρ/ρ vs N). Deliverable: an atlas table/figure mapping couplings → state
character (semiclassical / crumpled / hub-condensate). Mostly re-running
`roundness_analysis.py` + `distance_distribution.py` over the grid; new code:
one driver + the deficit-field correlation length (Part III's estimator).

### Part III — Level 2: curvature structure factor & constraint test (new light code; ~1–2 weeks)
1. New script `curvature_correlation.py`:
   - per seed: edge-deficit field δ(e) and vertex-aggregated deficit ρ(v);
     two-point C(r) binned by graph distance (edge and dual metrics), ensemble
     averaged over 32 reps (pattern after `distance_distribution.py`).
   - structure factor: project δρ on the low graph-Laplacian eigenmodes of the
     1-skeleton (`scipy.sparse.linalg.eigsh`, k ~ 100 lowest) and compute
     S(λ_k) = ⟨|δρ_k|²⟩.
2. Baselines for "generic": (i) the same estimator on degree-shuffled fields
   (white noise); (ii) a massive-scalar Gaussian field generated on the same
   graph. Discriminator: S(k→0) suppression relative to both baselines.
3. Sweep the atlas: does suppression appear anywhere in coupling space, and does
   it strengthen with N (constraint emerging in the continuum) or weaken
   (finite-size artifact)? **This is the flagship measurement of the program.**

### Part IV — Level 3 + spectroscopy of H (existing CSVs first, then instrumented runs; ~2–4 weeks)
1. **Variational spectroscopy from existing data:** build the cross-correlation
   matrix C_ij(t) over the logged observable basis from produce-run CSVs; solve
   the generalized eigenproblem C(t)v = λ(t)C(0)v (GEVP, as in lattice-QCD
   spectroscopy) to extract the 2–4 lowest gaps E_i and their observable
   content, per coupling and per N. Zero new simulation cost.
2. **Gap scaling / dynamical exponent:** E_gap vs N at fixed scaled couplings →
   z (prediction: z ≈ 2 if RK-generic; deviations locate the state on the
   Lifshitz axis, cf. the CDT–Hořava–Lifshitz connection).
3. **Mode-resolved runs for SO(4) multiplets (new instrumentation, small):** run
   chains at N = 1e3–3e3 dumping manifold snapshots every few sweeps
   (`--checkpoint-every` machinery or a lightweight snapshot flag); offline,
   project the deficit field of each snapshot onto the Laplacian eigenbasis of a
   reference (or per-snapshot) geometry; measure per-mode variance (equal-time)
   and per-mode τ (dynamics). Test: do eigenvalues/variances/τ's cluster into
   multiplets with S³-harmonic dimensions? Disk estimate: 1e4 snapshots ×
   ~100 KB ≈ 1 GB per run — fine.
4. **FDT check:** small objective perturbation (existing warmup-bias machinery
   applies one) vs correlation response; effective temperature must be exactly 1
   in equilibrium — deviations diagnose residual non-equilibrium.

### Part V — State design & extensions (open-ended)
1. **Inverse problem:** tune couplings to shape S(k) toward a target kernel over
   a k-window ("wavefunction engineering") using produce+certify as the inner
   loop. The certified library is the training set; the atlas (Part II) is the
   response surface.
2. **Tensor sector (the honest gap):** fixed edge lengths exclude spin-2 modes.
   The upgrade path is Regge-style length variables (or a dual-weight
   deformation) — a major D-core extension; scope it only if Parts III–IV give
   a positive signal in the scalar sector.
3. **Relational time:** with annealed matter (Stage 2 of the dynamics plan), a
   matter field can serve as a Page–Wootters clock; conditional equal-clock
   correlations define emergent evolution to compare against WdW-with-matter.
   Deferred until matter coupling exists.

### Priority order and why
1. **Part IV.1** (GEVP on existing CSVs) — zero simulation cost, immediately
   converts the seed campaign's byproduct into spectroscopy.
2. **Part III** (structure factor) — the flagship physics test; light new code.
3. **Part II** (atlas) — organizes everything; mostly reruns.
4. **Part IV.3** (multiplets) — the most novel result if it works; needs small
   instrumentation.
5. **Part V** — gated on III/IV outcomes.

## 6. Limitations and failure modes (stated up front)

1. **Scalar sector only** (fixed edge lengths): we test constraint suppression
   and background semiclassicality, not graviton physics. The Level-2 test is
   still meaningful — it is precisely the sector where gravity differs most from
   a generic field theory — but a positive result is "consistent with," not
   "demonstrates," a gravitational state.
2. **H is move-set dependent:** only ground-state (equal-time) statements are
   properties of Ψ; spectral statements are properties of our stochastic
   quantization scheme, except for universal exponents/multiplets.
3. **Euclidean, timeless:** no claim about Lorentzian dynamics; time is deferred
   to relational constructions.
4. **Inner product / measure ambiguity of WdW** is sidestepped, not solved: we
   commit to the normalizable counting measure on triangulations.
5. **Finite N:** all signatures must be tracked in N (the library's N-ladder is
   the tool); a "constraint suppression" that fades as N grows is a finite-size
   accident, one that strengthens is emergent physics.
6. **Negative results are informative:** if S(k→0) is flat-generic everywhere in
   coupling space, then log-local RK states on combinatorial superspace do not
   capture gravitational constraints — a real constraint on quantum-gravity
   model building (relevant to graphity/RK-type emergent-gravity proposals).

## 7. References

- B. S. DeWitt (1967), "Quantum theory of gravity. I. The canonical theory," Phys. Rev. 160, 1113. (superspace, WdW)
- J. A. Wheeler (1968), "Superspace and the nature of quantum geometrodynamics," in *Battelle Rencontres*.
- J. Hartle & S. Hawking (1983), "Wave function of the universe," PRD 28, 2960.
- J. J. Halliwell & S. Hawking (1985), "Origin of structure in the universe," PRD 31, 1777. (fluctuations about the HH state)
- G. Parisi & Y.-S. Wu (1981), Sci. Sin. 24, 483; P. Damgaard & H. Hüffel (1987), "Stochastic quantization," Phys. Rep. 152, 227.
- H. Risken, *The Fokker–Planck Equation* (Springer). (FP↔Schrödinger correspondence)
- D. Rokhsar & S. Kivelson (1988), PRL 61, 2376. (RK ground states from classical measures)
- E. Ardonne, P. Fendley, E. Fradkin (2004), "Topological order and conformal quantum critical points," Ann. Phys. 310, 493. (RK phenomenology, z=2, equal-time = classical critical)
- J. Ambjørn, A. Görlich, J. Jurkiewicz, R. Loll (2008), "The de Sitter universe from causal dynamical triangulations," PRL 100, 091304; PRD 78, 063544. (volume-fluctuation ↔ minisuperspace matching — closest precedent)
- J. Ambjørn et al. (2010), "CDT meets Hořava–Lifshitz gravity," PLB 690, 413.
- D. Page & W. Wootters (1983), PRD 27, 2885. (relational time)
- C. Isham (1992), "Canonical quantum gravity and the problem of time," gr-qc/9210011.
- M. Lüscher & U. Wolff (1990), Nucl. Phys. B 339, 222. (GEVP variational spectroscopy)
- T. Regge (1961), "General relativity without coordinates," Nuovo Cim. 19, 558. (length-variable extension)
