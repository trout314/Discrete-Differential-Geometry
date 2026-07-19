# Phased Pachner Moves: promoting the RK sampler toward quantum gravity

*Drafted 2026-07-19 (Aaron + Claude discussion). Companion to
`notes/combinatorial-lapse-plan.md` and `notes/ensemble-as-wavefunction-program.md`.
Status: speculation + concrete first measurements. See memory
`project-quantum-promotion` for the running log.*

## Starting point: what the RK construction already provides, and what it lacks

The equilibrium sampler defines a Hilbert space ℓ²(triangulations), a stoquastic
Hamiltonian H = −π^{-1/2} W π^{1/2} ≥ 0 (W = Markov generator, π = e^{−S}), and
an exact ground state Ψ₀ = √π with HΨ₀ = 0. This is a quantum theory whose
quantumness is inert: all-positive amplitudes (no interference), imaginary-time
dynamics only, and H is an evolution operator in external MC time rather than a
constraint. Promotion to a quantum-GR candidate needs three things:

1. **Phases** (sign structure → interference, Lorentzian character);
2. **Relational time** (H reinterpreted as a Wheeler–DeWitt constraint — note
   HΨ₀ = 0 already holds exactly — with time read from matter-clock
   correlations à la Page–Wootters; the defect complexes of the doping program
   are natural clock subsystems);
3. **A demonstrated gapless spin-2 branch in spec(H)** (the graviton), with the
   scalar channel gapped/projected — an empirical property of the H we already own.

This note is mostly about (1), plus the measurement plans for (3).

## Where phases live: gauge theory on the flip graph

Deform H_{T→T'} = −t e^{iθ(T→T')}. Constraints and structure:

- **Hermiticity**: θ(reverse move) = −θ(move) — orientation reversal is
  complex conjugation (time reversal built in).
- **Gauge freedom**: |T⟩ → e^{iφ(T)}|T⟩ shifts θ by coboundaries; physical
  content = **holonomies around closed loops of the flip graph** (Pachner
  cycles). The theory of phases is a lattice gauge theory on the space of
  geometries; its field strength lives on the loops.

## The central identification: move histories are 4D spacetimes, so θ is a Lagrangian

Each 3D Pachner move glues one 4-simplex onto the slice; an MC history is a
triangulated 4D cobordism (combinatorial-lapse observation); a flip-graph loop
bounds a closed 4D "bubble". Distant moves commute, so many orderings trace the
same 4D triangulation, and phase coherence across orderings holds **iff θ is
additive over 4-simplices, i.e. iff a local 4D Lagrangian density exists**.
Lagrangian locality is thus the anomaly-cancellation condition, not an input.
Whatever path-dependence remains lives on the genuine relations of the Pachner
groupoid — the natural home of the theory's quantum "curvature".

Two consistent Lagrangian families:

**(a) Dynamical — Lorentzian Regge / CDT class.** θ(move) = Lorentzian Regge
action of the glued 4-simplex, with CDT asymmetry α (spacelike slice edges vs
timelike interpolating edges). Each move type (1↔4, 2↔3, hinge) gets a definite
θ(α, couplings); inverses conjugate automatically. Combined path weight
e^{−S_E + iS_L}: the classical objective supplies a Euclidean regulator/measure,
the phase supplies interference — the "allowable complex metrics" structure;
Kontsevich–Segal-type admissibility bounds are the candidate selection
principle for which deformations define sensible theories.

**(b) Topological — gravitational θ-angles.** θ ∝ discrete topological density
(combinatorial Euler / Pontryagin-like terms) of the glued 4-simplex. Loop
holonomies are then topological invariants of the closed bubble: quantized,
automatically coherent, deformation-robust — discrete analogs of
Holst/Nieh–Yan/Pontryagin terms. Condensed-matter precedent (RK models +
plaquette flux → chiral liquids) predicts the qualitative physics: chiral,
time-reversal-breaking vacua with persistent probability currents —
impossible for any detailed-balance sampler, hence the unambiguous signature
of genuine quantumness.

**Lapse composition.** ADM weighting appears naturally:
θ(move) = N(x_move) · L_Regge(4-simplex), with N the measured lapse field
(local acceptance density). The lapse program and the phase program compose
into a single discrete N∫√g R.

## Data-gathering plans

Ordered by feasibility; 1–4 need no sign-problem technology.

1. **Chirality susceptibility from existing chains (linear response in θ).**
   First-order corrections to Ψ₀ under a θ-deformation are imaginary-time
   correlators of the phase operator, computable in the stoquastic ensemble.
   Measure orientation-odd / parity-odd move-current correlation functions in
   existing (or cheap new) runs → susceptibility of the vacuum to a
   gravitational θ-angle. Deliverable: does the geometry "want" chirality?

2. **Spin-2 spectroscopy (promotion criterion 3; independent of phases).**
   MC autocorrelators are imaginary-time Green's functions:
   C_O(τ) = Σ_n |⟨n|O|Ψ₀⟩|² e^{−E_n τ}. Channel-resolve by operator symmetry:
   scalar (curvature density / volume) vs transverse-traceless-like combinations
   of the deficit field on long-wavelength modes. Deliverable: gap structure by
   channel — a gapless spin-2 branch with gapped scalar = "the graviton is in
   the spectrum"; the single most important empirical target of the program.

3. **Exact diagonalization on enumerated micro-superspaces.** Enumerate the
   flip graph of small S³ triangulations; diagonalize the phased H exactly.
   Numerically exact interference between geometries (relative phase between
   two Pachner paths = action difference of their cobordisms); spectroscopy of
   the deformed vacuum vs θ and α; watch for level crossings / chiral ground
   states. A toy quantum gravity, honest at small N.

4. **Entanglement structure of Ψ₀.** Replica/SWAP estimators in the classical
   sampler give Rényi entropies across cuts: area law + corrections — the
   holographic sanity check on the vacuum.

5. **Sign-problem frontier (later).** Lefschetz thimbles / complex Langevin
   for intermediate sizes once 1–4 motivate specific θ values.

## Companion data plans: dopant-interaction characterization (matter sector)

Recorded here for completeness (details in memory `project-t3-library`):
complex-resolved statistics everywhere (complexes, not vertices, are the
particles); u_eff(r) = −ln g_cc(r) from equilibrated multi-replica ensembles;
mass-action spectrum ΔF_k (needs multiple equilibrated states — single-state
n₁ = 1 is insufficient); **host screening-halo profiles** around isolated
complexes (tests dressed-neutrality mechanism behind the hyperuniformity
result; halo range controls the flat-doping OCP verdict); identity-tracked
kinetics (birth/capture/escape rates, complex mobility D(size), detailed-balance
check = equilibrium vs slow precipitation); per-complex charge census
(measured: Q ≈ −0.65 rad/vertex uniform — charged droplets, host-screened).
