# Sampler Dynamics as a Dynamical System + Matter Coupling — Research Plan

*Program: treat the Metropolis/Pachner process as a dynamical system in its own
right; characterize its relaxation/transport; couple matter fields; test the
matter–geometry coupling against GR. Drafted 2026-07-11 in discussion with Aaron.*

## Framing

- **Fictitious time is respectable:** stochastic quantization (Parisi–Wu) —
  Euclidean QFT = equilibrium of a fictitious-time Langevin process. Our sampler
  is a discrete Langevin dynamics for 3D Euclidean geometry. Only *universal*
  dynamical content (exponents, scaling functions, FDT ratios) deserves physical
  interpretation, not raw sweep counts (move-set dependent).
- **Dimension-3 accident (statics):** 3D GR has no propagating gravitational DOF;
  Einstein eqs ⇒ flat away from matter; point mass = conical deficit ∝ mass. Our
  degree observables ARE deficit angles ⇒ the GR-coupling test is local, linear,
  finite-size-friendly.
- **Known danger:** 2D-DT lesson — matter with c>1 collapses geometry to branched
  polymers. First question of any annealed coupling: does the extended phase
  survive, and how does the phase boundary move with matter coupling g_m? The
  certified seed library is the uncoupled baseline.

## Stage 0 — Pure-geometry dynamics (no new infrastructure)

Characterize the dynamical system we already have:

1. **Relaxation spectrum.** Autocorrelation functions of vdv, edge_deg, hdv, and
   degree-distribution bin modes at several couplings inside the fluid window;
   extract slowest modes + τ's. (Fragments exist: gate τ estimates, τ∝1/g.)
2. **Dynamical exponent z.** τ vs N at fixed scaled couplings (τ ~ N^(z/3) if the
   relevant length is N^(1/3)). Diffusive vs activated/glassy relaxation.
   NOTE (ADM discussion): z also locates the emergent slice+time theory between
   isotropic (z=1) and Lifshitz-type anisotropic scaling — doubles as physics.
3. **Quench protocols.** Equilibrate at coupling A (use certified seeds), switch
   to B, record ⟨O(t)⟩ trajectories. Include Kibble–Zurek sweeps across the glass
   line g*≈0.01–0.03.
4. **Aging / FDT violation in the glass.** Two-time correlators C(t,t_w) and
   responses; effective-temperature ratio. Distinguishes true glass from
   slow-simple dynamics.
5. **Curvature transport / hydrodynamics.** Pin a small ball's edge degrees
   off-target, release, watch the disturbance spread: diffusive (r²~t)? Diffusion
   constant of curvature. Response-function half of FDT.

Starter list on current hardware/library: (1) and (3) at N=1e3–1e4 across the
existing β/N grid; (2) needs the N-ladder we already have (1e2→5.6e4 for k=0 and
ED5p0043_2); (5) at N=1e4, one family, ~days.

**Stage 0 TOOL IMPLEMENTED (2026-07-13): `scripts/quench.py`** — dedicated quench/relaxation driver. Each `--init` .mfd = one independent chain (use the 32 certified `_s0NN` replicas of a family as 32 ICs). Runs the sampler under a DIFFERENT (target) objective from t=0 (no burn-in skip, no warmup dispersal). Two-tier recording: Tier 1 dense (every `--thin` sweeps incl. t=0) observables + degree histograms → `chain_NNN.csv` (chain,sample,sweeps,num_facets,vdv,edge_deg,hdv,**energy**,accept) + `.hist.npz`; Tier 2 sparse ~`--snapshots` LOG-spaced full `.mfd` snapshots → `chain_NNN.snap_<sweeps>.mfd`; plus `manifest.json` (initial + target objectives, IC paths, snapshot schedule, git commit). ProcessPoolExecutor over chains (`--max-workers`). Snapshots round-trip (loadable full triangulations → any spatial observable post-hoc). This directly enables the quench protocols (3), aging/FDT (4 — from Tier 1), and curvature transport (5 — from Tier 2) measurements below.

## Stage 1 — Quenched matter probe (cheap; overlaps disordered-locality Plan B)

Gaussian scalar via graph Laplacian on frozen certified seeds: propagators, heat
kernels, spectral dimension. Baseline for backreaction. See
`notes/emergent-nonlocality-research-plans.md` Plan B.

## Stage 2 — Annealed coupling (engineering lift)

Joint measure exp(−S_geom − S_matter[φ,T]); alternate Pachner moves and field
updates; Pachner deltas must include the matter action (heat-bath resample fields
on created vertices). **Start with Ising per vertex** (cheap, 2D-DT benchmarks)
before continuous scalar. Measure: backreaction on ⟨vdv⟩/edge-degree
distribution; phase diagram in (edge target, β, g_m); survival of extended phase.
De-risk with a Python prototype at N≤1e3 before touching the D core.

## Stage 3 — The 3D GR test (the gem)

Localized source (mass term / pinned-field region at a marked site):
(a) deficit/degree anomaly localizes ON the source;
(b) linear in source strength;
(c) exponential return to background away from it (no long-range curvature — 3D
GR forbids propagating modes).
Pass all three = matching 3D GR matter–geometry coupling at long distance.
Long-range response would itself be a discovery (an emergent propagating mode).

## Stage 4 — Dynamical correspondence (horizon, not milestone)

FDT between matter perturbation and geometry response in fictitious time;
relational/matter-clock time (Page–Wootters); possible CDT-style foliation to
promote Euclidean statics to dynamics.

## ADM / spatial-slice reading (discussion 2026-07-11)

Aaron's alternative framing: the triangulation = a spatial slice in a discretized
ADM formulation of 4D GR, with sampler time as (some function of) coordinate
time. Assessment and literature anchors — see discussion notes:
- Canonical simplicial gravity (Dittrich–Höhn): Pachner moves on the spatial
  slice ARE the elementary time-evolution steps; tent moves. The rigorous
  scaffold for a genuine discrete-ADM dynamics (but constraint-driven, not
  Metropolis).
- Obstructions to literal ADM=MCMC: no momentum/extrinsic curvature in a Markov
  state (first-order vs second-order dynamics); stochastic relaxation to an
  attractor vs constrained Hamiltonian flow (H≈0, no equilibrium measure);
  detailed balance vs gauge dynamics.
- Defensible alternative: equilibrium ensemble = |Ψ[h]|² over 3-geometries
  (wavefunction-of-the-universe reading; sampler = stochastic quantization of a
  Hartle–Hawking-like state). Quantitative test: compare ensemble fluctuation
  spectrum to the graviton-vacuum Gaussian wavefunctional on TT modes.
- CDT connection: compare our equilibrium slice ensemble to the distribution of
  spatial slices inside 4D CDT; CDT's effective slice dynamics resembles
  Hořava–Lifshitz — Stage 0's z measurement locates us on that axis.

## Key references

- Parisi & Wu (1981), "Perturbation theory without gauge fixing," Sci. Sin. 24, 483. (stochastic quantization)
- Dittrich & Höhn (2012), "Canonical simplicial gravity," CQG 29, 115009 (arXiv:1108.1974); also Höhn (2014) on Pachner-move dynamics.
- Ambjørn, Görlich, Jurkiewicz, Loll (2012), "Nonperturbative quantum gravity," Phys. Rep. 519, 127. (CDT review; spatial-slice structure)
- Ambjørn et al. (2010), "CDT meets Hořava–Lifshitz gravity," Phys. Lett. B 690, 413.
- Hartle & Hawking (1983), "Wave function of the universe," PRD 28, 2960.
- Page & Wootters (1983), "Evolution without evolution," PRD 27, 2885. (relational time)
- Deser, Jackiw, 't Hooft (1984), "Three-dimensional Einstein gravity: dynamics of flat space," Ann. Phys. 152, 220. (3D point mass = conical deficit)
- Kazakov et al. / KPZ: Knizhnik–Polyakov–Zamolodchikov (1988); David (1988); Distler–Kawai (1989). (2D matter–gravity benchmark, c=1 barrier)
- Cugliandolo & Kurchan (1993), PRL 71, 173. (aging/FDT in glasses)
