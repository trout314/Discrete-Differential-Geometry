# A T³ seed library parallel to the S³ library

*2026-07-16. Status: adopted plan (Aaron). Infrastructure Stage T0 implemented
same day: `scripts/tcp_reference.py`. Companion context:
`notes/ensemble-as-wavefunction-program.md`, the QC/phason/TCP program
(memory: project-phason-tcp), and the certified S³ library workflow
(`scripts/seed_queue_runner.py`, `scripts/equilibrium_vdv.py`).*

## 0. Why a torus library

1. **Topology-controlled twin.** Every local observable we study (degree
   censuses, defect pair correlations, hyperuniformity estimators, lapse
   fields) should converge to the same thermodynamic limit on S³ and T³ if it
   is genuinely local. Deviations at matched (N, objective) measure
   finite-size topology effects — currently an uncontrolled systematic in
   everything we do. The exact combinatorial identities are topology-blind
   (χ = 0 and the link sum rule Σ(6−deg) = 12 hold for any closed 3-manifold),
   so the whole analysis stack runs unchanged.
2. **The flat reference without frustration.** S³ must carry net positive
   curvature, which is why our "near-flat" pin is offset (5.0043 at N=1e4).
   T³ admits globally flat metric structure: the flat pin 5.1043 is the honest
   target, and the TCP/QC physics (whose natural composition IS ≈5.10–5.11)
   lives there without the global-curvature tax.
3. **Crystalline and quasicrystalline endpoints are natively periodic.** The
   Frank–Kasper phases (A15, C14, C15, σ, ...) and QC approximants exist as
   *exact* T³ triangulations (Stage T0 below). That gives: perfect reference
   states at the EDV floor, ordered-side melting experiments (phase-boundary
   bracketing + hysteresis vs our disordered-side annealing), exact
   hyperuniformity baselines, and a σ-approximant ladder toward dodecagonal
   quasicrystal order.

## Stage T0 — TCP reference triangulations  [DONE 2026-07-16]

`scripts/tcp_reference.py`: periodic Delaunay of published Wyckoff positions
→ simplicial T³ quotient → validated .mfd.

- Structures: **A15** (Cr₃Si; q̄=5.1111; Z12₂Z14₆), **C15** (MgCu₂; 5.1000;
  Z12₁₆Z16₈), **C14** (MgZn₂; 5.1000; Z12₈Z16₄), **σ** (CrFe, 30 atoms,
  refined internal parameters; 5.1089; Z12₁₀Z14₁₆Z15₄ — the classical
  dodecagonal-QC approximant). All four validate EXACTLY at m=3: all edges
  {5,6}, literature Z-census reproduced per cell, χ=0, orientable, link sum
  rule clean.
- Method guards: sites perturbed ~1e-6 ONCE before tiling (consistent
  degeneracy resolution across periodic images — do NOT use per-point Qhull
  joggle); 2-cell margin; exactly one tet per translation class (centroid in
  fundamental domain); hard checks on f₃ count and facet uniqueness.
- Scale: N ≈ 10⁵ builds in ~2 s (A15 m=13); N ≈ 10⁶ ~ half a minute (m=28);
  10⁷ is minutes (m=60). Arbitrarily large flat crystals are free.
- Key identity for structure selection: q̄ = 6 − 12/CN. The classical phases
  sit at 5.1000 (C14/C15), 5.1089 (σ), 5.1111 (A15) — all within ±0.007 of
  the flat 5.1043. Our edge-pin dial resolves *between* real TCP phases.

## Stage T1 — pipeline plumbing (small; ~1 session)

The certify path is already topology-agnostic: `equilibrium_vdv.py` takes
`--topology` (pass-through to naming/metadata; default S3), gate features are
all local, Pachner moves preserve topology, and `build_seed_filename` prefixes
whatever topology string it is given (→ `T3_N1e4_...` naming for free).

Work items:
1. **Roots in `standard_triangulations/`**: commit small validated quotients
   (A15/C15/σ at m=2..3) as the T³ analogs of `standardSphere`. (m=2 needs a
   simplicial-quotient check; use m=3 if not.)
2. **`seed_queue_runner.py` / `tools/grid_sweep.py`**: add a `topology=` token
   to queue lines (default S3); pass `--topology` through to produce; restrict
   founding/nearest-family search to same-topology seeds (filename prefix
   match is enough); grow from a T³ root instead of the standard sphere when
   topology=T3.
3. **Sanity invariant**: χ=0 and orientability do NOT distinguish S³ from T³.
   Add a cheap occasional topology check to the T³ path (e.g. H₁ rank via
   GF(2) edge-cycle rank on a snapshot, or simply trust move correctness +
   provenance). Decision: trust moves + provenance tag; audit tool optional.
4. **`reburn`/analysis tooling**: nothing to do (all local); only
   `load_seed_metadata` consumers that assume the `S3_` prefix need a glance.

## Stage T2 — T³ roots and base families (~1 day of queue time)

1. Disorder a TCP crystal into a generic T³ triangulation: run the sampler at
   the base objective (edge pin only, no variance terms) from an A15 quotient
   until all crystal order is erased (track fFK → 0, Z-census → generic).
   This is ALSO the first melting experiment (see T4).
2. Grow/shrink to the standard N ladder (1e3, 1778, 3162, 5623, 1e4, …) with
   the existing grow machinery from the m=3 roots.
3. Produce base families (no variance penalty) at the three standard pins →
   `T3_N*_1e-1_ED*_2` names, certified by the standard gate.

## Stage T3 — the mirrored grid (queue-driven, ~week-scale background)

Mirror the S³ standard grid exactly, via queue lines with `topology=T3`:
- pins {5.0043, 5.1043, 5.2043} × k=2 × VDVs {1e-3, 2e-3, 4e-3, 8e-3, 1e-2}
  × N {1e3 … 1e4} to start (extend N as the S³ library does);
- the certified HDV combos that proved TCP-relevant (VDVs_2e-3 + HDVs_0p101 /
  0p228) at the three pins;
- priorities set so the flat pin (the physics target on T³) fills first.
Direct deliverable: topology-controlled twins for every measurement made this
week (census, defect g(r), hyperuniformity ladder, anneals).

## Stage T4 — the science the torus unlocks

1. **S³/T³ locality test**: matched-(N, objective) equal-time observables;
   any significant difference is a finite-size topology effect — quantifies
   systematics in the whole S³ program.
2. **Melting of real TCP phases**: heat A15/C15/σ at the matched couplings;
   locate loss of FK order from the ordered side; hysteresis loop against the
   disordered-side annealing results ⇒ is the TCP phase thermodynamically
   stable at finite g, or only the g→∞ limit? (The single most important
   open question from the annealing campaign.)
3. **Exact hyperuniformity baseline**: the window/spectral estimators on
   perfect A15/C15/σ (should → 0; validates the estimator stack), then on
   melted/annealed T³ states — the hyperuniform-network-vacuum hunt without
   S³'s global curvature tax.
4. **Approximant ladder → QC**: σ and larger dodecagonal approximants as
   competing initial conditions at the flat pin; does the equilibrium/annealed
   measure prefer approximant order, QC-like order, or the network glass?
5. **Energy ground truth**: evaluate the exact objective on A15 vs C14 vs C15
   vs σ as a function of pin — which real phase do our penalties select?
   (Cheap; do first. The pin resolves between phases at the 3rd decimal.)

## Costs / risks

- Queue time dominates (same as S³ library economics); T³ certification cost
  = S³ cost at matched N.
- σ/C14 use experimental internal parameters — combinatorics validated exact,
  but keep the parameter values in the .mfd comments for provenance.
- Naming collision risk: none (T3_ prefix).
- The only true topology-specific risk is a move-implementation bug silently
  changing topology; mitigated by provenance + optional H₁ audit (T1.3).
