---
name: fpkmc-design
description: FPKMC/heat-bath/event-chain design written (notes/FPKMC_DESIGN.md); S-matrix sweep KILLED pending D-side optimization; D-first performance policy; milestones M0-M6
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-26T03:06:42.732Z
---

2026-07-25 late: S-matrix sweep #3 KILLED by user decision — optimize
(likely D-side) before re-undertaking (~40-68 day exhaustive scope was the
trigger). Design doc written: **notes/FPKMC_DESIGN.md** — read it before
touching any of this.

Key decisions:
- Three samplers, one infrastructure: HB (segment heat-bath, exact
  equilibrium for the FULL model, build first), FP (first-passage KMC —
  the user's "carry until near a defect" idea completed with splitting
  probabilities + first-passage times so it is exact; frozen-background v1
  exact, dressed v2 approximate), EC (lifted event-chain, momentum sector,
  deferred).
- USER POLICY: code goes on the D side if there is even a small chance it
  needs to be performant. Planned C API: ddg_manifold_chain_walk,
  ddg_manifold_site_survey (washboard + clean-slot counts in one call),
  ddg_manifold_segment_scan, ddg_sampler_do_bistellar (targeted move WITH
  sampler bookkeeping). Linear algebra (tridiagonal birth-death chains,
  ~200x200) stays numpy.
- The sweep is a CONSUMER of the same D API (M5: rewrite sweep on
  chain_walk+site_survey+slide_at before relaunch).
- Critical invariant I1 (verify FIRST, V1): clean-slot symmetry
  n(j,j+4)_fwd == n(j,j+4)_bwd — licenses Boltzmann detailed balance on
  the washboard chain. R2: proposal constant nu must be measured (V2)
  against sampler.d internals before any kinetics claim.
- M0 requires re-deriving the splitting/FPT formulas symbolically (SymPy)
  before coding — the doc's phi-form is flagged as needing verification.

Uncommitted at time of writing: sweep fixes in knot_smatrix_sweep.py
(bidirectional walks, same-chain channels, transparent-by-theorem
prefilter, washboard cache, randomized class order), knot_collider.py
recs_out param, notes/FPKMC_DESIGN.md. Sweep partials in scratchpad
(sweep_pilot*.json, smatrix_full.json ~3 classes) — analysis-worthy:
V spectrum {0,+4.6,+7.8,+8.8,+10.4,+10.6,+11.2,+12.2,+12.6,+12.8,+16.0},
same-chain rows reproduce Phase 1. See [[smatrix-sweep-running]] (now
STOPPED), [[knot-collider-phase1]], [[crossing-collider-phase2]],
[[no-halo-verdict]].
