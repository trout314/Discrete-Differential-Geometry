# legacy/ — archived scripts

Superseded scripts kept for reference and reproducibility of earlier results.
They are **off the active surface**: nothing in `scripts/` or `tools/` imports
them. They may need path/environment fixups to run today (the `tools/`-based
ones have had a `sys.path` insert added so `seed_utils` still imports).

Archived 2026-07-10 (additions 2026-07-22; program archives 2026-08-16).
See git history for full provenance.

## Superseded seed-generation (→ `scripts/equilibrium_vdv.py --produce`)

- **`ramp_vdv.py`** — single-stage ramp of the VDV coefficient until acceptance
  hits zero. Earliest VDV-minimization approach.
- **`generate_seeds.py`** — multi-chain (K=4) R̂-gated seed generation without a
  VDV penalty. Its R̂-gating idea lives on inside `equilibrium_vdv.py`.
- **`equilibrate_seed.py` / `equilibrate_batch.py`** — March-era grow +
  fixed-sweep equilibration pipeline (batch → worker). Structural ancestor of
  the current `tools/reburn_batch.py → reburn_seed.py` pair, superseded by
  `grow_seed.py` + `equilibrium_vdv.py --produce`.
- **`generate_seed_params.py`** — emitted the `seed_params*.json` job tables that
  fed `equilibrate_batch.py`. (Those JSONs were removed from git; regenerate here
  if ever needed.)

## Retired C++ `manifold_sampler` convergence flow

- **`convergence_analysis.py`** — Gelman–Rubin R̂ / ESS CLI over CSV chain files
  written by the old C++ `manifold_sampler`. Thin wrapper over the still-live
  `discrete_differential_geometry.convergence` library.
- **`validate_seeds.py`** — seed-stationarity validation via fresh reference
  chains + R̂/Geweke; references the old `standard_triangulations/` layout.
- **`convergence_test.params` / `run_convergence_test.sh`** — launcher for
  parallel C++ sampler chains.

## Interactive viz

- **`live_histogram.py`** — live matplotlib degree-histogram window. Standalone;
  run with `PYTHONPATH=python`. Untouched since the earliest days of the project.

## One-shot migrations (already applied to the seed library)

- **`migrate_hdv_naming.py`** — renamed `_HDV_{coef}` seed families to scaled
  `_HDVs_{coef/N}` tokens (mirror of the VDV_→VDVs_ migration).
- **`migrate_history.py`** — back-filled the flattened `history` provenance
  field into pre-provenance seed `.mfd` headers.

## Superseded defect-dynamics diagnostics (2026-07 campaigns)

- **`decorr.py`** — run5h mobility/glassiness diagnostic. Both measurements live
  on: Jaccard-vs-lag in `scripts/defect_dynamics/mgas_analyze.py`, identity
  tracking (strictly deeper) in `run5h-passes/pass1_kinematics.py`.
- **`validate_centroids.py` / `validate2.py`** — one-off audits that certified
  the raw tree-lift centroid protocol (99.4 % exact-zero steps; glitch class
  0.09 %). Verdict recorded in `run5h-passes/pass1_kinematics.py`'s docstring; the protocol
  now also runs live inside `mobile_gas.py` (`cents` field).
- **`knotgas_sk.py`** — early "is the sub-Bragg plateau the frozen-gas Poisson
  floor?" probe; superseded by the shell-resolved, pooled `sl_verdict.py`.
- **`reconcile.py` / `final_fig.py`** — boundary-artifact control and headline
  figure of the finished curvature-length-scale campaign (results in memory /
  notes; regenerate from snapshots if needed).

## Superseded ADM/lapse prototype

- **`local_lapse_proto.py`** — Stage-1 lapse estimate from per-sweep facet
  diffs; superseded by `scripts/local_lapse.py` (Stage 2, exact per-move
  counters).

---

# Program archives (2026-08-16 cleanup, Phase 2)

Whole experimental programs whose scientific questions are ANSWERED, one
subdirectory per program. Each script keeps a working import bootstrap
(repo-relative paths; shared machinery still imports from `scripts/`,
`scripts/defect_dynamics/`, and the `discrete_differential_geometry`
package). Verdicts live in the memory files under `notes/memory/`.

## `catalysis/` — flicker catalysis of barrier crossings

**RETRACTED.** No evidence of catalysis: net ΔS is blind to it, and the
strict chord channel (U ≡ 0, per-move Metropolis) *cannot* lower a barrier —
a real mechanism would need move-bundling into one acceptance. See
`notes/memory/flicker-catalysis.md`.
`catalysis_search.py` / `catalysis_wide.py` (path search),
`catalysis_audit.py` (channel audit), `flicker_catalysis_barrier.py`
(barrier A/B).

## `percolation/` — chord-channel relaxation A/B

**NULL.** The chord channel does not accelerate relaxation: 1.6x *worse* on
total work (pure overhead), plain-work effect ~0. Gotcha for re-analysis:
chord episodes don't increment `totalTried`. See
`notes/memory/percolation-ab-test.md`. `perc_ab.py` (honest work
accounting), `perc_pilot.py` (density pilot).

## `ecmc/` — lifted ECMC flight of deg-3 chords

**LIFT LOSES.** Blob dispersal is transport-limited and nonlocal moves win
16-43x — but the unlifted diffusive slide beats the flight 6-9x at matched
budget. The wall stopping a flight at a defect complex is 100% the m² term
(pins + geometry contribute exactly 0). See `notes/memory/ecmc-blob-ab.md`,
`lifted-ecmc-transport.md`, `flight-contact-barrier.md`.
`ecmc_flight.py` (the assembled driver), `ecmc_ab.py` + `_analyze` +
`_figure` (the A/B; its blob-construction helpers moved to the active
`cimp_scan.py`), `contact_census.py` + `contact_census_merge.py` (what a
flight contact could hand off to).

## `colliders/` — knot-knot scattering

**S-MATRIX MEASURED, INTERACTION IS CONTACT-ONLY.** Phase 1 (same BC helix):
V(s) = 0 to 1e-14 until s = 4 contact; contact repulsive (+4.4..+12.8),
reversible; washboard 17-51 explains caging. Phase 2 (crossing chains): V = 0
at ALL angles until touch; contact chemistry is angular (chord-sharing
repels, angled touch is free — the decamer family). Plan C (harvested real
targets): pristine S-matrix transfers, free association only. The exhaustive
sweep (`knot_smatrix_sweep.py`) was killed after its relkey dedup proved
wrong both ways (exact Aut keys: 680 -> 271 classes). See
`notes/memory/knot-collider-phase1.md`, `crossing-collider-phase2.md`,
`harvest-collider-planC.md`, `smatrix-sweep-running.md`.
`knot_collider.py`, `crossing_collider.py`, `crossing_intrinsic_label.py`,
`harvest_collider.py`, `knot_smatrix_sweep.py`, `compound_lab.py`.

## `chord-bilocal/` — bilocal carriers over the strict chord channel

**BLOCKED by 5 structural impossibility results** (the master narrative is
`notes/memory/bilocal-program-saga.md`): both carriers conserve f₀ exactly,
cost factorizes iff supports are vertex-disjoint, transport is blocked by
the flat off-tube umbrella. The strict-chord *channel* itself survives (its
regression tests `twosided_chord.py` / `agg_knobs_test.py` stay on the
active surface). `bilocal_factorize.py`, `bilocal_range.py`,
`chord_tube.py` (umbrella from replayed episodes), `ep_cost.py`.

## `run5h-passes/` — the run5h defect-kinetics campaign analyses

**RESULTS RECORDED** (`notes/memory/defect-kinetics-run5h.md`,
`defect-travel.md`, `lifetime-vs-charge.md`): monomers blink, fused
multimers immortal, Q ≈ -5 per knot; ALL complexes caged (subdiffusive,
exponent 0.63), long-lived ones return home; no monotone charge->lifetime
law (the Σ_Z = 70 rung anomaly). `pass1_kinematics.py`,
`pass3_microdynamics.py` (`pass2_structure.py` stays ACTIVE — reused by the
amorphous census), `survivor_anatomy.py`, `release_run.py`,
`longlived_structures.py`, `lifetime_charge.py`, `defect_travel.py`.

## `diagnostics/` — resolved-bug investigations

All bugs RESOLVED; kept for the record of *how*.
`diag_state_divergence.py` / `diag_readmethod_mutation.py` (the DefectState
incremental-divergence bug of 2026-07-27: a read-side query method mutated
state), `leak_probe.py` (the embedded-GC rt_init leak,
`notes/memory/embedded-gc-rt-init.md`), `linkage_bias.py` (frame-overlap
linkage does NOT invent deaths).

## Single-file retirements

- **`f0-sector/wf_vs_planner.py`** — is the cheap proposal weight
  w(f) = [f has an edge of degree <= 4] as good as the planner? Result
  folded into the f0 program (`notes/memory/f0-surgery.md`).
- **`fpkmc/fp_dock_angles.py`** — registry-gauge dock angles, explicitly
  DEPRECATED in-file; superseded by the intrinsic (development-frame)
  `fp_dock_census_intrinsic.py`, which stays.
- **`hu-statics/defect_halo.py`** — the no-halo verdict, DECISIVE: the melt
  equals the pristine crystal beyond 0.4 cells; no screening cloud; static
  two-knot potential rigorously zero beyond contact
  (`notes/memory/no-halo-verdict.md`).
