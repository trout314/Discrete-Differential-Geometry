# legacy/ — archived scripts

Superseded scripts kept for reference and reproducibility of earlier results.
They are **off the active surface**: nothing in `scripts/` or `tools/` imports
them. They may need path/environment fixups to run today (the `tools/`-based
ones have had a `sys.path` insert added so `seed_utils` still imports).

Archived 2026-07-10 (additions 2026-07-22). See git history for full provenance.

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
  tracking (strictly deeper) in `pass1_kinematics.py`.
- **`validate_centroids.py` / `validate2.py`** — one-off audits that certified
  the raw tree-lift centroid protocol (99.4 % exact-zero steps; glitch class
  0.09 %). Verdict recorded in `pass1_kinematics.py`'s docstring; the protocol
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
