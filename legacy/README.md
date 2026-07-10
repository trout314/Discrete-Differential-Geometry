# legacy/ — archived scripts

Superseded scripts kept for reference and reproducibility of earlier results.
They are **off the active surface**: nothing in `scripts/` or `tools/` imports
them. They may need path/environment fixups to run today (the `tools/`-based
ones have had a `sys.path` insert added so `seed_utils` still imports).

Archived 2026-07-10. See git history for full provenance.

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
