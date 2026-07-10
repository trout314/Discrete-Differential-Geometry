# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A D language symbolic math package for discrete differential geometry, focused on efficient sampling of combinatorial n-manifold triangulations via Metropolis-Hastings with bistellar (Pachner) moves. The D core is exposed to Python via a C API and ctypes bindings.

## Build Commands

This project uses **Meson** as its build system with the **LDC** compiler. LDC must be on your PATH (set `DC=ldc2` if needed).

```bash
# Initial setup — debug build (one time)
DC=ldc2 meson setup builddir

# Build all targets
meson compile -C builddir

# Run tests
meson test -C builddir

# Build just the shared library (for Python)
meson compile -C builddir ddg_dlang

# Release build (optimized, asserts removed — use for production runs)
DC=ldc2 meson setup builddir-release --buildtype=release -Db_ndebug=true
meson compile -C builddir-release ddg_dlang

# Reconfigure after editing meson.build
meson setup builddir --reconfigure
```

### Python scripts auto-build the shared library

Scripts never load a stale `.so`. On `import discrete_differential_geometry`, the
loader (`python/discrete_differential_geometry/_dlang.py`) rebuilds the shared
library from the current source tree before loading it — a near-instant no-op via
`ninja` when nothing changed, an automatic rebuild when source changed, and a
hard error if the rebuild fails to compile (rather than silently running old
code). A directory lock serializes the build so the multi-worker drivers can't
race on it. **LDC and ninja (or meson) must be on `PATH`** wherever scripts run;
if no build tool is found the loader warns and loads the existing `.so`.

By default it loads the optimized **release** build (`builddir-release/`), so a
bare `python scripts/…` runs production-ready code with no extra flags. Env vars:

- `DDG_BUILD=debug` — load the debug build (`builddir/`, asserts on) instead of release.
- `DDG_LIBRARY=<path>` — override. A **file** is loaded as-is with **no** rebuild
  (pin an exact binary for reproducibility); a **directory** is treated as a build
  dir (kept fresh + loaded).
- `DDG_NO_AUTOBUILD=1` — skip the rebuild step (load whatever exists), e.g. on a
  box without the D toolchain, or to pin the current binary.

Both `builddir/` and `builddir-release/` must still be created once with
`meson setup` (see above); the loader compiles, it does not configure.

## Architecture

### D Core (source/)

- **`simplicial_complex.d`** — `SimplicialComplex` is the foundational type: a set of facets (maximal simplices) with vertex type as a template parameter. Supports insertion, removal, comparison, file I/O (`.sc` format), and iteration over simplices by dimension.

- **`manifold.d`** — `Manifold(dimension, Vertex)` wraps a `SimplicialComplex` and adds manifold-specific tracking: a `degreeMap` for n-simplices, `ridgeLinks` for codimension-1 faces, simplex counts per dimension, and total squared degrees. This bookkeeping enables fast Pachner moves. Supports `.mfd` file format for save/load.

- **`manifold_moves.d`** — `BistellarMove` type representing Pachner moves (center/co-center decomposition). Contains logic for determining valid moves and applying them to manifolds.

- **`sampler.d`** — MCMC step logic: `mcmcStep` function with speculative delta computation, hinge move support (dim=3), and per-move-type acceptance tracking.

- **`algorithms.d`** — Topological algorithms operating on manifolds: orientability testing, connected components, Euler characteristic, join of complexes.

- **`manifold_examples.d`** — Factory function `standardSphere` for generating standard sphere triangulations.

- **`ddg_capi.d`** — C API layer exposing D types/functions via `extern(C)` for Python ctypes bindings.

### Supporting Modules

- **`utility.d`** — `StackArray` (fixed-capacity array on the stack), range utilities, test assertion helpers (`shouldEqual`, `shouldBeSameSetAs`, `shouldBeEmpty`, `throwsWithMsg`), `dump` debug helper, and `flatDegreeInDim` lookup table for regular triangulation degrees.
- **`polygons.d`** — Polygon symmetry group computations (rotations, reflections) and triangulation enumeration.
- **`rational_number.d`** — Rational number arithmetic (templated on integer type, supports BigInt).
- **`rational_extension_vector.d`** — Vectors over rational extensions (for exact simplex point coordinates using square roots).
- **`factoring.d`** — Prime factorization utilities.

### Python Library (python/discrete_differential_geometry/)

Stable, reusable bindings to the D core. This is the public API for research scripts.

- **`_dlang.py`** — ctypes loader and declarations for the shared library.
- **`_manifold.py`** — `Manifold` class (owned handle, mutable).
- **`_manifold_view.py`** — `ManifoldView` class (borrowed handle, read-only).
- **`_sampler.py`** — `ManifoldSampler`, `SamplerParams`, `SamplerStats`.
- **`_simplicial_complex.py`** — `SimplicialComplex` class.

### Research Scripts (scripts/)

Active Python CLIs for experiments and analysis. They import the library but are not part of the installed package.

- **Seed generation:** `equilibrium_vdv.py` — the current production driver (fixed-β equilibrium sampling; `--produce` emits seed families with multi-chain R̂ gating). `grow_seed.py` — grows a small seed to large N to feed 1e6/1e7 chains. `anneal_vdv.py` — two-stage VDV annealing, kept for comparison against equilibrium.
- **Grid sweeps:** `explore_grid.py` / `produce_grid.py` — sweep the β/N × edge-target × N grid over `equilibrium_vdv.py --produce` (explore = short `--dry-run` map of the certifiable frontier, writes nothing; produce = production-length certify + copy into `seeds/`). Both are thin CLIs over the shared engine `tools/grid_sweep.py`; `--only-n/--only-edge/--only-bon` restrict to a subset or a single cell.
- **Analysis chain:** `distance_distribution.py` → `roundness_analysis.py` → `scale_curvature.py` (graph-distance distributions → roundness vs round S³ → scale-aware curvature). Each imports the previous.

Import bootstrap convention: compute `_ROOT` from `__file__`, then `sys.path.insert` `python/` (always), `tools/` (if using `seed_utils`), and the script's own dir (if importing a sibling script).

### Support Tools (tools/)

- **`seed_utils.py`** — shared helpers imported across `scripts/` and `tools/`: memory probes (`get_free_memory_gb`), seed filename encode/decode, `build_metadata_comments`, `load_seed_metadata`, git info.
- **`grid_sweep.py`** — shared engine for the grid-sweep drivers (`scripts/explore_grid.py`, `scripts/produce_grid.py`): standard-grid constants, nearest-family founding, per-cell `--produce` invocation + verdict parsing, and the `sweep()` loop with cell filters / prune predicate.
- **`reburn_batch.py` → `reburn_family.py` + `reburn_seed.py`** — the reburn pipeline: regenerate the whole seed library's metadata cheaply (orchestrator calibrates per-family burn-in, then spawns memory-capped `reburn_seed.py` workers).

### Legacy (legacy/)

Archived, superseded scripts kept for reference/reproducibility; nothing active imports them. See `legacy/README.md`.

### Data Files

- **`standard_triangulations/`** — Seed manifold files (`.mfd` format) used as starting points for sampling.
- **`data/`** — Test data files for unit tests (`.dat` format). Untracked experiment output also lands here; `data/` and `seeds_*/` are gitignored (committed unit-test fixtures stay tracked).

## Key Patterns

- Tests are co-located with source code using D's built-in `unittest` blocks. Test assertions use helpers in `utility.d`: `shouldBeSameSetAs`, `shouldEqual`, `shouldBeEmpty`, `throwsWithMsg`.
- Template-heavy D code: most types are parameterized on dimension and/or vertex type.
- The C API in `ddg_capi.d` uses opaque handles and dimension-dispatch (`switch (h.dim)`) to bridge D templates to C.
- Python library code goes in `python/discrete_differential_geometry/`; active research CLIs go in `scripts/`; shared script helpers and orchestration go in `tools/`; superseded scripts are archived in `legacy/`.
