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

One-off or evolving Python scripts for experiments and analysis. These import the library but are not part of the installed package. Free to experiment, duplicate, or delete.

### Data Files

- **`standard_triangulations/`** — Seed manifold files (`.mfd` format) used as starting points for sampling.
- **`data/`** — Test data files for unit tests (`.dat` format).

## Key Patterns

- Tests are co-located with source code using D's built-in `unittest` blocks. Test assertions use helpers in `utility.d`: `shouldBeSameSetAs`, `shouldEqual`, `shouldBeEmpty`, `throwsWithMsg`.
- Template-heavy D code: most types are parameterized on dimension and/or vertex type.
- The C API in `ddg_capi.d` uses opaque handles and dimension-dispatch (`switch (h.dim)`) to bridge D templates to C.
- Python library code goes in `python/discrete_differential_geometry/`; research scripts go in `scripts/`.
