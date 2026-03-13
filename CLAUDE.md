# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A D language symbolic math package for discrete differential geometry, focused on efficient sampling of combinatorial n-manifold triangulations via Metropolis-Hastings with bistellar (Pachner) moves.

## Build Commands

This project uses **Meson** as its build system with the **LDC** compiler. LDC must be on your PATH (set `DC=ldc2` if needed).

```bash
# Initial setup — debug build (one time)
DC=ldc2 meson setup builddir

# Build all targets
meson compile -C builddir

# Run tests
meson test -C builddir

# Build a specific target
meson compile -C builddir manifold_sampler

# Release build (optimized, asserts removed — use for production runs)
DC=ldc2 meson setup builddir-release --buildtype=release -Db_ndebug=true
meson compile -C builddir-release manifold_sampler

# Reconfigure after editing meson.build
meson setup builddir --reconfigure
```

## Architecture

### Core Types (source/)

- **`simplicial_complex.d`** — `SimplicialComplex` is the foundational type: a set of facets (maximal simplices) with vertex type as a template parameter. Supports insertion, removal, comparison, file I/O (`.sc` format), and iteration over simplices by dimension.

- **`manifold.d`** — `Manifold(dimension, Vertex)` wraps a `SimplicialComplex` and adds manifold-specific tracking: a `degreeMap` for n-simplices, `ridgeLinks` for codimension-1 faces, simplex counts per dimension, and total squared degrees. This bookkeeping enables fast Pachner moves. Supports `.mfd` file format for save/load.

- **`manifold_moves.d`** — `BistellarMove` type representing Pachner moves (center/co-center decomposition). Contains logic for determining valid moves and applying them to manifolds.

- **`algorithms.d`** — Topological algorithms operating on manifolds: orientability testing, connected components, Euler characteristic, join of complexes.

- **`manifold_examples.d`** — Factory function `standardSphere` for generating standard sphere triangulations.

### Supporting Modules

- **`utility.d`** — `StackArray` (fixed-capacity array on the stack), range utilities, test assertion helpers (`shouldEqual`, `shouldBeSameSetAs`, `shouldBeEmpty`, `throwsWithMsg`), `dump` debug helper, and `flatDegreeInDim` lookup table for regular triangulation degrees.
- **`polygons.d`** — Polygon symmetry group computations (rotations, reflections) and triangulation enumeration.
- **`rational_number.d`** — Rational number arithmetic (templated on integer type, supports BigInt).
- **`rational_extension_vector.d`** — Vectors over rational extensions (for exact simplex point coordinates using square roots).
- **`factoring.d`** — Prime factorization utilities.

### Applications (source/applications/)

Each application has its own `main()` guarded by `version (unittest) {} else`. Meson builds each as a separate executable target.

- **`manifold_sampler.d`** — The main application. Runs Metropolis-Hastings sampling of manifold triangulations. Reads parameters from a `.params` file (key-value format). Outputs `.mfd` manifold files and optional CSV data, edge graphs, and dual graphs.
- **`edge_graph.d`** — Extracts edge graph from a simplicial complex file.
- **`dual_graph.d`** — Extracts dual graph from a manifold file.
- **`surface_cross_S1.d`** — Computes the product of a surface with S^1.

### Data Files

- **`standard_triangulations/`** — Seed manifold files (`.mfd` format) used as starting points for sampling.
- **`data/`** — Test data files for unit tests (`.dat` format).

## Key Patterns

- Tests are co-located with source code using D's built-in `unittest` blocks. Test assertions use helpers in `utility.d`: `shouldBeSameSetAs`, `shouldEqual`, `shouldBeEmpty`, `throwsWithMsg`.
- Template-heavy D code: most types are parameterized on dimension and/or vertex type.
- Application `main()` functions are wrapped in `version (unittest) {} else` so they don't conflict with the test runner.
- The manifold sampler reads parameters from a plain-text `.params` file (see `source/applications/manifold_sampler.params` for an example).
