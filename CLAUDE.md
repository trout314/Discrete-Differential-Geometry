# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Reporting conventions (research sessions)

Notation, terminology, and report standards live in **`notes/CONVENTIONS.md`**
— use its symbols and terms in all responses, reports, and figures. Hard rules
(details in that file):

1. **State reports** on crystal-derived states include the standard block:
   provenance/lineage + #chains, action + couplings (λ vs λ_EDQ — never
   conflate), f_FK, N_gr + f_G1, e\* vs ⟨e⟩, n_ill, certification status.
2. **Figures** stamp action/couplings, box/host, and #chains; mark provisional
   (uncertified) data when it matters; print the full file path.
3. **Stationarity claims** use late-window slope ± block-bootstrap σ (never
   full-series OLS alone). **Transport claims** name the channel (tracer /
   population / collective). **"Certified"** only after a two-sided R̂_q pass.
4. **Run launches**: SNAP must be a multiple of TS; state instrumentation;
   resume only from verified-consistent snapshot pairs; wrap in `caffeinate -i`.

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
- **`symmetry.py`** — `CrystalSymmetry`: the EXACT automorphism group of a
  triangulation (as explicit vertex permutations), orbit maps for
  vertices/edges/faces/tets, stabilizers, and the full BC-chain enumeration
  with chain orbits, plus the orientation-preserving subgroup **Aut+**
  (`.orientation_preserving` is itself a `CrystalSymmetry`, so every orbit /
  stabilizer method works on it unchanged) and `.is_chiral(kind, obj)`. Use
  **Aut** for anything the sampler sees (action, rates, species — all
  Aut-invariant); use **Aut+** for quantities that exist only after an
  orientation is chosen (handedness, Wilson-line spinor signs, any
  pseudoscalar). Cached per `.mfd` in a `.sym.npz` sidecar. Use this
  instead of bucketing by WL colours / degree signatures / rounded coordinates
  — those give bounds, not orbits. Also home to `develop_partial` /
  `develop_total`, the shared covering-map traversal `crystal_grains.py`
  builds grains with. **Caveat (see the module docstring): Aut(K) is the
  symmetry of the abstract complex and may legitimately exceed the physical
  crystal's space group** — do not report it as a crystallographic order
  without checking the factorisation.
- **`crystal_grains.py`** — the covering-map crystalline-grain detector: the
  single authoritative crystalline/defect identifier (see its docstring for
  why local ball-matching is wrong). CLI shim at `scripts/crystal_grains.py`.
- **`fk_skeleton.py`** — FK census primitives: `edges_from_facets`,
  `vertex_classes`, `vertex_class_census` (pass `facets=` to split Z16 into
  `Z16_Td`/`Z16_D2` and get `FK_strict`), `skeleton_stats`.
- **`link_classes.py`** — exact classification of pure-{5,6} vertex links.
  These are the duals of the fullerenes C_{20+2·n6}, so **n6 does not
  determine lk(v)**: at n6 = 4 there are two classes, the T_d Friauf
  polyhedron (the FK Z16) and a D2 isomer that is edge-legal but not
  Frank-Kasper; the FK links are exactly those whose deg-6 edges are pairwise
  **non-cofacial** (the cheap per-vertex predicate the censuses use). Carries
  the exhaustive 5/6-sphere generator (`--enumerate KMAX` via the
  `tools/link_classes.py` shim). **Never bucket links by n6 alone.**
- **`tcp_reference.py`** — Wyckoff positions of the 10 TCP crystals
  (`STRUCTURES`) + `build_t3_triangulation` (periodic Delaunay → T³ quotient)
  + `reference_frac_positions`. The validate-and-save CLI stays at
  `scripts/tcp_reference.py`.
- **`chain_select.py`** — BC-chain class selection + provenance
  (`ChainClasses`, `chain_for_run`): makes the chain class an explicit,
  recorded input instead of an accident of tet 0.
- **`seed_utils.py`** — seed filename encode/decode, `.mfd` header metadata,
  memory probes, and the flattened provenance-history helpers
  (`build_metadata_comments` — lineage is mandatory).

### Research Package (python/ddg_lab/)

Shared research-grade machinery, one tier below the stable package: importable
without `sys.path` hacks, but still experiment-shaped (hard-coded couplings,
repo-relative data paths). 19 modules:

- defect core: `defect_state` (incremental defect bookkeeping over the
  accepted-move event stream), `event_replay`, `crystal_flicker`
- worm machinery: `worm_moves` (exact V-preserving move arithmetic),
  `worm_helix` (`bc_orbit`), `worm_slide` (the knot-slide), `worm_deg4_slide`,
  `tip_retract_search` (the Patch machinery), `dressed_generators`, `fk_moves`
- f0 sector: `f0_worm` (the extended-ensemble worm channel), `link_planner`
- crystals: `tcp_melt`, `dopant_pairs`, `graft_signature`
- geometry: `steiner_geodesic`, `heat_geodesic`, `sixhundred_cell`
- seed pipeline: `grid_sweep`

The old locations under `scripts/`, `scripts/defect_dynamics/`, and `tools/`
are **shims**: importing hands back the package module (`sys.modules` swap),
and running them as scripts re-executes the module under
`runpy.run_module(..., run_name="__main__")`, so both legacy `from X import Y`
(with the script dirs on `sys.path`) and every command line keep working.

### Research Scripts (scripts/)

Active Python CLIs for experiments and analysis. **`scripts/INDEX.md` is the
authoritative manifest**: every script classified by program and status
(shim / active / tool / validator / dormant), with docstring, last-commit
date, and import fan-in. It is generated by `tools/script_index.py` — edit
the CLASSIFICATION map there and re-run (`--check` fails on drift; the CLI
smoke tests derive their sweep list from it). Highlights:

- **Seed generation:** `equilibrium_vdv.py` — the current production driver (fixed-β equilibrium sampling; `--produce` emits seed families with multi-chain R̂ gating). `grow_seed.py` — grows a small seed to large N to feed 1e6/1e7 chains. `anneal_vdv.py` — two-stage VDV annealing, kept for comparison against equilibrium.
- **Grid sweeps:** `explore_grid.py` / `produce_grid.py` — sweep the β/N × edge-target × N grid over `equilibrium_vdv.py --produce` (explore = short `--dry-run` map of the certifiable frontier, writes nothing; produce = production-length certify + copy into `seeds/`). Both are thin CLIs over the shared engine `tools/grid_sweep.py`; `--only-n/--only-edge/--only-bon` restrict to a subset or a single cell.
- **Symmetry:** `move_site_census.py` — exact 2→3 move-site classes (= Aut face
  orbits) with curvature ledgers. `bc_chain_census.py` — every BC chain class of
  a crystal: length, screw, exact holonomy, winding orbit, tet/vertex orbit
  coverage, chirality under Aut+. Both over `symmetry.CrystalSymmetry`.
- **Analysis chain:** `distance_distribution.py` → `roundness_analysis.py` → `scale_curvature.py` (graph-distance distributions → roundness vs round S³ → scale-aware curvature). Each imports the previous.

Import convention: **new code imports the packages** —
`discrete_differential_geometry` and `ddg_lab` — with a single
`sys.path.insert` of `python/` computed from `__file__`. The historical
per-script bootstrap (inserting `scripts/`, `tools/`, and the script's own
dir to import sibling scripts) survives in existing scripts and is kept
working by the shims; don't add new script-to-script imports.

### Support Tools (tools/)

- **`script_index.py`** — generates `scripts/INDEX.md` (the classified
  manifest); the curated program/status map lives in its CLASSIFICATION dict.
- **`reburn_batch.py` → `reburn_family.py` + `reburn_seed.py`** — the reburn pipeline: regenerate the whole seed library's metadata cheaply (orchestrator calibrates per-family burn-in, then spawns memory-capped `reburn_seed.py` workers).
- `seed_utils.py`, `link_classes.py`, `chain_select.py`, `grid_sweep.py` —
  shims over the packages (see above).

### Legacy (legacy/)

Archived scripts kept for reference/reproducibility; nothing active imports
them (grep-enforced during the 2026-08 cleanup). Two layers, both documented
in `legacy/README.md`: flat files (early superseded scripts) and
**`legacy/<program>/` directories** — whole experimental programs whose
questions are answered (catalysis, percolation, ecmc, colliders,
chord-bilocal, run5h-passes, diagnostics, …), each with a verdict section
pointing at the `notes/memory/` file that records the result. Archived
scripts may import the packages and active scripts; never the reverse.

### Data Files

- **`standard_triangulations/`** — Seed manifold files (`.mfd` format) used as starting points for sampling.
- **`data/`** — Test data files for unit tests (`.dat` format). Untracked experiment output also lands here; `data/` and `seeds_*/` are gitignored (committed unit-test fixtures stay tracked).

## Key Patterns

- **Maximalist recording**: research scripts drive the sampler through
  `discrete_differential_geometry.recording.record_run` (or `Recorder.step`
  for custom loops), never bare `ManifoldSampler.run()`. Recording is
  unconditional and cheap (O(1) counters + D-side defect census per chunk,
  mid+final snapshots); equilibrium gating/certification is post-processing
  on the recorded `.rec.jsonl` series, deliberately decoupled.
- **Python tests** (`tests/`, pytest): the fast tier (`pytest -m "not slow"`,
  ~75 s) covers the packages, exact validator batteries, and a compile sweep
  of every manifest script; the slow tier (`pytest -m slow`) adds
  sampler-scale validator runs and a `--help`/import-rot sweep over every
  runnable CLI. Validator scripts are wrapped in subprocesses by
  `tests/test_validator_scripts.py`, so their campaign-scale CLI use still
  works. `tests/conftest.py` rebuilds missing reference crystals; tests
  pinned to unreproducible snapshots skip when the data is absent. **Run one
  pytest at a time per checkout** — concurrent suites interfere.
- Tests are co-located with source code using D's built-in `unittest` blocks. Test assertions use helpers in `utility.d`: `shouldBeSameSetAs`, `shouldEqual`, `shouldBeEmpty`, `throwsWithMsg`.
- Template-heavy D code: most types are parameterized on dimension and/or vertex type.
- The C API in `ddg_capi.d` uses opaque handles and dimension-dispatch (`switch (h.dim)`) to bridge D templates to C.
- Where code goes: stable, tested library code in
  `python/discrete_differential_geometry/`; shared research-grade machinery in
  `python/ddg_lab/`; research CLIs in `scripts/` (registered in
  `tools/script_index.py`); orchestration helpers in `tools/`; closed
  programs archived to `legacy/<program>/` with a verdict README entry.
