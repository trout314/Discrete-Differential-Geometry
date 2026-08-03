---
name: build-state
description: D core build status; the macOS unittest link failure is FIXED and meson test now passes
metadata: 
  node_type: memory
  type: project
  originSessionId: eaaf2503-7af1-42b2-a8c0-88f714d19b53
  modified: 2026-07-28T01:39:57.194Z
---

On this macOS/aarch64 box with LDC 1.32.2 + ninja:

- The shared library builds cleanly in both `builddir/` (debug) and
  `builddir-release/` (release). `python -c "import discrete_differential_geometry"`
  works — the auto-build loader rebuilds release on import as a no-op.

- **`meson test -C builddir` PASSES (1/1, ~17 s).** It BROKE 2026-07-27 and was
  FIXED same day. Cause: the always-on ridge guard commit ecddcc0 put an impure
  `fprintf`+`abort` lambda inside `insertFacet`, making the manifold-construction
  path impure, which the six explicit `pure @safe unittest` blocks (manifold.d
  L872, L1301, L1387, L1522, L1551, L1562) reject. (My earlier "pre-existing / LDC
  tightening" guess was WRONG -- it was our own guard.) FIX: laundered the guard's
  purity with the same cast trick as `reclaimCapacity` -- new static helpers
  `_abortInvariant` (`@trusted pure nothrow`) -> `_abortInvariantImpl` in
  manifold.d; guard still always-on + aborts, insertFacet stays pure/@safe.
  fprintf/abort are already nothrow @nogc, so only PURITY needed laundering.
  Superseding the old note that unittests couldn't link. Two fixes, in
  commit a6e4b7d:

  1. `meson.build` built rpath args as `'-L=-rpath=' + dir`, which reaches the
     linker as the single token `-rpath=PATH`. That joined spelling is a **GNU
     ld extension**; Apple's ld64 only takes `-rpath PATH` as two arguments and
     fails with `ld: unknown option: -rpath=…`. Fix: emit two `-L=` args,
     `'-L=-rpath', '-L=' + dir`. Works on both linkers. (The shared-library
     target never hit this — it uses `-link-defaultlib-shared=false` and
     doesn't take `d_link_args` at all.)
  2. `rational_number.d:43` asserted `"%.18f"` of pi equals the **x87 80-bit**
     expansion, which needs `real.mant_dig == 64`. On AArch64 D's `real` IS
     `double` (mant_dig 53), so it can never match. Now selects the expansion
     by `real.mant_dig`. Not a code bug — the exact-fraction assert passes.

  Also raised the meson `test()` timeout to 300 s: the suite is unoptimized and
  the 30 s default left no headroom.

- Gotcha when adding D unittests: `findProblems` rebuilds a whole
  SimplicialComplex and does connectivity + per-hinge link-circle checks — it
  cost ~13 s when called per-probe in a loop. Use `mfd.validateMaps` (cheap
  hash-map audit) inside loops and `findProblems` once at the end.

**ASSERTS ARE NOT IN THE PYTHON-LOADED LIBRARY (2026-07-27).** The
`shared_library('ddg_dlang', ...)` target in meson.build HARD-CODED
`d_args : [..., '-O3', '-release']` — so BOTH `builddir/` (buildtype
debug) and `builddir-release/` ship a `-O3 -release` lib with asserts +
bounds-checks STRIPPED. `DDG_BUILD=debug` therefore never gave asserts
in the loaded .so; asserts only ever ran in the `unittest` executable
(d_unittest:true, -O0) on tiny test manifolds. That's the "debug
tolerable at ~50 tets" recollection — the unittest suite, not a
production-scale asserts-on lib (which never existed).

MEASURED at N3=57996 (200-sweep pure-sampler bench, this box):
- release / "debug" lib (both -O3 -release): 6.17 / 5.96 sweeps/s (~3%
  = noise; the debug lib is byte-equivalent to release, asserts OFF).
- **checked (-O3, asserts+bounds ON): 5.18 sweeps/s = only ~16% slower.**
So asserts at production scale cost ~16%, NOT the feared ~20x. The slow
debug build was a mirage.

FIX ADDED (uncommitted): meson_options.txt option `lib_ndebug` (default
true = current -release behavior); meson.build now builds
`lib_dargs = ['-link-defaultlib-shared=false','-O3']` and appends
`-release` only if lib_ndebug. A CHECKED build:
`DC=ldc2 meson setup builddir-checked -Dlib_ndebug=false && meson compile
-C builddir-checked ddg_dlang`. Load it via
`DDG_LIBRARY=$(pwd)/builddir-checked/ddg_dlang.dylib` (still TODO: wire a
`DDG_BUILD=checked` branch into python/.../_dlang.py `_pick_build_dir`,
which currently only knows debug/release).

Also this session: added an ALWAYS-ON (release-safe) invariant guard in
manifold.d insertFacet (~line 447) replacing `assert(ptr.degree==2)` --
fires loudly (fprintf+abort) if a ridge/triangle reaches degree 3
(StackArray!(Vertex,2) overflow) or gets duplicate apexes ("two
triangles same link vertices"). Motivated by the reaction_census
edge-degree corruption hunt (tets stay correct, edeg drifts at ~sw8000);
prime suspect is a release-silenced D-core invariant. Pattern to adopt:
promote critical structural asserts to always-on guards.

See [[worm-sampler-program]] for the slide-rollback unittest this unblocked,
and [[regenerate-tcp-references]], [[tcp-r-c15-defect-state]],
[[reaction-census-campaign]].
