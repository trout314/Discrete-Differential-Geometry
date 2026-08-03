---
name: smatrix-sweep-running
description: "STOPPED 2026-07-25 (ran <1 day, ~3 of 680 classes). Its dedup key was WRONG both ways — rewritten onto exact Aut keys 2026-08-02; A-classes 680 -> 271 and --start-class numbering CHANGED, old partials not resumable"
metadata: 
  node_type: memory
  type: project
  originSessionId: c395a1df-09fb-4d44-bb80-ddc781dd7a5b
  modified: 2026-08-03T00:30:26.911Z
---

**STATUS: STOPPED.** Relaunch #3 went up 2026-07-25 late-late and was KILLED
by user decision the same night (see [[fpkmc-design]]): the exhaustive scope
came out at ~40-68 days, and the call was to optimize (likely D-side) before
re-undertaking. It reached roughly **3 of 680 A-classes**. Partials in the
scratchpad: `smatrix_full.json`, `sweep_pilot*.json`.

## Its dedup key was wrong in BOTH directions (found 2026-08-02)

`relkey` deduped configurations by a rounded relative-position hash of a BALL
of surrounding crystal vertices — keying on the ENVIRONMENT rather than on the
object. Two chords can sit in translation-congruent balls while oriented
differently within them. Measured on R m4, chainA (class 4, L=3252):

| | A-window classes |
|---|---|
| relkey, on the 810 stride-4 windows | 680 |
| EXACT, under the full Aut | **271** |

The two partitions **cross-cut** — neither refines the other:

- it missed every point-group equivalence (2.5x redundant work), and
- it wrongly MERGED: of the 130 windows it dropped, **87 are not
  Aut-equivalent** to the representative kept.

So the docstring's scope claim — "missed point-group symmetries only cause
harmless redundancy, never a gap" — was false in its second clause, and the
declared exhaustiveness with it. Because the run died at ~3 classes there is
almost no affected output; the pilot findings below are per-collision facts,
not class-census claims, so they stand.

## What changed 2026-08-02

All three uses (A-window classes, B-geometry dedup, washboard cache) now use
`CrystalSymmetry.config_key` — an exact equivalence certificate anchored on a
distinguished frame, resting on the free action of Aut on frames. `relkey` is
deleted. See [[crystal-symmetry-group]].

**RESUME GOTCHA: `--start-class N` means something different now.** The
partition is 271 classes, not 680, and the ordering is a fresh randomized
shuffle. Old partial JSONs cannot be merged with new ones by class index.
Start clean.

Runtime is genuinely uncertain, do NOT assume 2.5x: A-classes dropped 680 ->
271, but with A PINNED there is no residual symmetry (Aut acts freely on
frames, so any g fixing A's frame is the identity), which means distinct B
geometries the ball hash used to fuse are now all kept — more collisions per
class. One factor moved decisively in our favour, one against. Only a timed
run settles it, and the ~40-68 day estimate that triggered the kill should be
re-measured before deciding anything.

The sweep also now records its BC chain class (it walks class 4 of 14 by
default via the legacy `F[0]` seed) — see [[crystal-symmetry-group]] for why
that mattered.

## Pilot findings (still valid — per-collision, not census)

Same-chain rows reproduce Phase 1 exactly ((9,23,23,8) at +12.8, unfuses).
V_contact spectrum seen so far, NOT saturated: {0, +4.6, +7.8, +8.8, +10.4,
+10.6, +11.2, +12.2, +12.6, +12.8, +16.0}. Products include (9,24,25,9) — the
thermal 9-mer target of [[harvest-collider-planC]] IS a direct knot-knot
collision product — plus (9,20,18,6) and (9,21,19,6). The exact decamer
f=(10,22,18,5) has NOT appeared as a direct product in phases 1-2 or pilots.

History: #1 OOM-killed (Manifold load/del leaked ~70MB; fixed by the
Ledger-unwind rewrite, commit e2e1eb0). #2 killed for a COVERAGE HOLE the user
spotted — forward-only stretch walks silently dropped ~97% of chains through
the ball, including the same-chain head-on channel. #3 has bidirectional
stretches (walk_bidir), explicit same-chain geometries both orientations, and
a transparent-by-theorem prefilter (~57% of geometries dismissed by no-halo
without running).

Note the pattern across all three kills: every one was a COVERAGE or
CORRECTNESS defect in the enumeration, not a physics error. Treat any future
"exhaustive" claim here as a hypothesis to be tested, not an invariant.

## Also decided 2026-07-25

The SEGMENT-JUMP sampler (user's idea, sharpened): rejection-free heat-bath
jumps of a knot within the clear segment of its chain (segments bounded by
sites within contact range of any defect), exact for equilibrium by the
no-halo theorem + symmetric segment construction; washboard cache provides the
Boltzmann tables. Build as `scripts/defect_dynamics/segment_jump.py` driver
channel BEFORE any D-core work. Lifted/event-chain variant deferred to the
momentum program.

Analysis when resumed: fingerprint -> product/V table; saturation curve (new
fingerprints vs collisions); V spectrum by shared vertex count; compare the
product family against the thermal species table ([[flicker-background]],
[[crossing-collider-phase2]], [[knot-collider-phase1]] for washboard
conventions).
