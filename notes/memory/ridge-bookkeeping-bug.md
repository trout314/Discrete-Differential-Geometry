---
name: ridge-bookkeeping-bug
description: "RESOLVED 2026-07-25: 'ridge corruption' was unsorted center/coCenter through the targeted-move C API — core invariant is ascending vertex order and doMove's productUnion is a SORTED MERGE; fix = canonicalize at the C boundary + invariant asserts in BistellarMove"
metadata:
  type: project
---

RESOLVED (2026-07-25). The worm_slide crashes were NOT a core logic bug.
Root cause chain, fully proven:

- Core simplex invariant: vertices in ASCENDING order (assertValidSimplex
  asserts isSorted). All internal callers honor it (moves built from
  stored sorted keys).
- doMove's oldPiece/newPiece come from productUnion = cartesianProduct +
  MERGE — a sorted merge, garbage for unsorted inputs.
- ddg_manifold_do_bistellar_move VALIDATED with sorted copies but
  constructed BistellarMove from RAW caller-order slices. An unsorted
  coCenter (e.g. apex pair (976,68) from face_apexes order) therefore
  inserted facets and dimMap keys in NON-ASCENDING order: phantom keys
  invisible to sorted lookups (degree([39,976,68])=2 while
  degree([39,68,976]) crashed — the smoking-gun probe).
- Downstream symptoms: hashmap ArrayIndexError with index size_t.max,
  "center is not a simplex" for live triangles, all state-dependent.
- Debug-build asserts did NOT catch it despite b_ndebug=false —
  assertValidSimplex apparently not firing (unresolved side mystery:
  check LDC/meson assert flags someday).

FIX (commit after 4c4d10e): ddg_capi doIt sorts center/coCenter before
constructing the move; BistellarMove constructor now asserts sortedness
(documents the invariant at the type boundary). Query wrappers (degree,
face_apexes, edge_link(_cycle), illegal_edges, validate_maps) catch
Throwable so poisoned states report instead of aborting the process.

Validated after fix: the captured 29-move crasher replays with clean
per-move value audits; crossval_moves passes; worm_slide crystal test:
knot slid +50 steps and back, crystal restored to EXACT facet equality.

Diagnostic assets kept: HashMap.validate()/capacity(),
Manifold.validateMaps, ddg_manifold_validate_maps + python binding;
scratchpad crash_move_seq.json + replay_audit.py (session d562feea).
Lesson: canonicalize at every API boundary (same lesson as cocycle
save-path labels — see [[cocycle-detachment-bug]]).
