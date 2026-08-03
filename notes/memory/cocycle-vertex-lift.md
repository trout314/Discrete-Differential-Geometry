---
name: cocycle-vertex-lift
description: "Vertex positions are now maintained incrementally in D (O(1), fixed gauge); no more tree-lift rebuilds or gauge glitches"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T16:03:13.568Z
---

As of 2026-07-25 (commit 706b7e5) the cocycle's VERTEX LIFT is maintained in D,
not re-derived in Python.

**Why it works.** omega is a min-imaged real displacement, not a bare winding
indicator (`build_from_positions`: omega(u->v) = p(v) - p(u) - M k(u,v)), so
its integral is a genuine position. The invariant to keep is

    pos(v) - pos(u) == omega(u->v)  (mod boxM)   on every edge

and EVERY move preserves it with at most one assignment: 2->3 and the 4-4
diagonal telescope through the forced value; 3->2 and 4->1 only delete; 1->4 is
the sole vertex creation and the existing gauge choice omega(a->w)=0 puts w
exactly at a, so `pos[w] = pos[a]`.

**API** (opt-in; nothing existing changed): `s.enable_cocycle_positions(box)`,
`s.vertex_positions(verts)` -> (n,3), `s.check_cocycle_positions()`. D side:
`cocycleSeedPositions` (one BFS), `cocyclePosProblems`, and
`CocycleState.pos/boxM/posEnabled`.

**Payoff.** Per-sample at N3=57992: old `read_cocycle + tree_positions` 143 ms
vs new `vertex_positions(chord)` 0.008 ms (~19000x). More important than speed:
the gauge is FIXED for the run, so the "spanning-tree gauge glitch" that
pass1_kinematics filters (steps > 1 cell) CANNOT occur -- don't add glitch
filtering to new code that uses this path.

**GOTCHA -- vertex labels are RECYCLED.** A 4->1 frees a label and a later 1->4
reuses it elsewhere. Verified exactly against the event stream: under forced
churn (445 1->4, 763 4->1) all 75 apparent "position changes" were reused
labels, 0 drift. So *label identity is not object identity across a death* --
relevant to [[worm-sampler-program]]'s edge-key trackers (knot_worldlines,
deg4_moves).

**Testing trap.** An equilibrated chain never changes V, so a validation run at
the production volume pin exercises NONE of the 1->4/4->1 path. Force churn by
setting num_facets_target well off the current N3.

**Next bottleneck** (measured, not guessed): the O(N) species scan
`knot_chords_now` at 800 ms/sample, ~5x the sampler's 168 ms/sweep. Positions
are no longer the cost. See [[worm-sampler-program]].
