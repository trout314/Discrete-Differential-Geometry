---
name: crossing-collider-phase2
description: "Phase-2 DONE — V=0 to 1e-14 at ALL crossing angles (36-67 deg) until touch; contact chemistry IS angular — chord-sharing repels (+4.4), angled touch is FREE (V=0 exactly, 10-vertex compounds sumZ 135-143 = the decamer family); all reversible"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T23:49:09.553Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/crossing_collider.py` —
knot A static at a rung-70 chain-A site; candidate crossing chains found by
bounded BC-stretch walks (4 dropped-vertex choices per tet) from every tet
within 1 cell of A's chord; B's washboard measured on a clean manifold
(also validates the path), then B slid through the crossing with A present.
V washboard-corrected as in [[knot-collider-phase1]]. 6 crossings tested,
angles 36-67 deg, dmin 0.11-0.23 cells. JSON scratchpad `crossing_v1.json`;
figure `data/figs/crossing_collider_phase2.png`.

**V = 0 to 1e-14 at every angle and every approach distance until the
complexes touch.** No angular interaction at range — the no-halo additivity
holds in all directions, dynamically.

**But the CONTACT chemistry is angular:**
- Same-helix chord-sharing (Phase 1, angle 0): 9-vertex product
  f=(9,23,23,8), V = +4.4 to +12.8 (rung-dependent) — repulsive core. One
  67-deg crossing that also shared 2 chord vertices gave the same 9-vertex
  product at +4.4.
- Angled touch (36-52 deg, shared 0-3 vertices): 10-vertex two-chord
  products f=(10,21,18,6)/(10,24,22,7)/(10,23,21,7)/(10,22,19,6), sumZ
  135-143, at **V = 0.000000 EXACTLY** — association is FREE. Two knots in
  contact are iso-action with two separated knots; entropy alone decides
  the bound fraction.
- ALL products unfuse by a single clean slide — fully reversible.

**Decamer connection:** the immortal thermal decamer f=(10,22,18,5), sumZ
138 ([[flicker-background]], lives 1191 sw) sits inside this
free-association family (10 vertices, two chords, sumZ 135-143). Its exact
f-vector wasn't produced by these 6 geometries but plausibly by a nearby
one or a slide-step rearrangement of one. Its immortality is consistent
with V=0 association + washboard pinning, NOT with binding energy.

**Emerging interaction picture (the full S-matrix so far):** the knot gas
is an ideal gas with angular contact rules — zero interaction at range,
free association at angled touch, +4.4..+12.8 repulsive core at
chord-sharing. All contact outcomes reversible. Combined with the washboard
([[knot-collider-phase1]]): thermal defect matter = washboard-pinned free-
associating clusters. This IS the directional interaction of the ADM
program discussion, but it is a CONTACT rule, not action-at-a-distance.

Next candidates: Phase 3 thermal kinetics (slide channel on; check
detailed-balance rates against V=0/+4.4 predictions); hunt the exact
decamer f-vector over more crossing geometries; momentum/lifted sampler.
