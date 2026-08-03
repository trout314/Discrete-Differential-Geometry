---
name: crystal-grains-tool
description: scripts/crystal_grains.py — covering-map crystalline-grain detector (registry-aware)
metadata: 
  node_type: memory
  type: project
  originSessionId: eaaf2503-7af1-42b2-a8c0-88f714d19b53
  modified: 2026-07-21T02:52:58.198Z
---

2026-07-20: built `scripts/crystal_grains.py` to identify connected components
of the vertex graph consistent with a piece of a reference TCP crystal — genuine
crystalline grains, not merely local FK order. Motivation & design were worked
out with Aaron; see [[tcp-r-c15-defect-state]].

**Why:** `defect_census` (r=1 + weak heal) OVER-accepts (isolated liquid
icosahedron counts); `crystal_match` (r=2 exact WL cert, no heal) UNDER-accepts
(one perturbed edge in a ~60-vertex ball breaks it, so a 62%-ordered melt read
0.0% crystal). Crystallinity is translational REGISTRY (non-local), not local
coordination.

**Method — combinatorial covering map ("development"):** map sample tets to
reference tets preserving decorations (FK class + 6 hinge degrees); seed one tet
(fix a registry), BFS across shared faces FORCING neighbour images; a decoration
mismatch = grain boundary. A vertex is INTERIOR-crystalline iff its full star
develops consistently into one grain. Phase-exact via consistency (C14 can't
pass as C15 — stacking conflicts) yet tolerant (a covered vertex needs only its
own link clean, so grains grow up to defects).

**Key fix (subtle):** single-valuedness is checked MODULO the unit-cell site,
not the exact ref vertex. `tcp_reference` labels vertices `cell*ns + site`, so
`site = vertex % ns` (ns = atoms/cell, from `tcp_reference.STRUCTURES[name][1]`).
Without this, a sample grain wrapping its torus reaches a vertex at ref images
differing by a lattice translation and was wrongly cut (perfect crystal read
only 57%). With it: perfect crystal = 100%.

**Validation (all correct):** perfect C15/R → 100%, one grain, right phase,
excludes all other phases. Liquid null (fFK 0) → 0% even at min-size 1 (clean
floor ⇒ K can be small). Partial melts (where crystal_match r=2 = 0.0%): real
grains, monotone in melt depth — r_lam0.03/0.1/0.3 → 462 / 264 / (60+19, TWO
disconnected R grains) interior verts. Fast (seconds).

CLI: `python scripts/crystal_grains.py STATE.mfd [--ref c15 r] [--min-size K]`.
Reuses crystal_match.best_refs, fk_skeleton, dopant_pairs. Committed ba8e5b9.

**FIXED & VALIDATED 2026-07-20 (loose-heal landed, tests + real-data green;
uncommitted):** the tool had SPREAD defects by one ring — a single Pachner move
on a perfect crystal gave defects = touched UNION their 1-ring. Cause: interiority
required the FULL star to develop, but a defect vertex's degree/class change makes
_tet_ok fail for ALL its tets, so good neighbours had an incomplete star and got
flagged.

The fix has THREE parts (all in scripts/crystal_grains.py):
1. **FK-shell GATE** `_locally_fk(st,v)`: v interior only if every spoke edge (v,w)
   has hinge degree 5 or 6 with EXACTLY twelve 5s (the mandatory 12 fivefold
   disclinations of any all-5/6 triangulated 2-sphere; link-node-degree(w) ==
   hinge(v,w), so it's purely spoke degrees — no reference templates, and it
   deliberately OMITS neighbour-typ and link-edge degrees, which shift when a
   neighbour is defective and were what leaked the defect outward). A move breaks
   the shell for exactly the touched (subdivided → degree-3 spoke; flipped →
   loses a fivefold spoke, count5≠12), never their neighbours.
2. **Loose-heal registry** in `interior_vertices`: among v's ASSIGNED star tets,
   all must agree on ONE grain + ONE site; UNASSIGNED (defect-adjacent) tets are
   don't-cares. Signature UNCHANGED so tests call it as-is.
3. **`_drop_subsumed_grains`** (post-pass in find_grains): un-labels a grain whose
   vertex set ⊆ a single larger grain's — a redundant re-seed of the same crystal
   in a symmetry-related frame (happens when a defect locally blocks the flood, esp.
   on tiny tori). Needed so the seam of a 2-defect crystal isn't split across two
   frames. KEY: only subset-of-ONE-larger drops; holonomy-fragmented cross-phase
   grains (C15-vs-R = 56 similar-sized non-nesting patches) are left intact, so
   phase discrimination (which lives in the develop-time fragmentation) is preserved.

Diagnostics that pinned the design: C15-vs-R assigns ALL tets but into 56 grains
(no fused-frame site-conflict at m=2 → a global-frame conflict test does NOT
discriminate; the develop-time fragmentation does). The 2-defect R seam = one
1264-vertex grain + a 5-vertex re-seed subsumed by it.

Validation: all 13 tests in tests/test_crystal_grains.py PASS (perfect=∅,
1→4/2→3=touched only, two-move=two clusters, C15-vs-R<0.25V). Real-data melt of
T3_R_m2 (60 sweeps): perfect 100%/1 grain; lam 0.05/0.15/0.30/0.50 →
56.9/42.9/23.6/15.3% interior, largest grain 720/521/230/87 — monotone, tolerant
of partial order (crystal_match r=2 read these as 0%). Ready to commit
crystal_grains.py + tests/test_crystal_grains.py together.

**INTEGRATED into defect_census + K CALIBRATED 2026-07-20 (uncommitted):**
`scripts/defect_census.py` now has `--method {registry,local}` (default registry).
registry = `defect_cores_registry()` delegates to crystal_grains (build_struct /
find_grains / interior_vertices / grain_components), crystalline = vertices in
grains >= `--min-size`, defects = complement, then the existing r=1 machinery
reports Z-class composition + complex clustering on the cores. local = the old
one-sided accretion (kept: it uniquely catches wrong-Wyckoff-site NATIVE-class
dopants). Vertex ids align because both tools relabel via np.unique(F). CRYSTALS
(tcp_melt) and STRUCTURES (tcp_reference) share lowercase keys (r/c15/a15/sigma/c14);
ns from STRUCTURES[ref], guarded by V%ns==0 else exact-ref-vertex fallback.

**K calibration (R m2, V=1272):** melt-depth grain-size sweep is BIMODAL — one
dominant surviving domain (569/345/156/35 as lam 0.1→0.55) far above a tail of
COINCIDENTAL grains. Deep-liquid replicas (lam 0.7-1.0, fFK 0.32-0.43) cap the
coincidental grain at ~10-18. Real domains stay >=35 even half-melted. Ceiling is
extreme-value ~2.5*ln(V) (registry match per grain-extension is independent), so
it grows only ~ln V: ~30-36 at V~1e5-1e6. DEFAULT `--min-size 30` clears the
V~1e3 ceiling with margin and reads liquid null -> 0% / 0 grains, light/mid melt
-> single real domain. For V >> 1e5 raise K or calibrate per-run vs a matched
melted null (max_grain(null) x ~1.5-2).

Remaining follow-up: dislocation classification (infinite-lift / Burgers-vector,
deferred).
