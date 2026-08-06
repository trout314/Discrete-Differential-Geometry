---
name: crystal-grafting-program
description: "Grafting lumps of one TCP crystal into another via exact decorated-boundary-surface matching: L1/L2/L3 match hierarchy (L3 = zero-cost seam), graft_signature.py library validated, C36<->C14 slab control"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1fc2ced4-2271-4a60-ba17-bc68ad6d30ad
  modified: 2026-08-06T14:38:27.842Z
---

**Goal (Aaron, 2026-08-06):** graft a piece of one TCP crystal into another by
finding surfaces with isomorphic decorated graphs. Reframed: search over
LUMPS (tet sets) and match their decorated boundary surfaces.

**Match hierarchy (each exactly checkable, proofs in graft_signature.py docstring):**
- L1: simplicial iso of boundary 2-complexes → glued object is automatically a
  closed 3-manifold (disc ∪ disc = S² links).
- L2: + per-edge (total degree, lump-interior tet count d_in) preserved → ALL
  edges of the graft stay in {5,6}.
- L3: + per-vertex (cn, n6, n_interior_edges, n_interior_6edges) preserved →
  disclination web reconnects (no Z-class change, no broken lines; "lines
  can't end" is a theorem in the pure-{5,6} sector, see fk_skeleton.py).
  **An L3 graft has exactly zero local energy cost under any degree-only
  action — the seam is invisible to the sampler.**
- Vertex Z classes alone (the original idea) are NOT sufficient; edge splits +
  web continuity are the real conditions.

**Physics motivation:** the pin gap (e* − e_native) is extensive with opposite
signs across the library ([[crystal-library-gas-campaign]]); a mixed-sign L3
composite at the right volume fraction has ZERO net pin gap — no pure crystal
except R can do that.

**Infrastructure:**
- `scripts/graft_signature.py` — CrystalContext, lump_boundary (decorations +
  surface-manifoldness checks), EXACT canonical certificates by
  combinatorial-map BFS over all start darts (no WL hashing — exact, and the
  canonical vertex orders double as φ candidates), match_phis, graft surgery
  (fresh interior ids + compaction), validate_facets (link sum rule assert,
  {5,6} census, broken-line count, χ, orientability). `--selftest` = Z12
  star-swap on in-memory c15 m3: PASSES (12 decorated link symmetries,
  nontrivial φ, Z census preserved).
- `scripts/graft_c36_control.py` — positive control: C36 = hc Laves polytype
  interleaving C14-like/C15-like stacking. Enumerates basal slabs (cuts
  between atomic z-planes; c36 m4 has 16 planes/cell, c14 m4 has 8),
  certificates, within-crystal groups (same boundary + different interior =
  stacking-swap candidates), cross-crystal join, then performs + validates the
  best c14→c36 graft, saving to data/grafts/.
- Orientation gotcha: gluing along both slab boundary tori needs φ
  orientation-compatible on BOTH simultaneously or the result is
  non-orientable — handled by iterating match_phis candidates until
  validation passes.

**CONTROL RESULT (2026-08-06, m=4, kmin 3 kmax 12):** c36: 156 slabs → 6
distinct L3 boundary classes; c14: 78 slabs → the SAME 6 classes (identical
certs) — c36 and c14 share their entire basal boundary-class set, and
L1=L2=L3=78 cross pairs (for Laves basal slabs the surface complex already
determines the decorations). Every class contains slabs with different
interiors (stacking-swap candidates). Performed a THICKNESS-CHANGING graft:
12-plane c36 slab (1632 tets) out, 3-plane c14 slab (544 tets) in, first φ
validated: all {5,6}, 0 broken lines, χ=0, orientable, Z census pure
Z12:Z16 = 896:448 = 2:1 → a new 7616-tet mixed-stacking Laves polytope-box,
saved data/grafts/T3_C36m4_graftC14_0.8905-1.6405_f37616.mfd (+ summary JSON
c36_control_m4.json). ~1.5 min/crystal census, single core.

**Next steps:** C15 needs (111) cuts (layer coordinate u = (x+y+z) mod m, cubic
box); then cross-library graft compatibility matrix; then thermal stability of
a grafted seam under the minimal 3-term action; mixed-sign pin-gap composite.
