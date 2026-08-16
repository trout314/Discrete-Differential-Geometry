---
name: fk-vertex-link-classification
description: "n_6 determines lk(v) uniquely only up to n_6=3; n_6=4 has TWO legal links (Td Friauf + a D2 isomer), and non-cofacial six-edges is what caps FK at Z16"
metadata: 
  node_type: memory
  type: project
  originSessionId: a48cb7a2-9bd0-449c-81d5-b74518ecd7cb
  modified: 2026-08-13T15:23:33.682Z
---

In a closed 3-manifold triangulation with **all edge degrees in {5,6}**, lk(v) is a
triangulated S² with all degrees in {5,6} ⇒ n_5 = 12 always, deg v = 12 + n_6,
#tets at v = 20 + 2·n_6. These spheres are **exactly the duals of fullerenes
C_{20+2n_6}** (min degree 5 ⇒ 3-connected planar ⇒ Steinitz/Whitney, so the graph
determines the sphere).

Exhaustively enumerated (own backtracking generator, not quoted tables), agreeing
with fullerene isomer counts:

| n_6 | V  | classes | Aut orders | deg-6 induced subgraph of lk(v) |
|-----|----|---------|-----------|----------------------------------|
| 0   | 12 | 1       | 120       | — (Z12 icosahedron)              |
| 1   | 13 | **0**   | —         | **IMPOSSIBLE** (dual of the Grünbaum–Motzkin "no C22" gap) |
| 2   | 14 | 1       | 24        | 2 isolated (Z14, D6d)            |
| 3   | 15 | 1       | 12        | 3 isolated (Z15, D3h)            |
| 4   | 16 | **2**   | 24, 4     | isolated (Z16 Friauf, T_d) **vs perfect matching (D2 isomer)** |
| 5   | 17 | 3       | 20,4,4    | min 2 induced edges              |
| 6   | 18 | 6       | 4,2,12,2,12,6 | min 3 induced edges          |

**Key structural fact.** Adjacency of two degree-6 vertices in lk(v) ⟺ the two
six-edges at v are **cofacial** (span a triangle with v). In Z12/Z14/Z15/Z16-T_d
the six-edges at a vertex are pairwise NON-cofacial; from n_6 = 5 upward every
class has ≥2 cofacial pairs. So the FK list Z12/Z14/Z15/Z16 is exactly the set of
5/6-links with pairwise non-cofacial six-edges — the non-cofaciality condition is
what caps n_6 at 4 (verified by enumeration through n_6 = 8; classes
1,0,1,1,2,3,6,6,15 = the C_{20+2n_6} fullerene isomer counts). The minimum
cofacial-pair count rises monotonically past 4 — 2,3,5,6 at n_6 = 5,6,7,8 — so
the cap is not a near miss further up.

**DELIVERED as `ddg.link_classes (shim: tools/link_classes.py)`** (enumerator + fast detector + CLI), wired
into `fk_skeleton.vertex_class_census(eu, edeg, V, facets=…)` — optional arg, so
all ~10 existing callers are unchanged; with `facets` it adds `Z16_Td`, `Z16_D2`,
`FK_strict` — and into `defect_dynamics/amorph_census.py` as the `Z16*` column
plus `strict%`. Detector = per-vertex count of cofacial deg-6 pairs, one pass
over triangles.

**Validation:** 450 links re-derived from scratch on 3 crystals, 0 mismatches;
all 16 reference crystals give Z16_D2 = 0 (no false positives); on a flagged
amorphous state all 54 flagged links are graph-isomorphic to the enumerated D2
sphere (|Aut| = 4) and all 15 Z16 links to Friauf T_d (|Aut| = 24) — the cheap
predicate is exactly equivalent to the full isomorphism test.

**RESULT — the D2 link is NOT rare, it DOMINATES.** In the `fk_amorph` anneal
states non-FK D2 sites outnumber genuine Friauf Z16 by ~3–4× (a15_mu2_noratchet:
54 D2 vs 15 T_d; z_mu2_sl20: 73 vs 12). So ~78% of the n_6 = 4 population
previously counted as Z16 was never Frank–Kasper; the correction to f_FK is ~5
percentage points on top of the hub correction. Anything that scored those
anneals by Z16 content needs re-reading — see [[fk-amorphous-program]],
[[crystal-library-gas-campaign]].

**Consequence for legality bookkeeping:** "no illegal edges" is *strictly weaker*
than FK/TCP legality, and n_6 = 4 is the first place they diverge. See also
[[crystal-grains-tool]], [[six-web-gauge]] (the six-web is triangle-free and has
no cofacial pairs in any genuine FK state).
