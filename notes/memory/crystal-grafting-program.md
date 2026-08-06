---
name: crystal-grafting-program
description: "Grafting lumps of one TCP crystal into another via exact decorated-boundary-surface matching: L1/L2/L3 match hierarchy (L3 = zero-cost seam), graft_signature.py library validated, C36<->C14 slab control"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1fc2ced4-2271-4a60-ba17-bc68ad6d30ad
  modified: 2026-08-06T16:00:49.850Z
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

**C15 CROSS RESULT (graft_c15_cross.py):** epitaxy arithmetic: the (111) wrap
of a cubic m-box is a 2m x 2m hex torus of the Laves net, so **c15 m=2 is the
epitaxial partner of c14/c36 m=4**. c15 m2 (111): 12 slabs → 4 L3 classes.
Cross joins: c15×c14 L3=216 (L1=344), c15×c36 L3=432 (L1=688) — first case
where **L1 > L3**: some surfaces match topologically but NOT with decorations
(would be defective if glued; c36↔c14 basal had L1=L3). Grafted a cubic C15
(111) slab into BOTH hex hosts, first φ each: all {5,6}, 0 broken lines,
χ=0, orientable, Z12:Z16=2:1. Both grafts CHANGE total V by −192
(768→576 c14 host; 1536→1344 c36 host) with a zero-defect seam —
**vertex-count-changing grafts do NOT force non-FK seam vertices** (Aaron's
suspicion refuted within the Laves family; thickness-changing slab swaps
change V freely at L3). Saved: data/grafts/T3_C14m4_graftC15m2_f33264.mfd,
T3_C36m4_graftC15m2_f37616.mfd, c15_cross.json.

**BALL-GRAFT SEARCH (Aaron's ask: dV != 0 ball grafts, e.g. 3-vertex star
union → 2-interior filling).** Reuses the FK-move-search formalism
([[fk-move-search]]): ball graft = pair of FK-legal fillings of one FRAMED
boundary (= L2 decoration on S²) with different V_int. New
`scripts/graft_ball_search.py`: best-first refilling enumeration =
enumerate_fillings.py moves + the missing 1→4/4→1 channel (vertex-preserving
flips can NEVER change V_int) + exact framing bookkeeping (the old filler's
2→3 can silently change boundary-edge framing — my acceptance requires
frame_dev==0; NOTE this means enumerate_fillings.py itself may collect
framing-broken "fillings"). Lazy priorities: children scored by exact delta
bookkeeping, analyze only on pop. Acceptance = L2-exact; full L3/Z check by
validate_facets on the performed graft.
**c15 m3 results (max_bad 12, grow 12):** first flip out of a pristine
crystal filling costs ~8 (block-scale echo of "no single Pachner move
preserves FK"). Vertex stars (Z12, Z16) and edge cavities (edge1/2):
search EXHAUSTED → provably UNIQUE FK filling in radius. Face cavities
(V_int=3, 47-62 tets): cap-limited (100k pops); no dV hit yet, but best
V_int=2 filling reached total violation 3 (face2). Deep runs (400k
bad-first; 200k vfirst) launched 2026-08-06 evening.

**EDGE-CONTRACTION MOVE (Aaron's donor-free scheme, 2026-08-06 evening).**
"Delete two adjacent stars, cone the boundary from one new vertex" = EDGE
CONTRACTION (u,v)→w; validity = link condition = lump_boundary's pinch check.
`scripts/contract_relax.py`: exact local delta bookkeeping (verified vs full
recompute), apply, and A/B relaxation driver. **Measured on c15 m3:**
contracting a Z12-Z12 deg-5 edge costs remarkably little — CN(w)=17 with ALL
17 edges legal (m(w)=0, a pure-{5,6} "Z17", n6=5), just 3 ring edges go
deg-4; Δ(f0,f1,f3) = (−1,−6,−5), Σm² +10. Deg-6 ring edges drop to LEGAL 5.
**A/B verdicts:** (1) at c_imp≤0.5 the gas isn't f0-stuck (sampler sheds
inserted vertices in <100 sweeps; at 0.2 equilibrium f0 < crystal). (2) at
c_imp=3.0 f0 is FROZEN: prep-at-inflated-f3-target then quench → sampler
stuck at f0=650; +contraction channel reaches exactly 648 instantly (guard),
BUT its debris can't heal at 3.0 (S_B > S_A). (3) UNGUARDED contraction
ratchets BELOW crystal f0 (volume pin pays for it) and strands defect-rich —
the reverse (vertex-split) channel is needed for equilibrium; guard
(f0 > f0_ref) is the annealing stopgap. D-port design: split = choose vertex
+ splitting cycle in link; detailed balance needs per-vertex cycle counts
(bounded-length cycle enumeration in links, or restrict both channels to
deg(uv) ≤ L).

**Next steps:** cross-library graft compatibility matrix beyond the Laves
family (a15, sigma, z, mu, p, delta, r — expect the vertex-density/web
obstruction to appear as NON-EXISTENCE of shared L3 surfaces, degrading to
L2/L1 with quantifiable seam debt); thermal stability of a grafted seam under
the minimal 3-term action; mixed-sign pin-gap composite.
