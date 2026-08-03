---
name: defect-boundary-map
description: "2->3 defect is a LENS not a shortcut: exact boundary distance map of the ball, lambda=2*sqrt(6)/3 anisotropic scaling + double-image caustic"
metadata: 
  node_type: memory
  type: project
  originSessionId: 93ac5363-1537-4e86-9a6a-0dff5e3b1ea4
  modified: 2026-08-02T14:38:55.089Z
---

Built 2026-08-02. `python/discrete_differential_geometry/ball_boundary.py`
(library) + `scripts/defect_boundary_map.py` (CLI, figure, exact JSON table in
`data/figs/defect_boundary_map_n{6,8}.{png,json}`).

**The object.** A 2->3 replaces 2 tets by 3 on the SAME 5 vertices, so the two
complexes are identical outside a ball B and share dB = six unit equilateral
triangles. All metric consequences live in the pair d_B^R, d_B^D on dB x dB;
d_out (complement distance) is common mode and cancels in d_D - d_R. The map is
UNIVERSAL — depends only on the move's combinatorics, not the host — so compute
it once. Labels: A,B,C = shared face, P,Q = apexes. R: ABCP, ABCQ. D: ABPQ,
ACPQ, BCPQ.

**Headline: it is a lens, not a shortcut.** At order 8, 29% of boundary pairs
shorten, 40% LENGTHEN, 31% unchanged. My earlier "defect only shortens" was
about GLOBAL distances (which take the min with routes around B) — inside B
tangential paths get longer.

**Exact anisotropic scaling, lambda = 2*sqrt(6)/3 = sqrt(8/3) = 1.632993:**
- axial family (straddling the new edge PQ): d_D/d_R = 1/lambda = sqrt(6)/4
  EXACTLY, for the whole family (verified k=0..11 at order 12).
- equatorial family (straddling the destroyed face ABC): d_D/d_R = lambda
  EXACTLY, but only in the near-axis regime; beyond it a SECOND geodesic branch
  (around the back of the strut, past B and C) wins with rational lengths
  19/18, 10/9, 7/6 — an exact **double-image caustic**, crossover in
  t in (2/3, 3/4). This is the cosmic-string double image, exactly located.
- extremes: d(P,Q) 2sqrt(6)/3 -> 1; (A/3+2B/3, A/3+2C/3) 2/3 -> 4sqrt(6)/9.

**Mechanism.** 2->3 trades the SHEET ABC for the STRUT PQ. Through-paths
contract; anything that used to cut straight across ABC must now bend around PQ.

**Validation (all pass, orders 4/6/8).** 9 boundary edges = 1 exactly in both;
apex chord exactly sqrt(8/3) vs 1; symmetry; triangle inequality; d_B <= its own
boundary-surface distance; and an INDEPENDENT interior Steiner graph (which CAN
represent PQ-touching and dB-hugging paths that the strip enumeration excludes)
never beats the strip minimum, converging onto it O(h^2) from above.

See [[pl-geodesic-permeability]] for why the strip enumeration is complete here,
and [[no-halo-verdict]] — the exterior being literally identical is the PL
reason the static two-defect potential is exactly zero beyond contact.
