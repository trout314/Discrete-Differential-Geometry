---
name: pl-geodesic-permeability
description: Only degree>=6 edges can carry PL minimizers (the FK six-edge network); exactness limits and the certificate that closes them
metadata: 
  node_type: memory
  type: reference
  originSessionId: 93ac5363-1537-4e86-9a6a-0dff5e3b1ea4
  modified: 2026-08-02T14:39:11.938Z
---

Established 2026-08-02, in response to the user's (correct) objection that
minimizers between generic point pairs are NOT generic geodesics and can be
pinned to the 1-skeleton.

**Permeability threshold.** Transverse to a degree-k edge the geometry is a 2D
cone of angle Theta_k = k*theta, theta = arccos(1/3) = 70.5288 deg. A minimizer
can pass through the edge iff Theta_k > 2pi. Since 5*theta = 352.64 deg < 360 <
6*theta = 423.17 deg:

- degree 3, 4, 5 edges: geodesics AVOID them.
- degree 6: shadow wedge Theta-2pi = 63.17 deg = 14.93% of transverse directions.
- degree 7: 133.70 deg = 27.08%.

So the geodesic-capturing set is EXACTLY the FK six-edge network (10.0% of edges
in C15 m3, 432/4320) — the same object as the disclination skeleton. A 2->3
promotes 4 spokes 5->6 (new capture lines) and 2 spokes 6->7 (widening 63->134
deg). The 600-cell is entirely permeability-free (all edges degree 5, all links
positively curved), so it is the clean validation bed for pure strip unfolding.

**Vertices too.** Minimizer through a vertex iff diam(link) >= pi. Links are
spherical triangulations with all edges 60 deg, area (2k-4)(3*theta-pi): Z12
11.03 < 4pi, but Z14/Z15/Z16 = 13.23/14.33/15.44 > 4pi. Threshold k* = 13.3973,
which the TCP mean coordinations (C15 13.333, sigma 13.467, A15 13.500) straddle.
Z12 is not automatically safe either — its link carries its own negative cone
points wherever a degree-6 edge passes through. diam(link) per Z-class is still
UNCOMPUTED.

Identity: sum_v (4pi - area(link_v)) = 2 * S_Regge = 2(2pi f1 - 6 theta f3);
verified 45.7253 both sides on C15 m3.

**Exactness limits.** Steiner is UNHARMED by any of this (edge Steiner points are
shared across the whole tet fan, vertices are graph nodes), so bound-based work
stands. What breaks is "straight chord in an unfolded strip": valid only for
contact-free minimizers. Provably exact vertex-to-vertex distances are NOT
achievable in general — the analogous 3D problem (Euclidean shortest path among
polyhedral obstacles) is NP-hard (Canny-Reif 1987), and comparing sums of square
roots has no known polynomial algorithm. There is no 3D analogue of MMP /
Chen-Han.

**The certificate that works per-pair.** Upper bound = any realized path (exact
algebraic). Lower bound by Kantorovich duality: for PL u with ||grad u||^2 <= 1
on every tet (a RATIONAL, exactly checkable inequality), d(x,y) >= u(y)-u(x). If
U <= L the candidate is globally optimal and d = U exactly. Do NOT use a
cell-based lower-bound graph: cells on two faces of a tet touch along their
common edge, giving zero-weight edges, and since the 1-skeleton is connected the
bound collapses.

Contact-free lengths are single sqrt(rational) (comparisons trivial); one edge
contact gives ell = sqrt(dz^2 + (rho_A+rho_B)^2), degree <= 4. Minimizers touch
edges at ISOLATED points, never run along them (d/dt_in of the along-edge cost is
strictly negative). See [[defect-boundary-map]].
