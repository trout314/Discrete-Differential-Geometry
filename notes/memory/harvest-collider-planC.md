---
name: harvest-collider-planc
description: "Plan C DONE — knot vs REAL thermal targets (9-mer, curve-string, thermal knot) in frozen lam40x_snap15000; V=0.00 EXACTLY at range in the thermal background; free association at contact (all V=0, all unfuse); pristine S-matrix TRANSFERS"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-26T01:29:19.811Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/harvest_collider.py` —
knot C slid (deterministic slides, background frozen — sampler never runs)
into harvested long-lived complexes in `lam40x_snap15000` (census: 9-mer
(9,24,25,9) sumZ 127, curve-string (7,7,1,0) sumZ 92, 6 knots; the immortal
decamer is NOT in this snapshot — its worldline was born ~1300 sw later).
JSONs: scratchpad harvest_t0/t1/t3.json.

**Method that made it exact in a thermal state:** the washboard reference
comes from the PRISTINE crystal via the position->site registry (100.00%
match): find approach chains by mirroring seed windows to pristine sites and
walking the CANONICAL chain there (walking the thermal manifold near the
target routes through deformed tets and its mirror is not a pristine BC
chain — first version's bug); mirror the canonical chain back into thermal
labels (inverse site map). V(j) = [S(bg+C) − S(bg)] − W_pristine(site).

**Results (8 approaches, 3 targets):**
- V = 0.00e+00 EXACTLY at every range step — additivity holds in the
  thermal background with 8 other complexes present. The no-halo S-matrix
  transfers from pristine to thermal states.
- knot + 9-mer -> 14-mers (14,38,38,13)/(14,37,36,12), sumZ 197-199, V=0.
- knot + curve-string (7,7,1,0) -> 12-mers (12,22,15,4)/(12,19,11,3)/
  (12,20,10,3), sumZ 162 in ALL THREE geometries (a conserved compound
  ladder rung?), V=0.
- knot + thermal knot -> (10,23,20,6)/(10,22,19,6) sumZ 140/144, V=0 —
  IDENTICAL products to the pristine crossing collider
  ([[crossing-collider-phase2]]): direct transferability confirmation.
- Every fusion unfuses by one slide (reversible).

**Conclusion for the program:** free association + zero-range interaction +
reversibility hold for nature's own long-lived species in their native
background. The contact algebra measured on pristine crystal is THE
interaction physics of the thermal gas. Kinetic theory of the defect gas
can take the pristine S-matrix as exact input.

Note: no repulsive (+4.4-type) geometry appeared in these 8 approaches —
the picker takes smallest dmin first, which favors vertex-disjoint touches;
the sweep ([[smatrix-sweep-running]]) will quantify the branching.
