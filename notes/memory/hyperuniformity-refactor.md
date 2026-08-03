---
name: hyperuniformity-refactor
description: Structure-factor / hyperuniformity code unified into one package library; where the pieces live
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T19:15:58.137Z
---

Full refactor of the S(k)/hyperuniformity code (branch `refactor-hyperuniformity`,
2026-07-20/21). Goal: one shared S(k) implementation, one field library, no
script rolls its own. Done + validated (9 new tests + full suite 121 pass, all
7 scripts import-clean, bit-for-bit equivalence to the old inlined code).

Package modules (the single sources of truth), in `python/discrete_differential_geometry/`:
- **vertex_fields.py** — `FIELDS` registry + `edges_and_degrees(facets)->(eu,ecnt,deg,V)`
  in np.unique DENSE order. Fields: `n6`, `curvature_charge` (= ½Σδ_e, the Regge
  scalar-curvature charge = Hamiltonian-constraint quantity; θ=arccos(1/3)),
  `curvature_density` (charge/(deg/4)), `volume_charge`, `defect_indicator`.
  The naive (6−deg) "deficit" is DROPPED — it is identically 12 (link sum rule),
  a constant field with no S(k).
- **structure_factor.py** — `structure_factor(frac,basis,field,nmax)->(kmag,s_obs,s_null)`.
  The ONE real-k S(k); deterministic analytic permutation null. (was sk_torus steps 7-12)
- **graph_hyperuniformity.py** — graph-proxy estimators: `adjacency`, `bfs_order`,
  `window_ratio`, `spectral_ratio`, `lowpass_ratio` (matrix-free, robust on crystals).
  rng is now an explicit arg everywhere.
- **cocycle.torus_positions(facets,edges,omega)->(frac,basis)** — the coordinate
  builder (canonicalize→tree lift→lattice basis→harmonic gauge→frac).

Scripts are now thin drivers importing the package: sk_torus (`sk_one` adapter),
hyperuniformity (importable CLI, default field curvature_charge), adm_constraint
(fields_from_facets via FIELDS), curvature_hyperuniformity_g, w_liquid_ladder,
complex_analysis, reference_summary. Tests: `tests/test_hyperuniformity.py`.

Physics checks that must hold: perfect crystal → low-k S(k)≈0 (HU); melt →
clustered (>1); doped crystal → 0.01–0.06. See [[crystal-grains-tool]],
[[five-illegal-knot]].
