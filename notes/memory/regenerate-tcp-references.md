---
name: regenerate-tcp-references
description: Reference crystals and seed families are gitignored; how to rebuild them locally
metadata: 
  node_type: memory
  type: project
  originSessionId: eaaf2503-7af1-42b2-a8c0-88f714d19b53
---

`data/` and `seeds*/` are gitignored, so a fresh clone / new machine has NO
reference crystals or seed families — only the scripts and 3 tracked seeds in
`standard_triangulations/` (600_cell, T3_A15_m2_N368, T3_A15_m3_N1242).

**Reference TCP crystals** (`data/tcp_reference/*.mfd`) are regenerated exactly
and cheaply from published Wyckoff positions:
`python scripts/tcp_reference.py <struct> -m <mult>` (or `all`). Deterministic,
validates on build (all edges {5,6}, Z-census, χ=0, orientable). Sizes used by
`tcp_melt.CRYSTALS`: r m3=N24462 / rbig m4=N57984; c15 m4=N8704 / c15big m9=N99144;
a15 m6 / a15big m13; sigma m4 / sigmabig m7; plus c14/c36/z/mu/p/delta at m3–5.

**Seed families** must be re-sampled (expensive MCMC) via `equilibrium_vdv.py`
`--produce` / `grow_seed.py` — not regenerable for free. See [[tcp-r-c15-defect-state]].
