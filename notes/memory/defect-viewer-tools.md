---
name: defect-viewer-tools
description: interactive 3D defect viewer + sortable per-snapshot defect catalog (plotly HTML) — the visualization tools for defect complexes
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-29T23:06:41.713Z
---

`scripts/defect_dynamics/defect_viewer.py` — interactive 3D HTML scene per
defect complex (rotate/zoom/pan/hover/legend-toggle). Decoration per user
spec: illegal edges thick (orange deg-4, red deg-3, maroon deg-7+); legal
induced edges gray; anomalous-n6 nodes purple annotated n6−4; context =
closed star of the defect's vertices as translucent shell. Positions =
harmonic torus coords from `<snap>.cocycle.npz`, BFS-relative min-image
unwrap (single-ref min-image folds for big complexes). `--list/--complex
N/--all`; output `data/viz/`.

`scripts/defect_dynamics/defect_catalog.py` — sortable HTML table of ALL
complexes in a snapshot, rows linked to 3D scenes (shared plotly.min.js in
`data/viz/<base>/`). Columns incl. illegal-graph anatomy
(fragments/tips/branches/RINGS), Q rel, Rg (cells), decorated-species hash
(canonical_key; equal hash = same decorated species). `--note` for the
provenance stamp.

**Why:** built (2026-07-29, commit 68905e2) during the deg-4 worm program
([[deg4-worm-design]]) to inspect the long-lived bundle structures.

**How to apply:** needs `pip install plotly` + a cocycle-bearing snapshot.
lam35r_snap15000 catalog example: 24 complexes; big bundle = 14 fragments /
29 tips / 1 branch / 0 rings. Open `data/viz/<base>/index.html` LOCALLY
(links are relative). Harmonic chart smooths geometry near defects.
