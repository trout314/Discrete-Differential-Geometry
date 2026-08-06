---
name: a15-positive-gas
description: "A15 dilute net-positive-curvature defect gas: the window is c_imp 0.35-0.50 under the minimal flat-pin action; elementary species = deg-4 dimers (NOT the c15 (3,4,4) flicker); sign forced by e_nat > e*"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1fc2ced4-2271-4a60-ba17-bc68ad6d30ad
  modified: 2026-08-06T18:46:50.211Z
---

**Question (Aaron, 2026-08-06 late):** find a regime where A15 supports a
DILUTE population of NET-POSITIVE-curvature defects.

**Why A15 is the right host:** e_native = 5.1111 > e* = 5.1043, the largest
debt in the library (51.3 forced moves/1000 verts, [[crystal-library-gas-
campaign]]), with gap = f1 − 6f3/e* NEGATIVE — the flat pin pulls the mean
edge degree DOWN, so the debt currency is LOW-degree (positive-curvature)
structure. (Via f1 = f0 + f3 this equivalently means the pin wants MORE
vertices: f0* > f0_crystal — mirror of c15; the SPLIT side of the
contract/split channel would pay a15's debt.)

**Charge convention used:** per defect complex (connected illegal edges),
net_q = Σ_{illegal e}(5 − deg): deg-4 → +1, deg-3 → +2, deg-7 → −2.

**Census of the EXISTING campaign snapshots (a15 m5, f0=1000, final .mfd):**
- c_imp 0.00: melt, one percolated 5136-edge complex, net −33.
- 0.10–0.20: percolated but strongly net-positive (+1236, +1314).
- **0.35: 126 complexes, median size 2, max 66, 124 positive / 0 negative,
  total +643; degs {3:112, 4:465, 7:23}. Time series NOT stationary: n_ill
  climbed 247→723 over the span, just flattening at end; top1 40–80
  (6–12% of mass) — flirting with aggregation.**
- **0.50: STATIONARY dilute all-positive gas, n_ill 10–36 fluctuating
  (alive), sizes 2–5, only deg-3/4 — but sparse: debt mostly UNPAID (m²
  beats the pin).**
- 0.70/1.00: pristine.

**Elementary species: deg-4 DIMERS dominate (median complex size 2,
deg-4:deg-3 ≈ 4:1)** — a15's positive quantum differs from c15's (3,4,4)
flicker knot. Species census (valence/disclination instruments) not yet
run.

**Window: c_imp ∈ [0.35, 0.50] — INTERPOLATED (runs a15_c{0.35,0.40,0.45}L,
spans 6-9k sweeps; NOTE crystal_gas_scan has no --mcell arg, LIBRARY fixes
the cell):**
- 0.35 long (seed 43): AGGREGATES — n_ill ~800 but top1 median 177/max 385
  (≈half the mass in one cluster). Not dilute at equilibrium; the original
  campaign snapshot just hadn't coarsened yet.
- 0.40 (seed 41): stationary n_ill = 661±18, 135 complexes, 134 positive/1
  negative, totQ +625, size med/max 2/39 (top1 excursions to 97, ~5% mass).
  Boundary case: dispersed with mid-size clusters.
- **0.45 (seed 42): THE REGIME — 76 complexes, 75 positive/1 negative,
  totQ +205, size med/max 2/10, no aggregation; n_ill ~170-250 and still
  filling at span 6000 (rising from pristine).** 16k-sweep certification
  run launched (seed 45, a15_c0.45XL).
- 0.50: stationary but sparse (debt largely unpaid, n_ill ~20).
Trade-off axis: paid-but-coarsening (≤0.40) vs dilute-but-underpaid
(≥0.50); 0.45 pays ~a third of the quota as a dispersed gas.

**SPECIES LADDER (HTML catalogs, defect_catalog --narrow — REQUIRED on a15:
vertex-adjacency complex definitions percolate at 46% defect-vertex
density):** the elementary quantum is a SINGLE deg-4 edge (Q=+1); the
population is a deg-4 multimer ladder (0.45: 54/23/9/6 monomer/dimer/
trimer/tetramer of 126; 0.50-dilute: 21/12/5/2 of 43). Deg-3 species rare;
the lone deg-7 in the 0.50 box is CAGED inside a 3⁸4·7 cluster (screened
negative grain). Catalogs: data/viz/a15_c0.45_narrow/index.html and
data/viz/a15_c0.50_dilute/index.html (+ standalone_cx*.html self-contained
scenes). 0.45XL certification: late slope +0.21±0.10/obs — NOT yet
stationary at 16k sweeps (n_ill ~560 and creeping); 0.45 equilibrium
resembles 0.40's (mid-clusters); strictly-dilute-stationary is ~0.50.

Related: [[crystal-grafting-program]] (contract/split channel: the split
direction is a15's debt-payer — an f0-unpinned a15 gas run is the natural
follow-up), [[m2-only-gas]], [[edq-only-melting]].
