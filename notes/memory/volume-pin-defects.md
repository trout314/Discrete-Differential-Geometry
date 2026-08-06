---
name: volume-pin-defects
description: "Volume pin as a defect dial (channel on): +-5% strain on c15/a15 forces a dispersed ~500-650-edge gas, POSITIVE currency regardless of strain sign or host; f0 tracks the flat line f0*(f3); native-volume a15 does NOT refinance debt into vertices (interstitials too expensive)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1fc2ced4-2271-4a60-ba17-bc68ad6d30ad
  modified: 2026-08-06T21:18:38.028Z
---

**Experiment (Aaron, 2026-08-06 night):** volume pins ≠ native volume with
the contract/split channel on ([[crystal-grafting-program]]); c15 (e_nat <
e*) and a15 (e_nat > e*), ±5%. Instrument: crystal_gas_scan --vol-scale
--contract-split (committed a4eec1a+). All at c_imp 0.5, prob 0.02, 12k
sweeps, seeds 95-98; outputs data/crystal_gas/{c15,a15}_v{0.95,1.05}_cs.

**RESULTS (a15 pair final, c15 pair near-final):**
- Volume pin satisfied EXACTLY in all four; **f0 slides along the flat line
  to f0*(f3) − ~5** (c15 −5%: deleted 90 vertices, 1536→1446; +5%: created
  63; a15: 1000→953 / 1054). The channel + 1↔4 jointly do the work fast.
- **Strain forces a dispersed defect gas ~500-650 illegal edges** in
  140-190 complexes, top1 ≤ 19 (NO condensation), stationary (a15 −5%:
  506±12 no trend).
- **Currency is POSITIVE regardless of strain sign and host** (deg-4
  dominant ~320-370, deg-3 second, deg-7 only 3-12) — stretch does NOT
  summon negative structure (prediction falsified). m² prices deg-4 and
  deg-7 equally, so the asymmetry is structural/entropic — open question.
- **Host-independence**: c15 and a15 respond nearly identically once f0 is
  free; the native-curvature distinction washes out under strain.
- e_mean lands ~0.004 above e* (residual = the ~5-vertex f0 offset below
  the line — plausibly partly the bistellar 1↔4 uncorrected-asymmetry
  chemical potential ~ln2/vertex; see below).

**NATIVE-volume a15 channel runs (c_imp 0.45/0.60): NO refinancing** — f0
creeps only +2..+4 (f0* = 1009), defect gas stays; an a15 interstitial is
a fully-dressed defect cloud whose m² cost ≈ the gap² savings, unlike
c15's cheap legal-edged CN-17 contraction scar. The two crystals choose
OPPOSITE currencies for the same pin pressure — set by local geometry
cost, not debt sign. High two-way channel traffic with net absorbed by
bistellar 4→1s (loop current: the 1↔4 pair's known uncorrected proposal
asymmetry ≈ anti-vertex potential ln2 — the exact channel fights it;
measured f0 equilibria are the MIXED chain's).

**CONTROLS (channel OFF, both crystals −5%): ABSOLUTE RIGIDITY.** f3 frozen
at native (a15: 288, c15: 435 above target), f0 exact, ZERO illegal edges
all 12k sweeps, full volume penalty eaten (obj 8376 / 18997 vs 666 / 733
with channel). Mechanism: a pristine over-volume crystal has NO legal
first move downhill — 3→2 needs a deg-3 edge, 4→1 a CN-4 vertex, and any
path starts with a 2→3 costing +0.2·Δf3 in the volume term (+58 / +87),
so not even thermal flickers nucleate (they exist only at native volume);
the barrier GROWS with the strain it must relieve. Contraction tunnels it
in one composite step (Δf3 ≈ −5). ["c15 leaks via its thermal flickers"
hypothesis falsified.]

**FINAL 4-cell table (channel on):** c15 −5%: n_ill 609±19 top1 12; c15
+5%: 671±22, 15; a15 −5%: 506±12, 15; a15 +5%: 519±22, 24. Tension ~10%
denser than compression in both hosts (more deg-3); per-tet the a15 gas
is DENSER (0.09 vs 0.073/tet) despite the smaller absolute counts.

**UNIVERSAL QUANTUM (strained c15 catalog, data/viz/c15_v1.05_strained,
defect_catalog --layout mds — combinatorial layout since the cocycle lift
can't follow the channel yet):** the strain-born gas shows the SAME deg-4
multimer ladder as a15's native-debt gas (73/46/23/7 of 199) — the single
deg-4 edge is the universal positive quantum of the flat-pinned FK
ensemble, independent of host and driving force. All six deg-7-bearing
complexes show the CAGED-MONOPOLE motif (deg-7 wrapped in a deg-3/4
shell), same as a15's lone deg-7 — free negative curvature does not occur.

All single chains, uncertified. Related: [[a15-positive-gas]],
[[crystal-library-gas-campaign]].
