---
name: defect-kinetics-run5h
description: "run5h defect structure+dynamics — (3,4,4)-knot monomers blink, fused multimers are immortal; Q≈−5/knot; Rg~N^0.49"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-22T15:40:54.718Z
---

Three-pass analysis (2026-07-21) of the run5h R-m4 ensemble (8 chains, ~17.7k
sampled sweeps, frames every 150 sw in *.ts.jsonl, 4 snapshots+cocycles/chain,
.events.bin move stream). Scripts: scratchpad/pass{1,2,3}_*.py; figures
pass1_kinematics.png, pass2_structure.png, pass3_microdynamics.png.

PASS 1 (kinematics, 1997 tracks): two populations split sharply at size ~10.
Size 5-9: ephemeral blinkers — median life ≤150 sw (one frame), zero excursion,
~13 births/1000sw/chain, none survive. Size ≥10: immortal frozen cores — ≥15s
live the whole run (86% of founding survive), excursion ~0.3 cells/17.7k sw.
NO diffusing population at any size; MSD flat. Explains split R-hat/Jaccard
glassiness: frozen cores + fast blinker gas (n_illegal IAT ~257 sw).

PASS 2 (structure, 198 complexes/32 snapshots, vs R reference):
- State is ONE grain (98% crystalline) — cores are NOT grain boundaries.
- (3,4,4) knot = universal monomer: every size-5 has exactly 1 deg-3 + 2 deg-4
  edges; NO complex of size 1-4 ever exists. n3 counts knots in multimers
  (10-14: n3~2.2; ≥15: n3~3.7) with monomer n3:n4=1:2 ratio preserved →
  big cores are FUSED KNOT MULTIMERS (sizes ~11=dimer, ~21=tetramer).
- Universal charge quantum Q ≈ −1.0 rad per illegal vertex (Q≈−5/knot).
  NOTE (corrected): "above/below" chains both ran the SAME +8e-4 bump — they
  differ only in INITIAL state (two-sided certification: perfect vs
  over-coordinated start). So Q-consistency across them shows convergence of
  the two-sided design, NOT drive-independence.
- Rg ~ N^0.49: floppy string/random-walk geometry (knotted disclination line).
- Halo: monomers registry-bare (halo/N~1.6); multimers drag 3-5x off-registry.
- Blinker sites are spatially UNIFORM (KS p=0.52 vs null) — no templating by
  cores; earlier 82% recurrence stat was chance-level.

WHY IMMORTAL (hypothesis, confirmed shape via Pass 3): lone 5-knot has a local
create/annihilate channel (blinks). Fused multimers have none — partial
annihilation would need sub-5 illegal intermediates (never observed as states),
and monomers can't diffuse apart (no transport channel in the Pachner move set).
Kinetically immortal; count/arrangement frozen at quench value → each chain
its own basin.

PASS 3 (events.bin, 371k accepted moves): multimer cores are NOT dynamically
dead — they are the HOTTEST objects in the system: 17.8 ev/vtx/1000sw = 133x
bulk; 40% of ALL accepted moves happen at the ~0.7% core vertices (blinker
sites 76x; halo 3x; bulk vertex moves once per ~7700 sw). But the activity is
pure confined rattle: volume moves (1->4/4->1) are ~0% EVERYWHERE (volume pin
shuts the channel; facet count wobbles only +-4), so cores flicker exclusively
via balanced 2->3/3->2 (37%/37%) and 4-4 hinge moves — 4-4 is 26% at cores vs
10% bulk (2.6x enriched; the degree-redistributing move is the core's rattle
mode). Steady in time, no bursts. Dynamic heterogeneity: top 1% of vertices
carry 45% of all activity. Picture: dynamically trapped two-level rattling on
a closed micro-configuration set — glassy trap, not frozen stillness.

PIN-RELEASE TEST (bump 8e-4 -> 0, 4 chains x 5.1k sw from snap20050; scripts
scratchpad/release_run.py, data scratchpad/release/):
- MAJORITY THERMODYNAMIC: n_illegal ~70 -> ~21 within 600-1800 sw (blinker gas
  + most multimer mass dissolves when the drive is removed). below0 annealed
  COMPLETELY to the perfect crystal (n_illegal=0 by sweep 1200) — the stealthy
  vacuum is dynamically reachable.
- MINORITY METASTABLE: 0-2 survivors/box (sizes 9-13), stable >=4k sweeps
  against a strictly downhill drive (cimp=1.0/vertex). Barriers heterogeneous.
- SURVIVORS RESTRUCTURE: all have n3=0 (deg-3 edges healed) — pure OPEN CHAINS
  of deg-4 (+2 disclination) edges (checked: not closed loops; incidence mostly
  1), ~0.4 member overlap with ancestors, same location, Q ~ -7..-13. The
  pinned-state knot motif (with deg-3) is NOT what survives; the metastable
  object is a capped deg-4 segment. Protection mechanism ≠ loop topology;
  presumably a shrink barrier through deg-3/sub-5-illegal intermediates.
- KNOT-GAS S(k), CORRECTED 2026-07-22 (user caught the reference-convention
  error): complex charge MUST be referenced to the crystal CELL-MEAN q_R
  (~-0.002), NOT the median (= Z12 value +0.77, an icosahedral-vacuum ref).
  With the correct reference: knot monopole ~ -1.05 rad/knot (~-0.21/illegal
  vertex), bare Poisson floor = 0.0151, S_knot = 0.982, measured plateau
  0.0167 = 1.12x bare floor. So the frozen gas radiates like BARE monopoles:
  the earlier "93% halo screening" was an ARTIFACT of the Z12 reference —
  RETRACTED. All screening must come from the collective arrangement
  (mobility); none is local. (Q~-5/knot in pass2 tables = Z12-referenced;
  far-field monopole is ~-1.05.) ALSO: flat degree = 5.104299 (2pi/arccos(1/3));
  5.0043 was the POSITIVE-curvature pin of the old pos/flat/neg campaigns, NOT
  flat. R native 5.104225 is slightly below ideal flat; run5h pin native+8e-4
  is slightly ABOVE flat -> demands net negative curvature -> the knots supply
  it (dissolve when pin released). Neutrality sum-rule audit (does carrier
  monopole total match pin demand?) not yet done.

Centroid protocol (validated): ts.jsonl centroids are raw cocycle TREE-LIFT
coords — internally consistent in time (99.4% of unchanged-membership steps
exactly 0; 0.09% >1-cell spanning-tree gauge jumps — drop those), but ~0.8-cell
static offset vs harmonic coords. Use for RELATIVE motion only; recompute
harmonic (cocycle.torus_positions) centroids from snapshots for absolute
positions. Clusters are compact in min-image (no lift-splitting).

See [[five-illegal-knot]], [[curvature-length-scale]], [[crystal-grains-tool]].
