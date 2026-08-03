---
name: reaction-census-campaign
description: "8-hour campaign launched 2026-07-27 ~00:30: thermal merger/split reaction chemistry of defect complexes under FULL Pachner dynamics (new instrument reaction_census.py), plus high-statistics FP recombination-by-intrinsic-class; how to analyse when it lands"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-27T04:38:36.375Z
---

2026-07-27 (~00:30). Two campaigns launched, both deadline-driven, both
land ~07:00 the same morning. Machine: 8 cores / 16 GB, 7 workers, total
RSS ~2.5 GB (comfortable).

**THE GAP THIS CLOSES.** Every reaction measurement in the project so
far was in an artificial setting: the STATIC colliders (V(s) = 0 exactly
beyond contact) or the SLIDE-ONLY FP driver (P(rec|freed) = 0.70). Under
the full thermal dynamics the actual reactions are MERGES and SPLITS of
worldlines — which `defect_state.Worldlines` has always counted as two
integers and never characterised. Now it emits events.

**PHASE A — thermal reaction census** (5 workers, `data/reactions/`).
New instrument `scripts/defect_dynamics/reaction_census.py`; needed
`Worldlines` to gain `self.events` + `drain_events()` (birth / death /
merge / split, backward compatible — the 5 existing consumers are
untouched). Streams `<out>.events.jsonl`, `<out>.tracks.jsonl`,
`<out>.json`. Host R m4 N57984, **etarget 5.105025** (see gotcha), no
positions/cocycle at all so the run is frame-free.
Grid (name, lam, slide_prob, start snapshot):
  l40a 0.40 0.0  ab1_snap14000 | l40b 0.40 0.0  ab2_snap14000
  l40s 0.40 0.10 ab3_snap14000 | l35a 0.35 0.0  lam35_snap20000
  l35s 0.35 0.10 lam35_snap17000
Two crossed levers: DENSITY (lam 0.40 ~8-13 complexes vs 0.35 ~23) and
TRANSPORT (slide off/on — slide is DB-certified by V3b, so it changes
kinetics and NOT the stationary ensemble; a rate that moves means the
chemistry is diffusion-limited).

**EARLY NUMBERS at sweep 1000** (already informative): merge count ~=
split count in every run (steady state). l40a 26/23, l40b 58/65,
l40s(slide) 86/85 at ncomp 9, l35a 328/326 at ncomp 23. So slide roughly
TRIPLES the merge rate at LOWER density => strongly suggests
transport-limited. Throughput ~3 sweeps/s/worker under full load =>
expect ~70k sweeps and O(2000-20000) merge events per worker.

**PHASE B — FP recombination by exact intrinsic class** (2 workers,
`data/fpkmc/prodD/`, `run_recomb_campaign.sh`). Purpose: the
channel-dependent recombination (class 1/25 P(rec)=9/10 vs 64/225 1/5,
Fisher ~0.005) is underpowered. ~1 usable freed dock per 67 s per worker
=> expect ~700-750 vs the current 87. Host is **m2 deliberately**:
intrinsic classes are host-independent (measured) and m2 runs ~8x
faster. ONE EPISODE PER PROCESS — the D slide-graph scan leaks ~12
MB/scan (R6), so a long-lived process walks into GB RSS (observed 1.16
GB mid-episode).

**ANALYSIS** — `scripts/defect_dynamics/reaction_report.py` is written
and TESTED end-to-end on smoke data. Run it with
`--dir data/reactions --skip-sweeps 500`. It does: reduced rates
k2 = k_merge/<pairs> and k1 = k_split/<n> with BLOCK BOOTSTRAP errors
(intensive across the density ladder <=> the reaction picture is the
right language), the slide on/off transport verdict, K = k_merge/k_split
(an association free energy — the static potential is exactly zero
beyond contact, so K != ideal is entropic/kinetic in origin), compound
survival curves, and the species table (which induced shapes merge into
which). Smoke data already shows compounds (absorbed >=1 partner) live
median 35-100 sweeps vs solo 3-4 — the quantitative form of run5h's
"fused multimers are immortal".

**GOTCHA — e\* barely matters, and here is why.** I initially ran with
etarget = the crystal's NATIVE mean edge degree 5.104225 rather than the
convention e* = 5.105025 that every mgas ensemble used, and restarted.
Then measured: with the same seed the two targets give BYTE-IDENTICAL
trajectories over 100 sweeps. Reason: Sum_edges deg(e) is topologically
fixed at fixed N3, so Sum(deg - e*)^2 = Sum deg^2 + const(N3) — e*
enters ONLY as an additive constant plus the coefficient lam*e*/6, i.e.
a **0.016% change in stiffness**. Do not panic about small e* mismatches;
do keep 5.105025 for comparability.

**STILL UNCOMMITTED** at launch: reaction_census.py, reaction_report.py,
run_reaction_campaign.sh, run_recomb_campaign.sh, the Worldlines event
patch in defect_state.py. Offer to commit.
