---
name: knot-gas-campaign-lam035-m4
description: "7h defect-complex data campaign (2026-07-28) on the reproduced dilute gas -- m4 R, lam=0.35, FULL FK action zleg=0.6 cimp=1.0; reaction_census.py 8 chains; config recovery + how to analyse"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-28T23:02:59.195Z
---

**ABORTED then VINDICATED 2026-07-28.** I aborted the 8-chain campaign ~15 min in
when 2/8 chains "failed" reaction_census's self-audit (edge-degrees off by ~49.5k).
That was a FALSE ALARM: after a long bisection the real bug is in the AUDIT, not
the tracker. DefectState is CORRECT the whole time (probe: st.edeg ==
recompute(st.tets) == manifold at every step). ROOT CAUSE: audit built ref_edeg
from fk_skeleton.edges_from_facets(), which RELABELS vertices to dense 0..V-1, and
compared those dense-index edge keys against self.edeg (keyed by REAL vertex
labels) WITHOUT remapping. While labels stay contiguous they coincide (audit
passes); the instant vertex-label recycling (4->1) opens a GAP in the numbering
(~sw300 on the dense m4 start), np.unique compresses and dense!=label -> WHOLESALE
false diff (~49572, ~half the edges) that never recovers. The 2026-07-27
"duplicate-event" hypothesis in diag_state_divergence.py is DISPROVEN (apply-only
replay: 0 phantoms, clean 1500 sw).
FIX (committed pending): defect_state.py audit() lines ~609-618 -- remap
edges_from_facets' dense indices back through lab=np.unique(fac):
`ref_edeg = {tuple(sorted((int(lab[a]),int(lab[b]))): int(d) for (a,b),d in
zip(eu,edg)}`. Verified: real reaction_census past sw300 now clean.
IMPLICATION: the aborted campaign was almost certainly producing CORRECT data;
only the self-check cried wolf. Relaunch the 8-chain campaign -- instrument is
sound. (My earlier .get() edits to DefectState read methods were a WRONG
hypothesis, reverted.) State point + action below are VALID and reproduced.

(original launch notes below --- config still correct)

**7h campaign to collect defect-complex shape/lifetime/reaction data on the
reproduced dilute knot gas.**

**State point (recovered from reaction_census.py defaults + data/mgas snapshots):**
- Cell: `data/tcp_reference/T3_R_m4_N57984.mfd` (m4 R, big box -> dilute).
- lam=0.35, etarget=5.105025, `hinge_coef=lam*et/6`, num_facets_coef=0.1.
- `set_n6_potential(0.6*lam, 1.0*lam)` i.e. zleg=0.6, cimp=1.0 (STRONG legality;
  the 0.3/0.3 I first guessed was 3x too weak and melted). cimp=1.0 (strong m^2)
  gives diluteness; zleg=0.6 (n6) makes the dilute gas STABLE.
- Gas census: n_ill~195-244, ~27 components, sizes (illegal-edges/component) from
  2 up to ~29-35 -- broad variety. Dilute (~0.3%). NB two DIFFERENT census units
  in play: my scratch census counts VERTICES/component (max~9); ds.census in
  reaction_census counts ILLEGAL-EDGES/component (max~35). Same gas.

**System is FIRST-ORDER/BISTABLE (robustly established this session across m^2,
n6, full-action, m3+m4, all lam):** no two-sided single-phase dilute gas exists.
At lam=0.35 m4: pristine metastable branch ~48 (never nucleates) vs the GAS
branch ~195 (over-def from 198 AND ref-snap from 127 both converge to ~195, tight
+ stationary). So certify the gas WITHIN-BASIN (all chains start in the gas),
report as a stable stationary gas, NOT a two-sided cert. This is what the
historical lam=0.40/0.35 m4 campaigns implicitly did.

**Campaign = reaction_census.py x8 chains** (the validated instrument: streams
merge/split EVENTS, worldline TRACKS with lifetime+max_size+induced_shape(species
face-vector)+coordination, life_hist, census ts; tracking via DefectState/
Worldlines + sampler event log; FRAME none -- graph identity, no positions):
- Over-dispersed starts: 4 from lam35r_snap15000 (127, below) + 4 from
  scratch _m4_over.mfd (198, above); both flow to ~195.
- Slide control: 4 chains slide-prob=0.10 + 4 slide=0.0 (same stationary
  ensemble by DB; on/off separates transport- vs reaction-limited kinetics).
- seeds 101-108, --max-seconds 21600 (6h), out=`data/rxn_lam035_m4/c0..c7`.
- Launcher: scratchpad/launch_campaign.sh. Throughput ~3 sweeps/s/chain on m4
  (tracking overhead) -> ~65k sweeps/chain.

**ANALYSIS (after):** existing chain -- mgas_analyze.py, species_report.py,
reaction_report.py, sl_verdict.py, carrier_gr.py consume these outputs.
- Shapes: tracks.jsonl `shape` (induced_shape face-vector) + `coord`; census
  `sizes` for size distribution. (Geometry/anisotropy needs positions = offline
  cocycle pass; not captured here.)
- Lifetimes: `life_hist` + tracks `life_sw`,`max_size` -> survival curves;
  expect flicker peak (most births die in ~1 sweep) + long-lived tail.
- Reactions: events.jsonl merges/splits with participant sizes/ages/shapes;
  k_merge/k_split, K=k_merge/k_split (assoc. free energy), n^2 vs n scaling.
- Cert: within-basin qRhat over the 8 census ts (post-burn ~2000 sweeps) on
  n_components/n_illegal/max_size (use convergence.quantized_split_rhat).

GOTCHA: burn-in -- chains from 127 rise to ~195 over ~1500-2000 sweeps; discard.
Related: [[reaction-census-campaign]] [[flicker-background]] [[mobile-gas-liquid]]
[[edq-only-melting]] [[lifetime-vs-charge]] [[no-halo-verdict]].
