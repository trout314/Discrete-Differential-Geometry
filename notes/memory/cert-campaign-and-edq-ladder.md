---
name: cert-campaign-and-edq-ladder
description: "Crash-resume checkpoint (2026-07-22 eve): stiff EDQ rungs 1.2/1.6/2.4 done (no mobile-legal window); the 8-chain lam0.4 two-sided cert campaign was killed mid-flight and relaunched from snapshots"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-23T01:35:47.879Z
---

Session 1ed7f5d0 crashed 2026-07-22 ~18:36 local, killing 8 running
`mobile_gas.py` chains. Everything lives in that session's scratchpad — reuse it
verbatim, do NOT move: `/private/tmp/claude-501/-Users-atrout-Desktop-Discrete-Differential-Geometry/1ed7f5d0-df5c-4552-b65e-5060180eafe2/scratchpad`
(`r_m4.mfd`, `mgas/`, `run5h/`). Work through commit 5961d0b is pushed; nothing
uncommitted was lost.

STIFF EDQ-ONLY RUNGS COMPLETE (n6 OFF, aesthetic-action test; logs
`mgas/edq{12,16,24}.log`). Completes the ladder in [[six-web-gauge]]:
- lam_EDQ 0.4 -> crumple-melt; 0.8 -> dense disclination fluid.
- 1.2 -> defect fluid that DENSIFIES (nill 40 -> 134 over 12k sw): not stationary.
- 1.6 and 2.4 -> nill = 0 and <edeg> = 5.104225352112676 pinned to 16 digits for
  essentially all 12k sweeps (one 7-illegal excursion at the last record):
  FROZEN CRYSTAL, not a legal liquid.
VERDICT: the EDQ-only ladder has NO mobile-legal window — too soft crumples, too
stiff freezes into the crystal, and the intermediate rung drifts. R^2 alone can
never hold a *fluid* on the constraint manifold, so the six-web line-closure
(decoded n6) term is necessary. This strengthens (does not contradict) the
two-term-action verdict in notes/two_term_action_and_emergent_gauge.tex.

CERT CAMPAIGN (lam=0.40 two-sided R-hat gate; the ensemble is still uncertified,
qR-hat 1.062 / ESS 83 one-sided). 5 below-side extends reached sweep 15000/20000,
3 above-side chains (run5h lam=1 inits, relaxing down) reached 14000/22000 when
killed. Relaunched 2026-07-22 19:35 from the last snapshots with fresh seeds
(701-705, 801-803), outputs `mgas/{lam40,l40s20[1-4]}x2.*` (+5000 sw) and
`mgas/ab{1,2,3}b.*` (+8000 sw); runner script in this session's scratchpad as
`resume_cert.sh`. All 8 chains sit at nill ~ 40-100 with no side-dependent drift
— above and below overlap already, which is what the gate needs.
NEXT after the gate passes: pool with `quantized_split_rhat` (never hand-roll
R-hat — see [[five-illegal-knot]]), then the budgeted-illegality sampler from
[[six-web-gauge]].

FINE EDQ SWEEP DONE (2026-07-22 late eve; outputs `mgas/edq{125,130,135,140}.*`,
20k sw each, n6 OFF): NO stationary rung — all drift UP with significance
(slope nill/1000sw: 2.12 / 0.88 / 0.76 / 0.42±0.12 (3.5σ) at 1.25/1.30/1.35/
1.40). Slope falls monotonically, extrapolates to 0 near lam_EDQ ~ 1.45-1.5
(1.6 = frozen). Host is CLEAN everywhere: single r-phase grain, 97.8-98.4%
interior-crystalline, legalvert 0.995-0.998. Defects MOBILE at all rungs
(Jaccard -> 0 by 1.5-3k sw at 1.25-1.35, FASTER than the lam0.4 full-action
liquid; edq140 has a persistent ~10% core, J-plateau 0.11). Open: stationary
dilute gas near 1.45 vs slow upward relaxation everywhere. Discriminator =
two-sided bracket at 1.40/1.45 (above-side from dense edq12 states).

GOTCHA (my slip): mobile_gas snapshots fire only when (done-BURN) % SNAP == 0
with done stepping by TS — SNAP must be a MULTIPLE of TS (150). SNAP=4000 gave
1 snapshot (sw14000); SNAP=2500 with 5k span gave none. Use 3000/4500/6000.

mobile_gas.py NOW EMITS `cents` (commit 9a6645d): per-complex tree-lift
centroids each ts frame (validated 12/12 exact-zero steps on unchanged
membership) — pass1_kinematics-style MSD now possible for any new run; the
old runs (incl. the fine EDQ rungs) lack it. Legacy cleanup commit 5b44f3c
archived decorr/knotgas_sk/validate*/reconcile/final_fig/local_lapse_proto/
migrate_* to legacy/.

CERT GATE VERDICT (lam=0.40, 5 below + 3 above, 22k-sw common window):
qR-hat(nill) = 1.28 (1.51/1.36/1.28 as window grows — improving, FAR from
pass). Above chains disagree among themselves (1.47; means 51-89); ESS/chain
4.5-149. Needs SEVERAL-FOLD extension, not a top-up. Below chains now total
~40k sw each (orig 20k + x 15k + x2 5.1k); above ~22k (14k + b 8.1k). All
extension segments resume from `mgas/*x_snap15000` / `ab?_snap14000` lineage;
new segments will carry cents/MSD data.
