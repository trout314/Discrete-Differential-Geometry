---
name: fp-kinetics-findings
description: "First FP-driver physics (2026-07-26, R m2, lam=0.40 frozen background): tau(s0), isotropic dock angles, lam-INDEPENDENT D_slide via dS=0 flat directions, BC-chain precession period 8, recombination 0.70 with aligned-dock excess"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-26T23:28:51.006Z
---

Physical findings from the first FP production runs (M3 core validated;
host R m2 N7248, lam=0.40, slide channel, FROZEN background, time in
attempted moves; data/fpkmc/prod*, figures fp_production_report.png +
fp_dock_angles.png). Infra/how-to details in [[fpkmc-m1-status]].

**Transport (the big one):** D_slide = 6.6e-7 cells^2/attempt
(= 4.8e-3 cells^2/sweep), MSD exponent 0.88, and it is
**lambda-INDEPENDENT** (lam 1.0 gives 6.3e-7, exp 0.90). Mechanism: the
knot transports along dS = 0 flat directions (degenerate translates),
which no lambda can suppress — the washboard is nearly irrelevant to
mobility. Consequence: the thermal caging seen in [[defect-travel]]
(exp 0.63, range <= 1.26 cells) is ENTIRELY flicker/dressing, not
landscape. The bare landscape is a free highway.

**BC chains are coiled, not straight:** walk-direction autocorrelation
along the orbit oscillates with period 8 windows (+0.90 revival at lag
8), decaying to 0.16 by lag 48 — precession (helix-on-helix) plus slow
wander, persistence ~tens of windows (few cells; m2 box may truncate).
"Same orbit" does NOT imply parallel. Local rod axis estimator:
chain_walk +-8 windows (2 precession periods) end-to-end.

**Encounters:** tau(s0) KM medians 7.3e5 (s0=3) -> ~5e6 attempts
(s0=6-8); s0=10 drops (finite-box wrap, m2 = 2 cells — real tau(r)
needs m4). Docks ~83% OFF-CHAIN; dock-angle spectrum ~ isotropic
sin(theta) (mild 45-60 excess): a diffusing knot forgets direction
before docking. 251/253 docks additive (V = 0 to 1e-12); 2 contact
docks (V = +1.76, +4.48 — flights can absorb one step into the wall,
still exact; both at large angle).

**Recombination (pair correlation):** P(recombine|freed) = 0.70 pooled
(0.73 from sep 3, 0.63 from sep 2); t_unbind median 6.1e4 (contact is
repulsive, no bound state), t_return median 8.5e4 with tail to 1e7.
Partly finite-box. HINT: aligned docks (< 30 deg) recombine at 13/15 =
0.87 vs crossed 48/72 = 0.67 (~1.6 sigma) — possible directional
retention at the contact vertex; needs a few hundred more episodes.
