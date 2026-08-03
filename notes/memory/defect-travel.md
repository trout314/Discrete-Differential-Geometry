---
name: defect-travel
description: lam=0.40 no-slide transport — bigger complexes range farther (monotone in size) but ALL are caged subdiffusive; knots die at their birth site; long-lived complexes return home
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T21:59:53.942Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/defect_travel.py` —
per-worldline transport at accepted-move resolution, centroid via the D-side
incremental cocycle lift (complex-centroid tracer channel). 2500 sweeps,
seed 11, SAME trajectory as [[lifetime-vs-charge]] (verified by matching
worldline counts sweep-for-sweep), so travel and lifetime join per track.
JSON: scratchpad `travel_long.json`; figure `data/figs/defect_travel_lam40.png`.

**Do larger, longer-lived complexes travel farther? Yes in range, no in
transport.**

- Median max excursion is monotone in size: 0.000 cells at nv<=5 (2928 knots
  and everything smaller: EXACTLY zero — flicker dies where it is born, 87%
  of knots with net identically 0) rising to ~0.37 cells at nv=10.
- But motion is CAGED: largest excursion of all 3417 died tracks is 1.26 cells
  (box = 4 cells = 22 vertex spacings; 159 verts/cell so 1 spacing = 0.185
  cells). Median large-complex range is 1.0–2.5 vertex spacings. Nothing
  crosses even a third of the box. The 1191-sweep immortal decamer ranged only
  1.07 cells.
- maxexc^2 vs lifetime exponent: +0.64±0.03 (all), +0.46±0.05 (nv<=6),
  +0.63±0.07 (nv>=7); 1 = diffusive. Subdiffusive everywhere.
- net^2 vs lifetime exponent is NEGATIVE (−0.54±0.36): long-lived complexes
  end up BACK near their birth site. Range grows with life; net displacement
  does not. Pinned oscillation, not migration.

So at λ=0.40 with slides OFF the constrained liquid has NO defect transport
channel at these scales — motion is rattling in a cage of ~1 unit cell. This
sharpens the old "defects pinned over 100k sweeps" observation into a
per-track statement with a length scale. The obvious next experiment is the
same run with slides on (the slide channel was built precisely to unlock
this; see [[worm-sampler-program]]).

Protocol notes: centroid jitter from membership churn is handled by using net
and maxexc (jitter averages out) and NEVER path length (see run5h gotchas in
[[defect-kinetics-run5h]]); unwrap by minimal image between consecutive
per-move samples (steps << box/2 always). check_cocycle_positions audited
every 8 chunks — clean.
