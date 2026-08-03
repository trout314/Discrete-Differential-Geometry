---
name: defect-site-classes
description: "C15 has 4 inequivalent 2->3 sites, R has exactly 102 (6 ledgers) -- but the defect's OPTICS are identical at all of them, because the ball is convex"
metadata: 
  node_type: memory
  type: project
  originSessionId: 93ac5363-1537-4e86-9a6a-0dff5e3b1ea4
  modified: 2026-08-02T19:59:26.110Z
---

Measured 2026-08-02. Tools: `scripts/move_site_census.py` (classification) and
`scripts/defect_lens_field.py` (optics).

**Site classes.** Every face of a pristine TCP crystal is a legal 2->3, but the
sites are inequivalent. Exact counts via `symmetry.CrystalSymmetry` face orbits:

- C15 m3: |Aut| large, **4** classes, **4** distinct curvature ledgers (the
  ledger is a complete invariant there).
- R m3: |Aut| = 486, **102** classes, only **6** ledgers.

Ledger space is exactly 2 x 3: face edges in {(5,5,5), (5,5,6)} times spokes in
{six 5s, five 5s + one 6, four 5s + two 6s}. R realizes all six, C15 four. No
site has two degree-6 face edges or three degree-6 spokes -- the FK network is
too sparse. Worked example from earlier sessions, face (1,69,70) in C15, is the
17.6% MINORITY class, not a typical one.

**Universality of the optics -- the headline.** The boundary transit gain
A(p,q) = d_R(p,q) - d_D(p,q) on dB, measured in the FULL host so the option of
routing AROUND the ball is already taken, is IDENTICAL at all 102 inequivalent
R sites: max|A(site) - A(site_0)| = 5.6e-16 over the whole matrix.
A_max = 0.632993, A_mean = 0.002456, A_min = -0.488034, A_frac = 21.67%.

Reason: **B is convex in both configurations.** Every boundary edge has interior
dihedral theta (70.53 deg) or 2*theta (141.06 deg), both < 180 -- and the
pattern is EXCHANGED between R and D (R: face edges 2*theta, spokes theta; D:
the reverse), the sheet<->strut trade again. So geodesics between dB points
never leave the ball: max|d_full - d_ball_confined| = 2.2e-16, exterior strictly
helps on 0.000% of pairs. The around-route never binds, d|_dB = d_B, optics are
universal.

So site inequivalence is a CURVATURE/TOPOLOGY statement (illegal signature
(3,4,4) vs (3,4,4,4); how many new permeable degree-6 lines get grafted onto the
disclination network), NOT a metric one.

**GOTCHA that cost a measurement.** The global field (Delta = d_D - d_R from
random sources to all vertices) is dominated by source-placement noise: the
capture cone is narrow, so with a few dozen sources the statistic is set by
whether a source lands in it. Cross-sections ranged 0.01%-0.88% with the
WITHIN-class spread equal to the between-class spread -- and within-class spread
must be zero, since translation copies are geometrically identical. Never read
site dependence out of that statistic; use the canonical dB-sourced A instead.

See [[defect-boundary-map]] and [[pl-geodesic-permeability]].
