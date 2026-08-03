---
name: species-interactions
description: "lam=0.40 string-vacuum interaction is CONTACT-ONLY — g_bs(r) flat beyond 0.5 cells, directional coupling NULL at 0.5-1.5%, strings grow by absorbing adjacent flicker/fragments; string-string potential unmeasured (sessile = 1 config/run)"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T22:40:41.874Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/species_interactions.py` —
string-vacuum channel: 3544 worldlines, 17193 flicker-birth-vs-string pairs
with a MATCHED null (8 uniform points per birth scored against the same
instantaneous string set). String := age > 15 sw OR nv >= 7 (642 tracks).
2500 sweeps, seed 11 (same trajectory as the lifetime/travel runs), λ=0.40,
slide 0, R m4. JSON: scratchpad `si_long.json`; figure
`data/figs/species_interactions_lam40.png`. PROVISIONAL.

**Interaction is CONTACT-ONLY at this λ:**
- g_bs(r) = 1.00 ± 0.02 everywhere beyond ~0.5 cells, for every string
  species (incl. the immortal decamer (10,22,18,5) and curve-like (7,7,1,0)).
  NO long-range vacuum polarization at the few-% level. Short-r suppression
  (0.23 at r=0.12) is largely definitional — a 2->3 touching a string merges
  instead of registering a birth.
- **Directional coupling: NULL.** <P2> real−null of birth directions around
  the string principal axis (gyration tensor; elongation-gated) is zero at
  the 0.5–1.5% level in all 6 (elong, r) windows. Flicker births are
  isotropic around string axes. First direct test of the "directional
  interaction" hypothesis in the string-vacuum channel: not supported at
  λ=0.40 with slides off.
- **Growth pathway: absorption of adjacent births.** 120 absorption events;
  absorbed = knots (39) + small fragments (3,2,0,0)x33, (2,1,0,0)x9;
  median birth distance from absorber 0.35 cells (90th pct 0.88). Strings DO
  eat the vacuum, but only what lands on them.
- Knots born within 0.6 cells of a string die ~25% sooner (median 1.28 vs
  1.70 sw) — the absorption channel consuming them.
- String-string separations (2345 distinct pairs): consistent with uniform
  (0.88–1.05) but strings are sessile so a run gives ~1 independent
  configuration — NOT a measurement of the string-string potential.

**Reading with [[statics-hu-verdict]]**: crystal-grade HU coexists with
contact-only interactions — the screening that delivers HU is ultra-local
(the string's own halo + merge zone), not a long-range polarization cloud.

**Implication for the ADM program** (see conversation 2026-07-25): the
directional/long-range interaction hoped for in the emerging-gravity picture
does NOT arise spontaneously in the equilibrium no-slide gas. If it exists it
must be (a) at different λ / stiffer EDQ, (b) mediated by the slide/worm
channel (dynamic, not static), or (c) visible only in the string-string
sector — which needs the controlled two-knot experiment (frozen vertices +
directed slides) since natural configurations don't sample it.
