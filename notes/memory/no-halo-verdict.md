---
name: no-halo-verdict
description: DECISIVE — there is NO halo; melt = pristine crystal EXACTLY beyond ~0.4 cells; knot excess charge = exactly one flicker quantum; NO screening cloud exists (HU mechanism revised to dilution + volume pin); static interactions are rigorously contact-only
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T23:14:26.002Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/defect_halo.py` —
per-vertex differencing against the pristine reference via exact site
matching (cocycle lift positions are crystal-registered; gauge offset removed;
matching 100.00% on all 11 λ=0.40 snapshots). Figure
`data/figs/defect_halo_lam40.png`; JSON scratchpad `halo_lam40.json`.

**There is no halo.** P(class mismatch vs reference site) and Δq(v) are
EXACTLY ZERO (not small — zero, with up to 22k vertices/bin) for every legal
vertex beyond ~0.4 cells (≈2 vertex spacings) of any defect vertex. Knots
disturb NOTHING outside their own 5 vertices — even their contact shell is
pristine-classed. Strings have only a contact skin (P(mm) 0.137 at d=0.1,
0.007 at 0.3, 0 beyond), ~50% thicker along the string axis than equatorial
(0.133±0.020 vs 0.089±0.009, 2.1σ — the only anisotropy anywhere).

**Knot charge quantization:** dQ_cx (complex charge minus its sites' reference
charge) = −1.103 ± 0.000 for ALL 34 knots = exactly the flicker quantum
2(π−3θ) = −1.102571. Defect excess charge is volume-excess charge.

**NO screening cloud exists.** Legal-matrix excess charge: −0.10 rad/snapshot
vs defect excess ≈ −10 rad/snapshot (≈ ΔN3 × quantum, the global topological
excess of the extra tets). REVISES [[statics-hu-verdict]]'s mechanism claim:
crystal-grade HU at R≤2 is NOT active screening — it is (i) dilution (defect
point-process variance ~3 rad² vs lattice surface term ~180 at R=2) plus
(ii) at the largest scales, the VOLUME PIN: Q_defect_total ≈ quantum × ΔN3,
so c_N(N−N*)² directly suppresses global charge fluctuations. PREDICTION for
bigger boxes (m6/m8): an unscreened-gas R³ shoulder should appear between the
lattice regime and the pin-dominated regime (crossover est. R* ~ 30 cells at
this density).

**Rigorous consequence for interactions:** the action is local (sums over
edges/vertices), and two-defect configurations are EXACTLY crystal between
the skins ⇒ static energy is additive, and flicker fluctuations decorrelate
at contact ⇒ entropy additive too. The static two-knot potential is exactly
zero beyond skin overlap (~1 cell center-to-center). The planned
fixed-separation two-knot test at r > 2 spacings is predicted null — the
user's worry (2026-07-25) is confirmed and sharpened to a theorem-grade
statement. The Debye/λ-sweep idea is also moot: there is no field to screen.

**Where interaction physics actually lives now:** (a) contact/junction rules
(binding, absorption — 120 events measured in [[species-interactions]]);
(b) the DYNAMICAL web-mediated channel — slides moving knots along BC
helices, collision/scattering physics (the worm/momentum program is now the
ONLY route to action-at-a-distance, not just the momentum sector);
(c) near-melting criticality remains untested. For the ADM program: GR-like
structure must emerge from transport + conservation laws (hydrodynamics of
the defect gas), not static two-body forces — arguably MORE ADM-like, since
Einstein equations are local conservation statements.
