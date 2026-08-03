---
name: reaction-census-lam035-m4-results
description: 8-chain λ=0.35 R m4 reaction census RESULTS — clean chemistry on a non-stationary (melting) substrate
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-29T12:35:48.443Z
---

The 8-chain reaction-census campaign ([[reaction-census-campaign]]) finished
2026-07-29 (data/rxn_lam035_m4/, c0–c7, each 16k–30k sweeps, all 5h/18000s,
exit 0). All 8 chains **audit-clean for their full length** (0 audit_failures,
no truncation) — the [[knot-gas-campaign-lam035-m4]] audit-bug fix holds; every
event usable. Analyzed with scripts/defect_dynamics/reaction_report.py
--dir data/rxn_lam035_m4 --skip-sweeps 2000. Frame: none (graph identity).

Design: c0–c3 slide=0.10, c4–c7 slide=0.0; even chains start dilute
(lam35r_snap15000), odd start over-defected (_m4_over).

**LOAD-BEARING CAVEAT — not stationary.** Every chain's complex number drifts
UP: late-window (last 50%) slope of n_components +0.30 to +1.54 /ksw, 6/8 above
2σ (block-bootstrap), ALL central slopes positive. Both starts climb — they do
NOT bracket. This is the first-order melting drift: λ=0.35 m4 has no stable
dilute gas (confirms [[edq-only-melting]]). So all rates below are
QUASI-stationary on a drifting substrate, NOT certified equilibrium.

**Chemistry (quasi-stationary):**
- Association K = k_merge/k_split ≈ 0.98, tight across all 8 (0.96–1.00) —
  near ideal/non-interacting, consistent with static two-body potential = 0
  beyond contact ([[no-halo-verdict]]).
- Reduced rates roughly intensive: k2 (merge/⟨pairs⟩) ≈ 0.0016–0.0028,
  k1 (split/⟨n⟩) ≈ 0.027–0.040 (consistent to ~1.6× → bimolecular-merge /
  unimolecular-split picture holds).
- Transport test k2(slide-on)/k2(off) = 0.79 (<1): slides did NOT accelerate
  merging. Detailed-balance slides ⇒ diffusion-limited would push k2 UP; it
  went down. So reaction/contact-limited, not diffusion-limited (matches
  contact-only caged transport, [[defect-travel]]/[[species-interactions]]).

**Lifetimes/species:**
- Solo flickers median ~3.3 sw (the (1,0,0,0) monomer blink). Fused compounds
  (absorbed ≥1 partner) median ~15–24 sw, mean ~47–65 — 5–7× longer,
  heavy-tailed to ~10³ sw. Compound-survival curves COLLAPSE to one universal
  law across slide/start. Quantifies run5h "monomers blink, multimers
  long-lived" ([[defect-kinetics-run5h]]).
- Growth is monomer-accretion-dominated: top merge channel is absorbing a
  (1,0,0,0) monomer (products size 6–13); knot (5,10,9,3) absorption next.

Figure: data/rxn_lam035_m4/reaction_report.png.

NEXT (open): complex SHAPE / aspect ratio — user recalls complexes are long &
skinny (high aspect ratio); find the shape script and check this dataset.
Also open: deg-4 worm move ([[nonlocal-slide-move]] is the deg-3 analog).
