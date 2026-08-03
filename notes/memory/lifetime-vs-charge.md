---
name: lifetime-vs-charge
description: Defect lifetime is NOT monotone in charge; the ΣZ=70 rung is anomalously short-lived and reactive; Q_c is exactly conserved along a fixed-shape knot worldline
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T20:37:20.843Z
---

Measured 2026-07-25. Run: `scripts/defect_dynamics/lifetime_charge.py`, 2500 sweeps,
1 chain, seed 11, from `data/mgas/lam40x_snap15000`, R crystal m4 (N₃=57992),
λ=0.40 full action (EDQ e*=5.105025 + 0.6·U(n₆) + 1.0·m²), slide prob 0.
JSON: scratchpad `lc_long.json`; figure `data/figs/lifetime_vs_charge_lam40.png`.

**Answer: no monotone charge→lifetime law.** At fixed shape (the (3,4,4) knot,
f=(5,10,9,3), n=2375 died / 2 censored):
- Spearman ρ(ΣZ, lifetime) = −0.056, permutation p=0.007 — but drop the single
  dominant rung ΣZ=70 and it flips to **+0.038 (p=0.16)**. The monotone
  correlation is one rung, not a trend.
- Kruskal–Wallis across rungs 67–72: **H=48.8, 5 df, p≈2e-9** — the rungs really
  do differ, non-monotonically.
- Median lifetimes (sweeps): 67→1.007, 68→1.048, **69→1.360 (longest)**,
  **70→0.829 (shortest, n=1015 = 43% of sample)**, 71→1.125, 72→1.283.

**ΣZ=70 is the anomalous rung** (Q_c ≈ −1.661): most abundant, shortest-lived,
AND most likely to change shape mid-life (29% vs 8–18% at other rungs). Reads as
a shallow, easily-formed, easily-destroyed generic state.

**Charge is exactly conserved along a knot worldline.** All 2377 pure-shape
tracks had mode_frac = 1.0 on the joint (shape, ΣZ, n) key — ΣZ never changed
while the induced shape stayed the knot. Q_c changes only through shape change.

**Caveats:** 22% of knot tracks excluded for changing shape, and that cut is
itself rung-dependent (29% at 70), so the fixed-shape comparison is not cleanly
randomised. Narrow dynamic range (medians 0.83–1.36 sweeps) — this is the
transience regime, see [[mobile-gas-liquid]]. Single chain, uncertified.

**Why the confound matters:** raw lifetime-vs-charge just re-measures
lifetime-vs-size, since Q_c ≈ −0.33·n. Raw bins run 1.5 sweeps at n=5 to 117 at
n=12. Only the fixed-shape ladder tests charge on its own. See
[[defect-species-ladder]] for the ladder, [[five-illegal-knot]] for the knot.

**Two of my own errors on this measurement** (both caught by inspecting the
bins, not the correlation): a pilot version recorded shape and coordination as
independent modes, admitting 5-vertex "knots" with ΣZ=150 and yielding a
spurious r=+0.773; and the first long run never ran — `nohup … &` inside a
backgrounded wrapper died with the wrapper, while the task reported exit 0.
