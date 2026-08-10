---
name: quanta-strain-heal
description: "Few-quantum edge-degree strain + release on a15: the quantum is 1 tet = 1 sixfold edge; the two strain directions are KINETICALLY INEQUIVALENT (above-native = barrierless dilute flicker gas that heals to 100% pristine registry; below-native = melt or refuse, no dilute window)"
metadata:
  node_type: memory
  type: project
---

**Question (Aaron, 2026-08-10):** make a state with DILUTE defects and then
heal them (back to pristine, or to another edge-legal state). Proposed
instrument: a15 with the mean-edge-degree target ~3 quanta below native and
the f3 target moved to match, so no vertex change is needed; settle, then
restore native pins.

**Tools:** `scripts/defect_dynamics/quanta_heal.py` (two phases in ONE chain,
targets swapped mid-run via `set_num_facets_target` / `set_hinge_degree_target`
— both force an objective recompute, so this is safe), and
`quanta_heal_report.py` (verdict table + figure). Options that matter:
`--ramp-hold` (walk the target down one quantum at a time), `--heal-cimp`
(quench on release), `--contract-split` (the only edge-destroying composite
move). Data under `data/quanta_heal/{pilot,scan,lad,main,chan}`.

**THE QUANTUM (exact).** At fixed f0 a closed 3-manifold has f1 = f0 + f3, so
ē = 6f3/(f0+f3) is a function of f3 alone and its achievable values form a
lattice of spacing 6f0/f1² = 1.317e-4 for a15 m5 — **one tet**. This is the
same quantum the gas campaign calls one forced 2→3-equivalent move: e* sits
**51.72 quanta** below a15's native ē = 5.111111 (campaign: 51.3 forced
moves/1000 verts ✓). At dq quanta below native the EDGE-LEGAL answer is
n5 = 6000 **unchanged** and n6 = 750 − dq: straining by dq = asking the
six-web to shed dq units of line length with the fivefold web untouched.

**FLICKER SPECTRUM OF PRISTINE a15 (exact, all 11500 sites, combinatorial —
`scratchpad/flicker_spectrum_a15.py` logic):** exactly **4 species** of
one-2→3 excitation, (illegal degrees, Σm², share):
(3,4,4) 8 @13% | (3,4,4,4) 14 @9% | (3,4,4,7) 16 @26% | (3,4,4,4,7) 22 @52%.
Verified against the sampler objective to machine precision. This one table
predicts every freeze-out in the campaign: cost = c_eff(2dq+1) + c_imp·Σm².

**THE CENTRAL RESULT — THE TWO STRAIN DIRECTIONS ARE NOT EQUIVALENT.**
A 2→3 creates an edge and a tet; the crystal's elementary excitation therefore
carries **+1 quantum**. Destroying an edge (3→2) needs a degree-3 edge, which
pristine a15 has none of, and no single move makes one *net*. So:

- **ABOVE native (dq < 0) is barrierless** — the pin is paid one flicker per
  quantum, in the CHEAPEST species only.
- **BELOW native (dq > 0) has no dilute window** — c_imp ≥ 0.5 refuses
  (frozen at f3 = 5750 with the pin unpaid: at c_imp 0.6, A = 1.0 the strained
  cell is *perfectly pristine* and 0 accepted moves), c_imp ≤ 0.4 melts
  (n_ill 630–880, percolated at 0.3). Same wall the e*-target campaign hit
  ([[volume-pin-defects]]), now seen at only 3 quanta.

**ABOVE-NATIVE NUMBERS (c_imp 0.9, A = 1.0, B = 30, 4000+4000 sw, 3 seeds):**
dq = −3 → phase A ⟨f3⟩ = 5752.8 (target 5753), n_ill = 14.9 ± 0.3 in 5.9
complexes, top1 = 3.1, deg-7 ≈ 0.01 → **exactly 3 dilute (3,4,4)-type
flickers, expensive species priced out**. dq = −10 → 30.2 ± 1.3; dq = −25 →
71.5 ± 5.4 (n_ill per quantum falls 5.4 → 3.1 → 2.9: multimers get cheaper
per quantum). **Release heals completely at every dq**: first legal state
within 25–175 sweeps, late-window pristine fraction 0.84 ± 0.02 — the same as
the unstrained equilibrium, i.e. no memory of the strain.

**crystal_grains verdict (ref a15, min-size 10): the healed state is the
PRISTINE CRYSTAL, in registry** — 1000/1000 interior-crystalline, one grain,
identical to the reference, from both dq = −3 (strained: 985/1000, one grain)
and dq = −25 (strained: 933/1000, one grain). So the strained states are
genuine dilute-defects-in-crystal (contrast the ±5% strain alloy of
[[volume-pin-defects]], where the registry was destroyed).

**THE BELOW-NATIVE STATE IS REACHABLE WITH A BUDGET, NOT A PRICE.**
`set_illegal_budget(cap)` + c_imp 0.1 (cheap defects, hard cap) walks straight
down to the target: cap 12 → ⟨f3⟩ = 5747.4 ± 0.5, **f0 = 1000.000 ± 0.000**,
14 defect vertices in 3 complexes; cap 30 → 5747.2 ± 0.5, 42. Control (dq = 0,
cap 30) sits at 5750.35, so the strain response is 3.18 quanta. **crystal_grains:
ONE a15 grain, 97.2% / 95.1% interior-crystalline** — dilute defects in an
intact crystal. The SAME strain by price (c_imp 0.4) gives **0/1000
interior-crystalline, zero grains**: budget vs price is the whole difference
between a defected crystal and a melt.

**THERMODYNAMICS vs KINETICS, separated.** Lift the cap on the strained state
and keep the 3-quantum pins at c_imp 0.6 (`--start`, uncap_isostrain): it
HOLDS — ⟨f3⟩ = 5747.8, n_ill = 17.8 ± 3.9, slope −1.5 ± 1.0/ksw, **974/1000
interior-crystalline, one grain**. At that same c_imp a pristine crystal
refuses to strain at all (acceptance 7e-8). So the below-native dilute state is
not thermodynamically forbidden — it is only kinetically unreachable from
pristine. Deliverable state:
`data/quanta_heal/uncap/uncap_isostrain.A.final.mfd`.

**CAVEATS ON THE BUDGET GATE.** (i) A saturated cap is a wall in BOTH
directions: `below_dq3_cap12` phase B deadlocked completely (0 accepted in
4.6e7 proposals) because healing needs a 2->3 and every one would be the 13th
illegal edge. (ii) The capped chain never visits n_ill = 0 at c_imp 0.1 —
entropy fills the budget — so an isostrain quench (c_imp -> 0.4/0.6/0.9 with
the cap on) does NOT empty the reservoir either (n_ill stuck at 13). (iii)
Below-native release shows HYSTERESIS: uncapped at c_imp 0.6 the released
states sit at n_ill 21/64, flat slopes, vs 6.1 for a pristine start — no
return to baseline in 6000 sweeps. Contrast the above-native route, which
heals in 25 sweeps.

**THE CONTRACT/SPLIT CHANNEL DOES NOT TUNNEL IT.** Full 16k-sweep cells at
c_imp 0.7 and 0.9 (prob 0.02): **9.2e5 contract + 9.2e5 split proposals
reached Metropolis, 0 accepted** in each — the composite edge-destroying move
is itself priced out on a dilute crystal, so it is not the escape from the
below-native wall (`data/quanta_heal/chan`). Phase A stayed exactly pristine
(0 accepted moves of any kind); on release to native pins the same chains
immediately behave normally (first legal at sweep 25, pristine fraction 0.83
at c_imp 0.9), which is the internal control showing the phase-A freeze is
caused by the STRAIN, not by c_imp alone.

**EXACT INVARIANT (checked in every run to 2 decimals).**
n6 + 2n7 - n4 - 2n3 = 6f3 - 5f1 = f3 - 5f0, fixed by the pins alone: the
defects can only REDISTRIBUTE the net six-web excess, never change its total.
Measured 747.4 / 747.2 (strained) vs 750.3 (control); the c_imp 0.4 melt has
the same 747.2 with n4 = 412, n6 = 1254.

**PIN CALIBRATION (worth reusing).** At fixed f0 both pins are the SAME
quadratic in f3, joint stiffness c_eff = A + B(f0/f3_ref)²; f0 then follows
along the flat line f0 = f3(6/e_tar − 1), so **A (num_facets_coef) is the only
thing that pins f3** — at the campaign default A = 0.1, σ(f3) ≈ 2.2 quanta,
which swamps a 3-quantum strain. Use A ≈ 1–3. B = 30 pins the line.
Measured strain response in the branches that move: exactly **3.00 quanta**
(c_imp 0.4: 5750.50 → 5747.49; c_imp 0.3: 5750.37 → 5747.37).

Related: [[a15-positive-gas]], [[volume-pin-defects]],
[[crystal-library-gas-campaign]], [[flicker-background]],
[[crystal-grains-tool]], [[fk-amorphous-program]].
