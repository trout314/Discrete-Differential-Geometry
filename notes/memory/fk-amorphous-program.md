---
name: fk-amorphous-program
description: "Amorphous-FK program (2026-08-10): ebar = 6 - 12/Zbar is an IDENTITY once edges are legal, so mean edge degree = six-web density, not a free knob; PROOF + measurement that edge contraction also cannot map FK-legal -> FK-legal (it always mints one Z17+ hub); but on DISORDERED states contraction carries d_ill = -4/-5 moves that no Pachner move can"
metadata:
  node_type: memory
  type: project
---

**Question (Aaron, 2026-08-10):** revisit the old VDV-anneal experiment — make
AMORPHOUS triangulations with n_ill ≈ 0 and ⟨e⟩ near flat, now that the
contract/split (edge-collapse) channel exists. The old anneal produced all-{5,6}
edges but *not* the usual FK vertices.

## 1. The identity that reframes the problem

All edge degrees in {5,6} ⟹ the link sum rule Σ_{e∋v}(6−deg e) = 12 forces
**exactly 12 degree-5 edges at every vertex**, so

    Z(v) = 12 + n6(v),   E5 = 6 f0,   f3 = 5 f0 + E6,   f1 = 6 f0 + E6

and therefore

    ebar = 6 − 12/Zbar,   Zbar = 12 + n6bar,   n6bar = 2 E6/f0.

**⟨e⟩ is not an independent knob — it IS the six-web line density.**
Flat ⟨e⟩ = e* = 5.1042993 ⟺ Zbar* = **13.39733** ⟺ n6bar* = 1.39733, which sits
inside the FK window [13.3333 (C15), 13.5 (A15)]. "A range of ⟨e⟩ near flat" =
"a range of Z12:Z14:Z15:Z16 compositions". Also explains the old failure mode:
n_ill = 0 gives fullerene links automatically, but nothing stops n6 ≥ 5 —
Z17/Z18 **hubs** are edge-legal and non-FK. That is what the VDV anneal made.

## 2. Contraction cannot preserve FK-legality either (PROOF + measured)

Contract edge uv of degree r, ring x_1..x_r. Exactly three changes:
ring edges deg → deg−1; spoke pairs merge to deg(u,x)+deg(v,x)−4; uv deleted.
(Local rules **verified against literal contraction**, 5/5 sites, A15.)
Edge-legality of the result forces ring edges 6→5 and spokes 5+5→6 (5+6−4=7 and
6+6−4=8 are illegal). Then all r merged spokes at w have degree 6, so
n6(w) ≥ r = deg(uv) ≥ 5: **w is a Z17+ hub**. By reversibility no split works
either. So the FK-legal set is still dust under the *full* move menu
(2↔3, 1↔4, hinge, slide, contract/split) — [[fk-move-search]] extends verbatim.

**Measured (`scripts/defect_dynamics/fk_channel_census.py`, whole library):**
best contraction damage (Δ#illegal edges, Δ#non-FK vertices) from the pristine
crystal —

| a15 | c14 | c15 | c36 | delta | mu | p | r | sigma | z |
|---|---|---|---|---|---|---|---|---|---|
| (4,7) | (2,4) | (3,5) | (2,4) | (2,5) | **(0,1)** | (2,4) | (1,3) | (3,5) | **(0,1)** |

Only **μ and Z** have zero-illegal-edge contractions (1.15% / 2.13% of edges,
all r=6 → a Z20 hub); the other eight crystals cannot move at all without
minting illegal edges. Requirement is a 5- or 6-cycle *in the six-web* that is
the link of an edge — A15's web is straight lines, C15's diamond rings are not
edge-links. On μ/Z the one-hub state has **exactly 1** FK-restoring split (the
inverse), so a single-hub excursion on a crystal is rigid: crystals are traps.

## 3. What IS new: contraction is a cluster annihilator on disordered states

Damage is a *delta*, so on a defected state a move can be strictly downhill.
`fk_channel_census.py --state`:

| state | f_e | f_FK | free moves (Δ_ill≤0 ∧ Δ_nonFK≤0) | best Δ_ill |
|---|---|---|---|---|
| R crystal (pristine) | 100% | 100% | **0** (min damage (1,3)) | +1 |
| c15 +5% strained alloy | 94.9% | 46.4% | 105 (1.00%) | **−4** |
| a15 debt0 nh30 | 98.1% | 81.2% | 73 (1.11%) | **−5** |
| c15 debt0 nh30 | 99.0% | 87.2% | 26 (0.26%) | −2 |

A 2→3 can only nibble one defect quantum; a contraction removes 4–5 illegal
edges in ONE step. This is the lever the VDV anneal and `glass_quench` (which
arrested at f_e ≈ 0.52) never had.

## 4. Driver + protocol lesson

`scripts/fk_amorphous.py` — volume pin + flat pin (nh) + zleg hub penalty +
imp_lin illegal-edge fugacity, f0 free through contract/split, linear coupling
ramp, Recorder-based. **Lesson from the first pilot: do NOT ramp up from soft.**
Starting at zleg 0 / mu 0.5 with cs 0.25 melted the a15 nh30 state from 128 to
3189 illegal edges in 100 sweeps and it never recovered (f_FK 2.6% at the end).
The near-legal states are the precious resource — only ever tighten. Also
cs = 0.25 costs ~110 μs/move (≈40× a plain sweep); keep cs ≲ 0.05.

## 5. A/B verdict: the channel unfreezes the state, blind proposals then stall

Single chains from `a15_debt0_nh30.final.mfd` (f0 1000, f3* 5700, 128 illegal
edges, f_FK 81.2%), S = c_N(f3−f3*)² + 30(f1−6f3/e*)² + zleg·Σdist²(n6,FK) +
1.0Σm² + μ·Σm, seed 4001, 3000 sweeps. PROVISIONAL (1 chain/arm).

| arm | n_ill 0 → 3000 | sd | f_FK final |
|---|---|---|---|
| cs = 0, c_N 0.1 | 128 → 119 | **1.26** | 83.9% |
| cs = 0.05, c_N 0.1 | 128 → **106** | 4.78 | 84.2% |
| cs = 0.05, c_N 0.01 | 128 → ~118 (no gain) | — | 83.8% |

- **cs = 0 is exactly frozen**: sd 1.26 against √⟨n⟩ = 11.0 (8.7× sub-Poisson),
  hubs pinned at *exactly* 11 for 3000 sweeps. Not an equilibrium gas — a glass.
- The channel is what makes the population move at all, but its acceptance is
  **0.04% (contract) / 0.16% (split)** with blind uniform proposals, and accepts
  flatten after ~250 sweeps. The census says ~1% of edges are Δ_ill ≤ 0 sites,
  so the proposal — not the Metropolis factor — is the bottleneck. Softening the
  volume pin (which charges ≈ +2.4 for a Δf3 = −6 contraction) did NOT help,
  killing that diagnosis. **Fix = locally-informed proposals: restrict/weight
  contraction proposals to defect-adjacent or six-web-ring-rich edges with the
  exact Hastings correction.**
- **Stiffening ramp (zleg 2→8, μ 2→12 over 1950 sw, cs 0.05, seed 4002) is the
  only arm that moves — then arrests dead.** n_ill 128 → 75, f_FK 81.2 → 89.6%,
  **hubs rise 11 → 14** (the system pays the stiffening by converting illegal
  edges into Z17+ hubs — §2 in action, and exactly the old VDV endpoint: the hub
  sector is where the dynamics live). From sweep 2250 on, n_ill/hubs/imp/f0/f3
  are **bit-identical for 750 sweeps** — total arrest. The channel accepted only
  18 contractions + 44 splits in the whole run (0.006%), so this annealing was
  done by the ordinary 2↔3 sector, *not* the channel.
- ē drifts OFF flat as legality stiffens: 5.1052 → 5.1074 (Z̄ 13.415 → 13.444).
  Removing deg-4 (positive-currency) defects at fixed f3 lowers f1 and raises ē;
  the nh=30 flat pin loses to zleg/μ. The two dials genuinely compete.
- **Conclusion: fugacity is the wrong knob.** Defects are the transport medium,
  so pricing them out freezes the chain before it reaches n_ill = 0. What is
  needed is a *budget*, not a price: cap n_ill ≤ B (B ~ 5–20) so a fixed small
  reservoir keeps circulating, and harvest snapshots at n_ill = 0. Because
  S is **exactly degenerate** on the legal manifold at fixed (f0, f3) — every
  term is constant there — conditioning a budgeted chain on n_ill = 0 samples
  the *uniform* measure on legal states, in which crystals are measure-zero and
  amorphous states dominate. The target ensemble is maximum-entropy, so the
  entire problem is the move set, not the action.

Related: [[fk-move-search]], [[crystal-library-gas-campaign]],
[[volume-pin-defects]], [[six-web-gauge]], [[edq-only-melting]].
