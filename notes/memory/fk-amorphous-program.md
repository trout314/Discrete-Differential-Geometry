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

## 6. The illegality budget is BUILT (2026-08-10) — and it returns but does not travel

**Machinery.** `VertexPotState` now maintains `sumM` (= Σ_v m = 2·#illegal
edges, exact) and `lastDSumM` incrementally through all three potential-delta
functions; `IllegalBudget{cap, blocked}` gates the bistellar, hinge and
contract/split paths in `mcmcStep` with a plain `return false` (an infinite
energy — NOT a re-draw, which would reweight the proposal). C API:
`ddg_sampler_set_illegal_budget` / `..._illegal_budget_stats`; Python:
`ManifoldSampler.set_illegal_budget` / `.illegal_budget_stats`. Requires the n6
potential (it maintains the counter). D unittest asserts `sumM == 2·scan` and
`scan ≤ cap` at EVERY step over 4000 mixed steps at cs 0.25, plus a tight-cap
phase; **meson test 14/14 modules pass**. `--budget` on `fk_amorphous.py`;
new `scripts/fk_worm_harvest.py` harvests and fingerprints the n_ill = 0 returns.

**Clean slides are allowed alongside the budget** (dirty ones, nonlocal slide
and worm are refused, both directions): a clean slide preserves the multiset of
illegal degrees over its changed edges, hence n_ill globally, so it cannot
breach the cap — and it is the only transport move the reservoir has.

**MEASURED on the pristine R crystal (m2, f0 1272), 150–250 sweeps/cell:**

| B | μ_ill | returns to n_ill=0 | distinct legal states |
|---|---|---|---|
| 6 | 0.0 | 0 (fills to the cap) | — |
| 6 | 1.0 | 9 | **1** |
| 12 | 1.5 | 24–41 | **1** |
| 24 | 2.5 | 30 | **1** |
| 48 | 1.8 | 49 | **1** |
| 96 | 1.5 | 40 | **1** |
| 200 | 1.2 | 29 | **1** |
| 12 | 1.5 + clean slide 0.3 | 46 | **1** |

- **A flat reservoir (μ = 0) fills to the cap and never returns** — entropy
  owns the interior. A fugacity is still needed to make n_ill = 0 weighty.
- **Every return is to the SAME state.** Excursions retrace exactly. Fully
  consistent with [[defect-travel]] (all defects caged) and §2 (R has zero
  non-positive-damage contractions): the crystal is a deep trap.
- The clean-slide channel fired **2 times in 250 sweeps** — at μ 1.5 the chain
  is at n_ill = 0 for 92% of steps (mean 0.24, max 3), so there is almost never
  a degree-3 chord to slide. Fugacity and transport fight: cheap defects fill
  the reservoir, expensive ones never persist long enough to move.

**NEXT (clear from the above): a DIRECTED worm, not a thermal one.** Decouple
the worm's existence from its price — insert exactly one defect, slide it many
times while it is held open, then close it, with the insert/close pair carrying
the Hastings factors. That is the standard resolution of exactly this tension
and the machinery (clean slide + budget + maintained counter) is now all in
place. Also worth trying: start the harvest from a disordered legal state
rather than a crystal (none exists yet — chicken-and-egg, so grafting first).

## 7. MAXIMISING f_FK AT f_reg ≈ 0 (2026-08-11): hold μ at the threshold

**Goal (Aaron):** maximise the FK-legal VERTEX fraction (every incident edge
deg 5 or 6 — hubs explicitly allowed) while keeping f_reg ≈ 0. Any dispersed
defect population acceptable.

**Reframing that makes it tractable:** f_reg ≈ 0 is nearly FREE if you start
amorphous — melts do not recrystallise ([[quanta-strain-heal]] hysteresis) and
crystals are measure-zero. So this is a PREPARATION problem from a disordered
seed, not a sampling one. Start = the a15 c_imp 0.4 melt (f_FK 37.3%,
f_reg 0.00%). `fk_amorphous.py --ratchet-slack` added (cap = running-minimum
illegal-edge count + slack).

**RESULT (8 cells, 40k sweeps, zleg 0, cs 0, cn 0.1, nh 30, e* target):
f_FK 37.3% → 64.4% with ZERO crystalline grains** (1.59% interior-crystalline,
f_reg/f_FK = 0.025). Z̄ held at 13.397–13.401 against the FK-window target
Z̄* = 13.39733 for the whole run — the flat pin does exactly the job the
identity says it should. Hubs stayed at 78–87 (strict fk-coord 56.7%, so hubs
cost ~8 points). Z host final: f_FK 35.4% → **60.5%**, Z̄ 13.3977 / ē 5.10433 (targets
13.39733 / 5.10430 — the flat pin holds the composition essentially exactly),
0 grains — BUT 11.4% interior-crystalline in sub-threshold clumps, so
f_reg/f_FK = 0.19 against 0.025–0.049 on a15: the Z glass keeps markedly more
local registry than the a15 one. Prefer a15 as the amorphous seed.

**TWO NEGATIVE RESULTS THAT MATTER:**
- **The ratchet is NOT the active ingredient.** No-ratchet control 64.38% vs
  slack-5 64.09% — indistinguishable, and slack 5/20/50 span only 62.4–64.1%.
  At the working μ the fugacity already prevents the melt, so the cap never
  binds. (It IS needed at low μ: at μ = 0.1 this state blows up 491 → 923
  illegal edges in ONE sweep.) Keep it as a cheap safety net, not a driver.
- **RAMPING μ UP IS ACTIVELY HARMFUL** — the opposite of the driver's default
  design and of annealing intuition. μ = 2 const 64.4% | μ 2→4 56.3% |
  μ 2→6 54.7% | μ = 3 const 53.2%. Every stiffer cell has late slope
  **0.0 ± 0.1 /ksw (total arrest)** while the μ = 2 cells are still descending
  at −1.1 to −1.3 /ksw at 40k sweeps. **μ ≈ 2 is simultaneously the descent
  THRESHOLD (0.5 and 1.0 just pin at the cap) and the OPTIMUM.** Hold it;
  do not anneal it.
- Cost note: cs = 0.05 costs **11×** the wall time (17.1 vs 1.5 s/100 sweeps)
  for a marginal gain at ~500 illegal edges. Save the channel for the endgame
  where 2→3 can only nibble; the bulk descent is cheaper with cs = 0.

Figure: data/fk_amorph/ratchet/fk_ratchet.png; data under data/fk_amorph/ratchet.

## 8. THE CHANNEL IS THE WHOLE GAME (200k sweeps, both hosts, data/fk_amorph/long)

After the D optimisation made cs affordable (~3.1× on the channel: PairTable
1.11 × DeltaMap 1.93 × countCycles memo 1.43), six chains at 200k sweeps,
μ = 2 const, zleg = 0, **cs = 0.05 vs cs = 0 MATCHED CONTROLS** from the same
start and seed:

| | a15 ctrl | a15 cs.05 ×3 seeds | z ctrl | z cs.05 |
|---|---|---|---|---|
| f_FK | 69.44% | **86.11/86.41/86.61%** | 66.03% | **82.13%** |
| illegal edges | 253 | 93/95/103 | 425 | 204 |
| hubs | 84 | 68/72/78 | 120 | 97 |
| f_reg | 2.5%, 0 grains | 5.0–5.2%, **0 grains** | 15.2%, 0 grains | 19.8%, 0 grains |

**+17 points on a15, +16 on Z; three seeds agree to 0.5 points, so it is a
plateau not luck. f_reg stayed at ZERO GRAINS throughout — legality and
crystallinity did NOT rise together, which was the real risk.**

**MECHANISM, sharpest form: the Pachner-only control removed ZERO deg-3 edges
in 200k sweeps (62 → 63) while the channel cut them 5-fold (62 → 11–22).**
deg-4: 231 → 185 (control) vs 71–91 (channel). §3 confirmed at length on two
hosts — a 2→3 nibbles one quantum and cannot clear a deg-3 chord; a
contraction takes 4–5 illegal edges at once and eats the structures deg-3
lives in. Hubs FALL with the channel and RISE without it, so it improves
FK-coordination and edge-legality together (unlike the stiffening ramp).

Composition at the end: Z12-rich 485:150:93:73 at Z̄ 13.3988 vs the flat
target 13.39733 — an amorphous FK phase structurally unlike A15's 250:750:0:0.
CAVEAT: Z keeps 4× more local registry than a15 (f_reg/f_FK 0.19 vs 0.06);
prefer a15 as the amorphous seed.

**The states are now DILUTE**: a15_cs05_nr at sweep 160k has 37 complexes of
which 21 are the minimal 2-vertex object, largest 12 — so the dilute-gas
species instruments apply to an AMORPHOUS host for the first time.
Catalog: data/viz/a15_cs05_nr/index.html. Census tool:
`scripts/defect_dynamics/amorph_census.py`.

OPERATIONAL: multi-hour runs launched inside a harness background task get
KILLED when the task is reaped — use nohup + disown. `--snap-every` saved this
one. Resumes are NOT bit-continuations (fresh RNG, ratchet record re-inits),
and `--f3-target` must be passed explicitly or the driver re-pins to the
snapshot's own f3.

## 9. IS f_FK = 100% COMPATIBLE WITH BOTH PINS? YES, up to an irreducible bit

All edges in {5,6} ⟹ every link is a (5,6)-triangulation of S² (fullerene
dual), Euler forces exactly 12 deg-5 edges per vertex, so the ENTIRE legal
f-vector lattice is two integers (f0, E6):

    E5 = 6 f0,  f3 = 5 f0 + E6,  f1 = 6 f0 + E6,  Z = 12 + n6(v).

**e\* is IRRATIONAL and ē = 6f3/f1 is rational, so the flat pin is never
exactly zero for ANY triangulation, legal or not** — a property of the target,
not of legality. Shifting E6 by 1 moves ē by exactly one quantum 6f0/f1², so
E6 can always be chosen to land ē within HALF a quantum of e\* while f3 hits
the volume target EXACTLY. Best points: a15 m5 (1000, 699) → +0.33 q; a15
f0=1008 (1008, 704) → −0.26 q; z m6 (1512, 1056) → −0.38 q. Flat-pin penalty
at nh=30 is then 0.06–0.14, against ~4–8 for one flicker.

LOCAL CONSTRAINT: **n6(v) ≠ 1** (no fullerene has exactly one hexagon), so
admissible n6 is {0,2,3,4,5,…}; n̄6 ≈ 1.397 needs a mixture, which sits
comfortably between C15 (1.333) and A15 (1.5). VERIFIED IN DATA: n6 = 1 occurs
ZERO times in every state censused.

NOTE the 200k a15 chains inherited f3_ref = 5746, whose best legal point is
+1.74 q off e\*; **f3_ref = 5744 is the matched target** (costs ~0.9 vs 0.06 in
flat penalty — harmless, but worth fixing for a definitive run).

Related: [[quanta-strain-heal]], [[fk-move-search]], [[worm-sampler-program]], [[defect-travel]],
[[crystal-library-gas-campaign]],
[[volume-pin-defects]], [[six-web-gauge]], [[edq-only-melting]].
