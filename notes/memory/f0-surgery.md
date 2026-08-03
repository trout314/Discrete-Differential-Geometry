---
name: f0-surgery
description: "f0 is a hidden ensemble parameter (1<->4 kinetically dead everywhere); FIRST vertex removal achieved: 12-move collapse +23.4 then 4->1 refunds -35 => NET -11.6; edge_removal.py surgery primitive 32/32; dressed composites unfreeze the f-lattice"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-31T12:55:51.230Z
---

**Discovery chain (2026-07-30, m^2-only C15 gas + pin quench):**
1. Down-quench of the f1 pin parks at gap +10 NOT because target states
   don't exist (f0=1522 satisfies both pins, ~0.13 penalty) but because
   f0 is KINETICALLY FROZEN: 4->1 needs a Z=4 vertex (never forms;
   1->4 births 4 deg-3 edges, never accepts). Every ensemble ever run
   was implicitly microcanonical in f0. Fast sector (f3 via 2<->3)
   equilibrates to the constrained-at-fixed-f0 optimum EXACTLY
   (predicted f3=8722/gap+10.5; observed 8722-8726/+10).
2. Worm moves: ZERO effect on quench relaxation (f-neutral, ~2700
   accepted hops, both directions) -- sector orthogonality measured.
3. **Vertex removal EXISTS and is favored**: staged DFS (depth-4 stages,
   cage moves, exact rollback -- greedy dead-ends) collapsed a Z=10
   vertex to Z=4 in 12 moves (+23.40 = barrier), 4->1 refund -35.02,
   NET dS = -11.62. ~14 removals would close the whole +10 gap.
   Sampler capi FORBIDS 1<->4 through targeted moves -- complete the
   4->1 at Manifold level on a dup, reprice with fresh sampler.
   Script: scratchpad/vertex_removal_v2.py; state vrm_removed_657.mfd.
4. **edge_removal.py (committed de7cc4b)**: targeted edge surgery,
   32/32 across d0=3..6 on the quenched gas, 2-4 s each; d0=3 removal
   DOWNHILL (-5.93); deg-4 removal (old 0/97 impossibility) trivial
   without content-neutrality; costs crystallographically quantized.

**Design principle (user):** the move menu must not be curated by what
fires under our usual ensembles (hard proposal-starvation != soft dS
suppression). Complete bilocal/dressed menu = unconditional generators
for each simple Delta-f consistent with chi=0; on T^3 only Delta-f0 != 0
needs engineering (thermal bath supplies f3). Dressed generators are
single-region; bilocal PAIRS of complementary halves are f-neutral
(transport) -- pairing's role for f-movers is cost-sharing (uphill half
financed by downhill half, e.g. deg-3 removals at -5.9 in quenched
states). CORRECTED (user): f0 generators do NOT change the ensemble —
the Pachner alphabet is already ergodic across f0 (1→4 finite dS) and
the pins legislate f0 softly; generators change MIXING only (rare
fluctuation e^(−23) → one O(1)-acceptance proposal). Prior
certifications = correct metastable measurements conditional on their
frozen f0 slice; may relax to different f0 once the menu completes.

**LINK/COLLAR CALCULUS + TRANSMUTATION (2026-07-30, commit fcd6fc0):**
moves at v = 2D bistellar moves on lk(v); vertex removal = Steinitz
reduction of the link sphere. Every link flip has TWO 3D flavors
(2→3 expel / 3→2 absorb — surface dynamics is volume-mediated); chords
= blockers; planning state = (link, collar). FIRST Tier-1 deep
composite VALIDATED: 4-move octa-vertex ↔ deg-4-edge transmutation
(dressed_generators.edge_to_vertex/vertex_to_edge): precondition-free
at every deg-4 edge, exact round trips, dS antisym ±48.062 (m² gas,
orbit-identical). Docking architecture: collapse to octahedron →
transmute → hand deg-4 quantum to the gas. **PLANNER BUILT + VALIDATED (commit 290e765,
scripts/defect_dynamics/link_planner.py):** 2D best-first over
(link faces, ambient degrees, m-counters) with flavored flips +
deletions, lift to 3D. Planned dS == executed dS to MACHINE PRECISION
(3/3 targets) — collar calculus complete. 0.01 s vs ~9 min for 3D DFS.
v2 knobs: cost-optimal priority (v1 is Z-first: +20.6/+40.3 vs DFS
+16.4), insertions, lift-time repair. COST-OPTIMAL MODE DONE (optimize='cost': all targets +20.601 uniform, 0.4 s, exact match; gap to DFS +16.40 = missing insertions in v1 alphabet). **FORCE-MODE HARVEST DONE (2026-07-30, scratchpad/harvest_f0.py):**
14/14 planner-built removals on quench_down5q_wOFF, ZERO lift failures,
planned==executed dS exact every round, 20 s total. Gap +10.12 -> +0.49,
f0 1536 -> 1522 (exactly the predicted count). Pricing SELF-REGULATES:
net/removal decays -22.8 -> -0.4 as the gap closes (quadratic pin) --
exactly the acceptance profile an MH channel needs; harvest self-stops
at gap 0. THERMO TWIST (harvest_verdict.py): gap-closed sector sits at
late<S>=351+-19 vs quenched 261+-12; n_ill 63 -> 187 (~9 illegal
edges/vacancy debris; 14 vacancies at fixed volume pin). Instantaneous
surgery downhill but thermal re-dressing raises <S>; <S> does NOT decide
equilibrium (free energy does -- vacancy placement entropy is large).
Snapshot: harvest_f0_final.mfd. Pending: 3000-sweep anneal check
(n_ill debris vs equilibrium density for off-crystal <e>).
**MH CHANNEL v1 = CBMC: MEASURED NEGATIVE (2026-07-30, commit 63a6a14,
scripts/defect_dynamics/f0_channel.py):** stochastic reshaper
(flips/deletes/attaches/STOP), ONE shared rule
w=exp(-b*dS-eta*dZ-gamma*dPhi), per-path balance via retraced reverses
-- machinery VALIDATED (exact round-trips; bare-insertion alpha matches
analytic seed entropy +9.58). But shaping telescopes => acceptance pays
exp(2*(eta*dZ+gamma*dPhi)) at the boundary: completion-strength shaping
(2.5/3.0) costs e^-70 in alpha; slack-compatible shaping (<=+25 budget)
completes 0/15 walks even low-Z-seeded. completion*acceptance ~ 0
everywhere: stepwise CBMC CANNOT drive this composite (the +20 collapse
barrier + ~50-way branching beats any potential drift the acceptance
can afford). **SCHEME C BUILT (2026-07-31, scripts/defect_dynamics/f0_worm.py,
design doc 3.2):** extended-ensemble worm. Open sector (T, head vertex
v), weight zeta*e^{-S+U}; U = eta(Z0-Z)+gamma(PHI0-Phi) calibrated at
startup (C15 mean Z=13.4! Z12/Z16 mix); five moves with exact local
ratios (open/close-flag ~ zeta e^U f0, open-insert/close-41 ~ 1<->4 +
f3*pool factors, reshape = uniform over union-of-stars region
candidates with |A|/|A'| Hastings — membership move-symmetric by the
doc-3.2 lemma). Sampler-level targeted 1<->4 (capi commit 0991770,
label-pool bookkeeping) makes it one-sampler-in-place. Design
decisions: vertex head (tet = its closed shadow, transmutation at
sector crossing); region alphabet (star-only can't repair boundary
debris — DFS cage +16.4 vs star +20.6); pair family (2-head, V
association, created at coincidence, NO half-close, P(s) = correlator)
designed but not built; label-epoch semantics = doc 2.6. **v1 RUN AUDIT (2026-07-31, f0worm_quenched 400cyc DONE /
f0worm_harvested KILLED at 28cyc as INVALID):** balance formulas
verified consistent (round-trip episodes dS exactly 0; bail-outs dS 0;
quenched S stays at baseline 262-266 over 285 bailed episodes = no
debris heating; episode ledger exact). TWO DEFECTS FOUND: (1) reshape
acceptance ~0.1% (29/28752) — uniform-over-region proposals waste on
e^-10 candidates; heads FROZEN => removal walks never descend
(c4 0/1543), fresh insert stars can't anneal, and heads that accept
1-2 reshapes off Z=4 get STUCK (can't c4, can't afford flag e^-U).
(2) forced closure at MAXSTEP=3000 commits stuck-episode debris as
closed — 26/71 harvested episodes forced (+15-20 dS each, all
insert-episodes) => the observed f0 rise 1522->1527 was ARTIFACT
(forced insertions +26, partially cleaned by ordinary-sweep 4->1 of
the abandoned Z=4 vertices — the volume channel fires for the first
time ever on worm debris!). Quenched chain = valid null (f0 static,
no heating, kernel stable but non-mixing in f0). **v2 FIXES APPLIED (2026-07-31):** two-class reshape proposals (H =
support contains v, class split move-invariant, per-class count
Hastings, PH env); capped episodes EXACTLY UNDONE (full unwind,
recorded dS must be 0 — verified); episode tuning zeta=0.05 BCF=0.01
(episodes ~150+ steps, flag-exit alpha 0.65 at U~0). UMBRELLA STATUS:
convex spoke basis G=sum(6-d)^2 still insufficient — probe shows best
H-moves dS +3.4..6.2 vs dU +0.8..1.8 (re acc 0.5%). STAIRCASE
EXTRACTED (orbit-exact, both seeds identical): flips +9/+10.4, deletes
-4.8/-4.8/+9.2/-9.0/+8.1/-7.1, peak +29.4 at [4,4,4,4,4,5,5] star,
end +20.6 at [3,3,3,3]. Per-spoke linear fit residuals +-7 (delete
costs context-dependent via side-spoke crossings) — star-descriptor
potentials CANNOT flatten this. v3 options: planner-exact cost-to-go U
(deterministic w/ fixed nodecap; S-U exactly flat on corridors; needs
speed check) or frozen Wang-Landau on visited-star table. **v3 UMBRELLA CAMPAIGN VERDICT (2026-07-31): PYTHON KINETIC CEILING.**
Iterations measured (reshape acc): linear basis 0.1% -> convex 0.5% ->
corridor multiset table (planner-harvested, min-agg; mean-agg inflates
fresh-star credit -> spurious inserts) 1.9% -> +patience(8000 steps) no
closures -> FULL WL (90 cyc, delta 1.0->0.016, |UTAB| 4312 multisets,
frozen+saved .utab.json, UTAB_LOAD reload support) + 250 production
cycles: ZERO removals (218 of-episodes ALL capped+exactly undone),
15 insert round-trips (oi->c4 identity, balance exact), f0 static.
ROOT CAUSE: corridor traversal needs ~100k steps (corridor-move
proposal rate ~4% of H-draws x context-noise acc ~5% -> ~0.1%/step
progress; 10-deep corridor diffusion ~ 100k) vs MAXSTEP 3000; multiset
key too coarse (+-5-7 collar-context residue) and richer keys explode
WL fill time; planner-guided proposals dead (100ms/call, unsolvable at
small nodecap from Z>=10). Machinery fully VALIDATED (balance ratios,
exact undos dS=0, sector crossings, table freeze/reload). **D PORT DONE + FIRST BALANCED f0 TRANSITION (2026-07-31):**
sampler.d wormF0Episode (~550 lines): 33k steps/s (300x Python), exact
undo verified (dS=0), objective tracked exactly. Key design: H-kernel
(head-in-support moves, exact |H|/|H'| Hastings, O(star) enum, cache
snapshot/restore on reject) + GLOBAL repair kernel (chooseRandomMove,
head/1<->4 auto-reject; dU==0 identically off-head since spoke degs
change only when head in support) + sector crossings with label pool.
Debug journey (all diagnosed via zmin/nZ4/umax probes): WL table =
corridor-entrance damage; linear fallback = monster-star runaway (fixed
w/ caps, then confining coeffs); fitted fallback n3=+9.15 = deg-3
farming; min-agg multi-orbit table = context-blind TRAPS (one move to
foreign +30.6 entry). WINNING DESIGN = SINGLE-ORBIT TUBE: replay ONE
cost-optimal corridor (11 multisets, context-exact by orbit
quantization) + flat off-tube OFFPEN=-3 + seed bias e^{-mu Z} (Z=2+D/2
from tet-degree, O(1)) mu=3, zeta=50, bcf=1e-4, bc4=0.1, caps 31/-15.
RESULT: head 1401 removal COMMITTED through full MH (8032 steps, umax
+29.4 = corridor peak, f0 1536->1535, S -7.7) + insert round-trips
climbing the same corridor. KEY THEOREM: closed-sector measure is
U-independent => tube may be REBUILT between episodes (each episode
balanced under its own frozen U) => per-commit retube = exact adaptive
umbrella (f-drift solved). RUNNING: two-sided production
(f0worm_prodq seed505 / f0worm_prodh seed606, 500cyc x 6ep,
maxstep 1e6, per-commit retube) = the f0 free-energy measurement.
UNCOMMITTED: D episode engine + capi + bindings + driver. NEXT after
runs: commit; brackets/label epochs (doc 2.6); pair family.

**Phase-1 reversibility lab (2026-07-30,
scratchpad/phase1_reversibility_lab.py):** reverse-check pass rate =
6/6 self-mirror sector (deg-3 edges, dS antisymmetric ±5.93), 0/24
deeper (naive mirrors walk to different attractors). VERDICT:
deterministic+reverse-check = fast path for self-mirror sectors only;
general deep sector needs CBMC/Rosenbluth stochastic constructors
(retraced reverse paths; alphabet symmetry ⇒ nonzero reverse support).
Bonuses: general engine (dressed_generators.py) cut vertex-removal to
net −22.82 (barrier +16.40, was −11.62/+23.40); costs are
CRYSTAL-ORBIT-QUANTIZED (same-orbit vertices → bit-identical
trajectories) — weight tables tiny in ordered states.

**f-adaptive tube (2026-07-31, BUILT + VALIDATED):** tube values
decompose EXACTLY as cum_k = M_k (local m^2 part, f-independent,
orbit-quantized) + [g(f0+d_k) - g(f0)] with g(f1,f3) = 0.1(f3-F3T)^2 +
(f1 - 6 f3/e*)^2 the only f-NONLINEAR action terms, and d_k = the
corridor state's exact integer (df1, df3) offset ('23' flip +1,+1;
delete -1,-1). D stores the skeleton (cum, df1, df3, build-f) and
RECOMPILES the frozen scalar table at each episode open from the live
f-vector (sampler.wf0Compile; U frozen within episode => dU==0 lemma
for global repair moves survives; per-episode recompile exactly
unbiased). One tube now serves the whole f0 range; retube demoted to
every RETUBE_EVERY=10 commits (m^2 refresh only), UCAP_HI default
raised 31->45 (compiled peak grows with gap). VALIDATION: (1) adaptive
table at build-f => bit-identical same-seed episode trajectories vs
scalar table; (2) same corridor re-measured after 15 far 2->3 flips
(f +15,+15): predicted staircase matches to 1.4e-13 (pin swings up to
8 units, e.g. [3,3,3,3] d=(-2,-2) 20.60->12.77). Files: sampler.d
(WF0Skel, wf0PinPart, wf0Compile), ddg_capi.d (config +df1/df3/tube_f
args, per-episode compile), _sampler.py set_worm_f0(f0_ref=,
triple values), f0_worm.py build_orbit_tube returns (tab, f_ref).
PRODUCTION (old per-commit-retube code, still valid, running since
2026-07-31 08:55): prodq f0 1536->1533 gap +10.3->+8.2 (162 commits,
walking DOWN toward gap-closed); prodh HOLDING f0 1522 gap +0.14 (79
commits) -- early two-sided evidence the defect-financed gap-closed
state is the free-energy minimum.

**AUTO-ZETA (2026-07-31, the kinetic blocker, FOUND + FIXED):** the
harvested chain committed ZERO f0 moves in 582 episodes (207/207 insert
episodes reached Z=4, sat there, and abandoned at maxstep) -- NOT
physics, a zeta miscalibration. Exact (sampler.d:5211):
  log a_c4 = -dS41 - U - log z + log(aoi/bc4) - log f3 - log lmax
  log a_oi = mirror image (detailed balance is exact either way)
The -log(f3*lmax) ~ -17.4 proposal-density penalty must be paid by the
insert cost dS14, which is PIN-DOMINATED and so state-dependent:
measured +39.8 at gap +10 (quenched) vs +18.9 at gap 0 (harvested).
=> log a_c4 = -0.4 quenched (a=0.6, matches observed 77% round-trip
rate) but -13.6 harvested (a=1e-6, dead). A fixed zeta prices the
sector over a NARROW GAP BAND only; outside it one leg saturates at
a~1 and its inverse dies, so the walk is pumped into the open sector
and undone at maxstep. WORSE: log z* falls steeply as the gap closes
(+3.4 at gap +10 -> -11.4 at gap +0.5), so a quenched chain WALKS INTO
the trap -- prodq2's plateau at f0 1532-1533 (gap +7.7, z* ~ 0.03 vs
running 50) was ~1000x suppression, not thermodynamics. FIX (driver
only, zeta is a free ensemble param -- closed measure is
zeta-independent, same theorem as retubing):
  log z* = dS14 - U([3,3,3,3]) - log f3 - log lmax - log(bc4/aoi)
with dS14 = mean over 24 uniformly-random-tet 1->4 insert/undo probes
(matches the engine's own randomFacetOfDim(3) proposal). f0_worm.py:
calib_zeta() + pin_part() + tube_u(); recalibrated EVERY commit (cheap)
vs retube every RETUBE_EVERY=10 (expensive); AUTOZETA=0 restores fixed
ZETA_D. RESIDUAL: zeta fixes the MEAN only; log a_c4 = dS14_site -
dS14_mean, so cheap sites insert easily and resist closing (site spread
+12.6..+37.8) -- escape is via cf after corridor growth, watch it.
RESULT: harvested side went 0 commits/582 eps -> f0 1522->1524 in
cycle 0; quenched 4 removals/300 cyc -> 1 removal in cycle 0.
RUNNING (2026-07-31, auto-zeta, seeds 909/1010, 500cyc x 6ep,
maxstep 1e6): f0worm_prodq3 (quenched, z 204->7.9) and f0worm_prodh3
(harvested, z 1.8e-6->1.1e-3). NOTE harvested f3=8699 is BELOW the
volume target 8704, so its insertions move f3 toward target -- the
gap-closed state may genuinely want to grow. Prior logs kept as
.v1/.v2.chan.jsonl (all exactly balanced, usable as extra chains).

**TWO-SIDED f0 RESULT (2026-07-31, auto-zeta chains, PROVISIONAL):**
BOTH CHAINS COMPLETE (500 cyc each). quenched prodq3 1536 -> 1533
(gap +10.1 -> +8.1); harvested prodh3 1522 -> 1526 (gap +0.5 -> +2.8).
Both sides moved TOWARD each other and
stalled => BRACKET f0_eq in [1526, 1533], gap in [+2.8, +8.1].
=> NEITHER endpoint hypothesis survives: not gap-open f0=1536/n_ill 63,
not fully defect-financed f0=1522/n_ill 189 -- the ensemble wants
PARTIAL financing, somewhere in between. NOT converged: separation 7.0
vertices, no overlap (78 pooled sd), so this is a bracket, not a
measurement. CAUTION on the quenched side: late-window slope is
0.00 +- 0.00, which looks like a clean stationarity pass but is
STUCK-NESS -- 491 insert round-trips vs only 3 removals in the WHOLE
500-cycle run (all three early, none after) and 236 abandoned episodes. Its 1533 is an upper bound,
not equilibrium; the corridor barrier still gates the removal leg.
Harvested side is genuinely mixing (6 removals AND 12 insertions, slope
-0.06 +- 0.12) and its 1526 is the more trustworthy end.
**zeta instability found (fix pending):** log10 zeta sd 1.4 (quenched) /
2.2 (harvested) => alpha swings e^+-3 around 1, alternating between
insert-favouring and remove-favouring regimes. NOT the OFFPEN fallback
(U4 == -3 never fired, 0/919). Cause: U([3,3,3,3]) varies up to 10 units
between REBUILT corridors (different seed vertex/orbit each retube), and
log zeta* = dS14 - U4 - const inherits it. FIX: anchor every rebuilt
tube to a canonical U([3,3,3,3]) -- a pure gauge choice, since a uniform
shift of U trades exactly against zeta (anchor_utab already exists in
f0_worm.py). Then zeta* tracks only dS14 (sd 1.7).

Related: [[m2-only-gas]], [[deg4-worm-design]] (bilocal v2 doc
notes/bilocal-worm-design.md).
