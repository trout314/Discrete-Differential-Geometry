---
name: edq-only-melting
description: EDQ-only R crystal has NO dilute defect phase (first-order melt); the m^2 CLUSTERING penalty does NOT cure it at lam=0.35 -- bigger m^2 WIDENS the coexistence (empties crystal side, liquid branch untouched)
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-28T07:58:11.229Z
---

**CORRECTED (2026-07-28): the m^2 clustering penalty does NOT give a stable
dilute equilibrium at lam=0.35 -- the "stable dilute ~75" was METASTABLE
supercooling, exposed by the formal over-dispersed gate.** Action = EDQ +
volume pin + V(m)=CIMP*lam*m^2 per vertex (m = #incident illegal edges), NO n6 U:
`s.set_n6_potential(0.0, CIMP*lam, tilt=[0]*5)`.

FORMAL GATE (8 chains, 4 pristine + 4 over-defected, quantized_split_rhat +
complex census as observables, scripts/certify_edq_clust.py): at CIMP=1 the two
directions provably DON'T meet -- pooled n_ill climbs 142->178 monotone, qRhat
5-7, never certifies. First-order coexistence, not equilibrium.

m^2-COEF LADDER at lam=0.35, PER-CHAIN long runs to 6000 sweeps (m2_moderate.py,
1 pristine + 1 over-defected chain per CIMP). The m^2 coef (imp=CIMP*lam) has a
sharp qualitative structure -- CIMP=1 was in the MELTING regime, not metastable:
  CIMP=1.00 (imp0.35): prist MELTS (48->73->136 @6000, still climbing, no
            plateau); over runs up to ~300. Crystal has NO durable dilute state.
  CIMP=1.25 (imp0.44): SWEET SPOT -- BOTH branches stable+flat over 6000 sweeps:
            prist ~58, over ~120 (liquid stops running away). Bistable
            coexistence, nonzero dilute gas WITH complexes (max~5-9, ~75 comps).
  CIMP=1.50 (imp0.53): prist stable ~15 (small, few complexes).
  CIMP=1.75 (imp0.61): prist stable ~13.
  CIMP=2.00 (imp0.70): prist ~3;  CIMP>=3: prist EXACTLY 0 (frozen, overshoot).
Mechanism: m^2 = NUCLEATION BARRIER. Too small (CIMP=1) -> crystal slowly melts +
liquid runs away. Moderate (1.25-1.5) -> STABILIZES both branches (crystal stops
melting, liquid capped). Too big (>=2) -> crystal frozen to ~0. The transition
stays FIRST-ORDER/bistable across the whole moderate window (two stationary
branches that don't merge) -> over-dispersed-from-both-directions can't certify a
SINGLE phase. Trajectory fig: scratchpad/m2_moderate_trajectories.png.
WITHIN-BASIN CERT of the CIMP=1.25 liquid FAILED (certify_liquid_cimp125.py, 8
chains from lam=0.30 melt starts spanning n_ill 98-192, gate on census): the
branch is NOT stationary. Over gates 14->46 (2500+560..1840 sweeps) pooled n_ill
DRIFTS UP monotonically 119->124, objective 230->240 (roughly linear, not
plateauing), qRhat only creeps 3.4->2.9 (never near 1.01), tau_grow stuck ~2-2.4,
volume-pin qRhat stuck ~1.4. Signature of a SLOWLY-LEAKING metastable branch --
the "flat ~120 over 6000 sweeps" from a single chain was too short; with
over-dispersed starts + the census as observable it's clearly still melting.
STRUCTURAL CONCLUSION: in EDQ+m^2 the moderate-density complex-rich state is
INTRINSICALLY a metastable coexistence branch (first-order); it does not certify
as a single phase. A CERTIFIED complex-rich ensemble needs either (a) the DENSE
liquid at a lam below the crystal spinodal (unique phase, both directions
converge -- but dense) or (b) the FULL FK action (EDQ + n6 legality + m^2), which
historically DID give certified dilute complex-rich states (lam35 snapshots
n_ill~50-155) -- the n6 term is likely what makes a sparse gas a true single
phase (see [[edq-only-melting]] note "FK term is what makes a sparse gas").
EARLIER (WRONG) READ, now corrected: called CIMP=1 "metastable coexistence" from
a 3000-sweep pooled gate; the pristine side is actually MELTING at CIMP=1.

ISOLATING n6 vs m^2 (2026-07-28, m3 R lam=0.35, convergence test pristine vs
over-defected). set_n6_potential(zleg_coef, imp_coef): U(n6)=zleg*dist^2(n6,
{0,2,3,4}) per vertex (n6=#incident deg>=6 edges; FK-legal = Z12/14/15/16) +
V(m)=imp*m^2. Campaign convention: zleg_coef=zleg_scale*lam, imp_coef=cimp_scale
*lam; production FRESH_DEFAULTS zleg_scale=cimp_scale=0.3 (at lam=1.0!).
  - n6 ALONE (zleg 0.3*lam=0.105, m^2 OFF): pristine INSTANTLY melts 0->~9500,
    ONE giant percolating cluster (~3500 verts, comps=1); over-def also ->~9500.
    SINGLE-PHASE (converges both ways) but DENSE. => n6 does NOT provide
    diluteness; the m^2 term does. Refines "FK term makes sparse gas": it's m^2.
  - FULL action lam=0.35, zleg=cimp=0.3*lam=0.105 (production scale): MELTS to
    ~3990 dense, max_size~3350 (one giant cluster), all 4 starts converge.
    Single-phase but dense. The 0.3-scale legality is 3x weaker than the
    cimp=1.0*lam=0.35 m^2 that gave the (metastable) dilute branch -> too weak
    at lam=0.35. NO dilute phase here.
KEY: at lam=0.35 NONE of {m^2-strong, n6-alone, full-action-0.3} gives a dilute
SINGLE phase. Diluteness needs stronger effective coupling. Both knobs scale with
lam (EDQ stiffness = lam*e*/6; legality = 0.3*lam), so pushing lam UP should
cross the first-order transition into a DILUTE unique phase. lam ladder
{0.50,0.65,0.80} full-action running to locate it (lambda_ladder.py). The
historical certified dilute states were lam=0.40 on m4 (different crystal) or the
lam=1.0 baseline -- NOT lam=0.35 m3. Campaign target = find dilute single-phase
point on m3 R, then run 7h defect-complex data campaign there.

The original EDQ-only finding below (no clustering penalty):

**Measured 2026-07-28 on m3 R (T3_R_m3_N24462, N3=24462), EDQ + volume pin
ONLY (no FK n6/m^2 term, no variance penalties), target e*=5.105025.**
Action = 0.1(N3-N*)^2 + lambda*(e*/6)*sum(deg-e*)^2.

KEY NEGATIVE RESULT: the EDQ-only R crystal has NO stable DILUTE defect gas.
Equilibrium is either PRISTINE (stiff, lambda above ~1.2-1.5) or a DENSE LIQUID
(~5000-11000 illegal edges), with a FIRST-ORDER melting transition between them
and STRONG metastable supercooling. The dilute defect gas seen in the FULL
action (lam35 snapshots, n_ill~50-155) is created by the FK term (the m^2
illegality penalty + n6 potential); strip FK and it's all-or-nothing.

n_ill (illegal-edge count) vs sweeps from a pristine start, count equilibrates
in ~200 sweeps at low lambda but the metastable pristine survives LONGER the
closer lambda is to the transition, then melts:
  lam 0.35: ~11200 (heavy liquid)   lam 0.50: ~8800
  lam 0.80: climbs to plateau ~5300 by ~2000 sweeps
  lam 0.85: still climbing ~5100 @4000   lam 0.90: climbing ~4360 @4000
  lam 0.95: looks dilute (~135) to ~1500 sweeps then ACCELERATES away
            (1102 @4000, still rising) -- pure metastability
  lam 1.00: 72 @1000 sweeps was METASTABLE (would melt if run longer)
  lam 1.50/2.0/3.0: 0 defects (frozen pristine) at 1000 sweeps
So the earlier "dilute at lam=1.0" reading was a supercooled-pristine artifact,
not equilibrium. Transition is somewhere lam ~1.0-1.3; NOT yet bisected.

The liquid is DEGREE-4 (disclination-line) DOMINATED at all lambda (deg4 >> deg3
>> deg>=7 low), because EDQ cost (4-e*)^2=1.2 << (3-e*)^2=4.4 so degree-4 is the
cheapest defect. NB: the non-local slide channel [[nonlocal-slide-move]] moves
DEGREE-3 chords, of which there are ~0-100 here -- so it is NOT the mixing tool
for this deg-4 population.

Timing: m3 R runs ~450-500k moves/s EDQ-only (no pot); 1 sweep = 24462 moves,
so ~1000 sweeps/min. Easy to generate samples once a regime is chosen.
NEXT (user to pick): (a) analyze the dense liquid (e.g. lam0.80 ~5300, burn
>=2500 sweeps); (b) bisect the first-order transition (long runs lam 1.0-1.3);
(c) add back JUST the m^2 term for a stable dilute gas. See
[[tcp-r-c15-defect-state]], [[statics-hu-verdict]], [[no-halo-verdict]].
