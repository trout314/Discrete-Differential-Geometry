---
name: mobile-gas-liquid
description: "The constrained knot LIQUID at lam=0.40 (m4) — mobile carriers, walker transport, measured layered dielectric screening"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-23T04:29:12.534Z
---

Mobility campaign results (2026-07-22; scripts scripts/defect_dynamics/
{mobility_sweep,mobile_gas,mgas_analyze,cross_spectrum}.py; data scratchpad/
mgas/). Coupling scale lam multiplies the run5h objective (edge pin + n6 FK
push; lam=1 = run5h frozen regime).

PHASE MAP (m4, bump +8e-4): lam>=0.5 defect-free crystal; lam=0.40 STATIONARY
CONSTRAINED KNOT LIQUID (the state the ADM/OCP program needs): nill~59,
legaledge 0.9995, Jaccard->0 by 6k sw (full turnover). lam=0.35 gas but slowly
densifying (non-stationary over 20k sw). lam=0.30 RUNAWAY into percolated web
(nill->3.9k, edeg overshoots). Spinodal between 0.30-0.35. (m2 window was
~0.35; shifts with box.)

CARRIER INVERSION at lam0.4 vs lam1: multimers (>=10) frozen-immortal at lam1
-> EPHEMERAL here (med 150 sw); sub-5 minimal units (sizes 2-4, deg-4-edge
scale) NEVER exist at lam1 -> the LONGEST-LIVED species here (med 1200 sw, max
6300). 5-knots blink (<=150 sw) in both. Birth rate 45/1000sw. Carrier identity
is lam-dependent.

WALKERS (revised 2026-07-23; original claim was membership-based only):
long-lived size-3/4 carriers migrate — visit 24-33 distinct vertices over 1-4k
sw, some ending ZERO overlap with birth site (>=2-shared-member tracking). The
"~1 vertex-spacing per ~100 sw" rate was a VISITATION rate, not displacement.
FIRST DIRECT MSD (cents-instrumented ab?b/x2 segments, 214 tracks): rms ~1.3
vertex spacings per ~1500 sw, SUBDIFFUSIVE in-window (MSD/lag falls 0.63->0.19
by 1800 sw) — order of magnitude slower than the visitation figure; caged
rattling + rare hops, walker-as-fast-transport UNCONFIRMED. Caveats both ways:
survivor/blinker bias understates the walker tail; visitation overstates
translation. Needs long cents runs (stage-1 cert extension). Transport channel
vs lam=1 frozen still qualitatively real (full turnover, migration).

SCREENING MECHANISM MEASURED (cross_spectrum.py, carrier/halo/bulk partition of
dq): S_cc = 0.0069 = point-monopole floor (carriers radiate BARE — form-factor
alternative dead). ALL cross terms negative: 2S_ch=-0.006 (dielectric cloud,
directly imaged), 2S_hb=-0.016 (halo in turn screened by bulk — LAYERED onion
carrier<-halo<-bulk). Autos 0.030 cancel to 0.0035 total (~88% destructive).
Liquid plateau/floor = 0.33 vs FROZEN 1.12 (frozen control: expect no negative
cross terms — verify xspec_frozen.log).

FINAL VERDICTS (30 snapshots, 5 chains, jackknife; commit dcfaf6f):
- S_N flat at 1.00 all shells; S_Q within 1-2sig of 1: NO Stillinger-Lovett
  onset at k>=1 — an EQUILIBRIUM statement now.
- g_cc(r) ~ 1 everywhere (|u_eff| <~0.1-0.2 kT, r>=0.5 cells): carriers are an
  IDEAL GAS — no measurable constraint-induced pair interaction.
- plateau/floor ~0.4 FLAT in k: medium = dielectric-CONSTANT response (bound
  polarization), no Debye downturn in-window. Frozen control: no screening,
  even ANTI-screening at low k (+cb). Mobility switches on the dielectric.
- TT-CHANNEL HIERARCHY (the big one): liquid inverts into the GR pattern —
  scalar 0.0023+/-0.0004 vs TT/V/L ~0.011 (5x gap; frozen state had all ~equal
  0.02, scalar LEAST suppressed). Vector NOT gapped -> model has H-like
  constraint, no momentum sector — measured independently, matches
  construction. Host = single crystal grain (98.8% registry) = superionic
  analog: solid host + liquid carrier subsystem. lam30 endpoint: registry 1.7%
  but 64% FK vertices (closest-yet FK-disordered state; hysteresis vs melt-up).
- k->0 (class-I HU) still OPEN: scalar suppression flat in k; needs modes
  below k=1 => bigger box (m6/m7).
Charge budget lam40: carriers supply ~ALL realized curvature demand
(-12.1 vs -10.3 rad; background +1.8); pin only ~18% satisfied (soft).

LAM=0.35 INSTRUMENTED (2026-07-23, chains lam35c from snap11000 + lam35r from
snap20000-with-regenerated-cocycle, 20k sw each, cents+guard clean): still
DENSIFYING through nill~230 (two seeds agree; +3-5 nill/1000sw; NOT the lam0.30
runaway scale). MSD out to 9600-sw lags shows CAGE-ESCAPE CROSSOVER: local
exponent alpha rises 0.44 (2.4-4.8k sw) -> 0.75 (4.8-9.6k sw), rms 2.8 vertex
spacings and climbing; one complex lived the full 20.1k sw. The lam0.40
"subdiffusive/caged" verdict = same curve truncated at its plateau. Species at
lam0.35 vs EDQ-only: full action makes (3,4,4) knots (91% d4 + 9% d3 edges,
q=-1.19+/-0.18); EDQ-only makes pure deg-4 open WORM segments (107/107 illegal
edges d4, no d3, unquantized sizes/charges) — n6 is the knot-binding term.

Corrected conventions (user-caught): complex monopole MUST use CELL-MEAN qR
reference (not median/Z12): knot ~ -1.05..-1.5 rad. Flat degree = 5.104299;
"5.0043" = the old POSITIVE-curvature pin. See [[defect-kinetics-run5h]],
[[curvature-length-scale]], [[cocycle-detachment-bug]].
