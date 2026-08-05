---
name: ecmc-blob-ab
description: "RESULT: blob dispersal IS a transport-limited state where nonlocal transport beats plain dynamics 16-43x, but the LIFTED ECMC flight gets only 3.7-6.8x -- the unlifted diffusive slide beats it 6-9x at matched transport budget"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0f4e228e-f94b-48d8-abdc-8231608d3f78
  modified: 2026-08-05T13:31:42.469Z
---

## VERDICT (2026-08-05)

Question: is there an out-of-equilibrium state the lifted ECMC flight
channel ([[lifted-ecmc-transport]], ecmc_flight.Flight) relaxes better than
basic dynamics? **YES for "nonlocal transport", NO for "lifting".**

C15 m3 N3672, k=12 fliers packed in a radius-4 ball, 6 chains/arm,
20000 sweeps, S = 30(f3-f3*)^2 + 1.0(f1-6f3/e*)^2 + 0.5*sum m^2.
Paired on identical initial blobs (same seed -> same blob across arms).
Speedup = first-passage WORK vs arm A, paired log-ratio.

    arm  transport channel                       spread thresholds 3.3/3.5/3.7/3.8
    A    plain (local Pachner + 4-4 hinge)        1x  (2.9M / 8.3M / 11.4M / 14.9M work)
    B    lifted ECMC flight, ~8 steps/sw          3.7x  4.8x  5.7x  6.8x   (z +3.0..+8.6)
    M    diffusive nonlocal slide, ~8/sw MATCHED  16.2x 42.9x 36.8x 41.2x  (z +4.3..+12.9)
    N    diffusive nonlocal slide, ~184/sw        17.0x 42.9x 64.5x 46.6x  (z +4.0..+9.6)

**M is the decisive control**: same transport-attempt budget as B (~8
steps/sweep), same f-preserving composite move, no lifting -- and it beats
the flight by 6-9x. Not a budget confound (M ~= N despite 23x less
transport), not a work-accounting artifact (B's episodes are only 4.4% of
its total work, so even zeroing the flight's overhead moves 6.8x -> ~7.1x).
Figure: data/ecmc_ab/ecmc_ab_four_arm.png

## Why the lift loses -- HYPOTHESES, not measured

1. **Selection, not speed.** An episode drives ONE chord for 40 steps; the
   diffusive channel picks a fresh uniform chord every attempt. Dispersal
   needs ALL 12 fliers moved, so uniform coverage beats depth on one.
   Fundamental tension: the lift NEEDS persistence on one object to be
   ballistic, but this relaxation is limited by moving MANY objects.
2. **Rung confinement.** Q is conserved for a whole episode (ecmc_flight
   docstring: "no Q-refresh yet, so one run samples one rung sector"); free
   hops reach only same-rung sites. The diffusive slide picks any slot/step
   and pays dS, so it changes rung freely.
3. Per-hop displacement of the flight (median k) was NOT recorded -- add it
   before trusting any reach-based explanation.

## Building the state: the two corrections that mattered

**Barrier = nfc MINUS the m^2 relief, not nfc.** Pilots at nfc 0.1 and 5.0
had ZERO time resolution (every threshold crossed in the first 20-sweep
chunk; B/A = 1.047 = exactly the neutral overhead 1/(1-0.044)). At nfc=5 the
f-vector is already exactly pinned (bistellar accepts balance 29<->29) yet
kill/rebirth teleport ran wide open, because killing a CROWDED flier
relieves ~10-25 units of c_imp*m^2 and pays the pin. nfc>=30 closes it
COMPLETELY (zero 3->2 accepts ever; nfc 30/60/120 bit-identical). Then only
the f-preserving 4-4 hinge move remains for arm A -- an honest local
competitor, since both the flight composite and the nonlocal slide are
exactly f-preserving (net df3 = 0) and so pay ZERO pin cost.

**spread_mean is the order parameter, not m^2.** sum m^2 saturates in ~10
sweeps (375 -> ~178, local relief is cheap), then dips and RISES again at
late work in every arm (the equilibrium gas populating) -- non-monotone, so
useless as progress. spread_mean = mean graph distance of defect vertices
from the FIXED blob center; packed = 3.13, fully dispersed = 3.98 (= mean
graph distance from a fixed vertex on C15 m3, V=648, eccentricity 7).

## FALSE STATIONARITY -- a warning for the certification standard

Arm A passes the late-window slope test on ALL FOUR observables (spread,
m2, n_ill, n3; |slope| < 2 sigma, chain-level bootstrap) while sitting at
spread 3.57 -- when the other arms show equilibrium is ~4.0. The
late-window-slope test cannot distinguish EQUILIBRATED from KINETICALLY
ARRESTED. A stationarity pass on a single dynamics is not evidence of
equilibrium; it needs a faster channel as a witness. Cf the glassy
bistability in [[m2-only-gas]] and [[edq-only-melting]].

## Balance status: NOT sharply tested

At late work all four arms overlap in a common band on both observables
(spread 3.5-4.2, m2 120-180), consistent with a shared stationary law --
but arm A never equilibrates, so this is not a sharp bias test. The
flight's audits (exact reverse proposal after every accepted gated move,
dS_sum == tracked objective, no AUDIT_FAIL) establish local reversibility
only. A real test needs long runs from an ALREADY-dispersed state comparing
stationary distributions with and without the channel.

## Code

- scripts/defect_dynamics/ecmc_ab.py -- driver (arms A/B/N, blob builder,
  self-consistent pins, incremental .ab.jsonl sidecar, work accounting:
  B's work = d_tried + 2*nonlocal_slide_at calls since totalTried never
  sees the kernel -- the [[percolation-ab-test]] lesson).
- scripts/defect_dynamics/ecmc_ab_analyze.py -- paired log-ratio
  first-passage + late-window stationarity table.
- scripts/defect_dynamics/ecmc_ab_figure.py -- the four-arm figure.
- ecmc_flight.Flight gained `expect_free=` (the free => dS=0 assert is
  EDQ-specific; under pins+m^2 same-rung hops were still dS=0 on C15, the
  free_not_free counter never fired).
- GOTCHA fixed mid-analysis: a CENTERED smoother on the first-passage
  series leaks later values backwards and reports crossings before they
  happen, biasing the FASTEST arm most (N's 3.3-threshold speedup fell
  47x -> 17x on the fix). Use a trailing window; skip index 0 (t=0, zero
  work, shared by all arms, and log(0) poisons the ratio).

Data: data/ecmc_ab/{main_nfc30,matchedN}/ (gitignored).
