---
name: percolation-ab-test
description: "RESULT: the chord channel does NOT accelerate percolation relaxation -- 1.6x WORSE on total work, matching pure overhead; plain-work effect consistent with zero"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-08-02T14:22:27.127Z
---

## VERDICT (2026-08-01): the chord channel does not accelerate this relaxation

A/B at cimp=0.25, C15 m4 N8704, 6 chains/arm, rmax=0 (certified config),
CHUNK=10 sweeps, NEP=2400 episodes/chunk.

    observable                    armA        armB    B/A       z
    fp P>=0.25  TOTAL work   3,186,494   5,244,785  1.646   +6.58
    fp P>=0.50  TOTAL work   4,270,527   6,858,760  1.606   +5.63
    fp P>=0.25  PLAIN only   3,186,494   2,940,196  0.923   -1.06
    fp P>=0.50  PLAIN only   4,270,527   3,977,927  0.931   -0.85

Episode work-fraction 40.3%. If episodes contributed EXACTLY NOTHING the
total-work ratio would be 1/(1-0.403) = 1.674. Measured 1.646 / 1.606 --
within 2-4% of pure overhead. The plain-work ratio isolates the same
thing: 0.92-0.93 at |z| < 1.1, NOT significant, 95% CI 0.79 to 1.10.

Even at the optimistic end the episodes buy <= 27% of plain work while
costing 67% in overhead: the channel cannot pay for itself at this work
fraction wherever in that interval the truth lies.

NOT a "channel never fired" artifact: 6132 commits/chain, 4.48% commit
rate. It worked as designed; the design does not help here.

Internal consistency check PASSED (measured ratio ~= neutral prediction),
which validates the work accounting -- the thing most likely to have
produced a false positive.

## Arm C (rmax=2, beta=2) 2026-08-02: SAME VERDICT

    arm   ep_frac  attempts/ep  commit rate
    B     40.3%       24.5        4.48%
    C     12.2%       12.9        1.79%

    arm observable        armA         arm   ratio      z  neutral
    C  fp>=0.25 TOTAL  3,186,494   3,744,407  1.175   1.97   1.139
    C  fp>=0.25 PLAIN  3,186,494   3,228,525  1.013   0.17   1.000
    C  fp>=0.50 PLAIN  4,270,527   3,840,332  0.899  -1.86   1.000

95% CIs on C's plain-work effect: 0.87-1.18 (fp25), 0.80-1.01 (fp50).
Both contain 1.0. The fp50 z=-1.86 is NOT claimable: below threshold,
one of eight comparisons, and fp25 on the SAME chains is 1.013 / z=0.17.
Two thresholds disagreeing on one dataset = noise.

C's smaller total-work penalty is not the channel doing better, it is
the channel doing LESS, with the penalty tracking its reduced overhead.

## KEY TRANSFER FAILURE: the knob sweep did not generalise

On the QUENCHED state (n3~17) rmax=2/beta=2 gave the MOST commits
(182 vs 62 baseline, ~3x throughput). At PERCOLATION densities it gives
LESS THAN HALF arm B's commit rate (1.79% vs 4.48%) and its episodes
decline so early they make only 12.9 attempts each. Throughput measured
on a dilute quench does NOT predict behaviour on a melting substrate.
This was the main reason to expect arm C to win.

## Episode cost scales with DEFECT DENSITY, not just N

wf0ChainSites phase 2 walks one BC chain per start window and the number
of start windows scales with n3 (site count at k=20: 628 at n3=11 ->
2194 at n3=63). So low-density smokes badly underestimate both runtime
and NEP calibration -- this cost a 3x runtime miss (25 min/chain
predicted, 76 actual) AND the ep_frac miss above. ALWAYS calibrate at
the density the run will actually reach.

DRIVER FLAW: perc_ab.py writes its JSON only after ALL seeds finish, so
killing a long run discards completed chains (the .rec.jsonl files are
flushed live but lack the illegal-edge decomposition needed for P).
Write per-seed results incrementally.

Drivers: scratchpad/perc_ab.py; results ab25_{thermal,chord}.json,
ab25r2_chord22.json.

---

The test design that survived scrutiny for proving the chord channel
does something, after two earlier framings failed.

## Why the earlier framings failed

**Equilibrium population A/B is structurally void.** Measuring n_deg4
with/without the channel at equilibrium tests DETAILED BALANCE, not
usefulness: a correct channel MUST leave the stationary population
unchanged. Measured 2026-08-01: pristine C15 m4 at cimp=0.7, e*=5.1
(pristine objective EXACTLY 0.000000, so the crystal is the true ground
state) is flat over 87.1M attempted moves — n3 4.45+-1.60, n4
14.83+-4.44, n_ill 19.29+-5.12, all late-window slopes < 2 sigma
(block bootstrap). No accumulation whatsoever.

**Percolation from pristine IS valid** because it is a RELAXATION: both
arms share the same endpoint, so a balanced channel is free to change
the RATE. That is the loophole the population test lacked.

## Calibration (2026-08-01, C15 m4 N8704, m2-only, edge pin ON, e*=5.1)

Order parameter P = S_max / n_ill on the illegal-edge graph;
chi = sum' s^2 n_s / sum' s n_s (primed = excluding largest) peaks AT
the transition and is threshold-free.

    cimp   fp P>=.25   fp P>=.50   chi_pk   chunks to percolate
    0.35     43.8M      NEVER       55.75   240 (capped, P->0.249)
    0.30     13.7M      17.9M       83.05    95   <-- USE THIS
    0.25      3.70M      5.23M      35.63    28   (fallback, 3x faster)
    0.20      1.09M      1.31M      15.98     8   (no resolution)
    0.15      0.44M      0.65M      27.18     5   (no resolution)

cimp=0.30: ~63 s/chain, sharpest chi peak, and the chi peak (18.07M)
agrees with the P>=0.50 crossing (17.85M) to 1%. Climb is gradual with
a resolved chi peak => assembly, not a single nucleation event (but ONE
SEED only; chain-to-chain scatter is still unmeasured and sets the
required chain count).

## Arm A baseline + power (6 chains, cimp=0.30, 2026-08-01)

    observable        mean (work)         sd    CV   sd(log)
    fp P>=0.25         12,810,346  1,273,020  0.10     0.100
    fp P>=0.50         15,857,258  1,614,443  0.10     0.102
    chi peak           18,758,929  3,104,039  0.17     0.157

CV = 0.10, NOT ~1 => first passage is NOT exponential => gradual
assembly, not single-event nucleation. The first-order hazard does not
bite at cimp=0.30. 89-120 chunks/chain, ~63 s/chain.

Power for a two-sample log-mean test, n = 2*(2*sd_log/ln(ratio))^2:
10 chains/arm resolves 1.1x at 2 sigma; 3 chains resolves 1.2x.
~11 min per arm. Driver: scratchpad/perc_ab.py.

## BLOCKER: the chord episode is O(N) and ~1200x too slow

Measured 2026-08-01. Per ATTEMPTED move: plain sampler 2.8 us; chord
episode 3294 us (rmax=0) / 3716 us (rmax=2,beta=2) => 1171x / 1321x.
At a 50% episode work-fraction that is ~720-810 s/chunk = ~20 h/chain.
Infeasible.

N-scaling (C15 m3 N3672 vs m4 N8704): 47.1 -> 115.1 ms/episode, ratio
2.44 vs N ratio 2.37 => LINEAR in N, attempts/episode ~constant. Fit:
cost ~ 0.0135 ms * N, negligible constant -- essentially ALL the cost is
O(N) sweeps.

CULPRIT (corrected 2026-08-01 -- my first diagnosis was WRONG). It is
NOT the wf0ChainSites facet sweep: measured directly via
ddg_sampler_chain_sites, that costs only ~1.5 ms/call (debug build),
~7.5 ms of the 115 ms episode.

The real hot spot is **`Manifold.star()`**, manifold.d:879:

    auto star(S)(S simplex) const {
        return facets.filter!(f => simplex.isSubsetOf(f));
    }

`facets()` (manifold.d:845) allocates a dup'd slice for every one of the
N facets AND SORTS the result; filter cannot save it because facets() is
eager. So EVERY star() call is O(N log N) with N allocations -- measured
~1.4 ms at N=8704.

Worse, FIVE call sites in sampler.d call `mfd.star(v.only)` purely to
grab ONE seed facet and then `break`: 4839, 5029, 5533, 6653
(wf0RegionClean, fires at open and every close attempt), 6850 (in the
walk). Each pays full O(N log N) for one tet. wf0ChainSites' own comment
already warned "mfd.star() is avoided entirely -- it allocates".

REAL FIX: an O(1) "some facet containing vertex v" accessor on Manifold
(per-vertex witness facet, maintained through facet insert/remove, same
pattern as Deg3Set's witness). Then all five sites become O(1).

FAILED ATTEMPT (do not repeat): rebuilding wf0ChainSites phase 1 from
the degree-3 edges outward via wf0EdgeLink. Counts came out IDENTICAL
across all states and all kmax (so the set logic is right), but
wf0EdgeLink calls star(), turning one O(N) sweep into n3 x O(N log N):
1524 -> 10571 us/call at n3=11, 1696 -> 88326 us at n3=63. Reverted.
Any "build outward from the edges" approach needs an O(1) witness, which
is exactly what Deg3Set stores and what a bare edge does not have.

Validation harness that DID work: scratchpad/chainsites_validate.py --
melt at a fixed seed (plain sampler only, so states are bit-identical
across builds), call ddg_sampler_chain_sites for several kmax, diff the
counts between old and new builds. The COUNT is the balance-relevant
quantity since `pick` is uniform over [0, n).

Same bug class as the planner's O(|edeg|) CollarState.__init__ and its
link_edges()-inside-the-loop: global scans where a local walk would do.

## THE control that is not automatic

`totalTried` is incremented in exactly two places in ddg_capi.d, BOTH
inside the main run() loop. **The worm episodes never touch it.** So
Recorder's `d_tried` counts arm A's moves only and arm B gets free work.
True arm-B work = d_tried + sum over episodes of (nH + nG), both
exposed in the `worm_chord_strict_episode()` dict. Without this any
"speedup" is an artifact. See [[strict-chord-channel]].

## Calibration gotchas found the hard way

- A lone flicker is NOT one component. Its deg-3 apex chord shares no
  vertex with its three deg-4 face edges, so one 2->3 on pristine gives
  components {1, 3, 1}. Compound threshold is size >= 4, not > 5.
- One 2->3 on pristine makes 3 deg-4 edges (spectrum {3:1, 4:3, 7:1}),
  so RAW n_deg4 measures flicker density, not species formation.
- lam is the EDQ coupling (hinge_degree_target_coef = lam*e*/6). The
  m2-only gas has lam = 0 EXACTLY. reaction_census.py defaults are
  --zleg 0.6 --cimp 1.0 BOTH scaled by lam, num_hinges_coef=0; so
  lam=0.35 means zleg=0.21, cimp=0.35. Setting zleg=cimp=0 on a lam row
  deletes the n6 potential and percolates instantly — that is a
  mis-specification, not a property of lam=0.35.
- The census "~22 complexes at lam=0.35" was on R m4 N57984, 6.6x the
  C15 m4 box. Complex counts do not transfer between cells.

## Channel-specific caveat

The chord channel harvests deg-3 edges, so on pristine crystal it is
IDLE until thermal flickers appear. It cannot accelerate nucleation of
the first flicker, only subsequent assembly => the effect should be ~0
early and grow with density. A "speedup" in the earliest window is a
BUG signal, not a result. Also the strict closure needs the vacated
region clean, which gets harder as density rises — run the A/B at both
regionMax=0 and regionMax=2 (see [[strict-chord-channel]]).

Scripts: scratchpad/perc_pilot.py (density ladder),
deg4_observable.py (flicker calibration), deg4_density_scan.py.
