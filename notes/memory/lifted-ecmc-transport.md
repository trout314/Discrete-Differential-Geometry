---
name: lifted-ecmc-transport
description: "Lifted/ECMC chord transport: same-rung hops are EXACTLY free (validated vs the sampler); the hop does not traverse so there is no barrier; mean free path to contact 55-83 steps at lam=0.35, so ell must be ~50-150 steps to see collisions"
metadata:
  type: project
---

Started 2026-08-02 after [[bilocal-program-saga]] concluded that plain bilocal
moves buy nothing (factorization + the invariance theorem make a reversible
disjoint-support bilocal move exactly two independent local moves). The two
escapes are BUNDLING and LIFTING; this is the lifting thread. Library:
`python/discrete_differential_geometry/ecmc.py`.

## The framing error worth not repeating

I first built a flight sampler that integrates the uphill action along the ray,
as if the chord must TRAVERSE the chain. **It does not.** The non-local slide
is a 3->2 at the source and a 2->3 at the target, so the manifold changes only
at the endpoints and

    dS = c * (Q_target - Q_source)

depends on the endpoints ALONE. Whatever lies between is never visited and
costs nothing. There is no barrier to climb. `sample_flight` is still correct
for the LOCAL +-4 slide (which really does step through) and for physical
kinetics; it is the wrong model for the non-local move. The user spotted this.

## Same-rung hops: exactly free, validated against the D core

`scripts/defect_dynamics/samerung_validate.py`, R m2, pure edge action:

    468 legal (slot, steps) targets   max |dS - c*dQ| = 0.000e+00  EXACT
    124 same-rung targets             max |dS|        = 0.000e+00
                                      acceptance      = 1.000000000000
                                      displacement    = 1..15 chain steps
    mean acceptance, uniform proposal    0.307
    mean acceptance, same-rung proposal  1.000

That 0.307 is the entire 4-22%-on-real-melts story from
[[nonlocal-slide-move]]: the channel proposed a uniform STEP, landing off-rung
about two thirds of the time. `hop_target`/`hop_index` are inverse by
construction, so the proposal is symmetric and acceptance is exactly 1 -- no
Hastings factor, no rejection.

GOTCHA: zero the vdv coefficients explicitly (num_hinges_coef,
hinge_degree_variance_coef, codim3_degree_variance_coef) or the residual
returns and hops stop being free. See [[vdv-hdv-conflicts-with-tcp]].

## The rung landscape along a chain is SHORT-CORRELATED

    free-step fraction     0.36-0.37   (matches the recorded 37-45%)
    same-rung run length   mean 1.5-1.6, max 5-8
    stride sweep 1..12     P(same rung) 0.25-0.39 at EVERY stride

No stride has long runs; the apparent period-4 pattern in the rung sequence is
a local coincidence. So a scheme that must STEP along the chain has persistence
~1.6 and gains only a constant factor. Hops do not care -- they skip.

## Mean free path to contact (lam=0.35, R m4, ~2.2% defect vertices)

    steps to contact   mean 74-83   median 55-56   p10 7-8   p90 154-179
    rays on short closed chains with NO defect at all: ~5% (period 105-132)

**ell guidance: order 50-150 chain steps.** Set ell ~ 10 hops and you terminate
before nearly every collision and never see a handoff; the ADM-relevant
momentum transfer only happens if ell reaches the mean free path. This is the
opposite of what I guessed before measuring.

## Two structural facts

**The chain through a defect's own site does not exist.** For a degree-3 chord
(p,q) with link triangle (x,y,z), the face (x,y,z) is exactly the one the 2->3
destroyed, so `face_apexes` throws there. Any "walk the chain from the defect"
code must start from the post-annihilation state.

**Passing over another defect is ALLOWED, stopping is a modelling choice.**
Since dS depends only on the endpoints, a hop may legitimately jump OVER
another defect provided neither endpoint's support touches it. So:
  * FORCED: the target support must be vertex-disjoint from other defects,
    else factorization fails and dS is not 0.
  * OPTIONAL: whether a hop may pass over a defect. Allowing it gives the
    fastest sampler and NO collisions, hence no momentum transfer; forbidding
    it creates the collision structure the ADM program wants. The maths does
    not decide this -- the physics question does.

## Why "go as far as possible until contact" is wrong

Deterministic maximal flights collapse to a period-2 orbit between the two
blockers bounding a segment; the interior is never sampled. This is the known
non-ergodicity of the Bouncy Particle Sampler without refreshment. Hard-sphere
ECMC avoids it with a fixed total displacement ell that cuts the last flight
off mid-stride; BPS uses a Poisson refresh rate. Either works; "maximal" does
not. And ~5% of rays have no contact at all, so a cutoff is needed regardless.

## MEASUREMENT-ARTIFACT BUG (cost two wrong answers)

My first two mean-free-path answers -- "100% of rays unblocked at 400 steps",
then "at 4000 steps" -- were pure instrumentation failure. The walk died on its
FIRST step (the broken-chain fact above), the generator returned empty, and the
caller read `hit is None` as "never blocked". Both reported nulls were noise.

Caught only because a later script printed the chain period and it came out 1,
which is impossible. FIX AND LESSON: a dead walk must be an ERROR, not a value.
Same shape as every other bug this session -- a silent degradation that
produced a plausible number.

## TRANSMIT/ABSORB CENSUS (2026-08-03) -- the blocker is almost always dead

Rays walked through the MELT state's own face->apex map from crystalline seeds;
first support touching a defective vertex is the CONTACT, the cluster it
touches is the BLOCKER. Defective := imp>0 (edge degree <5 or >6) OR
n6 not in {0,2,3,4}; verified 0 defective on the pristine R m4 reference.
A blocker TRANSMITS iff it holds a deg-3 edge with >=1 free (dS=0) slide.

                     lam=0.35 (2.2-2.7% def)   lam=0.40 (0.61% def)
    transmit                14.5%                    28.1%
    absorb                  81.7%                    47.4%
    no contact               3.8%                    24.5%
    mean free path         87 steps                 214 steps

15000 / 7500 rays over 6 / 3 snapshots. Steps-to-contact at lam=0.35
reproduces the independent earlier measurement (mean 87 vs 74-83, median 59
vs 55-56, p10 9 vs 7-8) -- instrument cross-checked. NOTE lam is STIFFNESS:
higher lam = FEWER defects (lam0.40 0.61%, lam0.35 2.2%, lam0.30 26%).

Among rays that contact anything, transmit is 15% (lam.35) / 37% (lam.40), so
expected handoffs before the current dies = 0.18 / 0.59: **the chain absorbs on
its first collision either way.** Momentum persistence ~103 / ~340 chain steps.

WHY -- a population imbalance: lam=0.35 R m4 has only **3-6 degree-3 edges**
against ~230 defective vertices in ~24 clusters. Nearly every cluster is built
purely of deg-4 edges (no deg-3, no deg>=7 anywhere). The recurring
size-5/deg3=1/deg4=2 cluster is the (3,4,4) flicker knot; everything bigger is
deg-4 scaffolding. The mobile species is a TRACE IMPURITY in an immobile deg-4
matrix.

Holding a deg-3 edge is NOT sufficient: the size-48 cluster at snap15000 has
two and gives **0 free targets in 2880 trials**; at a 120-step scan still 0 in
2880 while free chords scale linearly (22%->23%). Burial is complete, not a
scan-window artifact. **Compounds are walls, not conductors** -- which kills the
"through-compound handoff" (activity threading a cluster chord-by-chord, the
way ECMC traverses dense phases): the internal relay needs each chord to have
some free hop, and buried chords have none.

Also measured: EVERY non-local slide is legal (360/360 per chord over 12 slots
x 30 steps). There is no "illegal" blocking -- blocking is entirely dS != 0.

CONSEQUENCE: the bottleneck is not the handoff rule, it is why the deg-4 matrix
is immobile. See [[deg4-worm-design]] -- the deg-3-CATALYZED worm move is the
one known deg-4 transport channel, and the colliding chord is itself a deg-3 in
range, so a nominal "absorb" event may really be a worm event.

## DEG-4 WORM AS THE ECMC RECIPIENT (2026-08-03 measurements)

The worm channel ([[deg4-worm-design]]) is fully built: oracle
`worm_deg4_slide.py`, walker `worm_walk.py`, D core `wormEnumerate`/
`tryWormMove` + mcmcStep channel, Python `set_worm_prob` / `worm_stats` /
`worm_enum` / `worm_at`. NOT cocycle-safe (gated off when a cocycle is
attached).

**EVERY worm candidate is dS = 0 EXACTLY** (472 on lam35r_snap15000, 506 on
lam40_snap14000; all |dS| < 1e-9, none downhill). This is forced, not luck:
the class preserves the illegal edge->degree MAP and is f-vector neutral, so
n_edges and Sum(deg) = 6*n_tets are both fixed; with n5+n6 and 5n5+6n6 fixed,
n5 and n6 are individually determined, so the WHOLE degree multiset is
preserved and any multiset-function action is invariant. CAVEAT: per-vertex
terms (n6 potential, m^2 clustering) are NOT multiset functions and DO break
it -- worm_walk saw dS in {-0.42, 0, +0.42} under a richer action.

Availability is gated two ways (lam35r_snap15000 / lam40_snap14000):

    interior deg-4     0/15, 0/3     0% live      <-- NEVER moves
    tip                21/78, 12/20  27% / 60%
    lone               3/28, 3/10    11% / 30%
    overall            24/121, 15/33 20% / 45%

and by graph distance to the nearest deg-3 catalyst: dist 1 -> 73%/100% live,
dist >= 8 -> 0% live. Matches the recorded "~74% of firings catalyst-gated".

**The worm is strictly a TIP move.** Interior deg-4 edges are 0% live in both
states -- consistent with the tip motif (1x deg-4 + 10x deg-5, charge 12).

Causal test (plant one 2->3 on a face touching a dead anchor at catalyst
distance >= 5, re-enumerate): opened only **4/56**. UNDER-POWERED -- it tried
12 faces per anchor with one 2->3 each and did not restrict to tips or search
catalyst placement. Treat as "a randomly placed catalyst rarely suffices",
NOT as a refutation of catalysed handoff.

## State

Built and tested (235 tests): `face_rung`, `chain_rungs`, `uphill_staircase`,
`sample_flight` (+FlightResult), `same_rung_positions`, `hop_target`,
`hop_index`, `resolve_event` + EVENT_POLICIES (reflect / handoff_sigma /
handoff_chain / rotate_slot / refresh).

Next: ell wiring, then the skew-DB two-sided test (a plain DB test FAILS by
design here), then the ballistic diagnostic -- MSD ~ t^2 below ell is the
go/no-go, since no reversible channel can produce t^2 at any step size.

Related: [[nonlocal-slide-move]], [[bc-washboard-not-free-spirals]],
[[bilocal-program-saga]], [[crystal-symmetry-group]], [[no-halo-verdict]].
