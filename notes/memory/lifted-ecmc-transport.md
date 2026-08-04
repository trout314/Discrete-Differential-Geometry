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

## CARRIER SWITCH + CHAIN-FOLLOWING WORM (2026-08-03)

The deg-3 chord is the WRONG carrier. It is a trace impurity (3-6 edges vs ~230
defective vertices) whose free moves are conditional and rare, and the census
above shows its momentum dies on the first collision. The deg-4 TIP is the
right one: ~78 tips, and every worm candidate is dS = 0 unconditionally. The
deg-3 chord is not a rival carrier but the CATALYST -- the quasiparticle is a
deg-4 tip dressed with a deg-3 catalyst, which is why lifting the bare chord
was doomed (it tracked half the object).

**wormEnumerate is exhaustive local search, not construction.** Anchor
octahedron (edge + link 4-cycle) -> stage-1 star of tets touching the anchor ->
every legal 2->3 / 3->2 whose support touches the anchor -> commit m1 -> find
disturbed vertices -> stage-2 star -> m2 of the opposite kind -> degree-
bookkeeping class check (gone4/new4 <= 2 each, catalyst <= 1 each, everything
else legal->legal) -> escape test. The ONLY spatial constraint anywhere is the
escape test (leave the cavity). Nothing prefers near or far landings.

**GREEDY MAX-DISPLACEMENT ANTI-SELECTS CHAIN ALIGNMENT.** Chain-aligned
candidates have median displacement 0.237 cells, off-chain 0.357; at one anchor
all 6 off-chain candidates ranked in the top 13 of 60 by displacement and the
single longest was off-chain. So a greedy march leaves every chain at step one
(chain-run 0 on 6/6 tips), while an undirected random march stays on 1-2 steps
and goes nowhere (straightness 0.10-0.34).

I concluded from that "chain-following and transport are in tension". **WRONG,
and the user caught it** -- the random control was UNDIRECTED, so it tested
nothing about small DIRECTED steps. Only a few tets of progress are needed if
the pattern repeats.

**CHAIN-FOLLOWING MARCH WORKS** (carry the frame; pick the candidate with the
SMALLEST positive offset along that same frame and direction; never reverse):

    offsets used across 46 steps:  {1: 42, 3: 3, 4: 1}   <- one tet at a time
    steps sustained:  12,9,12,4,3,6 of 12   median 8   2/6 reach full length
    net per flight:   0.34-1.36 cells       straightness 0.35-0.75
    greedy for comparison: 12/12, net 4.34-4.93, straightness 0.89-0.96

Reading: flights are FINITE and that is the EVENT condition, not a failure --
"no candidate advances along my chain" is geometrically well defined, and
median ~8 steps is the lifted walker's natural flight length. Sub-unity
straightness is expected: a BC chain is a helix, so net is axial drift while
path is arc length. One flight (~1 cell) already matches the caged thermal max
of 1.26 cells per defect LIFETIME.

TRADE: greedy gets ~0.36 cells/step net but is many-to-one with NO valid
kernel and no sigma structure; chain-following gets ~0.1 cells/step but sigma
is a discrete slot (frame + direction), reversal is exact (walk the chain
backwards), and events are well defined. **Build the driver on chain-following**
-- it also removes the cone, the chart and the cocycle-gating question.

OPEN: (1) do flights compose across events (unbounded transport) or does the
helix close a turn and circulate? 12 steps cannot distinguish ballistic from
circulating -- needs NSTEP ~ 40 on the tips that sustain. (2) can the walker
re-acquire a chain at an event and continue with a related direction?

## WHY THE LIFT NEEDS SIGMA -- proved, not assumed (2026-08-04)

On pristine C15 the worm head is a rigid species (95/95 candidates preserve its
decorated patch key) and admits a fully LOCAL move template: anchor (a,b) plus
its link 4-cycle c0..c3, with the deg-3 catalyst on one link edge; the other
three link edges name BC rays that reach a displaced destination at offset
k = 1, both directions. No chart, no distance, no chain id.

**Every such rule is a period-2 orbit.** All three link-edge choices x both
directions run 40/40 steps with composition preserved and graph distance 0
throughout; adding a vertex-disjointness filter does not help, since a 2-cycle
A->B->A is disjoint at each step. Forced by symmetry: |Stab| = 2 and that
element swaps the head's two ends, so the head is locally reversal-symmetric
and "forward" is not a function of local data.

    template -> WHICH MOVE        sigma -> WHICH WAY

So carrying sigma is not a convenience, it is the only way to break a symmetry
the local data cannot. Minimal sigma is ~one bit (which end is the front). This
also explains why the chain-following march moved: its no-backtrack rule was a
1-step memory acting as a crude sigma. Details in [[deg4-worm-design]].

## CARRIER SPLITS BY SPECIES (2026-08-04 correction)

I switched the carrier wholesale from deg-3 chord to deg-4 tip above. That was
too broad. Anything a single 2->3 creates has a deg-3 AXIS edge whose 3->2
undoes the entire excitation, so `nonlocal_slide_at` already annihilates the
whole object and rebuilds it arbitrarily far down a BC chain -- and does it
better than the worm (arbitrary range vs one tet; 12 slots with an exact
inverse; acceptance exactly 1 on same rung). Measured: all four C15 site
classes are 480/480 legal with 54-186 free slides, relocating the species
intact and vertex-disjoint. See [[deg4-worm-design]].

    flickers / anything one 2->3 makes  ->  deg-3 non-local slide (built)
    deg-3-FREE deg-4 bundles            ->  worm (the only option)

The transmit/absorb census's 82% absorb was exactly the second class, so the
worm is genuinely needed -- but the pristine heads I marched were the FIRST
class, which already had a better channel. Do not test the worm on flickers.

## THREE-MOVE FLIGHT DESIGN + THE HOP IS BALANCED (2026-08-04)

Design converged with the user to three move types on a deg-3 chord:

    1 EXIT   dS-gated: leave a complex onto the crystal washboard at rung Q
    2 TRAVEL dS = 0 same-rung hops at CONSTANT Q -- free, no attrition
    3 ENTER  dS-gated: into the next complex; bounce (sigma -> -sigma) if rejected

The rung is CONSERVED during free flight, so Q acts as a carried kinetic
energy: exit onto a high rung and you arrive at the next complex still high,
where more entries are downhill. Only the endpoints are gated.

**User caught a reversibility BUG in an earlier version:** the inverse of every
uphill ENTRY is a downhill EXIT, and a travel rule restricted to dS = 0 can
never propose a downhill move. Entries would have had no reverse and the chain
would pump flickers into complexes -- a broken stationary distribution, not
slow mixing. Fixed by keeping EXIT as an explicit move type (1 <-> 3 pair,
2 pairs with itself). Symptom the user spotted first: "it never stops in the
crystal".

Why travel must stay dS = 0: under the general factorized-Metropolis rule every
uphill step charges, Lambda accumulates ~0.32/step (c = 0.34034, Q in
{46,48,50,52}), so against E ~ Exp(1) a flight dies after **~3 steps** against
a mean free path of 87. Correct and useless. dS = 0 travel has no attrition and
the Q46/Q50 free webs are crystal-spanning (Q48/Q52 fragment) -- so transport
is RUNG-DEPENDENT, two conducting energies and two localised ones.

Also measured: restricting the entry 2->3 to the sigma-ray vs allowing any of
the 4 faces containing the target deg-4 edge. Along the ray, 1 of 480 slides
converts, at dS = +10. Over all containing faces: 24 legal, 20 convert, dS on
an exact ladder {0,2,4,6,8,10}, **2 of them FREE**. The BC walk was landing on
the worst face -- the cost was an artifact of the restriction. CAVEAT: that
chord was planted ADJACENT to the bundle (its link triangle shared bundle
vertices, so supports are not disjoint and the rung formula does not apply);
an arriving chord may behave differently.

**SLOT STRUCTURE AND BALANCE (pristine R m2, steps = 1, 4 sites):**

    12 slots -> 6 distinct 1-step targets, multiplicity EXACTLY 2 each
    24 ordered pairs (C0 -> C1): fwd slots 2, rev slots 2, mismatches 0
    dS exactly antisymmetric: max |dS_fwd + dS_rev| = 0.000e+00
    dS ladder: exact even integers {0, +-2, +-4, +-6}

The 2-fold degeneracy is the TRAILING-LINK TRANSPOSITION: swapping u2 <-> u3 in
the window [u0,u1,u2,u3] leaves face (u1,u2,u3) unchanged as a set, so the
1-step target is identical -- but the window differs, so the rays diverge from
step 2 on.

CONSEQUENCE, and it splits two ways:
  * VALIDITY -- the degeneracy is symmetric (2 vs 2), cancels in the proposal
    ratio, no Hastings factor. The atomic hop satisfies detailed balance as is.
  * PERSISTENCE -- sigma must be the 4-vertex FRAME, not the slot. Holding a
    slot does NOT hold a ray; a flight carrying only the slot forks every hop
    and drifts onto a different chain. nextChainWindow is a bijection on
    frames, so continuation is unique there. (I first stated this as a validity
    requirement; it is only a persistence one.)

**ALL THREE CLOSED (2026-08-04 late):**

FLIP INVOLUTION, concrete:  F(e, p0 p1 p2) = (e', p2 p1 p0) -- swap chord
endpoint, reverse link order. Derived by conjugating the sliding-window map;
verified through the D core: M(F(M(C,w))) = F(C,w) **18/18 PASS** (6 sites x 3
frames, hop lengths k = 1..20, every dS exactly 0, every committed arrival
matching the frame walk). This IS the lifted skew-DB check for the
deterministic travel kernel; with the 2-vs-2 proposal balance the travel phase
is validity-complete. Rung antisymmetry is why the return hop stops exactly at
the start: intermediate positions are off-rung both ways.

ENTRY/EXIT AT A REAL COMPLEX, ARRIVING CHORD (pristine R m2 + CREATE_SEQ
bundle; launches created >= 3 from the bundle, flights 20-43 steps into
contact): 4/4 BALANCED. Every entry has EXACTLY ONE reverse (slot, k)
recreating the source chord, same k, dS exactly antisymmetric (max |dS_e +
dS_x| = 0.00e+00). No Hastings factor survives. The pristine 2-vs-2 slot
degeneracy splits to 1-vs-1 near the defect -- but symmetrically.

TWO STRUCTURAL FINDINGS:
  * FIRST CONTACT IS DOCKING, NOT CONVERSION. All four entries landed the
    chord ADJACENT to the bundle (arrival shares one bundle vertex), dS in
    {-4,-2,+4}; no bundle deg-4 changed degree. Faces containing a bundle
    deg-4 lie DEEPER than first contact, so entry is "dock", with "convert" a
    subsequent move from the docked position.
  * RUNG-DEPENDENCE IS DRAMATIC: one launch had ZERO free sites in 23 steps
    (isolated high rung -- travel cannot move that chord along that frame at
    all); another had 24/43 free. The Q48/Q52-fragment story visible in a
    single flight profile. Exit-rung selection will strongly shape transport.

## DRIVER BUILT + THE HANDOFF IS THE ONLY COUPLING (2026-08-04/05)

`scripts/defect_dynamics/ecmc_flight.py` (commit eeb261e): ONE uniform kernel
-- scan the frame walk for the first site that is clean-and-on-rung-Q or
touching a complex; propose the slide there; Metropolis; flip on rejection.
Travel/entry/exit/deeper are all this one rule. Validated 6 seeds x 1500
steps: drift exactly 0, 62/62 live reverse audits, full event repertoire.
Rung-dependence at driver level: Q46 seed flew 1491/1500 with 5 bounces; the
Q52 seed 75/1500 with 1414. Assembly bug fixed: rung-Q proposals from a DOCKED
chord are the EXIT and legitimately carry dS != 0; the dS = 0 assertion
applies only to clean sources.

**BETA SCAN (1.0 / 0.5 / 0.25, x3 seeds): temperature opens penetration
(deeper accepts 7 -> 145 -> 200; bounces 1414 -> 888 -> 221) but
bundle_changed = 0 ALWAYS -- and that is a THEOREM, not a statistic.** The
background (state minus flier) is invariant under every slide (that is also
the driver's caching optimisation), so under the default handoff (activity
follows the axis chord) a conversion is transient: the freed edge returns to
deg-4 the moment the flier slides on. The heavy sector is EXACTLY frozen at
any temperature. The freed-edge handoff is not an enhancement -- it is the
ONLY mechanism by which this channel can move bundles. (The user's original
"convert the next deg-4 and transfer momentum to it" was the unique coupling
all along.)

## THE HANDOFF RULE (2026-08-05): equivariant bijection, chosen not discovered

Post-conversion INTERLOCK, verified 20/20 sites on the CREATE_SEQ bundle:
freed edge e4 = (u,v) has link {a1,b1,n} containing BOTH endpoints of the
axis chord C1 = (a1,b1); C1's link is Ft = (u,v,x) containing e4. Mutual
local recognizability -> the kernel can branch (slide vs handoff) reversibly
with no memory.

Corridor geometry UNDERDETERMINES the rule: 1-3 of 12 freed-edge frames
continue into the outgoing corridor, unique for only 8-11 of 12, stable at
depths 3/5/8; greedy tie-breaks break bijectivity AND flip-equivariance at
every site. Unlike the flip, the handoff cannot be read off the geometry.

The interlock's matching distinguished structure ({u,v} pair + x singleton
vs {a1,b1} pair + n singleton) gives a canonical algebraic rule:

    h(e, pi) = (first of {u,v} in pi,  sigma)
    sigma: n at the position x held in pi; e at the lower remaining
           position, e' at the higher            (2 x 3 x 2 = 12 -> 12)

**Verified 20/20 sites: bijective, flip-equivariant h(F(w1)) = F(h(w1)).**
Reverse handoff g = F h^-1 F, so conversion+handoff is skew-balanced with
deterministic proposals. Physics quality: ZERO backscatter (0/240 h-images
meet the incoming corridor -- not automatic, backscattering frames exist) but
only 14% direct forward alignment (mostly lateral). The equivariant bijection
is NOT unique -- composing with any equivariant permutation of the freed
frames stays valid, so alignment is optimisable later without touching
validity.

REMAINING: wire h into ecmc_flight as the branch at conversions (branch prob
p, audits extended to handoffs), then the first run where bundle_changed CAN
fire: does the bundle move under flicker bombardment with the coupling on?

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
