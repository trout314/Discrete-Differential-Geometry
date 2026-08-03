---
name: strict-chord-channel
description: "Strict-closure chord (flicker) channel — balanced (two-sided cert RE-VERIFIED 2026-08-02) and now 15.7x faster, but MEASURED NOT to accelerate relaxation (see percolation-ab-test)"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-08-02T14:22:52.364Z
---

The bilocal chord channel whose discrimination lives in the CLOSURE
CRITERION rather than in U. Both marks are pure flags (no move at open
or close); the close fires when one mark is a degree-3 chord, the other
ABSENT, and the absent one's region carries no degree-3 edge.

**Why not an umbrella.** The open-sector weight is zeta e^{-S+U}, so a
high U DRAWS the walk to a state while alpha_close ~ e^{-dS-U} BLOCKS
closing from it — same parameter, same sign, both effects. Measured: big
positive tube entries sent abandonment to 101/200; clamping traded that
for lost opens; refining the key to stop aliasing made the tube exactly
inert. No U makes a state unattractive AND hard to close from.

**Why selectivity cannot be a closure criterion either.** A criterion is
a predicate on the open state; the reverse of the close is the open, so
both gates must be the SAME predicate — hence the close condition holds
already at the open and a null closure is always available. "f changed
since the open" is not such a predicate, and storing the opening f does
not help (the reverse open would set the reference to the CLOSING f).
Selectivity therefore belongs in the MEASUREMENT: WormF0Result.df is a
per-episode f census, and post-selection on it is exactly unbiased.

**Auto-calibration (both derived together, frozen per episode).** The
open makes no move, so there is no dS term and
    zeta2* = -(log n3 + log nSite + log p_close)
in closed form, NO probe (unlike calib_zeta for the vertex carrier).
p_close is not independent — it sits inside that formula — so it is
derived too: p_close = pclTarget/maxstep, giving mean episode length
maxstep/pclTarget and abandonment ~e^-pclTarget. At zeta2 = 0 the open
saturates at +9.1 and the close dies at e^-9.1: NOTHING closes.
Deriving both took commits 2 -> 33 per 200 episodes, 76% reactions.

**BALANCE TEST (2026-07-31, 4000 episodes each). The BALANCE half is
valid; the CONVERGENCE half was ill-posed — see the design error.**
  start f3=8724 n3=17 -> late-half n3 16.3 +- 1.0, slope -0.19/1000ep
  census: 34 down / 33 up / 24 null, net df3 = -1
  start f3=8699 n3=44 -> late-half n3 45.3 +- 1.4, slope -0.18/1000ep
  census: 89 down / 82 up / 59 null / 5 (+2) / 1 (-2), net df3 = +1
Up and down reaction rates agree to ~5% within EACH run and the net
drift is ~0 — a genuine within-sector balance signature, and it retires
the earlier worry (n3 44->49 over 200 episodes was small-sample noise).

**VALID TWO-SIDED TEST (same f0 sector, 5000 episodes each) — PASSED.**
Both starts built from quench_down5q (f0=1536); the HI one prepared with
60 forced 2->3 flips, which conserve f0 and raise n3/f3 directly.
  LO  start n3=17,  f3=8724 -> late-30%% n3 17.8 +- 0.9, f3 8724.8, S 268
      slope +1.37 n3/1000ep    census 48 up / 47 down / 36 null
  HI  start n3=119, f3=8784 -> late-30%% n3 26.0 +- 1.7, f3 8735.3, S 601
      slope -3.24 n3/1000ep    census 48 up / 69 down / 24 null
The CENSUS ASYMMETRY is the cleanest evidence and it is exactly right:
at equilibrium (LO) up/down = 48/47, ratio 1.02, no bias; displaced (HI)
up/down = 48/69, ratio 0.70, a clear restoring bias toward removal. A
balanced sampler shows no bias at equilibrium and a restoring bias when
displaced -- both observed, with the SAME up-rate (48) in each, i.e.
only the removal channel responds to the displacement.
HI relaxed n3 119 -> 26 and S 4671 -> 601, closing 92%% of the initial
102-chord gap. It is still descending (slope -3.24) while LO is flat, so
it has NOT stalled -- an earlier reading of three identical checkpoints
as a stall was a plateau inside the noise. Residual gap ~8 in n3 after
5000 episodes; full equilibration would need a longer run.

DESIGN ERROR, do not repeat: the two starts were the f0-campaign's
quenched (f0=1536) and harvested (f0=1522) states, which differ by
exactly the 14 forced vertex removals of the harvest — and THIS CHANNEL
CONSERVES f0 EXACTLY (pure-flag marks, 2<->3 walk moves, global kernel
skips 1<->4; every df census entry has first component 0). The two runs
therefore sat in DISJOINT f0 sectors and could never converge. Any
conclusion drawn from their non-convergence — about basins, mixing, or
n3 being a flat direction — is void. A valid two-sided test needs two
states with the SAME f0 and different n3, e.g. built by forced 2->3
flips, which conserve f0 and raise n3 directly.

**RE-VERIFIED 2026-08-02** after the wRev>0 fix AND the _vertexWitness
change (which shifted the RNG stream), same script, 5000 episodes each:
  LO  n3 17  -> late 17.9 +- 2.2   slope +1.30 +- 0.74 /1000ep
      census 50 up / 51 down / 29 null   ratio 0.98
  HI  n3 119 -> late 26.1 +- 1.6   slope +0.84 +- 0.45 /1000ep
      census 68 up / 99 down / 52 null   ratio 0.69
      S 4670.7 -> 548.8
Late means match the original to 0.1 (17.8/26.0), ratios to 0.01-0.04,
and the ~8 residual reproduces (26.1 - 17.9 = 8.2). f0 conserved in
EVERY commit, both arms (asserted df0 == 0 across the whole census).

RETIRE ONE CLAIM: "the up-rate is IDENTICAL (48) in both arms" does NOT
reproduce (50 vs 68), and total commits differ between arms (130 vs
219), so up-counts scale with overall activity and there is no reason
they should match. The 48/48 was coincidence. The ROBUST signature is
the up/down RATIO (0.98 at equilibrium vs 0.69 displaced). Reading the
certificate as "up-rate must be identical" would false-alarm on a
correct channel.

Also note HI's late slope is window-dependent: over the second half it
is ~0 (already relaxed), not the -3.24 recorded from a window that still
contained the descent. Compare late-window MEANS and RATIOS, not slopes
across differently-defined windows.

**BUT IT DOES NOT HELP.** [[percolation-ab-test]]: measured 2026-08-01/02,
neither rmax=0 nor rmax=2/beta=2 accelerates percolation relaxation; the
plain-work effect is consistent with zero and the overhead is a 1.6x
penalty at 40% work-fraction. A balanced, fast, correct channel that
does not do the job it was built for.

See [[flicker-catalysis]] (the wide head class), [[bilocal-factorization]].
