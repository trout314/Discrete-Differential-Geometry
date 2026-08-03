---
name: nonlocal-slide-move
description: "Non-local slide (undo 2->3 + redo 2->3 down the spiral) ported to D + validated; scripted entry live, mcmcStep channel + DB proposal still to wire"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-28T04:01:51.327Z
---

**New sampling move, prototyped + validated 2026-07-27**
(scratchpad/nonlocal_slide_proto.py). Builds on [[bc-washboard-not-free-spirals]].

**PERF FIX -- O(N) hint (2026-07-28), commit 782e9c8:** the channel built its
slide hint with `mfd.link(chord)` = O(N) (star scans ALL facets via StarRange),
once per proposal -> ~690x slower than baseline at N3=58k (cost scaled 8x from
N3=7248 to 58k, matching size). FIX: Deg3Set stores a WITNESS facet per chord;
draw validates it O(1) via writeFaceApexes, falls back to O(N) mfd.link (then
caches) only when stale; the nonlocal move supplies the witness for the chord it
creates (moved defect always O(1)); rebuildDeg3 iterates FACETS so initial
chords get witnesses too. Result 58k: 542 -> 107534 moves/s (~200x); ON is 3.5x
of OFF (real per-move work). GOTCHA: mfd.link / mfd.star / mfd.contains are ALL
O(N) (facets.any / StarRange) -- never call them per-move in the sampler; use
the ridge maps (writeFaceApexes O(1)) or a cached witness.
TRIED on real melts (full action + n6 pot): lam40 (4->11 chords) 22% accept,
lam30 (36->44) 4% accept, set stays CONSISTENT vs recount; ~1.4x spread of
degree-3 activity ON vs OFF. Modest because the FULL action makes long jumps
mostly uphill (Metropolis-rejected) -- the washboard is real energy, not a
mechanical barrier. Biggest wins expected in the dilute high-lambda regime and
with the lifted/ECMC variant (walk-to-contact, no reject).

**NATIVE mcmcStep CHANNEL WIRED + VALIDATED (2026-07-28):** the 1/n_3 channel
now runs inside the D sampler. Added to sampler.d: NonlocalSlideConfig
(prob, maxStep) + Deg3Set!Vertex (swap-remove array + edge->index map, O(1) draw)
+ reconcileDeg3/rebuildDeg3; mcmcStep gained a `nlSlide`/`deg3Set` channel that
picks a deg-3 chord uniformly (deg3Set), draws slot x steps(1..maxStep),
tryNonlocalSlide(metropolis, n3Before=|deg3Set|). The set is GATED (only live
when nlSlide.prob>0), rebuilt at each run start, and maintained incrementally by
reconcileDeg3 at EVERY commit site: hinge, bistellar (over `ball`), the old
slide (deg3Set threaded into trySlideMove), and the nonlocal move (over its
annSup+redoSup). Null pointers => existing runs untouched. C API:
ddg_sampler_set_nonlocal_slide_prob(h,prob,max_step),
ddg_sampler_nonlocal_slide_stats, ddg_sampler_deg3_count (diagnostic). Python:
ManifoldSampler.set_nonlocal_slide_prob / nonlocal_slide_stats / deg3_count.
VALIDATED (scratchpad/native_channel_validate.py): (1) prob=1 single defect,
60k moves via sampler.run, samples crystal-wide Boltzmann (66% accept);
(2) thermal+nonlocal prob=0.3, ONE 60k-move run, deg3 set size == full recount
(3==3) => reconcile stays EXACT across all move types, no drift. meson test
1/1 OK. So the channel is production-usable in equilibrium runs now.

**MECHANISM EXTENDED for the equilibrium channel (2026-07-28):** tryNonlocalSlide
now also returns dn3 (exact change in the degree-3 edge count, counted over the
two moves' O(1) supports) + the arrival chord, and applies the 1/n_3 Hastings
factor n3/(n3+dn3) in metropolis mode (param `n3Before`, -1 to omit). C API
ddg_sampler_nonlocal_slide_at gained out_dn3/out_ta/out_tb; python
nonlocal_slide_at now returns (dS, dn3, arrival_chord). Validated
(scratchpad/nonlocal_channel_validate.py): (A) dn3 == full recount exact over
1004 moves incl. a dn3!=0 case; (B) **k_f==k_r for k>1**, 300 checks 0
mismatches (the closure test the user asked for); (C) 1/n_3 + n_3-factor
Metropolis samples crystal-wide Boltzmann to ~0.01. meson test 1/1 OK.
DESIGN AGREED: proposal = pick deg-3 edge sigma uniformly (1/n_3) x fixed 12
slots x symmetric-k; direction lives in slot orientation; accept
min(1, e^{-dS} n3/(n3+dn3)). Dilute-isolated-(3,4,4) gas => dn3=0 => factor=1 =>
plain Metropolis exact. STILL TO DO for the D-NATIVE mcmcStep channel: a
MAINTAINED degree-3 chord set in SamplerState (swap-remove array + edge->index
map, updated at each commit site over the move support; GATE on channel enabled
to keep existing runs untouched) so the 1/n_3 draw + n3 are O(1). No such set
exists yet (illegal_edges scans O(E)); that's the next piece.

**D PORT DONE + VALIDATED (same day):** source/sampler.d `tryNonlocalSlide`
(annihilate 3->2 + walk `steps` nextChainWindow + redo 2->3; exact dS via
speculativeBistellarDelta; exact rollback); source/ddg_capi.d
`ddg_sampler_nonlocal_slide_at(h,a,b,slot,steps,mode,out_dS)` (0 trial/1 force);
python ManifoldSampler.nonlocal_slide_at(a,b,slot,steps,commit). Validation
(scratchpad/nonlocal_slide_d_validate.py): dS exact to 4e-14, trial rollback
exact, samples crystal-wide Boltzmann to ~0.005. STILL TO DO: wire the sampler
CHANNEL in mcmcStep with a SYMMETRIC proposal (scripted entry only so far, like
slide_at) -- that's where detailed balance lives. Pre-existing manifold.d
unittest purity break (NOT this work) blocks `meson test`; release lib fine.

THE MOVE: instead of the 4-move local slide, move a 2->3 defect by (1) a 3->2
that UNDOES the creating 2->3 (annihilate -> pristine locally), then (2) a 2->3
that RE-CREATES it ANY number of tets down the BC spiral. Net dS = c*(Q(m)-Q(n))
= difference of the two 2->3 creation costs. FREE <=> same rung Q(m)=Q(n).
Only 2 Pachner ops, no deriveSlideFrame / degree-4 collision machinery -- it
steps THROUGH the pristine intermediate, so it can't get frame-blocked.
The local +-4 slide is the max-step=4 special case (M1=undo, M3=redo-4-down,
M2/M4 only bridge the endpoint OVERLAP; a distant redo needs no bridge).

VALIDATION (R m2, one spiral L=105, single (3,4,4), pure edge-degree action):
- Real do_bistellar ops: max|action_real - (S0 + c*Q + offset)| = 5.7e-14
  (machine 0) over 8000 moves; N3 conserved; all moves mechanically valid.
- Samples correct Boltzmann pi ∝ exp(-c*Q): rung occupancy obs/theory
  46:.530/.540 48:.187/.182 50:.253/.246 52:.029/.031.
- MIXING: tau_int(position) = 158.4 (local +-4) vs 2.4 (non-local) = ~66x, and
  since local is diffusive (tau~L^2) the win GROWS with L (762,1626 -> 100s-1000s).

USE IT FOR EQUILIBRIUM SAMPLING ONLY. For real dynamics/transport (FPKMC,
worldlines, non-equilibrium) keep it LOCAL (max-step=4): the 3->2 annihilates
the defect to vacuum, so a distant redo is re-nucleation, not motion -- breaks
the worldline, gives infinite D, no caging. Plan: max-step knob interpolates
(4 = dynamics limit, large = fast sampler). User wants equilibrium now,
dynamics later.

TODO / caveats: detailed balance needs a SYMMETRIC proposal (pin the spiral
deterministically from chord+orientation so fwd/reverse walk the same loop) and,
multi-defect, reject targets colliding with another defect. Hastings-clean
variant: propose uniformly among SAME-RUNG windows -> free moves accepted at
alpha=1. Same rung != same species on rungs 48/50/52 (add [5,5,6] filter to
preserve (3,4,4) identity). NEXT (user to pick): (1) port to D sampler.d + C API
w/ max-step knob; (2) L^2 scaling test + same-rung-uniform proposal;
(3) full lam=0.40 thermal trial vs current sampler. See [[worm-sampler-program]],
[[fpkmc-design]], [[fp-kinetics-findings]].
