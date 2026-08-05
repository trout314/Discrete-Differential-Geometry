---
name: flight-contact-barrier
description: "DECISIVE: the wall that stops a flicker entering a defect complex is 100% the m^2 anti-clustering term we chose -- pins+geometry contribute EXACTLY ZERO; without it 100% of contacts would cost < 3, in every rung sector"
metadata:
  type: project
---

Measured 2026-08-05. Tool: `scripts/defect_dynamics/contact_census.py`
(+ `_merge.py`). Data: `data/ecmc_ab/cc_dec_*.json`.

## THE RESULT

Each rejected contact re-priced with the m^2 coefficient temporarily zeroed
(`set_n6_potential(0.0, 0.0)` around a `commit=False` re-evaluation of the
IDENTICAL proposal, then restored). 8,084 rejected contacts, C15 m3 blob:

    rung   n      dS median   pins+geom   m^2     m^2 share   frac<3 w/o m^2
    all    8,084     16.00       0.00     16.00      100%         100.0%
    Q=38     360     22.00       0.00     22.00      100%         100.0%
    Q=40   1,629     17.00       0.00     17.00      100%         100.0%
    Q=42   2,939     16.00       0.00     16.00      100%         100.0%
    Q=44   1,871     15.00       0.00     15.00      100%         100.0%
    Q=46   1,029     12.00       0.00     12.00      100%         100.0%

**Pins and geometry contribute EXACTLY ZERO to the contact barrier, in every
sector.** The entire 12-22 wall is `V(m) = cimp * m^2` from
`set_n6_potential(zleg, imp)`, m = #incident edges of degree outside {5,6},
cimp = 0.5 in the A/B. We wrote a quadratic penalty on defects having
defective neighbours, and it is the only thing stopping defects touching.

The apparent RUNG ordering of the barrier (22 at Q=38 down to 12 at Q=46) is
not rung physics: low-Q flights sit near BIGGER blockers (Q=40 -> 12.3
vertices, Q=42 -> 13.3, Q=46 -> 7.2) and m^2 grows with neighbour count.

## cimp=0 IS NOT A CONTROL

Setting cimp=0 for a whole run gives ZERO contacts and rungs running to -30:
the state degenerates without the term (it is what holds the mobile gas --
[[m2-only-gas]]). The barrier must be decomposed per-proposal, as above, not
by turning the term off globally.

## THE HANDOFF RECEIVER -- the mechanism I got wrong twice

`partner()` does NOT look for a pre-existing deg-3 chord. `minus_flier` ADDS
1 to each of the three link edges when it removes the flier, so an edge that
is deg-3 LIVE is deg-4 in the background -- and `partner()` searches exactly
the link-vertex pairs for LIVE degree 3. **The receiver is a background
deg-4 edge that the flier's own 2->3 pushed down to deg-3** (the user's
description, verified in code).

So the abundance question is deg-4, not deg-3, and fused complexes ARE deg-4
structures: mean 4.78 (size 7-9), 7.14 (10-14), 13.75 (>=15) deg-4 edges per
complex over 74 thermal snapshots. Receivers are everywhere.

RETRACTED along the way: "the heavy sector cannot receive a lift" (wrong --
it is made of receivers); "62-94% of blockers hold a carrier" from the blob
(an artifact of `build_blob` packing 12 fliers, each of which IS a deg-3
edge -- but it measured the wrong edge degree anyway); and the whole
Q-refresh line, which assumed a missing MECHANISM where there is only a
barrier.

## WHAT SURVIVES

* ~50% of flight episodes have zero free (same-rung) continuation anywhere
  within a 60-step walk, identically early and late. Geometric, unaffected by
  m^2. But with contacts at dS < 3 a flight would not NEED a free
  continuation -- it could enter and hand off.
* `self.Q` is set once (ecmc_flight.py:112) and never updated -- not after a
  handoff, not after a dock. A real bug: post-handoff the flight hunts the
  OLD carrier's rung. One-line fix, still unfixed.
* The A/B verdict ([[ecmc-blob-ab]]) was measured with p_hand=0.0 (zero
  handoffs) and refresh hardcoded every 10 steps. Both now exposed/known.

## THE OPEN MODELLING QUESTION

We want defects to collide and exchange momentum, but the m^2 term is what
keeps the gas mobile and unclustered ([[m2-only-gas]]) and lam=0.30 shows the
runaway when the constraint is too weak. Either there is a cimp window, or
the term needs a different SHAPE -- a short-range core plus a flat outer
region rather than a quadratic that grows with every neighbour. Not yet
decided; needs the cimp scan (collision + handoff rate vs clustering
runaway).

Related: [[lifted-ecmc-transport]], [[ecmc-blob-ab]], [[m2-only-gas]],
[[bc-washboard-not-free-spirals]], [[no-halo-verdict]].
