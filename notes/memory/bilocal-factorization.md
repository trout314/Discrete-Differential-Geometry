---
name: bilocal-factorization
description: Bilocal step-1/2 MEASURED — exact factorization iff supports are vertex-disjoint; conserving pairs are pin-free and often exactly free; BC-chain index is NOT a separation coordinate
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-31T15:31:21.471Z
---

Measured 2026-07-31 on quench_down5q_wOFF (C15 m4, f=[1536,10260,17448,8724],
gap +10.12) with the worm action (F.fresh: c_v=0.1, c_h=1.0, e*=5.1067907,
cimp=0.7). Scripts (scratchpad): `bilocal_factorize.py` (BC-chain scan,
T1 non-conserving + T2 conserving), `bilocal_range.py` (general sites by
BFS shell). JSON: `bilocal_factorize.json`, `bilocal_range.json`.

**1. The analytic pin cross term is EXACT.** For two regions with f-offsets
dA, dB the global part never factorizes; it carries an s-INDEPENDENT
X(dA,dB) = 2 c_v dA3 dB3 + 2 c_h dwA dwB, dw = df1 - 6 df3/e*. For two
knots X = +0.2612 exactly. Subtracting it leaves residuals at 1e-14 at
every separation. This independently validates the same quadratic
identity `wf0Compile` uses (see [[f0-surgery]]).

**2. s_min is NOT a distance — it is SUPPORT DISJOINTNESS.** A 2->3 on
face (a,b,c) with apexes (x,y) changes only edges inside {a,b,c,x,y}, so
the m^2/edge support IS those 5 vertices. Measured: graph distance 0
(shared vertex) => residual up to +29.4; distance >= 1 => 1e-13 EXACTLY,
at every d from 1 to 8, both along a chain and over general sites. So the
criterion is necessary and sufficient: two head regions factorize iff
their vertex supports are disjoint. For deep corridors the support is the
union over the corridor's moves -- bigger, still explicit and cheap. Do
NOT impose a minimum separation; impose disjointness.

**3. BC-chain index is NOT a separation coordinate.** The length-240 chain
re-approaches itself: graph distance runs 0,1,1,1,2,2,2,...,6 then back to
0 at s=48 and s=96, where the residual jumps to +5.6/-8.4. A chain-index
guard would silently permit OVERLAPPING regions. Any bilocal guard must
use graph distance / explicit support disjointness. (Consistent with the
site-survey finding that the clean slide network branches off-chain.)

**4. CONSERVING pairs are pin-free, exactly.** Create knot at A +
annihilate knot at B (dA+dB = 0): measured pin content = 0.0e+00 at every
separation, and dS_AB = washboard(A) - washboard(B) exactly -- quantized
in multiples of 1.4, values {0,-1.4,-2.8,-4.2,-5.6,-9.8,...,-22.4}, and
EXACTLY ZERO for 25-62% of pairs at every distance out to the box
diameter. Compare the single-head f0 move, which pays +22.75 in pin alone
at this gap. THIS IS THE DESIGN CLAIM CONFIRMED: the paired conserving
move lives on the pins' constraint surface and should fire trivially where
the single-head channel needed the whole tube/auto-zeta apparatus
([[f0-surgery]]).

**5. Box geometry.** BFS reach = 8 hops covers all 1536 vertices, so d=8
IS the whole manifold on m4 -- "long range" here means d<=8, and every
d>=1 is equally exact.

Consequences for the pair design ([[deg4-worm-design]],
notes/bilocal-worm-design.md 3.2): U factorizes as U(h1)+U(h2) using the
EXISTING single-head tube (no new umbrella at range); V(s) === 0 beyond
support contact, so there is no long-range V table to build; the pair
proposal picks B anywhere with disjoint support; zeta_2 should be far less
state-sensitive than the single-head zeta since the pair cost is pin-free.
CAVEAT: measured with the ELEMENTARY head op (one 2->3). Re-check
disjointness at full corridor depth before relying on it there.
