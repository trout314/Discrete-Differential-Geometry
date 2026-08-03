---
name: flicker-background
description: "The (3,4,4) knot IS the shape of one 2->3 on pristine crystal; ~96% of defect births at lam=0.40 are that flicker; long-lived defects are sparse/curve-like, not dense"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T21:17:17.308Z
---

Measured 2026-07-25. Tools: `scripts/defect_dynamics/crystal_flicker.py`
(exhaustive enumeration of single 2->3 on a pristine crystal) and
`flicker_fraction.py` (attributes observed worldline births against that list).

**A single 2->3 on pristine R always makes f = (5,10,9,3) — the knot shape.**
Proof, then verified on all 115968 triangles of T3_R_m4_N57984: the move changes
exactly nine edge degrees — abc by −1, the six side edges by +1 — and creates xy
at degree 3; ALL nine have both endpoints in {a,b,c,x,y}, so no outside vertex
can change legality. 100.00% of sites give one connected component with
f=(5,10,9,3). **The (3,4,4) knot is not an excitation of the melt — it is the
shape of one Pachner move against a perfect crystal.**

47 distinct decorated species; 6 illegal-edge sigs, always exactly one degree-3
edge: (3,4,4,4,7) 41.7%, (3,4,4) 23.8%, (3,4,4,4) 13.6%, (3,4,4,7) 9.9%,
(3,4,4,4,7,7) 8.0%, (3,4,4,7,7) 3.0%. #4s = degree-5 edges in abc (2 or 3),
#7s = degree-6 side edges (0–2). ΣZ ladder 62–76, mode 70. m2 and m4 give
identical species sets (counts exactly 8x) — crystal-local, not a box artifact.

**At lam=0.40 the per-move census is ~95% flicker.** 800 sweeps, 1021 births:
96.5% born by 2->3 with a key in the pristine list; 94.8% are 2->3 in / 3->2 out
with median life 1.91 sweeps; 97.7% have f=(5,10,9,3). Only 3.5% born otherwise
(28 by 4-4, 8 pre-existing).

**But birth share != population share.** Time-integrated (life x n_verts), flicker
is only 71% of the standing defect population; the 3.5% non-flicker births carry
29%. And 49 of 57 tracks living >20 sweeps were flicker-born — being born as
flicker does not doom a complex. The 4-4 move is the channel that makes
non-knot shapes, and they live far longer (99th pct 589 vs 56 sweeps).

**GOTCHA:** the "2->3 into disturbed background" class came out 0.00% — NOT
because it doesn't happen, but because a 2->3 adjacent to an existing defect
merges into that complex instead of creating a new worldline, so it is not a
birth by construction.

**Long-lived defects are the structural opposite of the flicker** (from the
2500-sweep run, see [[lifetime-vs-charge]]): 166/3544 exceed 20 sweeps. They are
big (survival rate rises monotonically 0.8% at n=5 to 100% at n=12), SPARSE
(median 1-skeleton density 0.71 vs 1.00 — the knot is K5, maximally dense; the
two longest have f=(7,7,1,0), zero tets, genuinely curve-like), and constantly
rearranging (5/166 pure-shape; mean mode_frac 0.383 vs 0.868). Q_c per vertex
−0.194 vs −0.332: same total charge spread thinner.

Figures: `data/figs/flicker_rungs_R.png` (pristine rung distribution vs observed
— melt enhances ΣZ 67–70, suppresses >=71, and rungs 74/75/76 never fire despite
being 5.9% of pristine sites).

Implication for any species table: report the flicker attribution alongside, or
the census is mostly measuring the move set, not the state. See
[[defect-species-ladder]] and [[reporting-conventions]].
