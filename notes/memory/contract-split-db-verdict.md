---
name: contract-split-db-verdict
description: "FIXED: the circulation was a LABEL-BOOKKEEPING DB bug (min(u,v) survivor rule + an over-eager pool push). Survivor coin + hole-only push -> DB residual 6e-18 by exact enumeration, and the measured circulation 1.36/sweep -> 0.01"
metadata: 
  node_type: memory
  type: project
  originSessionId: e3632925-132f-4abe-992c-644edc57ae77
  modified: 2026-08-13T18:52:09.939Z
---

Re-opened 2026-08-13 the worry (from commit e8051df) that the edge
contract/split channel breaks detailed balance. **The channel passes every
statistics-free test; the previous attribution of the residual circulation to
it is RETRACTED.**

## What was measured

1. **Reproducer, 20 s**: `validate_contract_split.py --circulation`. Net f0
   flux per family in a mixed stationary chain (churned sphere, f3 240,
   nh 0.1, cs 0.3): channel **+1.08/sweep**, bistellar **-1.08/sweep**,
   z ~ 30. Not a transient (stable 800 -> 3000 sweeps). Ring-cap scan
   3,4,5,6,7,8 -> 1.34, 1.49, 1.31, 1.08, 0.85, 0.72: **present at maxRing=3**,
   so NOT the long-cycle DFS / FK-catalog cycle machinery.
2. **dS is exact**: tracked objective vs the closed form
   FC(f3-T)^2 + NH(f1 - 6f3/hdt)^2 (using f1 = f0 + f3, Euler) -- max drift
   **0.000e+00** over 10^6 moves at ring 3/6/8. Not an action bug.
3. **The formula is EXACTLY reversible.** Full enumeration from a real state
   (f0 27, f3 133, ring 3): every split path (w, gamma, side) grouped by
   target, every contract path back, exact K(x->y)/K(y->x) vs exp(-dS) ->
   **1.0000 for all 865 targets**. Forward/reverse path multiplicities match
   1<->1 (696), 2<->2 (36), 4<->4 (133), so the sub-kernel pairing is a real
   bijection.
4. **The implementation matches the formula**: branch-resolved rate audit,
   500k single-move trials from a fixed state, measured/predicted
   ring 3 split **0.992 +- 0.010**, contract **0.999 +- 0.010**;
   ring 6 split **0.997 +- 0.004**, contract **0.991 +- 0.012**.
   (An earlier 60k-trial version had only +-2% and was NOT sensitive enough --
   the imbalance to detect is only ~4%.)
5. Debug build (asserts live) fires nothing.

## The methodological point that matters

**Circulation does not localise.** A biased kernel that CONSERVES f0 (4-4
hinge, 2<->3) tilts pi just as effectively, and then *both* f0-moving families
show equal and opposite currents. So "circulation != 0" proves only
pi != exp(-S) -- it does not convict either f0 mover. This is what the
e8051df reasoning got wrong.

## Bistellar side, re-derived

The unified proposal loop is a rejection loop, so q(m|x) = c(m)/C(x) with
c = 1 for a facet center (1->4) and 2 otherwise. The applied correction
`(c_rev/c_fwd) * (f3/f3_after)` covers c and f3 but NOT the loop
normalisation, so the default chain is stationary for exp(-S) * C(x)/f3(x).
**Measured `d ln(C/f3) / d f0 |_f3 = +0.0003 +- 0.0004`** (C computed exactly
in Python: 1->4 + 2->3 + 3->2 + 4->1 + hinge). So that residual is f0-flat and
CANNOT be a vertex fugacity. The `f3new = f3 + 2*center.length - 2 - dim`
and c-assignment in the patch are algebraically correct.

## The f0-conserving sector is CLEAN too (audited 2026-08-13)

Same fixed-state rate-audit method, 1.2M single-move trials. `bistellar_tries`
/ `hinge_tries` are bumped AFTER the validity check, so q (proposal) and
q*alpha (acceptance) are audited SEPARATELY. Loop model: q(m|x) = c(m)/C(x),
c = (facet routes) * min(1, 2/d) keep-coin:

    1->4  center facet,    d=1 -> c=1     2->3  center triangle, d=2 -> c=2
    3->2  center edge,     d=3 -> c=2     4->1  center vertex,   d=4 -> c=2
    hinge deg-4 edge -- NO keep-coin at all: 4 routes * 1/2 diagonal = c=2/move

    move    q ratio             q*alpha ratio
    1->4    0.9997 +- 0.0018    0.9968 +- 0.0041
    2->3    0.9993 +- 0.0013    1.0006 +- 0.0022
    3->2    0.9977 +- 0.0027    0.9977 +- 0.0027
    4->1    1.0074 +- 0.0085    1.0074 +- 0.0085
    hinge   1.0038 +- 0.0024    1.0038 +- 0.0024

Every kernel matches its design to <=0.5%. The hinge preserves the f-vector,
so with the variance coefs at 0 its dS is exactly 0 and it is always
accepted; its reverse is a hinge at the added edge (also degree 4), so
c_rev = c_fwd and taking no Hastings factor is correct.

## THE CONTRADICTION (open, and it is sharp)

Established, each independently:
* the channel is reversible for **exp(-S)** (exact enumeration, 865 targets);
* the default bistellar+hinge sector is reversible for **exp(-S) * C(x)/f3(x)**
  (derivation checked term by term; C's composition validated to 0.1%);
* so the two families disagree by exactly the potential ln(C/f3), and
  **<delta ln(C/f3)> per 1->4 move = -0.0034 +- 0.0003** (240 real moves,
  f3 257, C ~ 993);
* yet the circulation is a genuine STEADY STATE -- flat over 30k sweeps
  (1.03 -> 1.01, ratio 1.011), and **size-independent at ~14% per accepted
  channel move** (targets 60/120/240/480 -> 14.3, 13.5, 14.5, 15.6%).

+14% imbalance needs an acceptance offset of ~0.28 nats. The only structural
difference between the two targets supplies 0.0034. **Short by a factor 83.**

Note the 14% holds at f3 ~ 134 -- the very state where the rate audits found
<=1% error. Those audits are AGGREGATE per branch, so compensating per-target
errors (e.g. a wrong cycle count N on a subclass of links, over-accepting some
splits and under-accepting others) would survive them.

### Next test
**Per-target rate audit**: from a fixed x with cs = 1, key each accepted move
by the resulting labeled state and compare the per-target measured frequency
against the per-target prediction (`apply_split` already generates both).
Aggregate agreement with per-target disagreement is exactly the signature a
compensating-error bug would leave.

## n_6 specialisation: the FK catalog (checked 2026-08-13)

There is exactly ONE n_6/FK-specialised piece in the channel:
`link_cycles.matchCatalog`, the Z12/Z14/Z15/Z16 fast path -- and the SPLIT
branch uses it while the CONTRACT branch always DFS-counts the merged link.
Two different code paths for the same quantity, in opposite directions.

The n_6 = 4 two-class discovery is exactly the hazard, and the stakes are real:
the catalog's Z16 is the **T_d Friauf**, the **D2 isomer has the same 16
vertices / 28 faces** (so matchCatalog's cheap guard cannot separate them),
and their cycle counts agree through length 5 then diverge at **length 6**:

    Z16 T_d (catalog)  28 42 96 282 960 3237
    Z16 D2  (isomer)   28 42 96 296 1008 3212

i.e. a false match would be SILENT at maxRing <= 5 and wrong by 3% in N at the
production maxRing = 6. Real D2 links exist in production states:
`data/crystal_gas/r_c0.50.final.mfd` has 5 Z16_D2 beside 170 Z16 T_d.

**The guard holds.** matchCatalog keys on exact triangulation isomorphism
(propagateIso walks the dual graph, so every face is verified, and it ends
requiring a complete injective vertex map) -- not on n_6 or (V, F). Fed a real
Z16_D2 link from the R crystal it returns **-1** and falls back to DFS with the
correct counts (verified by compiling a probe against source/link_cycles.d;
`data/extract_d2.py` regenerates it).

### But this exposed the real coverage gap
**Every circulation and audit number in this file came from a host where the
catalog NEVER fires**: churned spheres are 100% "impure" links (0 genuine
catalog matches), and the a15 crystal MELTS to impure under the pure geometric
action (FK fraction -> 0.000, circulation still +15.3% / +10.8% per accepted
move at ring 3 / 6). So the channel has never been DB-tested in the regime it
is actually used in, where matchCatalog fires on ~100% of vertices.

To hold FK structure the n_6 potential must be on (cimp > 0) -- which is how
production runs. Note dS is then no longer a closed form of (f0, f3), so the
rate audit needs the potential modelled; **circulation needs no theory and
still works**. Also note the existing `version(ExpensiveAsserts)` guard checks
only the catalog COUNT, not the transported cycle LIST, so a bad `perm` would
still slip through.

## RESOLVED FOR PRODUCTION: circulation is ZERO with the n_6 potential on

Requested run, host `data/crystal_gas/r_c0.50.final.mfd` (f0 1272, f3 7258),
hdt = 6 f3/f1, cs 0.6, ring 6, circulation per ACCEPTED channel move:

    zleg=cimp   FK frac   Z16_D2   chan acc/sweep   circ per accepted move
      0.6        0.746      42          1.05         -0.07% +- 0.62%
      0.4        0.367      77          4.00         +0.26% +- 0.28%
      0.0        0.000       0        236.21         +9.59% +- 0.10%

**With the potential on the bias is gone** (14% would have been 22-50 sigma).
The differentiator is the POTENTIAL, not the host, not size, not the catalog,
and NOT the chemical potential -- mu is ~-1.5 in both the zleg 0.6 crystal run
and the churned-sphere runs. Note the zleg 0.4 row is 68% impure links, close
to the sphere's composition, and still shows zero.

**So production runs (n_6 potential on) are unaffected.** The bias lives in
the potential-off regime, which is a test artefact, not a production setting.

### Mechanism candidate: the hard capacity caps are ASYMMETRIC
`tryContractSplit` has STAR_CAP = 96 and `link_cycles.maxLinkVerts` = 32:

    split at w       needs |star(w)| < 96
    reverse contract needs |star(w)| + |gamma| < 96

so a split at |star(w)| in [96 - |gamma|, 96) is accepted while its reverse
contraction is rejected as `noValid` -- a ONE-WAY move pumping f0 UP, which is
the observed sign. Measured |star| / |lk| distributions:

    churned sphere (14% circ) : |star| max 58,  |lk| max 31 -- 0.00% in window
    R melted (+9.6% circ)     : |star| max 278, |lk| max 141 -- 7.94% in window
    R pristine FK (0% circ)   : |star| max 36,  |lk| max 20  -- 0.00% in window

7.94% in the one-way window vs +9.59% measured is a striking match for melted
R, and FK states never come near either cap -- consistent with the zero. But
the churned sphere has 0% in the window and still shows 14%, so **the caps are
a real asymmetry worth fixing but are NOT the whole story**.

### Remaining lead for the potential-off sphere case
The potential forbids exactly what the bias needs: it drives every edge degree
into {5,6}, so gamma of length 3-4 (and the degree-3/4 edges a short-ring
split creates) never occur. Circulation on the sphere was LARGEST at
maxRing = 3 (1.34/sweep). So the bias looks concentrated in SHORT-RING
splits/contracts -- precisely the sub-kernel the potential suppresses. Next:
resolve the rate audit per ring length |gamma| rather than per branch.

### Caveat on the "exact DB enumeration"
`K(x->y)/K(y->x) = exp(-dS)` is an ALGEBRAIC IDENTITY of Metropolis-Hastings
whenever forward path i and reverse path i use the same (q, q') pair. What
that test really established is (a) the forward/reverse path multiplicities
are a genuine bijection and (b) the contract branch's ingredients computed on
y agree with the split branch's on x. It does NOT verify that q is the true
proposal probability -- that is what the rate audits do, and they are
AGGREGATE per branch.

## Per-|gamma| audit: CLEAN (2026-08-13) -- that lead is dead

df3 identifies both direction and ring length exactly (split df3 = +L,
contract df3 = -L), so the branch audit can be resolved per L. 1.5M
single-move trials, maxRing 6, measured/predicted:

    state drawn from a BISTELLAR-ONLY chain (f0 28, f3 128)
      L      split               contract
      3   0.9914 +- 0.0045    0.9989 +- 0.0128
      4   1.0036 +- 0.0036    1.0382 +- 0.0236
      5   0.9980 +- 0.0028    0.9760 +- 0.0222
      6   0.9942 +- 0.0036    0.9473 +- 0.0260

    state drawn from the MIXED CIRCULATING chain (f0 29, f3 131)
      3   1.0068 +- 0.0049    1.0155 +- 0.0154
      4   1.0058 +- 0.0044    1.0056 +- 0.0076
      5   0.9965 +- 0.0061    0.9927 +- 0.0087
      6   1.0101 +- 0.0091    1.0040 +- 0.0087

No compensating error by ring length, and no state-class effect. The
bistellar+hinge audit re-run on a mixed-chain state is also clean (all <=1%:
1->4 1.0009, 2->3 1.0010, 3->2 0.9999, 4->1 1.0087, hinge 0.9960). The caps
skipped 0 vertices and 0 edges at these states.

**Both kernels are locally correct on the states the circulating chain
actually visits.** Local rate correctness is NECESSARY but NOT SUFFICIENT for
reversibility: it constrains the total rate out of x, not the DISTRIBUTION
over targets. Per-L is a coarse per-target test and it passes; a full
per-target test is statistically infeasible (per-target probability ~3.6e-5
against a 14% effect).

### The one feasible refinement left
Bin by a RICHER measurable key than L -- e.g. (df3, change in the edge-degree
histogram, change in the vertex-degree histogram). Computable from the
resulting facets for the measurement and from each enumerated move for the
prediction, and coarse enough to keep per-bin statistics. That is the only
remaining way to catch redistribution among targets.

## `run(exact=True)` is NOT exact -- separate bug

`ddg_sampler_run_exact` is its own loop: it excludes the contract/split
channel, hinge moves and the potential, ignores `set_bistellar_hastings`, and
uses `log(V_before) - log(V_after)` -- which is the right Hastings only for a
proposal uniform over valid moves. `chooseRandomMove` gives 1->4 HALF weight
(the min(1, 2/d) clip at d = 1), so run_exact carries the same ln-2
anti-vertex bias it is supposed to remove. Correct weights: q(m) = c(m)/C,
C = 2V - f3.

## Both fixes APPLIED and verified (2026-08-13)

### run_exact -- REAL BUG, FIXED, verified
`ddg_sampler_run_exact` used `log(V_before) - log(V_after)`, right only for a
proposal uniform over valid moves. `chooseRandomMove` gives the 1->(dim+1)
move HALF weight (min(1,2/d) clips at d=1), so q(m) = c(m)/C with
C = 2V - f_dim, and the factor is `(c_rev/c_fwd) * (C_before/C_after)`.
Now implemented; `set_bistellar_hastings(False)` reproduces the old formula.
Verified against the independently-corrected mcmcStep, meter <mu> = dS/df0:

    run_exact OLD    <mu> -1.0127 +- 0.0234
    run_exact FIXED  <mu> -0.6951 +- 0.0228
    mcmcStep ref     <mu> -0.7163 +- 0.0108
    FIXED - ref = +0.0212 +- 0.0252  (0.8 sigma -- agreement)
    OLD   - ref = -0.2964 +- 0.0258  (11.5 sigma -- the bug)

**Every past `run(exact=True)` chain carried this ln-2 anti-vertex bias.**

### STAR_CAP -- real asymmetry, closed, but PROVABLY UNREACHABLE (no effect)
Split gated on |star(w)| < 96 while the reverse contraction needs
|star(w)| + |gamma| < 96; now gated on nT + clen so both directions see the
identical condition. But it changes nothing measurable, and there is a proof:
a vertex link is a triangulated 2-sphere, so **F = 2V - 4** (confirmed
numerically: |star| max 278 vs |lk| max 141). Hence
`maxLinkVerts = 32` <=> |star| <= 60, while the STAR_CAP window is [90, 96) --
**the symmetric linkOverflow gate always fires first, so the asymmetric gate
was never reachable.** The earlier "7.94% of vertices in the one-way window"
was misleading: those vertices were already excluded by linkOverflow.

Measured after the fix (unchanged, as predicted):
    melted R, potential off : +9.27% +- 0.12% (was +9.59% +- 0.10%)
    churned sphere ring 3/6 : +1.3648 / +1.0730 per sweep (was 1.34 / 1.08)

So the circulation cause is STILL OPEN. Note also that linkOverflow excludes
11.4% of vertices of a melted state from the channel entirely (symmetric, so
not a DB bug, but a real coverage limit).

## Flux spectrum by ring length (2026-08-14) -- new ENSEMBLE instrument

Added per-ring-length accept counters to ContractSplitConfig, exposed as
`ddg_sampler_contract_split_by_len` / `ManifoldSampler.contract_split_by_len()`.
Splits with |gamma| = L and contractions of a ring-L edge are each other's
reverses, so each L is a sub-family closed under reversal and must show zero
net flux separately. Churned sphere, ring cap 6, cs 0.3, 6000 sweeps:

    L   splits  contracts    net    net/(s+c)        share of total net
    3     8531       6176  +2355  +16.01% +- 0.70%        50.3%
    4     7413       6035  +1378  +10.25% +- 0.75%        29.4%
    5     5333       4761   +572   +5.67% +- 0.82%        12.2%
    6     3494       3117   +377   +5.70% +- 0.88%         8.1%
    ALL  24771      20089  +4682  +10.44%

Present at EVERY L, decaying monotonically with L -- not localised to one
sub-family. That is the signature of an external potential tilting pi with the
channel merely responding, not of an L-specific bug.

## Everything now ruled out, quantitatively

* **The tilt is not ln(C/f3).** Measured along the CHANNEL's own splits (not
  just 1->4): **+0.0013 / +0.0006 / +0.0022 / +0.0059** for L = 3/4/5/6. A flux
  fraction phi needs delta = 2*atanh(phi), i.e. 0.323 / 0.206 / 0.113 / 0.114.
  **Short by 20x to 250x.**
* **The proposal loop is a proper rejection loop** -- `continue` only on the
  degree gate / keep-coin / frozen / invalid, and `return false` on Metropolis
  rejection (checked in source). So pi_B ∝ exp(-S)*C/f3 stands; it is NOT the
  rejection-free variant, which would have given pi ∝ exp(-S)*Z(x) and a far
  bigger tilt.
* **The path bijection holds at ring 6 too**: 21668 distinct split targets,
  (fwd, rev) multiplicities exactly (1,1) 19950, (2,2) 1322, (3,3) 264,
  (4,4) 132.
* **C has no missing category**: all five proposal-weight ratios came out at
  1.000 with no common offset, and a missed landable class would have scaled
  all five by the same factor.

### The gap, stated precisely
For a mixture of two reversible kernels, flux/throughput ~ epsilon, the
log-target difference per move. Measured flux/throughput = 0.104 => epsilon
~ 0.21 nats. Every measurement of the actual target difference says
epsilon <= 0.006. **A factor ~35 that nothing so far explains.**

### The decisive next experiment
Stop probing and CONSTRUCT: on a small manifold (f3 pinned to ~12-20), BFS the
reachable labeled state space, build the full mixed transition matrix by
exhaustively enumerating both kernels, and solve for the stationary vector.
That either shows pi = exp(-S) (so the fault is scale- or state-class
dependent) or hands over the offending (x, y) pair outright. No statistics,
no theory -- it terminates the search either way.

## SOLVED (2026-08-14): it is the LABEL BOOKKEEPING, not the geometry

Built the mixed chain EXACTLY -- both kernels from the validated replicas, on
the labelled state space (facets + the unusedVertices stack), f3 window
[5,10], 32 states -- and solved for the stationary vector. It reproduces the
phenomenon: channel flux **-19.19% per accepted channel move**. Then the
offending pair fell straight out:

    K(31->14) = 8.571e-03      K(14->31) = 0
    state 31: f3=5 f0=5 pool=(0,)  live=[1,2,3,4,5]
    state 14: f3=8 f0=6 pool=()    live=[0,1,2,3,4,5]

State 31 has HOLEY labels: 0 is retired into the pool. A split takes
`fresh = pool.back = 0`, minting a vertex labelled 0. The reverse would be
contracting (0, w) -- but `tryContractSplit` keeps `min(u,v)` = 0 and retires
w instead. **A split that creates a label smaller than its split vertex has no
reverse.** One-way move, f0 pumped, done.

Two components, both label discipline:
 (a) the contraction keeps min(u,v), so only ONE of the two survivor choices
     is ever proposed while being weighted as the whole move -- worth exactly
     a factor 2, which is the long-suspected ln 2;
 (b) `fresh` comes from a LIFO stack and a split from an EMPTY pool allocates
     a brand-new label leaving the pool empty, while its reverse contraction
     always PUSHES -- so the pool state cannot return either.

Fix, verified in exact arithmetic (total |pi K - pi K| over all pairs):

    FIX=0  current code                   chan 8.872e-03   bist 5.866e-03
    FIX=1  + survivor coin (a)            chan 3.632e-03   bist 5.866e-03
    FIX=2  + canonical fresh label (b)    chan 4.879e-18   bist 1.301e-17

i.e. **both families become exactly reversible for their derived targets**:
the channel for exp(-S) and the bistellar+hinge sector for exp(-S)*C/f3. So
both derivations were right all along; the label discipline was the bug. Note
(b) fixes the BISTELLAR family too -- 1<->4 shares the same pool.

The fix in the sampler:
  * `fresh` = a canonical choice (smallest free label), not a stack pop;
  * the contraction picks its survivor by a fair coin, both endpoints allowed;
  * hence `-log 2` on the contract branch's q_fwd and on the split branch's
    q_rev (splits become 2x less favoured relative to contractions).

### Why every earlier test missed it
Every fixed-state rate audit measured the AGGREGATE rate out of a state; this
bug is about WHICH LABELLED STATE you land in, which those cannot see. Worse,
the exact-enumeration replica (`db_exact.py`) used `fresh = max(label) + 1`,
so fresh > w always -- **precisely the benign case**. The original pair test
was designed around it too, its docstring saying labels are "compacted to
0..n-1 so a split creates fresh = n ... and the reverse contraction of (w, n)
keeps w < n, restoring x bit-for-bit". The blind spot was written into the
test suite from the start.

### Why the n_6 potential hid it
The irreversible sub-case needs the pool to hold a label SMALLER than a
typical split vertex. Pool entries come from channel contractions (which
retire max(u,v), a large label) and from bistellar 4->1 (which retires
whatever degree-4 vertex it finds -- often small). With the potential on,
degree-4 vertices are illegal, 4->1 is suppressed, and the pool holds mostly
large labels -- so the bad case rarely fires. Hypothesis, consistent with
every regime measured; worth confirming with a counter.

## THE FIX, validated and shipped (2026-08-14)

Two changes, both label bookkeeping, each validated to MACHINE ZERO by exact
enumeration of a 1404-state chain (f3 window [5,12], where holes coexist):

1. **Survivor coin.** The contraction picks which endpoint survives by a fair
   coin instead of always keeping `min(u, v)`. Hastings: `+log 2` on the
   contract branch's q_fwd and `-log 2` on the split branch's q_rev.
2. **Hole-only pool push** (`retireLabel`). Push the retired label onto
   `unusedVertices` UNLESS the pool is empty and the label equals `fVector[0]`
   after the move -- the allocator would hand that one back anyway, so pushing
   creates a second name for one state. Applies to the channel's contraction
   AND to the bistellar 4->1 (mcmcStep and run_exact both).

Incremental validation, total |pi K - pi K| over all pairs, 1404 states:

    variant                          channel      bistellar
    current code                     8.320e-03    4.952e-03
    survivor coin only               3.077e-03    4.952e-03
    hole-only push only              5.242e-03    5.692e-18   <- fixes bistellar
    BOTH                             5.870e-18    5.692e-18   <- both exact

Note the push fix alone makes the BISTELLAR family exactly reversible (4->1's
retired vertex is geometrically determined, so it needs no coin); the channel
needs both. An earlier attempt using `fresh = smallest free label` reached
zero only on the 32-state window where holes never coexist -- at [5,12] it
left 2.3e-03. The push rule must test `pool empty`, not just the label value:
a NON-empty pool returns its own back, so a retired label must be recorded.

`ddg_sampler_set_label_fix(bool)` / default ON; OFF reproduces old chains.

### MEASURED after the fix
Circulation diagnostic, churned sphere, f3 240, cs 0.3 -- the statistic that
sat at 30 sigma through this whole investigation:

    ring   before (per sweep)   after (per sweep)        z
      3         +1.3648          -0.0146 +- 0.0357     -0.41
      4         +1.4878          +0.0029 +- 0.0417     +0.07
      6         +1.0816          +0.0672 +- 0.0395     +1.70
      8         +0.7203          -0.0658 +- 0.0295     -2.23   -> PASS

Per-ring-length flux, per accepted channel move:

    L      before        after
    3     +16.01%    +0.65% +- 0.72%
    4     +10.25%    -1.19% +- 0.76%
    5      +5.67%    -1.34% +- 0.83%
    6      +5.70%    +0.84% +- 0.97%
    ALL   +10.44%            -0.29%

**This CHANGES SAMPLED RESULTS for every chain that ran the channel**, and the
push fix also changes bistellar-only chains (4->1 pool bookkeeping). Chains in
the production regime (n_6 potential on) already showed zero circulation, so
they should shift little.

## Post-fix confirmations (2026-08-14)

**The flag is a clean reproducibility switch.** Same host/params, churned
sphere f3 240, cs 0.3, ring 6, 3000 sweeps, per accepted channel move:

    set_label_fix(False)   +11.10% +- 0.38%     (reproduces the old bias)
    set_label_fix(True)     -0.17% +- 0.46%     (fixed)

**Production is unaffected**, as expected -- it never had the bias to begin
with. R crystal, n_6 potential zleg=cimp=0.6, cs 0.6, ring 6, 1800 sweeps:

    before the fix   -0.07% +- 0.62%   (FK 0.746, 42 Z16_D2)
    after  the fix   +0.00% +- 0.69%   (FK 0.777, 29 Z16_D2)

So chains run with the potential on need no reinterpretation; chains run
WITHOUT it (bare geometric action, any churned/melted state) carried the
~10-16% current and should be treated as suspect.

`ManifoldSampler.set_label_fix(bool)` now wraps the C entry point.

## Tools

`scripts/validate_contract_split.py` gained `--circulation` (with
`--max-ring` bisection) and `--fugacity` (reads the gap in nats). The
enumeration/audit scripts live in the job tmp dir and should be promoted to
`tools/` if this is picked up again. See [[bilocal-program-saga]],
[[volume-pin-defects]].
