---
name: flicker-catalysis
description: "NO evidence of flicker catalysis — net dS is structurally blind to it, and the strict channel has U≡0 with per-move Metropolis so it cannot lower any barrier; the real mechanism needs BUNDLING"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-08-01T13:35:10.802Z
---

**This memory previously claimed catalysis EXISTS and is common. That
claim does not survive scrutiny and is RETRACTED.** What survives, and
what replaces it:

## What is still true

The original "flicker is catalytically inert" verdict WAS an artifact of
an over-narrow move class. `wf0ChordEnumH` required the support to
contain BOTH chord endpoints, but a chord's degree counts tets holding
both, so that rule guarantees the chord's degree changes and the closure
condition breaks. Measured on harvest_f0_final, depth-1:

    both endpoints:  0 / 125   left the chord annihilable (0.0%)
    >= 1 endpoint: 2186 / 2286 left it annihilable (95.6%)

The widened class is correct and is what the D kernel now uses.

## Why the catalysis claim was wrong

**95.6% is a CLOSABILITY statistic, not a catalysis statistic.** It says
the widened class is compatible with closing. It says nothing about
whether the flicker enabled anything.

**Net dS cannot detect catalysis, structurally.** In
`catalysis_wide.py` the probe flicker is created AND destroyed inside
the path, so it contributes exactly zero to the net. Every bit of the
"best net dS −8.33" comes from the OTHER moves. Catalysis is a
statement about the BARRIER,

    B(P) = max_k (S_k − S_0),

and the search optimized the endpoint instead. It could not have
detected catalysis even if present.

**The representative path is a spectator.** In
`0: (27,42,47)->(26,32); 1: (33,47)->(27,32,54); 2: inverse of 0`,
move 0's support is {27,42,47,26,32} — it contains 47 but NOT 33, so no
tet holding both 33 and 47 is touched, deg(33,47) is unchanged, and
move 1 was already legal in T0.

## The decisive structural fact

In `wormChordStrictEpisode` the umbrella is **U ≡ 0** (dropped when
discrimination moved into the closure criterion) and every walk move is
individually Metropolis:

    head:   lah = -dh + log(nH0) - log(nH)      // Hastings counts only
    global: la  = -(dB + dP)                     // plain Metropolis

Sequential Metropolis moves **cannot transfer action between them**.
Annihilating the flicker releases its dS at its own step; the next
uphill move still pays e^-X in full. So the channel cannot lower a
barrier, and no measurement will find catalysis in it.

What the channel DOES give is real but different: ~40 localized attempts
per episode at two marked sites (vs ~1 from an ordinary sweep), plus a
bilocal relocation the local sampler would take enormous time to find.
Accelerated, balanced sampling — not catalysis.

## The mechanism that would actually work: BUNDLING

Price the annihilation and the barrier in ONE acceptance:

    sequential: e^-Q * e^-X, each accepted alone  ->  e^-X
    composite:  e^-(X-Q)

so the flicker's Q genuinely pays for the barrier. **MEASURED
2026-08-01: the flicker annihilation quantum is dS = -5.930 with sd
EXACTLY 0.000 over all 17 flickers in quench_down5q_wOFF** — perfectly
quantized, and independently equal to the design doc's self-mirror
sector value (d0 = 3 edges, 6/6, "dS exactly antisymmetric +-5.93").
Creation on random pristine faces costs +14.20 +- 6.18 (min +6.19),
i.e. the surviving population sits at the cheapest site class. (Do NOT
confuse this with the rung ladder Q in {46,48,50,52} of
[[bc-washboard-not-free-spirals]] — that is a CHARGE ladder with
dS = c*dQ, not an action.) Both halves are self-mirror (2->3
inverted by 3->2), so the composite's inverse is the same kernel
mirrored: **Tier 1** in notes/bilocal-worm-design.md §2.5 (the 6/6
region, knot slide as existence proof), NOT the Tier-2 / scheme-B region
already dead at e^-70 (see [[pair-carrier-calibration]]). Open design
question: the rule pairing the B-move with the A-annihilation must be
computable from either endpoint so the frame probability cancels.

## Note on screening

A "support-disjointness screen" is meaningful ONLY for
`catalysis_wide.py`'s internal probe flicker. It does NOT apply to the
strict channel, where the two marks are REQUIRED to sit at different
locations — that is the move's purpose, not a defect. See
[[bilocal-factorization]], [[strict-chord-channel]].

GOTCHA that stands: `m.illegal_edges()` returns a TUPLE (pairs, degs) —
`len()` of it is 2, the arity; use `len(pairs)`.
