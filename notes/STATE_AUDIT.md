# Duplicated sampler state — audit

Every entry below is information that is *derivable from the manifold's
facet set* but is also stored somewhere else and maintained incrementally.
Each such copy is a place where a bug can hide silently, because nothing
crashes when a cache drifts — the run just samples the wrong thing.

Two bugs of exactly this class are already on record: the ridge
bookkeeping bug (unsorted center/coCenter through the C API vs
sorted-merge `productUnion`) and the cocycle label skew. This audit is an
attempt to enumerate the remaining surface before it produces a third.

Status: **first pass, 2026-08-01.** Inventory + risk ranking done; no
remediation yet.

---

## A. The same incidence relation, cached four ways

All four answer some version of "which facets touch this thing", with
four different owners, lifetimes, and staleness policies.

| # | Cache | Location | Lifetime | Staleness policy |
|---|---|---|---|---|
| A1 | `_vertexWitness` — vertex → some facet | `manifold.d:69` | always | **validated on read**, rescans if stale |
| A2 | `RidgeInfo.link` — ridge → its ≤2 apexes | `manifold.d:37` | always | assumed correct; maintained in insert/remove |
| A3 | `Deg3Set` — edge → witness facet, degree-filtered | `sampler.d:2033`, instances `ddg_capi.d:1640,1646` | **conditional** | assumed correct; **no validator** |
| A4 | `Live.v2t` / `Live.edeg` | `scripts/defect_dynamics/worm_deg4_slide.py` | per-script | mirror-on-write via `L.do()` |

Only **A1** checks itself. It is also the newest (added 2026-08-01), which
is backwards — the oldest and least-guarded caches carry the most risk.

**A3 is the highest-risk item in the whole audit:**

- Its lifetime is gated on `nlSlideCfg.prob > 0` or `wormCfg.prob > 0`
  (`ddg_capi.d:1923-1925`) — flags that have nothing to do with whether a
  caller needs the set. Any new consumer must remember to extend that
  condition or silently gets a stale/empty set.
- `add()` stores a `[-1,…]` sentinel when no witness facet is supplied
  (`sampler.d:2060`), so "witness unknown" is a *legal* state every
  consumer has to handle. `reconcileDeg3`'s witness argument only covers
  edges inside the facet passed, so sentinels are normal, not exceptional.
- A missing entry is undetectable from the set itself. If `wf0ChainSites`
  is ever driven from it, a missed degree-3 edge means fewer enumerated
  sites, and since `pick` is uniform over `[0, n)` that **changes the
  proposal density with no error and no crash** — it would quietly
  invalidate the balance certificate.

**A4** diverges the moment anything mutates the manifold without going
through `L.do()`. Native `sampler.run()` does exactly that, so any script
mixing `L.do()` with `s.run()` has a stale `Live` — this is a live
footgun, not a hypothetical.

---

## B. Derived scalars, incrementally updated

| # | Value | Where | Repair path |
|---|---|---|---|
| B1 | `currentObjective` | `ddg_capi.d:22` | `recomputeObjective`, called defensively at ~13 boundaries |
| B2 | `numSimplices` (f-vector), `totSqrDegrees`, `_vertexDegrees` | `manifold.d:105-106,45` | none at runtime; checked only in unittests (`manifold.d:2453`) |
| B3 | `VertexPotState.n6`, `.mImp` | `sampler.d:511-512` | rebuild exists (`sampler.d:~535`) but is **init-only** |
| B4 | `VertexPotState.total` | `sampler.d:513` | none; only ever `+= dS` (`sampler.d:640,715`) |

**B1 is the model the others should copy.** It has a cheap incremental
path *and* a full recompute invoked at every state-changing boundary, so
drift is bounded by one move. Nothing else in the codebase has that.

**B4 is two derivations from ground truth** — `total` is derived from
`n6`/`mImp`, which are derived from edge degrees — and it only ever moves
by accumulation. A single missed update persists for the whole run and
biases every acceptance thereafter. Given `real` accumulation over ~10⁸
moves, this deserves a periodic re-derive-and-compare even if the logic
is correct.

---

## C. Worm-channel caches

| # | Item | Note |
|---|---|---|
| C1 | `WormF0Params.skel` + `tubeF1`/`tubeF3` | **Good pattern**: stores the f-vector it was built against, i.e. an explicit staleness stamp. The only cache here that records its own provenance. |
| C2 | `hingeTries`/`hingeAccepts`/`bistellarTries`/`bistellarAccepts` local copies in the run loop | `ddg_capi.d:1930-1945, 1974-1976`. Copied to locals, written back around callbacks and after the loop. Any early-return path added between those points loses counter updates. Hazard class, not a confirmed bug. |
| C3 | `CocycleState` | Already produced one bug (label skew on the save path). Not re-audited here. |
| C4 | `MoveCounters` | Opt-in; low blast radius (diagnostic only). |

---

## Findings

1. **Four caches of one relation, one validator.** The invariants differ
   per cache and are documented nowhere central. A1's validate-on-read is
   the right pattern and should be the template.
2. **`Deg3Set` can silently corrupt a proposal density.** It is the only
   duplicated state whose drift changes sampled *physics* rather than
   just a diagnostic, and it is the one with no validator. It should not
   acquire new consumers (e.g. `wf0ChainSites`) without a debug-build
   cross-check against a ground-truth sweep.
3. **`validateMaps()` is not what it sounds like.** It checks *HashMap
   internal integrity* (`manifold.d:702` → `hashmap.d:142`), not semantic
   agreement between the maps and the facet set. It is exposed to Python
   as `validate_maps()` and called from exactly one place in the tree
   (`scripts/defect_dynamics/worm_slide.py:922`) — no active driver runs
   it.
4. **`--audit-every` audits the tracker, not the sampler.**
   `reaction_census.py:282` calls `st.audit(v)`, which is the species
   tracker's own consistency check. Nothing periodically audits sampler
   or manifold state during a production run.
5. **The newest cache is the best-guarded.** Guarding has been added
   reactively, per-incident, rather than as a policy.

## Proposed remediation (not yet done)

- **One `auditState()`** that re-derives every cache above from the facet
  set and diffs, returning a description of the first mismatch. Exposed
  to Python and wired into the existing `--audit-every` cadence so
  production runs check themselves. Cost is O(N) per call, which is
  nothing at an audit interval of tens of chunks.
- **Debug-build assertions** on the two caches whose drift changes
  physics (A3, B4): re-derive and assert on every use under
  `version(assert)`, so the debug build is a correctness oracle and the
  release build stays fast. This is the pattern that would have caught
  the ridge bug immediately.
- **Retire A4** (`Live`) in favour of reading the D side, or at minimum
  give it an explicit `assert_fresh()` that scripts call after any native
  `run()`.
- **Adopt C1's provenance stamp** more widely: a cache that records the
  f-vector (or a move counter) it was built against can detect its own
  staleness cheaply, without a full re-derive.
