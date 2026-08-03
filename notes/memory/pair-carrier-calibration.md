---
name: pair-carrier-calibration
description: "Bilocal VERTEX pair channel: conserves f0 exactly; role-signed umbrella + calib_zeta2 got roundtrips working, but TRANSPORT is blocked by the single-orbit tube being flat off-tube"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-08-01T12:47:20.192Z
---

**Both bilocal carriers conserve f₀ exactly.** The vertex pair episode
opens with one 1→4 (create at ball A) + one flag (adopt at ball B) and
closes with one 4→1 + a drop, so exactly one vertex is created and one
destroyed — `closedHow 4` (transport) and `6` (roundtrip) are both
f-neutral. Measured net df = (0,0,0,0). The chord channel conserves f₀
structurally. So **neither can move the f₀ read-out**; the only
f₀-changing channel is the single-ball create/destroy. The pure-f₀
composite is E = vertex − 3·chord, Δf = (+1, 0). See [[f0-surgery]].

**Invariance result (kills a whole class of tuning).** With
laOpen = log ζ₂ − bracket and la_close = bracket_rev − log ζ₂,
P(open)·P(close) ≤ e^(−Δ), Δ = bracket − bracket_rev. **No ζ₂ can fix a
bracket mismatch** — only reducing Δ helps. Shifting ζ₂ trades one leg
against the other at constant product.

**Bugs found and fixed (all uncommitted as of 2026-08-01):**

1. `set_worm_pair(zeta2=NaN)` silently killed the channel — auto stores
   ζ₂ = 0 and the pair open takes `log(cfg.zeta2)` = −∞. 0 opens in
   2000 episodes, reported as a clean zero-commit run. Now a loud error;
   auto is chord-only (the pair density has a per-draw log(W/wv) term,
   not just state counts).
2. **Role-signed umbrella.** The tube is a one-way ramp (U +3.85 at
   z=11 → +20.60 at (3,3,3,3)). Applied to both balls it drove both to
   z=4, so the close subtracted 20.60+20.60 where the open added
   20.60+0.85 — a ~19.8 deficit, 100% cap-undone, zero closes. Fix:
   created ball carries −U, adopted +U; the close must price under the
   **reverse** episode's roles (deleted↔created, kept↔adopted), else
   balance breaks by twice the umbrella.
3. `calib_zeta2` modelled the seed weight with the neighbour count;
   the engine uses `_vertexDegrees` = **facets in the star** (mean
   4f₃/f₀ = 22.8, not 2f₁/f₀ = 13.4). `len(L.v2t[v])` is the right
   proxy.
4. `cfg.bcp` was **inert** for the pair — p_close was `1 − ph − pg`
   (=0.10, mean episode 10 steps). The in-code claim that "the config
   validates ph + pg + bcp == 1" is false. Same class as the chord
   p_close bug.
5. Diagnostics aliased `res.umax`, which the pair open already seeds
   with `ball[cre].u + ball[adp].u` (= −19.75) — readouts were floored
   there. `zmin`/`nZ4` tracked `ball[1]`, an arbitrary index that was
   the created ball half the time, making the counter look healthy.

**Working:** at μ = 0, 94/200 roundtrips, median close log α = −0.26,
6 abandoned. The seed bias was the last barrier: μ=3 on the tet-degree
gives e^(−1.5Δd), so the open drew a very probable vertex (log(W/wv)≈4)
while transport deposits at an ordinary site the reverse open would
draw at log(Wc/wc)=24.

**STILL BLOCKED — transport never fires**, and the UFB fallback was
necessary but NOT sufficient. Chain of measurements:

- Coverage. The crystal has only **16 distinct spoke multisets** over
  1533 vertices (913 share one). The single-orbit tube covers
  **1/1533** vertices — the one it was built from. `build_umbrella(16
  seeds)` covers **967/1533 = 63%**, and 48 seeds gives 63.3%, i.e.
  saturated.
- The linear fallback cannot substitute. Residuals against the exact
  staircase: sd 3.29, range [−6.13, +6.33], and the per-STEP residual
  jump reaches **7.52** — an e^−7.5 barrier partway along a corridor
  that needs ~10 consecutive accepted moves. (The code comment already
  said a per-spoke sum can't do this, ±7; now measured on this fit.)
  Wiring UFB in did unstick the ball from a hard 12/12 to occasionally
  reaching 10–11, but never 4.
- Even on the 63%-coverage harvested table the adopted ball stalls
  (zmin median 12, best 11, z=4 visits 0 in 200 episodes).

**Root cause, deeper than coverage:** U is a function of the spoke
MULTISET, but a collapse is a specific PATH. A random walk over head
moves with a staircase potential still has to find that path; being
covered at the start says nothing about the intermediate states it
actually wanders into. The single-ball channel never had to solve this
— its open-insert head is a fresh z=4 vertex already AT the corridor
end, so it only has to not leave. **Confirmed by prodq3's own numbers:
491 round-trips against 3 genuine removals.** Collapsing an arbitrary
existing vertex has never worked in this machinery.

**DO NOT propose "put the planner in D as a guided proposal" — that is
scheme B and it is already REFUTED.** `scripts/defect_dynamics/
f0_channel.py` (commit `63a6a14`, "CBMC removal/insertion — balance
validated, measured negative") is exactly that design, built and
balance-validated in Python. Verdict from its commit message: shaping
telescopes, so acceptance pays exp(2(η·dZ + γ·dΦ)) at the boundary —
shaping strong enough to complete collapse walks (η≈2.5, γ≈3) costs
**~e^−70 in α**, while shaping inside the acceptance slack (~+25)
completes **0/15** walks. completion × acceptance ≈ 0 in every corner.
Scheme C (the umbrella worm, `f0_worm.py`) was the documented response
to that failure — and the measurements above are scheme C failing on
the same sub-problem from the other side.

**The line is the two-tier synthesis in notes/bilocal-worm-design.md
§2.5.** Tier 1 = self-mirror alphabets (flicker create/absorb, 2-move
composites, deg-3 edge surgery), measured **6/6** with exact dS
antisymmetry, "the PRIMARY bilocal transport implementation", knot
slide as existence proof. Tier 2 = deep/dressed halves, measured
**0/24**, "no frame fixes that". The **chord carrier is Tier 1** —
which is why it works and certifies. **Vertex transport is Tier 2**,
because it needs a deep collapse of an arbitrary crystal vertex; the
pair channel's roundtrip closure works only because it falls back to
Tier 1 (create then immediately absorb).

Realistic options: accept the vertex carrier as roundtrip-only (real
fixed-f corridor transport, just not vertex relocation), or look for a
self-mirror vertex FRAME in the spirit of the deg-3 surgery result.
`harvest_f0.py` is force-mode (no Hastings) and stays valid for
PREPARING states (14/14 removals, gap +10.12 → +0.49 in 20 s), never
for measuring. See [[bilocal-factorization]], [[strict-chord-channel]].
