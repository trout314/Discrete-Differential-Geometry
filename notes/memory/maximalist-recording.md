---
name: maximalist-recording
description: ALWAYS record all cheap observables when driving the sampler — use/build the record_run flight-recorder; bare s.run() in research scripts is a smell
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-30T15:28:08.723Z
---

The user's maximalist data-collection philosophy: any script that drives
`ManifoldSampler` must record everything cheap, every run, without
per-script discretion. Recording is UNCONDITIONAL; equilibrium
gating/certification is post-processing on the recorded files — do not
couple the two (that coupling is why ad-hoc runs kept losing data).

**Why:** repeated data losses from hand-rolled run loops — reaction_census
lacked f3(t)/⟨e⟩ series and final .mfd saves (blocked the e*-vs-⟨e⟩ report
and post-hoc certification); killed chains lost their in-memory
timeseries; ad-hoc scan scripts each record a different subset. Also the
2026-07-30 GC leak would have been caught immediately if gc_stats were in
every run log.

**How to apply:**
- Use `discrete_differential_geometry.recording.record_run(s, sweeps=,
  out=, chunk=)` instead of bare `s.run()` in every research script (build
  it if it doesn't exist yet — design agreed 2026-07-30, session
  0122FV6jn8urXWnuVD1gAmTo).
- O(1) tier every chunk: sweep, wall, objective, f0–f3, ⟨e⟩=6f3/f1,
  per-move-type acceptance deltas, channel stats (slide/nonlocal/worm),
  potential total, D-heap used/free (leak sentinel).
- Header once: SamplerParams, couplings (λ vs λ_EDQ distinct), host file
  + hash, seed, git commit, argv, build flavor.
- Census tier (n_ill, degree hist, components) every K chunks; snapshots
  (.mfd + cocycle sidecar) mid + FINAL, always.
- Queued enhancement: capi `ddg_manifold_illegal_edges` O(E) dump to make
  the census tier ~free.
- Related: [[reaction-census-campaign]], [[embedded-gc-rt-init]].
