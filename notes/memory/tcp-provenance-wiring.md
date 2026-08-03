---
name: tcp-provenance-wiring
description: TCP crystal lineages wired into the provenance/history machinery (2026-07-20)
metadata: 
  node_type: memory
  type: project
  originSessionId: eaaf2503-7af1-42b2-a8c0-88f714d19b53
---

2026-07-20: wired the TCP-crystal side into the flattened-provenance history
machinery (`tools/seed_utils.py`), which previously only knew sphere-rooted S³
lineages. Motivation: `tcp_melt.py` (and the whole TCP-experiment family) saved
plain historyless `.mfd` via `v.save()`, silently dropping lineage.

Changes:
- `seed_utils.py`: new `is_root_leg_from(f)` — a lineage ROOT is `"sphere"` OR
  `"crystal:<phase>"`. `verify_history`/`history_fields` now use it (were
  hard-coded to `"sphere"`); `history_root` header reports the actual root token.
  Backward compatible (sphere paths unchanged).
- `tcp_reference.py`: stamps a `from_="crystal:<name>@wyckoff"` build root leg
  (op="build", sweeps=0) into every reference `.mfd`. Also fixed a pre-existing
  quirk where `git_commit` printed the whole (commit,dirty) tuple.
- `tcp_melt.py` `run_point`: reads source root via `read_history`, appends a
  `"melt"` leg via `build_metadata_comments` (records tried/accepted + full
  objective incl. n6 potential). Falls back to a synthesized crystal root if the
  source predates stamping.

**Still NOT wired (follow-up):** `dope_hold`, `quench`, `fk_anneal`, `anneal_vdv`,
`w_liquid_ladder` all still use plain `.save()`. Same pattern applies if lineage
matters for their outputs. See [[tcp-r-c15-defect-state]].

Test artifacts built: small references `T3_C15_m3_N3672`, `T3_R_m2_N7248` (rooted),
and partial melts under `data/melt_test/` (c15_lam*/r_lam*). fFK~0.6 melts have a
PERCOLATING defect set (crystal is minority residue); lighter lam (0.03) leaves
crystal the majority — pick level per the disconnected-domain test.
