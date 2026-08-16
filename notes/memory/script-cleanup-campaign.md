---
name: script-cleanup-campaign
description: "2026-08-15 script-library cleanup — plan in notes/CLEANUP_PLAN.md, manifest scripts/INDEX.md via tools/script_index.py; two-tier lib split + legacy/<program>/ archival + full tests; awaiting Aaron's Phase-0 review"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2ba2630d-5dd2-4f4e-a2ae-d6d1091a81b0
  modified: 2026-08-16T14:52:00.840Z
---

Cleanup campaign for the ~213 research scripts, started 2026-08-15.

**Decisions (Aaron, 2026-08-15):** two-tier library split (mature hubs →
`python/discrete_differential_geometry/`, research-shared →
new `python/ddg_lab/`); closed programs → `legacy/<program>/` dirs with
verdict READMEs; full test depth (tier-1 units + 13 latent validators →
pytest + `--help` smoke sweep over all active CLIs); Claude classifies,
Aaron reviews before anything moves.

**Artifacts:** `notes/CLEANUP_PLAN.md` (the plan, phases 0–4),
`scripts/INDEX.md` (classified manifest), `tools/script_index.py`
(generator; curated program/status map lives in its CLASSIFICATION dict —
edit there and re-run; `--check` flags unclassified files).

**Key survey facts:** import hubs living in scripts/: defect_state (31
importers), fk_skeleton (34), cocycle_check (23 — scripts import a
*validator* for its helpers), dopant_pairs (17), f0_worm (16), worm_helix
(12); 190 files carry the sys.path.insert bootstrap. Census: 23 lib,
59 active, 14 tool, 13 validator, 67 dormant, 37 closed.

**Status:** Phase 0 DONE (2026-08-16): Aaron approved the classification
with one change — geometry+optics merged into one `geometry` program,
`defect_lens_field.py` closed → tool (follow-up: generalize it beyond the
2→3 lens study when next touched; TODO in CLASSIFICATION map). Final
census: 23 lib / 59 active / 15 tool / 13 validator / 67 dormant / 36
closed. Next: Phase 1a — tier-1 promotions + tests (leaf-most first:
link_classes, chain_select, seed_utils, then crystal_grains, fk_skeleton,
cocycle helpers, tcp_reference). Then archive closed → ddg_lab → remaining
tests → docs. Archive before ddg_lab (don't promote code whose only
consumers are closed).

Related: [[crystal-grains-tool]], [[crystal-symmetry-group]],
[[memory-lives-in-repo]] (memories are tracked in the repo — commit this
with the plan).
