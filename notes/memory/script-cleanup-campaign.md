---
name: script-cleanup-campaign
description: "2026-08-15 script-library cleanup — plan in notes/CLEANUP_PLAN.md, manifest scripts/INDEX.md via tools/script_index.py; two-tier lib split + legacy/<program>/ archival + full tests; awaiting Aaron's Phase-0 review"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2ba2630d-5dd2-4f4e-a2ae-d6d1091a81b0
  modified: 2026-08-16T15:10:33.103Z
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
closed. Phase 1a DONE (2026-08-16): fk_skeleton, link_classes,
tcp_reference (+reference_frac_positions from cocycle_check, +vertex_classes
from dopant_pairs), chain_select, crystal_grains, seed_utils all live in
`python/discrete_differential_geometry/`; old locations are sys.modules-swap
shims or thin CLIs re-exporting the core, so legacy `from X import Y`
sys.path imports keep working. 5 new test files (incl. the cocycle ladder as
pytest and a fast-vs-exact link-oracle crossval); suite 311 fast + 2 slow
(`-m slow`), all green. GOTCHA: `dd/pass2_structure.py` reads sys.argv[1] at
module import (pre-existing). Next: Phase 2 — archive the 36 closed scripts
to legacy/<program>/ (before ddg_lab, so we don't promote code whose only
consumers are closed); then Phase 1b ddg_lab, remaining tests, docs.

Related: [[crystal-grains-tool]], [[crystal-symmetry-group]],
[[memory-lives-in-repo]] (memories are tracked in the repo — commit this
with the plan).
