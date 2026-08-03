---
name: tcp-r-c15-defect-state
description: State of the R/C15 TCP crystal-defect campaign at the 2026-07-20 remote-crash resume
metadata: 
  node_type: memory
  type: project
  originSessionId: eaaf2503-7af1-42b2-a8c0-88f714d19b53
  modified: 2026-07-21T03:15:15.820Z
---

Context at 2026-07-20 (resuming after a remote Claude session died; remote
memories did NOT sync — reconstructed from repo). Work = crystal-defect analysis
of R and C15 Frank-Kasper TCP phases in the T³ seed-library program
(see `notes/t3-library-plan.md`, `notes/physics-review-2026-07.md`).

**The crystals** (`tcp_melt.CRYSTALS`, files in gitignored `data/tcp_reference/`):
- **R**: q̄=5.1042, per-cell census Z12:81 Z14:36 Z15:18 Z16:24 (fFK=1). small=r m3 N24462, large=rbig m4 N57984.
- **C15** (MgCu₂): q̄=5.1000, census Z12:16 Z16:8 only (no Z14/Z15). small=c15 m4 N8704, large=c15big m9 N99144.

**Where the campaign stood** (HEAD c2af7ad). All 4 crystals (r/a15/sigma/c15)
calibrated in `reference_campaign.py`: FLAT_VAC_MU={r:1.1, c15:3.0},
MELT_LAM={r:0.4, a15:0.5, sigma:0.5, c15:0.4}. Recent threads:
- `r_solubility.py` — R dopant-solubility isotherms (species Z14/Z15/Z16, µ sweep,
  --edge-target sweep, Z12 pin-scan queue `r_pinscan_queue.sh`); analyze via `r_isotherm.py`.
- `defect_census.py` — as of 2026-07-20 a REGISTRY-ONLY reporting front-end over
  `crystal_grains` (the old r=1 grow/heal `local` method was removed). Defects =
  vertices not in a crystalline grain >= --min-size (default 30); reports Z-class
  breakdown + complex clustering. See [[crystal-grains-tool]].
- Supporting: `dope_hold.py` (segment/lineage dopant holds), `tcp_melt.py` (melting).
  NOTE: `crystal_match.py` (the brittle exact ball-certificate phase matcher) was
  DELETED 2026-07-20 — superseded by crystal_grains; its `best_refs`/`REF_GLOB`
  helpers moved into crystal_grains.py. crystal_grains is now the single
  authoritative crystalline/defect identifier.

**To run any of this locally** first regenerate the needed reference .mfd
(see [[regenerate-tcp-references]]); seed families need re-sampling. Build is fine
except the unittest link ([[build-state]]).
