---
name: cocycle-detachment-bug
description: "RESOLVED as phantom: 'cocycle detachment' was a SAVE-PATH label-canonicalization skew in mobile_gas, not a D-core tracking failure; all snapshot pairs repaired; fix + guard in mobile_gas"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-24T05:38:39.111Z
---

CORRECTED DIAGNOSIS (2026-07-24; supersedes the 2026-07-23 version of this
memory). There was NEVER a D-core cocycle-tracking failure. The 2026-07-23
"detachment" (lam35 resume failing with "cocycle missing edge"; snap14k-20k
pairs mismatching) was a file-format skew: `Manifold.save` renumbers vertex
labels to rank order (canonicalize), but `mobile_gas.py` passed RAW live
sampler labels (holes from 1<->4 churn) straight to `coc.save_cocycle`, which
does not canonicalize. Whenever the live label set had holes at snapshot time,
the saved (.mfd, .cocycle.npz) pair skewed — while the IN-MEMORY state stayed
perfectly consistent (that's why the edge-set guard added to mobile_gas never
fired: it compares in-memory, where nothing was ever wrong). `dope_hold.py`
always canonicalized before saving (line ~377) — mobile_gas just missed it.

Consequences, all resolved:
- ALL affected pairs repair losslessly via `coc.canonicalize_labels` on the
  npz: 24 files repaired in place across scratchpad mgas/ + run5h/ on
  2026-07-24 (133 were already consistent, 0 unrepairable), including the
  original lam35_snap{14,17,20}000 previously declared unusable.
- `mobile_gas.py` now canonicalizes before saving (fix next to the guard; the
  in-memory guard stays as a true-detachment tripwire).
- cents/MSD streams were NEVER affected (computed in-memory, self-consistent).
  The lam35c/lam35r MSD results stand.
- `regen_cocycle.py` (registry-based chart regeneration via crystal_grains
  covering map) was built for the phantom but remains valid + validated — the
  right tool if a chart is ever GENUINELY lost. Note lam35r runs on a
  regenerated chart (fresh gauge, winding not continuous with lam35 history);
  lam35c runs on the original chart lineage.
- `web_transport.py` canonicalizes cocycle labels on load and refuses truly
  mismatched pairs.

Lesson: any new driver persisting cocycles MUST canonicalize
(`coc.canonicalize_labels`) before `save_cocycle` — grep dope_hold for the
pattern. See [[mobile-gas-liquid]], [[cert-campaign-and-edq-ladder]].
