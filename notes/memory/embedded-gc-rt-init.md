---
name: embedded-gc-rt-init
description: "RESOLVED: dlopen'd D lib never ran rt_init -> GC.collect was a silent no-op; EVERY D allocation under Python leaked permanently. Fixed via ddg_runtime_init at load. Keep D hot paths allocation-free anyway."
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-30T04:32:42.555Z
---

**The bug (found 2026-07-30, killed the first worms-on campaign):** a
dlopen'd D shared library does NOT get `rt_init()` (only D executables do,
from main). Without it the calling thread is never registered with druntime,
so the stop-the-world collector silently bails: `GC.collect()` frees 0 bytes
FOREVER (the capi's `try/catch(Exception){}` wrappers swallowed any error),
and every GC allocation grows the heap permanently. Invisible for years
because all hot paths are allocation-free by design; the worm channel's
~3 MB/enumeration garbage turned it into ~13 GB/chain/hour (8 chains
swamped the 16 GB box: 27 GB swap, load 24).

**Fix (two independent layers):**
1. `ddg_runtime_init()` in ddg_capi.d (Runtime.initialize()), called by
   `_dlang.py` immediately after CDLL load, hard error if it fails.
2. Worm channel made allocation-free anyway: `wormEnumerate`'s per-m1 AA
   degree snapshot -> fixed parallel arrays (`origE[21]/origD[21]`, exact
   capacity bound 10+1+10); `tryWormMove`'s three per-call
   `new WormCand[](512)` -> lazily-initialized static TLS scratch buffers.

**How it was cornered:** RSS/Dheap probe (scratchpad leak_probe.py) with
`ddg.gc_stats()` before/after `ddg.gc_collect()`; standalone-D leaktest ran
the same wormEnumerate with a STABLE heap (proving code clean, embedding
broken); 5x load+del Manifold showed 344 MB of provably-dead objects
surviving collect (proving it global, not worm-specific).

**Verification after fix:** worm probe (worm 1e-3, 40 sweeps, N=58k):
+1.4 MB RSS total (was +1738 MB), Dheap oscillating 40-70 MB; collect
reclaims load+del garbage (+123 MB freed); meson test 0 fail; D-vs-oracle
crossval unchanged (D 472 >= oracle 454, dS multiset inclusion, all 121
anchors).

**Rules going forward:**
- Any new D-side channel that allocates per-proposal WILL show up in
  long campaigns; prefer fixed/static scratch buffers (codebase style).
- `Manifold.memoryDiag()` (added) prints container len/cap for leak hunts;
  a capacity that ratchets under move churn is a bug.
- Python-side leak probe pattern: `ddg.gc_stats()` + `ddg.gc_collect()`;
  "used" that collect can't shrink in a quiescent state = real retention.
- Related gotcha context: [[deg4-worm-design]], [[build-state]].
