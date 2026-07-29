"""Verify + characterize the accidental 2-move 'fusion' on the D core:
  M1: 2->3 on face (7874,7997,8021), apexes (7841,7881)
  M2: 3->2 on edge (7874,7997), link (7841,7865,7881)
Before/after: illegal edges near the site, deg-4 fragment structure, complex
membership (DefectState components), f-vector, and the local degree ledger."""
import os, sys, itertools
from collections import Counter, defaultdict
import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import discrete_differential_geometry as ddg
import defect_state as dsm

SNAP = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")
SITE = {7841, 8021, 7874, 7997, 7881, 7865}

def snapshot_info(m, tag):
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    d4 = [e for e, d in ill.items() if d == 4]
    adjv = defaultdict(list)
    for a, b in d4:
        adjv[a].append(b); adjv[b].append(a)
    # fragments touching the site
    near = [e for e in ill if set(e) & SITE or
            any(v in SITE for x in e for v in adjv.get(x, []))]
    print(f"\n[{tag}] illegal total {len(ill)}  "
          f"{dict(Counter(ill.values()))}")
    print(f"  illegal edges at/adjacent to the site:")
    for e in sorted(near):
        va = len(adjv.get(e[0], [])); vb = len(adjv.get(e[1], []))
        print(f"    {e}: deg {ill[e]}  frag-valence ({va},{vb})")
    st = dsm.DefectState(m)
    comps = [c for c in st.components() if set(c.verts) & SITE]
    for c in comps:
        print(f"  complex touching site: n={len(c.verts)} sig={c.sig} "
              f"nodes={c.nodes}")
    return ill

m = ddg.Manifold.load(SNAP, 3)
f0 = list(m.f_vector)
ill0 = snapshot_info(m, "BEFORE")

assert m.has_bistellar_move([7874, 7997, 8021], [7841, 7881])
m.do_bistellar_move([7874, 7997, 8021], [7841, 7881])
assert m.has_bistellar_move([7874, 7997], [7841, 7865, 7881])
m.do_bistellar_move([7874, 7997], [7841, 7865, 7881])

f1 = list(m.f_vector)
ill1 = snapshot_info(m, "AFTER")
print(f"\nf-vector: {list(f0)} -> {list(f1)}  "
      f"({'PRESERVED' if list(f0) == list(f1) else 'CHANGED'})")
ch = {e: (ill0.get(e), ill1.get(e)) for e in set(ill0) | set(ill1)
      if ill0.get(e) != ill1.get(e)}
print(f"illegal-set delta: {ch}")
