"""The flicker annihilation quantum: how much action does killing one
degree-3 chord actually release?

This is the financing unit for a bundled (annihilate-here, pay-there)
composite, so it has to be a DIRECT measurement, not inferred from an
accepted-dS tail. For every degree-3 edge in the state: apply the 3->2
that annihilates it, record dS, undo exactly.

Split by ISOLATION, because a bare flicker on pristine crystal (the
elementary (3,4,4) knot, one 2->3) and a degree-3 edge embedded in a
larger complex are different objects. Isolated == no OTHER degree-3
edge anywhere in the union of the two endpoint stars, which is exactly
the strict channel's region-clean test.

Also measures the CREATION side on pristine faces, which should mirror
the isolated annihilation if the quantum is well defined.
"""
import os, sys
from collections import Counter
import numpy as np
_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg

m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
S_ref = s.current_objective
fv = [int(x) for x in s.manifold.f_vector]
print(f"state f = {tuple(fv)}   S = {S_ref:.3f}")


def region_edges(e):
    """every edge in the union of the two endpoint stars"""
    tets = set()
    for x in e:
        for t in L.v2t[x]:
            tets.add(tuple(sorted(t)))
    return {(t[i], t[j]) for t in tets
            for i in range(4) for j in range(i + 1, 4)}


def link_of(e):
    ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
    lk = tuple(sorted({x for q in ts for x in q} - set(e)))
    return lk if len(lk) == 3 else None


deg3 = sorted(e for e, d in L.edeg.items() if d == 3)
print(f"degree-3 edges present: {len(deg3)}")

iso, emb, skipped = [], [], 0
for e in deg3:
    lk = link_of(e)
    if lk is None:
        skipped += 1
        continue
    others = sum(1 for q in region_edges(e)
                 if q != e and L.edeg.get(q, 0) == 3)
    S0 = s.current_objective
    try:
        L.do(e, lk)
    except Exception:
        skipped += 1
        continue
    dS = s.current_objective - S0
    L.do(lk, e)
    assert abs(s.current_objective - S0) < 1e-9, "annihilation not undone"
    (iso if others == 0 else emb).append(dS)

for nm, v in (("ISOLATED (region clean)", iso), ("EMBEDDED in a complex", emb)):
    a = np.array(v) if v else np.zeros(0)
    if not len(a):
        print(f"{nm:24s}  n=0")
        continue
    print(f"{nm:24s}  n={len(a):4d}   dS(3->2) = {a.mean():+7.3f} "
          f"+- {a.std():.3f}   range [{a.min():+7.3f}, {a.max():+7.3f}]")
    print(f"{'':24s}  => RELEASE -dS = {-a.mean():+7.3f}   "
          f"quantized values: "
          f"{sorted(Counter(np.round(a, 3)).items())[:4]}")
print(f"skipped (no valid 3->2): {skipped}")

# creation side: 2->3 on faces whose edges are all pristine (deg >= 5)
def apexes(face):
    fs = set(face)
    ts = [t for t in L.v2t[face[0]] if fs <= set(t)]
    if len(ts) != 2:
        return None
    ap = tuple(sorted({x for t in ts for x in t} - fs))
    return ap if len(ap) == 2 else None


seen, cre = set(), []
for v, ts in L.v2t.items():
    if len(cre) >= 400:
        break
    for t in ts:
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if face in seen:
                continue
            seen.add(face)
            if any(L.edeg.get((min(face[a], face[b]), max(face[a], face[b])),
                              0) < 5 for a in range(3)
                   for b in range(a + 1, 3)):
                continue
            ap = apexes(face)
            if ap is None or ap in L.edeg:
                continue
            S0 = s.current_objective
            try:
                L.do(face, ap)
            except Exception:
                continue
            cre.append(s.current_objective - S0)
            L.do(ap, face)
            assert abs(s.current_objective - S0) < 1e-9
            break
a = np.array(cre) if cre else np.zeros(0)
if len(a):
    print(f"\nCREATION on pristine faces  n={len(a):4d}   "
          f"dS(2->3) = {a.mean():+7.3f} +- {a.std():.3f}   "
          f"range [{a.min():+7.3f}, {a.max():+7.3f}]")
assert abs(s.current_objective - S_ref) < 1e-6, "global drift"
print("state restored exactly")
