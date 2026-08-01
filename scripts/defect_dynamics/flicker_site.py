"""Is the -5.930 flicker sitting on PRISTINE background?

For each existing degree-3 edge: annihilate it, then look at the site it
leaves behind -- the link face lk and the degrees of lk's three edges,
plus whether any degree-4 edge survives in the region. If the vacated
site is pristine (all 5s and 6s, no deg-4), then -5.930 IS the
pristine-crystal flicker quantum and its mirror is +5.930 at that face.
If the site carries deg-4 edges, the flicker lives on a defected
background and 'pristine flicker' is a different number.
"""
import os, sys
from collections import Counter
_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg

m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
S_ref = s.current_objective


def link_of(e):
    ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
    lk = tuple(sorted({x for q in ts for x in q} - set(e)))
    return lk if len(lk) == 3 else None


def region_degs(vs):
    tets = set()
    for x in vs:
        for t in L.v2t[x]:
            tets.add(tuple(sorted(t)))
    ed = {(t[i], t[j]) for t in tets
          for i in range(4) for j in range(i + 1, 4)}
    return Counter(L.edeg.get(q, 0) for q in ed)


face_sig, site_sig, recreate = Counter(), Counter(), Counter()
deg3 = sorted(e for e, d in L.edeg.items() if d == 3)
for e in deg3:
    lk = link_of(e)
    if lk is None:
        continue
    S0 = s.current_objective
    L.do(e, lk)                      # annihilate
    # the vacated site: the link face's own three edges
    fe = tuple(sorted(L.edeg.get((min(a, b), max(a, b)), 0)
                      for a, b in ((lk[0], lk[1]), (lk[0], lk[2]),
                                   (lk[1], lk[2]))))
    face_sig[fe] += 1
    rd = region_degs(lk)
    site_sig[tuple(sorted(rd.items()))] += 1
    S1 = s.current_objective
    L.do(lk, e)                      # recreate -- exact mirror
    recreate[round(s.current_objective - S1, 3)] += 1
    assert abs(s.current_objective - S0) < 1e-9

print(f"{len(deg3)} flickers")
print(f"\nlink-face edge degrees after annihilation (the site's own face):")
for k, n in face_sig.most_common():
    print(f"   {k}  x{n}")
print(f"\nrecreation dS at the SAME site (exact mirror): "
      f"{dict(recreate)}")
print(f"\nregion degree spectra around the vacated site:")
for k, n in site_sig.most_common(4):
    print(f"   x{n}: " + " ".join(f"d{d}:{c}" for d, c in k if d))
assert abs(s.current_objective - S_ref) < 1e-6
print("\nstate restored exactly")
