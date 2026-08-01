"""Within (5,5,6): is there a sharp sub-class at 5.930, and how many
faces are in it? That count, not the signature count, is the N_s the
entropy closure needs."""
import math, os, sys
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
NS = int(os.environ.get("NS", "600"))


def sig_of(f):
    return tuple(sorted(L.edeg.get((min(a, b), max(a, b)), 0)
                        for a, b in ((f[0], f[1]), (f[0], f[2]),
                                     (f[1], f[2]))))


def apexes(f):
    fs = set(f)
    ts = [t for t in L.v2t[f[0]] if fs <= set(t)]
    if len(ts) != 2:
        return None
    ap = tuple(sorted({x for t in ts for x in t} - fs))
    return ap if len(ap) == 2 else None


faces, seen = [], set()
for v, ts in L.v2t.items():
    for t in ts:
        for i in range(4):
            f = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if f in seen:
                continue
            seen.add(f)
            if sig_of(f) == (5, 5, 6):
                faces.append(f)
N556 = len(faces)
step = max(1, N556 // NS)
ds = []
for f in faces[::step][:NS]:
    ap = apexes(f)
    if ap is None or ap in L.edeg:
        continue
    S0 = s.current_objective
    try:
        L.do(f, ap)
    except Exception:
        continue
    ds.append(round(s.current_objective - S0, 3))
    L.do(ap, f)
    assert abs(s.current_objective - S0) < 1e-9
a = np.array(ds)
c = Counter(ds)
print(f"(5,5,6) faces total {N556}; sampled {len(a)}")
print(f"dS distribution, most common values:")
for v, n in c.most_common(8):
    print(f"   {v:>8.3f}  x{n:4d}  ({100.0*n/len(a):5.1f}%)")
frac = sum(n for v, n in c.items() if abs(v - 5.930) < 1e-6) / len(a)
print(f"\nfraction at exactly 5.930 = {100*frac:.1f}%")
Nstar = N556 * frac
print(f"=> N* (true site count)   = {Nstar:.0f}")
print(f"   predicted n3 = N* e^-5.930 = {Nstar * math.exp(-5.930):.1f}"
      f"   (observed 17)")
assert abs(s.current_objective - S_ref) < 1e-6
print("state restored exactly")
