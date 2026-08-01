"""The flicker formation spectrum, resolved by SITE CLASS.

A 2->3 on a face drops its three edge degrees by one, so the site class
(the face's edge-degree signature BEFORE the flip) is a state function
readable from either side: an existing flicker's link face sits at
(4,4,5) and its class is that +1 = (5,5,6). That symmetry is what makes
it legal to condition a move class on flicker TYPE.

Per class: how many such faces exist, what a flicker there costs, the
equilibrium occupancy N*exp(-dS), and the stored energy that class
actually contributes to the reservoir, occupancy * dS. Battery capacity
and availability are governed by the SAME exponential, so the question
is whether the energetic classes are populated at all.
"""
import math, os, sys
from collections import Counter, defaultdict
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
KMAX = int(os.environ.get("KMAX", "40"))


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


byclass = defaultdict(list)
seen = set()
for v, ts in L.v2t.items():
    for t in ts:
        for i in range(4):
            f = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if f in seen:
                continue
            seen.add(f)
            byclass[sig_of(f)].append(f)

print(f"{'class':>12} {'N_faces':>8} {'n':>4} {'dS':>9} {'sd':>6} "
      f"{'occupancy':>10} {'stored':>8}")
rows = []
for cls, fl in sorted(byclass.items(), key=lambda kv: -len(kv[1])):
    if len(fl) < 8:
        continue
    ds = []
    for f in fl[:: max(1, len(fl) // KMAX)][:KMAX]:
        ap = apexes(f)
        if ap is None or ap in L.edeg:
            continue
        S0 = s.current_objective
        try:
            L.do(f, ap)
        except Exception:
            continue
        ds.append(s.current_objective - S0)
        L.do(ap, f)
        assert abs(s.current_objective - S0) < 1e-9
    if not ds:
        continue
    a = np.array(ds)
    occ = len(fl) * math.exp(-a.mean())
    rows.append((cls, len(fl), len(a), a.mean(), a.std(), occ,
                 occ * a.mean()))
    print(f"{str(cls):>12} {len(fl):>8} {len(a):>4} {a.mean():>9.3f} "
          f"{a.std():>6.3f} {occ:>10.3f} {occ * a.mean():>8.3f}")

tot = sum(r[5] for r in rows)
print(f"\ntotal predicted flicker occupancy = {tot:.2f}   "
      f"(observed n3 = {sum(1 for e,d in L.edeg.items() if d==3)})")
print(f"total stored energy = {sum(r[6] for r in rows):.2f}")
assert abs(s.current_objective - S_ref) < 1e-6
print("state restored exactly")
