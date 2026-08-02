"""Is n_deg4 a usable observable, or is it just a flicker by-product?

On pristine crystal a 2->3 on a (5,6,6) face drops all three face edges
by one -> (4,5,5). So EVERY live flicker manufactures deg-4 edges for
free. If so, counting deg-4 edges measures the transient flicker
population, not the formation of stable deg-4 species, and the A/B test
needs a flicker-excluded observable.

Measures, on the pristine crystal: the degree spectrum, then the same
after placing one 2->3, and the deg-4 count split by whether the edge
touches a live deg-3 edge.
"""
import os, sys
from collections import Counter
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg

START = os.environ["START"]
m = ddg.Manifold.load(START, 3)
s, L = F.fresh(m)


def spectrum():
    return Counter(L.edeg.values())


def d3():
    return [e for e, d in L.edeg.items() if d == 3]


def d4():
    return [e for e, d in L.edeg.items() if d == 4]


def support(e):
    """vertices of the tets on edge e"""
    ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
    return {x for q in ts for x in q}


def split_d4():
    """deg-4 edges that DO / DON'T touch a live deg-3 edge"""
    d3v = set()
    for e in d3():
        d3v |= support(e)
    near, far = 0, 0
    for e in d4():
        if set(e) & d3v:
            near += 1
        else:
            far += 1
    return near, far


print(f"start: {os.path.basename(START)}")
print(f"f = {tuple(int(x) for x in s.manifold.f_vector)}")
sp = spectrum()
print(f"pristine degree spectrum: {dict(sorted(sp.items()))}")
print(f"  n_deg3 = {len(d3())}   n_deg4 = {len(d4())}")

# place one 2->3 on a legal (5,6,6) face and re-measure
placed = None
for t in list(L.v2t[sorted(L.v2t)[0]])[:1]:
    pass
cnt = 0
for v in sorted(L.v2t):
    for t in list(L.v2t[v]):
        for i in range(4):
            f = tuple(sorted(x for j, x in enumerate(t) if j != i))
            ts = [q for q in L.v2t[f[0]] if set(f) <= set(q)]
            if len(ts) != 2:
                continue
            ap = tuple(sorted({x for q in ts for x in q} - set(f)))
            if len(ap) != 2 or ap in L.edeg:
                continue
            sig = tuple(sorted(L.edeg[(min(a, b), max(a, b))]
                               for a, b in ((f[0], f[1]), (f[0], f[2]),
                                            (f[1], f[2]))))
            S0 = s.current_objective
            L.do(f, ap)
            placed = (f, ap, sig, s.current_objective - S0)
            break
        if placed:
            break
    if placed:
        break

f_, ap_, sig_, dS_ = placed
print(f"\nplaced one 2->3 on face {f_} sig={sig_}  apexes={ap_}  "
      f"dS={dS_:+.3f}")
sp2 = spectrum()
print(f"after one flicker:        {dict(sorted(sp2.items()))}")
print(f"  n_deg3 = {len(d3())}   n_deg4 = {len(d4())}")
near, far = split_d4()
print(f"  deg-4 touching a live deg-3 edge: {near}")
print(f"  deg-4 NOT touching one:           {far}")
print("\ndelta vs pristine: "
      f"{ {k: sp2[k] - sp[k] for k in set(sp) | set(sp2) if sp2[k] != sp[k]} }")

L.do(ap_, f_)
print("\nstate restored")
