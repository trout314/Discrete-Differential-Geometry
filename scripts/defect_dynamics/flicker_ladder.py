"""Where does the -5.930-per-flicker ladder break?

Greedily add flickers around a collapse target, one at a time, each
chosen to minimise the total barrier
    total_k = sum_i (dS_place_i - 5.930) + max_j (cumulative dS of the
                                                 re-planned corridor)
and record the marginal gain total_{k-1} - total_k.

PREDICTION: the ladder holds at one quantum per flicker while the
flickers' supports stay VERTEX-DISJOINT -- that is exactly the
condition under which bilocal factorization was measured exact to
1e-13 -- and must bend once they overlap, because the contributions
stop being additive. So the disjointness column should flip at the
same k where the marginal gain departs from 5.930.

Support of a 2->3 on face (a,b,c) with apexes (x,y) is {a,b,c,x,y}.
"""
import os, sys
import numpy as np
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg
import link_planner as LP

HOME = 5.930
KMAX = int(os.environ.get("KMAX", "6"))
NF = int(os.environ.get("NF", "24"))
NODECAP = int(os.environ.get("NODECAP", "200000"))

m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
S_ref = s.current_objective


def replay(v, ops):
    applied, cum = [], []
    S0 = s.current_objective
    for op in ops:
        if op[0] == "delete":
            u = op[1]
            tt = [t for t in L.v2t[v] if u in t]
            lk = sorted({x for t in tt for x in t} - {v, u})
            cen, coc = F._e(v, u), tuple(lk)
        else:
            _, ab, xy, flavor = op
            cen, coc = ((tuple(sorted((v,) + ab)), xy) if flavor == "23"
                        else (ab, tuple(sorted((v,) + xy))))
        L.do(cen, coc)
        applied.append((cen, coc))
        cum.append(s.current_objective - S0)
    return applied, cum


def undo(applied):
    for cen, coc in reversed(applied):
        L.do(coc, cen)


def plan(v):
    fvn = [int(x) for x in s.manifold.f_vector]
    ops, c = LP.plan_collapse(L, v, F.CIMP, 1.0, F.ETARGET, fvn[3],
                              fvn[1], F.F3T, nodecap=NODECAP,
                              optimize="cost")
    return ops


def apexes(f):
    fs = set(f)
    ts = [t for t in L.v2t[f[0]] if fs <= set(t)]
    if len(ts) != 2:
        return None
    ap = tuple(sorted({x for t in ts for x in t} - fs))
    return ap if len(ap) == 2 else None


def candidates(v):
    out, seen = [], set()
    for t in list(L.v2t[v]):
        for i in range(4):
            f = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if f in seen:
                continue
            seen.add(f)
            a = apexes(f)
            if a is None or a in L.edeg:
                continue
            out.append((f, a))
    return out


v = int(os.environ.get("V", "0")) or sorted(
    {x for e in L.edeg for x in e},
    key=lambda u: (len(F.nb_of(L, u)), u))[0]
ops = plan(v)
ap0, cum0 = replay(v, ops)
bar0 = max(cum0)
undo(ap0)
print(f"target v={v}  z={len(F.nb_of(L, v))}  "
      f"unassisted barrier {bar0:.2f}  endpoint {cum0[-1]:.2f}")
print("B_tot = barrier from the ORIGINAL state (sets how fast z=4 is "
      "reached)")
print("E_tot = z=4 state cost from the ORIGINAL state (sets how OFTEN "
      "it is there)\n")
ops0, cum0b = None, None
print(f"{'k':>2} {'reloc_tot':>10} {'barrier':>9} {'endpoint':>9} "
      f"{'B_tot':>8} {'E_tot':>8} {'gain':>7}  site")

placed, supports, reloc_tot, prev = [], [], 0.0, bar0
for k in range(1, KMAX + 1):
    best = None
    for f, a in candidates(v)[:NF]:
        S0 = s.current_objective
        try:
            L.do(f, a)
        except Exception:
            continue
        d = s.current_objective - S0
        o2 = plan(v)
        if o2:
            ap2, cum2 = replay(v, o2)
            bar = max(max(cum2) if cum2 else 0.0, 0.0)
            end = cum2[-1] if cum2 else 0.0
            tot = reloc_tot + (d - HOME) + bar
            if best is None or tot < best[0]:
                best = (tot, d - HOME, bar, f, a, end)
            undo(ap2)
        L.do(a, f)
        assert abs(s.current_objective - S0) < 1e-9
    if best is None:
        print(f"{k:>2}   (no placement leaves a plannable target)")
        break
    tot, rel, bar, f, a, end = best
    sup = set(f) | set(a)
    disj = all(not (sup & p) for p in supports)
    print(f"{k:>2} {reloc_tot + rel:>10.2f} {bar:>9.2f} {end:>9.2f} "
          f"{tot:>8.2f} {reloc_tot + rel + end:>8.2f} "
          f"{prev - tot:>7.2f}  {f}+{a}")
    L.do(f, a)                      # keep it for the next level
    placed.append((f, a))
    supports.append(sup)
    reloc_tot += rel
    prev = tot

for f, a in reversed(placed):
    L.do(a, f)
assert abs(s.current_objective - S_ref) < 1e-6
print("\nstate restored exactly")
