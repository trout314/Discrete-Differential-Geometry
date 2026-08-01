"""Does a flicker placed next to a target vertex LOWER its collapse
barrier? The direct test of the relocation/catalysis picture.

For each target vertex v:
  UNASSISTED -- plan_collapse(v), replay, record the running cumulative
                dS. barrier = max_k cum, cost = final cum.
  ASSISTED   -- place a flicker on a face near v (2->3), re-plan, replay,
                and record the barrier measured from the ORIGINAL state.

Charging the placement: relocating an EXISTING flicker is one bilocal
move, so it costs dS(create at the new site) - 5.930 (the home-site
annihilation gain) and the intermediate is never visited. That net is
what enters the barrier, which is the whole point of doing it
bilocally -- a diffusive walk would have to cross every state on the
way, a bilocal move pays only the endpoints.

Catalysis iff  reloc + barrier_assisted  <  barrier_unassisted.
"""
import os, sys
import numpy as np
_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg
import link_planner as LP

HOME = 5.930                      # measured home-site annihilation gain
NV = int(os.environ.get("NV", "6"))
NF = int(os.environ.get("NF", "24"))
NODECAP = int(os.environ.get("NODECAP", "200000"))

m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
S_ref = s.current_objective
fv = [int(x) for x in s.manifold.f_vector]


def replay(v, ops):
    """apply a plan step by step, returning (applied, cum profile)"""
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


def sig(f):
    return tuple(sorted(L.edeg.get((min(a, b), max(a, b)), 0)
                        for a, b in ((f[0], f[1]), (f[0], f[2]),
                                     (f[1], f[2]))))


verts = sorted({x for e in L.edeg for x in e},
               key=lambda u: (len(F.nb_of(L, u)), u))[:NV]
print(f"f = {tuple(fv)}   targets: {verts}\n")
print(f"{'v':>5} {'z':>3} {'unassisted':>22} {'best assisted':>34}")
print(f"{'':>5} {'':>3} {'barrier':>10} {'cost':>11} "
      f"{'reloc':>8} {'barr_a':>8} {'total':>8} {'site':>10}")

for v in verts:
    z = len(F.nb_of(L, v))
    ops = plan(v)
    if not ops:
        print(f"{v:>5} {z:>3}   (no plan)")
        continue
    ap_, cum = replay(v, ops)
    bar_u, cost_u = max(cum), cum[-1]
    undo(ap_)
    assert abs(s.current_objective - S_ref) < 1e-6

    # candidate placement faces: faces of tets containing v
    cand, seen = [], set()
    for t in L.v2t[v]:
        for i in range(4):
            f = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if f in seen:
                continue
            seen.add(f)
            a = apexes(f)
            if a is None or a in L.edeg:
                continue
            cand.append((f, a))
    best = None
    for f, a in cand[:NF]:
        S0 = s.current_objective
        try:
            L.do(f, a)
        except Exception:
            continue
        d_place = s.current_objective - S0
        ops2 = plan(v)
        if ops2:
            ap2, cum2 = replay(v, ops2)
            bar_a = max(cum2) if cum2 else 0.0
            undo(ap2)
            reloc = d_place - HOME
            total = reloc + max(bar_a, 0.0)
            if best is None or total < best[0]:
                best = (total, reloc, bar_a, sig(f), (f, a))
        L.do(a, f)
        assert abs(s.current_objective - S0) < 1e-9
    # GREEDY SECOND FLICKER: fix the best single placement, try again
    best2 = None
    if best is not None and best[4] is not None:
        f1, a1 = best[4]
        S0 = s.current_objective
        L.do(f1, a1)
        d1 = s.current_objective - S0
        cand2, seen2 = [], set()
        for t in list(L.v2t[v]):
            for i in range(4):
                f = tuple(sorted(x for j, x in enumerate(t) if j != i))
                if f in seen2:
                    continue
                seen2.add(f)
                a = apexes(f)
                if a is None or a in L.edeg:
                    continue
                cand2.append((f, a))
        for f, a in cand2[:NF]:
            S1 = s.current_objective
            try:
                L.do(f, a)
            except Exception:
                continue
            d2 = s.current_objective - S1
            ops3 = plan(v)
            if ops3:
                ap3, cum3 = replay(v, ops3)
                bar2 = max(cum3) if cum3 else 0.0
                undo(ap3)
                reloc2 = (d1 - HOME) + (d2 - HOME)
                tot2 = reloc2 + max(bar2, 0.0)
                if best2 is None or tot2 < best2[0]:
                    best2 = (tot2, reloc2, bar2, sig(f))
            L.do(a, f)
            assert abs(s.current_objective - S1) < 1e-9
        L.do(a1, f1)
        assert abs(s.current_objective - S0) < 1e-9
    if best is None:
        print(f"{v:>5} {z:>3} {bar_u:>10.2f} {cost_u:>11.2f}   (no assisted plan)")
    else:
        tot, rel, ba, sg = best[:4]
        mark = "  <-- CATALYSIS" if tot < bar_u - 1e-9 else ""
        print(f"{v:>5} {z:>3} {bar_u:>10.2f} {cost_u:>11.2f} "
              f"{rel:>8.2f} {ba:>8.2f} {tot:>8.2f} {str(sg):>10}{mark}")
        if best2:
            t2, r2, b2, s2 = best2
            m2 = "  <-- BETTER" if t2 < tot - 1e-9 else ""
            print(f"{'':>5} {'2fl':>3} {'':>10} {'':>11} "
                  f"{r2:>8.2f} {b2:>8.2f} {t2:>8.2f} {str(s2):>10}{m2}")
    assert abs(s.current_objective - S_ref) < 1e-6

print("\nstate restored exactly")
