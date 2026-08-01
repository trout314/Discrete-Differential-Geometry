"""Does the cheap proposal weight w(f) = [f has an edge of degree <= 4]
actually concentrate on the faces that lower the collapse barrier?

For every candidate face near a target vertex, measure the TRUE value --
place the flicker, re-plan the collapse, record
    total = (dS_place - 5.930) + max_k(cumulative dS)
-- and compare the ranking that induces against the one w(f) induces.
w is only useful if the faces it flags are the ones that pay.

Reported: enrichment (how much more likely a w=1 draw is to be good
than a uniform draw), where the true best sits in the w ranking, and
the barrier split between the two groups.
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

HOME = 5.930
NV = int(os.environ.get("NV", "3"))
NF = int(os.environ.get("NF", "40"))
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


def sig(f):
    return tuple(sorted(L.edeg.get((min(a, b), max(a, b)), 0)
                        for a, b in ((f[0], f[1]), (f[0], f[2]),
                                     (f[1], f[2]))))


verts = sorted({x for e in L.edeg for x in e},
               key=lambda u: (len(F.nb_of(L, u)), u))[:NV]
allrows = []
for v in verts:
    ops = plan(v)
    if not ops:
        continue
    ap_, cum = replay(v, ops)
    bar_u = max(cum)
    undo(ap_)
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
    rows = []
    for f, a in cand[:NF]:
        sg = sig(f)
        w = 1 if min(sg) <= 4 else 0
        S0 = s.current_objective
        try:
            L.do(f, a)
        except Exception:
            continue
        d = s.current_objective - S0
        o2 = plan(v)
        if o2:
            ap2, cum2 = replay(v, o2)
            tot = (d - HOME) + max(max(cum2) if cum2 else 0.0, 0.0)
            rows.append((tot, w, sg, d - HOME))
            undo(ap2)
        L.do(a, f)
        assert abs(s.current_objective - S0) < 1e-9
    if not rows:
        continue
    rows.sort()
    allrows.append((v, bar_u, rows))
    good = [r for r in rows if r[0] < bar_u - 1e-9]
    nw = sum(1 for r in rows if r[1])
    gw = sum(1 for r in good if r[1])
    print(f"\nv={v}  unassisted barrier {bar_u:.2f}   "
          f"{len(rows)} candidate faces, {nw} with w=1")
    print(f"  {'total':>8} {'reloc':>8}  w  sig")
    for tot, w, sg, rel in rows[:6]:
        print(f"  {tot:>8.2f} {rel:>8.2f}  {w}  {sg}")
    print(f"  ... best-{len(rows)}: worst {rows[-1][0]:.2f}")
    print(f"  catalytic faces (total < unassisted): {len(good)}/{len(rows)}"
          f";  of those, w=1: {gw}/{len(good)}")
    if len(good) and nw:
        pu = len(good) / len(rows)
        pw = gw / nw
        print(f"  P(good | uniform) = {pu:.2f}   P(good | w=1) = {pw:.2f}"
              f"   enrichment {pw / max(pu, 1e-9):.2f}x")
    bw = [r[0] for r in rows if r[1]]
    b0 = [r[0] for r in rows if not r[1]]
    if bw and b0:
        print(f"  mean total  w=1: {np.mean(bw):7.2f}   "
              f"w=0: {np.mean(b0):7.2f}")
    print(f"  true best is w={rows[0][1]}, sig {rows[0][2]}")
    tot = np.array([r[0] for r in rows])
    rel = np.array([r[3] for r in rows])
    print(f"  corr(reloc, total) = {np.corrcoef(rel, tot)[0,1]:+.3f}"
          f"   reloc range [{rel.min():+.2f}, {rel.max():+.2f}]"
          f"   total range [{tot.min():.2f}, {tot.max():.2f}]")
    # what does plain Metropolis acceptance (w ~ e^-dS_place) deliver?
    pw = np.exp(-(rel - rel.min()))
    pw /= pw.sum()
    print(f"  E[total] under uniform proposal      = {tot.mean():7.2f}")
    print(f"  E[total] under w ~ exp(-dS_place)    = "
          f"{float((pw * tot).sum()):7.2f}")
    print(f"  E[total] under an oracle (best face) = {tot.min():7.2f}")
    bysig = {}
    for t, w, sg, r in rows:
        bysig.setdefault(sg, []).append(t)
    print("  mean total by PRE-placement signature:")
    for sg, v2 in sorted(bysig.items(), key=lambda kv: np.mean(kv[1])):
        print(f"     {str(sg):>12}  n={len(v2):2d}  {np.mean(v2):7.2f}")
assert abs(s.current_objective - S_ref) < 1e-6
print("\nstate restored exactly")
