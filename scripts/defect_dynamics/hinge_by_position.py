"""The ONE validated deg-4 move (the 4->4 hinge flip) applied to real deg-4
lines, with FULL-action dS resolved by line position: tip (valence-1 end),
interior (valence-2 middle), branch-adjacent, isolated single edge.

dS is the full shape action in lam=1 units (EDQ with estar coupling + zleg*U(n6)
+ cimp*m^2); the acceptance weight is exp(-lam*dS). This is the length-0 local
flip, NOT the (unbuilt) transport worm -- but it is a real deg-4 move and shows
whether ends are cheaper to move than middles."""
import os, sys, itertools
from collections import Counter, defaultdict
import numpy as np
_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import discrete_differential_geometry as ddg
import worm_moves as wm

LAM, ESTAR, ZLEG, CIMP = 0.35, 5.105025, 0.6, 1.0
SNAP = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")

def edeg_of(tets):
    d = Counter()
    for t in tets:
        for e in itertools.combinations(sorted(t), 2): d[e] += 1
    return d

def link_cycle(tets_uv, u, v):
    pairs = [tuple(sorted(t - {u, v})) for t in tets_uv]
    if len(pairs) != 4 or any(len(p) != 2 for p in pairs): return None
    adj = defaultdict(list)
    for x, y in pairs: adj[x].append(y); adj[y].append(x)
    if any(len(a) != 2 for a in adj.values()) or len(adj) != 4: return None
    start = min(adj); cyc = [start]; prev = None; cur = start
    while len(cyc) < 4:
        nx = [w for w in adj[cur] if w != prev]
        if not nx: return None
        prev, cur = cur, nx[0]; cyc.append(cur)
    return cyc

def flip_tets(u, v, cyc, diag):
    p, q = (cyc[0], cyc[2]) if diag == 0 else (cyc[1], cyc[3])
    other = [x for x in cyc if x not in (p, q)]
    r, s = other
    rem = [frozenset((u, v, cyc[i], cyc[(i + 1) % 4])) for i in range(4)]
    add = [frozenset((p, q, u, r)), frozenset((p, q, r, v)),
           frozenset((p, q, v, s)), frozenset((p, q, s, u))]
    return rem, add, (p, q)

def dS_full(d0, d1):
    """full shape action diff, lam=1 units (as worm_slide.dS_between)."""
    x = ESTAR - int(ESTAR)
    dS = 0.0; changed = set()
    for k in set(d0) | set(d1):
        a, b = d0.get(k), d1.get(k)
        if a == b: continue
        changed |= set(k)
        if a is not None: dS -= (ESTAR / 6.0) * ((a - ESTAR) ** 2 - x * (1 - x))
        if b is not None: dS += (ESTAR / 6.0) * ((b - ESTAR) ** 2 - x * (1 - x))
    def counters(dd):
        n6 = defaultdict(int); mm = defaultdict(int)
        for (a, b), d in dd.items():
            if a not in changed and b not in changed: continue
            if d >= 6: n6[a] += 1; n6[b] += 1
            if d not in (5, 6): mm[a] += 1; mm[b] += 1
        return n6, mm
    n60, mm0 = counters(d0); n61, mm1 = counters(d1)
    for v in changed:
        dS += ZLEG * (wm.U_zleg(n61.get(v, 0)) - wm.U_zleg(n60.get(v, 0)))
        dS += CIMP * (mm1.get(v, 0) ** 2 - mm0.get(v, 0) ** 2)
    return dS

def main():
    m = ddg.Manifold.load(SNAP, 3)
    tets = {frozenset(int(x) for x in f) for f in m.facets()}
    d = edeg_of(tets)
    e2t = defaultdict(list)
    for t in tets:
        for e in itertools.combinations(sorted(t), 2): e2t[e].append(t)
    ill = {e: k for e, k in d.items() if k not in (5, 6)}
    d4 = [e for e, k in ill.items() if k == 4]
    # deg-4 line valence
    adjv = defaultdict(list)
    for a, b in d4: adjv[a].append(b); adjv[b].append(a)
    val = {v: len(nb) for v, nb in adjv.items()}
    def classify(e):
        a, b = e; va, vb = val[a], val[b]
        if va == 1 and vb == 1: return "isolated"
        if max(va, vb) >= 3: return "branch"
        if 1 in (va, vb): return "tip"
        return "interior"
    # D-core validation of one flip
    for e in d4:
        cyc = link_cycle(list(e2t[e]), *e)
        if cyc and not d.get(tuple(sorted((cyc[0], cyc[2]))), 0):
            m2 = ddg.Manifold.load(SNAP, 3); m2.do_hinge_move(list(e), cyc, 0)
            got = {frozenset(int(x) for x in f) for f in m2.facets()}
            rem, add, _ = flip_tets(e[0], e[1], cyc, 0)
            print(f"D-core validation: abstract==do_hinge_move: "
                  f"{(tets - set(rem) | set(add)) == got}")
            break
    # census by class
    res = defaultdict(list); dn = defaultdict(list)
    base_ill = sum(1 for k in d.values() if k not in (5, 6))
    for e in d4:
        u, v = e
        cyc = link_cycle(list(e2t[e]), u, v)
        if cyc is None: continue
        for diag in (0, 1):
            p, q = (cyc[0], cyc[2]) if diag == 0 else (cyc[1], cyc[3])
            if d.get(tuple(sorted((p, q))), 0): continue
            rem, add, _ = flip_tets(u, v, cyc, diag)
            if not all(r in tets for r in rem): continue
            nt = tets - set(rem) | set(add)
            nd = edeg_of(nt)
            cls = classify(e)
            res[cls].append(LAM * dS_full(d, nd))     # acceptance dS = lam*dS_shape
            dn[cls].append(sum(1 for k in nd.values() if k not in (5, 6)) - base_ill)
            break
    print(f"\n{SNAP.split('/')[-1]}: deg-4 edges by line position, 4->4 hinge flip")
    print(f"full-action acceptance dS = lam*dS_shape (lam={LAM}); "
          f"accept ~ exp(-dS)\n")
    print(f"{'class':10s} {'n':>4s} {'dS med':>7s} {'dS mean':>8s} "
          f"{'dS min':>7s} {'<accept>':>9s} {'dn_ill med':>10s}")
    for cls in ("tip", "interior", "branch", "isolated"):
        a = np.array(res[cls]);
        if len(a) == 0: continue
        acc = np.exp(-np.clip(a, 0, None)).mean()
        print(f"{cls:10s} {len(a):4d} {np.median(a):7.3f} {a.mean():8.3f} "
              f"{a.min():7.3f} {acc:9.3f} {np.median(dn[cls]):10.1f}")

if __name__ == "__main__":
    main()
