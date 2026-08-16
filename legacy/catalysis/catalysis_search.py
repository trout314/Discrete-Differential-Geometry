"""Find a CATALYSED path: a flicker that enables a net-downhill
rearrangement it could not do alone.

Why a targeted seed. A 2->3 on face (a,b,c) drops the degrees of its
three face edges by one. In a pristine region every edge sits at 5 or 6,
so they land on 4 or 5 -- never 3 -- and no downhill 3->2 is unlocked.
That is exactly why the head kernel measured accH = 0 everywhere. So the
flicker can only catalyse where a face edge is ALREADY at degree 4: the
flip takes it to 3 and a 3->2 becomes available.

Search: seed on faces carrying a degree-4 edge, create the flicker, then
best-first over the head-move class (support contains BOTH chord
endpoints -- the same class wf0ChordEnumH enumerates) to depth DEPTH.
A node is a HIT when the flicker can be annihilated and the net action,
measured against the flicker-free start, is negative -- i.e. the whole
create/work/annihilate cycle paid for itself.
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts", "defect_dynamics"))
sys.argv = ["t", "unused", "unused"]

import f0_worm as F  # noqa: E402
import discrete_differential_geometry as ddg  # noqa: E402

MFD = os.environ.get("MFD", os.path.join(SCR,
                                         "quench_down5q_wOFF.final.mfd"))
DEPTH = int(os.environ.get("DEPTH", "4"))
NSEED = int(os.environ.get("NSEED", "40"))
NODECAP = int(os.environ.get("NODECAP", "4000"))


def apexes(L, face):
    fs = set(face)
    ts = [t for t in L.v2t[face[0]] if fs <= set(t)]
    if len(ts) != 2:
        return None
    return tuple(sorted({x for t in ts for x in t} - fs))


def chord_tets(L, c):
    return [t for t in L.v2t[c[0]] if c[0] in t and c[1] in t]


def head_moves(L, c):
    """Moves whose support contains BOTH chord endpoints (the class the
    D head kernel enumerates). Excludes the 3->2 on the chord itself."""
    out, seenf, seene = [], set(), set()
    for t in chord_tets(L, c):
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if face in seenf:
                continue
            seenf.add(face)
            ap = apexes(L, face)
            if ap is None or len(ap) != 2 or ap in L.edeg:
                continue
            if not (set(c) <= set(face) | set(ap)):
                continue
            out.append(("23", face, ap))
        vs = sorted(t)
        for i in range(4):
            for j in range(i + 1, 4):
                e = (vs[i], vs[j])
                if e in seene or e == tuple(sorted(c)):
                    continue
                seene.add(e)
                if L.edeg.get(e, 0) != 3:
                    continue
                ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
                lk = tuple(sorted({x for q in ts for x in q} - set(e)))
                if len(lk) != 3:
                    continue
                if not (set(c) <= set(lk) | set(e)):
                    continue
                out.append(("32", e, lk))
    return out


def can_annihilate(L, c):
    e = tuple(sorted(c))
    if L.edeg.get(e, 0) != 3:
        return None
    ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
    lk = tuple(sorted({x for q in ts for x in q} - set(e)))
    return lk if len(lk) == 3 else None


m = ddg.Manifold.load(MFD, 3)
s, L = F.fresh(m)
S_ref = s.current_objective
pairs, degs = m.illegal_edges()   # a TUPLE (pairs, degs) -- len() is 2
from collections import Counter as _C
print(f"start S={S_ref:.4f}  n_ill={len(pairs)}  "
      f"{dict(sorted(_C(int(d) for d in degs).items()))}")

# seeds: valid 2->3 sites whose face carries a degree-4 edge
seeds = []
seen = set()
for v, ts in L.v2t.items():
    for t in ts:
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if face in seen:
                continue
            seen.add(face)
            d4 = any(L.edeg.get(tuple(sorted((face[a], face[b]))), 0) == 4
                     for a in range(3) for b in range(a + 1, 3))
            if not d4:
                continue
            ap = apexes(L, face)
            if ap is None or len(ap) != 2 or ap in L.edeg:
                continue
            seeds.append((face, ap))
print(f"seed faces carrying a degree-4 edge: {len(seeds)}")

hits = []
depth1 = []
for si, (face, ap) in enumerate(sorted(seeds)[:NSEED]):
    S0 = s.current_objective
    try:
        L.do(face, ap)
    except Exception:
        continue
    chord = tuple(sorted(ap))
    # best-first over head moves, tracking the applied path
    best = None
    nodes = 0

    def rec(path, depth):
        global nodes, best
        nodes += 1
        if nodes > NODECAP:
            return
        lk = can_annihilate(L, chord)
        if lk is not None:
            Sa = s.current_objective
            L.do(chord, lk)
            net = s.current_objective - S0
            L.do(lk, chord)
            assert abs(s.current_objective - Sa) < 1e-9
            if len(path) >= 1 and (best is None or net < best[1]):
                best = (list(path), net)   # NON-TRIVIAL paths only
        if depth >= DEPTH:
            return
        for kind, x, y in head_moves(L, chord):
            cen, coc = (x, y) if kind == "23" else (x, y)
            try:
                L.do(cen, coc)
            except Exception:
                continue
            path.append((cen, coc))
            rec(path, depth + 1)
            path.pop()
            L.do(coc, cen)

    # DIAGNOSTIC: does any single head move leave the chord at degree 3?
    keep3 = 0
    hm = head_moves(L, chord)
    for kind, x, y in hm:
        try:
            L.do(x, y)
        except Exception:
            continue
        if can_annihilate(L, chord) is not None:
            keep3 += 1
        L.do(y, x)
    depth1.append((len(hm), keep3))

    rec([], 0)
    # undo the flicker
    lk = can_annihilate(L, chord)
    if lk is not None:
        L.do(chord, lk)
    assert abs(s.current_objective - S0) < 1e-6, "seed not restored"
    if best is not None:
        hits.append((face, ap, best[0], best[1], nodes))
        if si < 12 or best[1] < 0:
            tag = "HIT " if best[1] < -1e-9 else "best"
            print(f"  seed {si:3d} {tag} net dS {best[1]:+8.4f} "
                  f"in {len(best[0])} moves ({nodes} nodes)")

# does the best path reach a DIFFERENT configuration, or is it a
# move/inverse cycle that returns to the start?
def statehash():
    return hash(frozenset(tuple(sorted(t)) for ts in L.v2t.values()
                          for t in ts))
distinct = 0
for face, ap, path, net, _ in hits[:12]:
    h0 = statehash()
    L.do(face, ap)
    chord = tuple(sorted(ap))
    for cen, coc in path:
        L.do(cen, coc)
    lk = can_annihilate(L, chord)
    changed = None
    if lk is not None:
        L.do(chord, lk)
        changed = (statehash() != h0)
        L.do(lk, chord)
    for cen, coc in reversed(path):
        L.do(coc, cen)
    L.do(chord, face) if False else None
    lk2 = can_annihilate(L, chord)
    if lk2 is not None:
        L.do(chord, lk2)
    if changed:
        distinct += 1
print(f"\nof the first {min(12,len(hits))} best paths, "
      f"{distinct} reach a DIFFERENT configuration "
      f"(the rest are move/inverse cycles)")
if depth1:
    tot = sum(a for a, _ in depth1)
    k3 = sum(b for _, b in depth1)
    print(f"\ndepth-1 diagnostic: {tot} head moves tried, "
          f"{k3} left the chord annihilable "
          f"({100.0 * k3 / max(tot, 1):.1f}%)")
print(f"\nseeds tried {min(len(seeds), NSEED)}, "
      f"catalysed paths found {len(hits)}")
if hits:
    hits.sort(key=lambda h: h[3])
    f2, a2, path, net, _ = hits[0]
    print(f"best: net dS {net:+.4f}, {len(path)} intermediate moves")
    for k, (cen, coc) in enumerate(path):
        print(f"   {k}: {cen} -> {coc}")
assert abs(s.current_objective - S_ref) < 1e-6, "global drift"
print("state restored exactly")
