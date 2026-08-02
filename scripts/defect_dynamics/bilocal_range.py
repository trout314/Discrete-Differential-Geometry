"""Bilocal step-1/2, part B: factorization and teleport cost vs TRUE
GRAPH DISTANCE, with sites chosen anywhere (not restricted to one BC
chain).

Part A showed the local interaction switches off the moment the two
regions stop sharing a vertex -- but also that BC-chain index is NOT a
proxy for spatial separation (the chain re-approaches itself at s=48
and s=96, where the residual jumps back to several units). So here B is
chosen by BFS shell from A over ALL 2->3 sites in the manifold.

Measures, per graph distance d:
  - max |residual| after subtracting the analytic cross term X (should
    be machine zero for d >= 1),
  - the CONSERVING pair cost dS_AB (create knot at A, annihilate knot at
    B; df = 0 exactly, hence pin-free), and the fraction of pairs that
    are exactly FREE (same washboard rung).
"""
import json
import os
import sys
import time
from collections import defaultdict

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, os.path.join(_ROOT, "scripts", "defect_dynamics"))
sys.argv = ["bilocal", "unused", "unused"]

import f0_worm as F  # noqa: E402
import discrete_differential_geometry as ddg  # noqa: E402

C_V, C_H = 0.1, 1.0
MFD = os.environ.get("MFD", os.path.join(SCR,
                                         "quench_down5q_wOFF.final.mfd"))
PER_D = int(os.environ.get("PER_D", "12"))     # samples per distance
OUT = os.environ.get("OUT", os.path.join(SCR, "bilocal_range.json"))


def g(f1, f3):
    return C_V * (f3 - F.F3T) ** 2 + C_H * (f1 - 6.0 * f3 / F.ETARGET) ** 2


def cross(dA, dB):
    (a1, a3), (b1, b3) = dA, dB
    wa, wb = a1 - 6.0 * a3 / F.ETARGET, b1 - 6.0 * b3 / F.ETARGET
    return 2.0 * C_V * a3 * b3 + 2.0 * C_H * wa * wb


def flip_sites(L_):
    """All valid 2->3 moves: (face, apex pair) with exactly two tets on
    the face and the apex edge absent."""
    out = {}
    for v, ts in L_.v2t.items():
        for t in ts:
            for i in range(4):
                face = tuple(sorted(x for j, x in enumerate(t) if j != i))
                if face in out:
                    continue
                fs = set(face)
                sh = [q for q in L_.v2t[face[0]] if fs <= set(q)]
                if len(sh) != 2:
                    continue
                ap = tuple(sorted({x for q in sh for x in q} - fs))
                if len(ap) != 2 or ap in L_.edeg:
                    continue
                out[face] = ap
    return out


def bfs_from(L_, srcs, cap=40):
    from collections import deque
    dist = {v: 0 for v in srcs}
    dq = deque(srcs)
    while dq:
        v = dq.popleft()
        d = dist[v]
        if d >= cap:
            continue
        for t in L_.v2t[v]:
            for u in t:
                if u not in dist:
                    dist[u] = d + 1
                    dq.append(u)
    return dist


def main():
    t0 = time.time()
    m = ddg.Manifold.load(MFD, 3)
    s, L = F.fresh(m)
    fv = [int(x) for x in s.manifold.f_vector]
    print(f"manifold {os.path.basename(MFD)} f={fv} "
          f"gap {fv[1] - 6*fv[3]/F.ETARGET:+.3f}")

    sites = flip_sites(L)
    print(f"2->3 sites: {len(sites)}")

    faceA, apA = next(iter(sorted(sites.items())))
    opA = (faceA, apA)
    dist = bfs_from(L, set(faceA) | set(apA))
    print(f"head A face {faceA} apex {apA}; BFS reach "
          f"{max(dist.values())} hops over {len(dist)} vertices")

    buckets = defaultdict(list)
    for face, ap in sites.items():
        d = min(dist.get(v, 999) for v in set(face) | set(ap))
        buckets[d].append((face, ap))

    S0 = s.current_objective
    rows = []
    print(f"\n{'d':>3} {'n':>4} {'max|resid|':>11} {'dS_AB values':>34} "
          f"{'free%':>6}")
    for d in sorted(buckets):
        if d > 24:
            continue
        cands = sorted(buckets[d])[:PER_D]
        worst, costs = 0.0, []
        for faceB, apB in cands:
            if set(faceB) & set(faceA) or set(apB) & set(apA):
                pass  # overlapping is exactly what d==0 means; keep it
            Spre = s.current_objective
            try:
                L.do(faceB, apB)               # base: knot at B
            except Exception:
                continue
            # dS_A: create knot at A | knot at B
            SA0 = s.current_objective
            try:
                L.do(*opA)
            except Exception:
                L.do(apB, faceB)
                continue
            dA = s.current_objective - SA0
            fA = [int(x) for x in s.manifold.f_vector]
            L.do(opA[1], opA[0])
            # dS_B: annihilate knot at B
            L.do(apB, faceB)
            dB = s.current_objective - SA0
            L.do(faceB, apB)
            # dS_AB: both
            try:
                L.do(*opA)
                L.do(apB, faceB)
            except Exception:
                L.do(apB, faceB)
                continue
            dAB = s.current_objective - SA0
            fAB = [int(x) for x in s.manifold.f_vector]
            pin = g(fAB[1], fAB[3]) - g(
                *[int(x) for x in (fA[1] - 1, fA[3] - 1)])
            # restore
            L.do(faceB, apB)
            L.do(opA[1], opA[0])
            L.do(apB, faceB)
            assert abs(s.current_objective - Spre) < 1e-9, "undo drift"
            X = cross((1, 1), (-1, -1))
            resid = dAB - dA - dB - X
            worst = max(worst, abs(resid))
            costs.append(round(dAB, 4))
            rows.append({"d": d, "faceB": list(faceB), "dS_A": dA,
                         "dS_B": dB, "dS_AB": dAB, "resid": resid})
        if not costs:
            continue
        uniq = sorted(set(costs))
        free = 100.0 * sum(1 for c in costs if abs(c) < 1e-9) / len(costs)
        us = ",".join(f"{u:g}" for u in uniq[:6])
        print(f"{d:>3} {len(costs):>4} {worst:>11.2e} {us:>34} "
              f"{free:>5.0f}%")
    json.dump({"mfd": MFD, "rows": rows}, open(OUT, "w"), indent=1)
    assert abs(s.current_objective - S0) < 1e-9
    print(f"\nwrote {OUT} ({time.time() - t0:.0f}s)")


main()
