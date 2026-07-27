#!/usr/bin/env python3
"""Is the path -> parallel-transport map injective? (exact check)

Claim under test: the transport rotation essentially determines the
path (from a fixed base tet, mod backtracking). Equivalent statements:
  (R,t)-injective  <=>  no nontrivial developed-flat loop
                        (identity rotation AND zero translation);
  R-injective      <=>  no nontrivial loop with identity rotation at
                        all (translations allowed).

Method: enumerate ALL non-backtracking dual-graph paths from a base
tet of the pristine R m2 crystal up to length --maxlen, developing
incrementally (exact rational). Each path is keyed two ways:
  K1 = (end tet, canonical Wilson quaternion)            [R only]
  K2 = (end tet, full end placement -- 12 rationals)     [R and t]
Distinct paths sharing K2 are flat loops; sharing K1 only are
rotation-identity loops with net translation. Also censuses the
one-step quat norms (the 3-adic generators) and checks all closed
non-backtracking loops for identity rotation.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import development as dev
from discrete_differential_geometry.development import Quat, _SEED

REF = "data/tcp_reference/T3_R_m2_N7248.mfd"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxlen", type=int, default=7)
    ap.add_argument("--base", type=int, default=0)
    args = ap.parse_args()

    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    ctx = dev.TransportContext([tuple(int(x) for x in f) for f in F])
    t0 = frozenset(int(x) for x in F[args.base])

    # one-step quat census
    norms = defaultdict(int)
    for nb in ctx.dual[t0]:
        q = ctx.wilson_line([t0, nb])
        norms[q.norm] += 1
    print(f"one-step Wilson quat norms from base: {dict(norms)}")

    canonical = [_SEED[i] for i in range(4)]

    def wq(placement, tet):
        order = ctx.canon_order(tet)
        placed = [placement[v] for v in order]
        Rg = dev._rotation_between(canonical, placed)
        RgT = [tuple(Rg[i][j] for i in range(3)) for j in range(3)]
        q = Quat.lift(RgT)
        return q, tuple(x for v in order for x in placement[v])

    seen_R = {}          # (end, quat) -> path
    seen_Rt = {}         # (end, placement12) -> path
    coll_R, coll_Rt = [], []
    loops_checked = 0
    id_loops = []
    npaths = 0

    def dfs(path, placement, prev):
        nonlocal npaths, loops_checked
        cur = path[-1]
        if len(path) > 1:
            npaths += 1
            q, pl = wq(placement, cur)
            kR = (cur, q.canonical())
            kRt = (cur, pl)
            if kRt in seen_Rt:
                coll_Rt.append((seen_Rt[kRt], list(path)))
            else:
                seen_Rt[kRt] = list(path)
                if kR in seen_R:
                    coll_R.append((seen_R[kR], list(path)))
                else:
                    seen_R[kR] = list(path)
            if cur == t0:
                loops_checked += 1
                if q.canonical() == Quat(1, 0, 0, 0):
                    id_loops.append(list(path))
        if len(path) > args.maxlen:
            return
        for nb in ctx.dual[cur]:
            if nb == prev:
                continue
            shared = cur & nb
            drop = next(v for v in cur if v not in shared)
            new = next(v for v in nb if v not in shared)
            p2 = {v: placement[v] for v in shared}
            p2[new] = dev.reflect(placement[drop],
                                  *(p2[v] for v in shared))
            merged = dict(placement)
            merged.update(p2)
            dfs(path + [nb], merged, cur)

    dfs([t0], ctx.canon_placement(t0), None)
    print(f"{npaths} non-backtracking paths (len <= {args.maxlen}), "
          f"{loops_checked} closed loops among them")
    print(f"(R,t) collisions (flat loops): {len(coll_Rt)}")
    print(f"R-only collisions (rotation-identity loops with net "
          f"translation): {len(coll_R)}")
    print(f"closed loops with identity rotation: {len(id_loops)}")
    for a, b in coll_Rt[:3]:
        print(f"  FLAT: {[sorted(t) for t in a]}  vs  "
              f"{[sorted(t) for t in b]}")
    for a, b in coll_R[:3]:
        print(f"  R-COLL: len {len(a)-1} vs len {len(b)-1}")
    if not coll_R and not coll_Rt and not id_loops:
        print("\nINJECTIVE at this depth: the Wilson rotation alone "
              "determines the path (from this base, mod backtracking).")


if __name__ == "__main__":
    main()
