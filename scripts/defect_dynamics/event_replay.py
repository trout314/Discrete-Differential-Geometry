#!/usr/bin/env python3
"""Offline replay of the accepted-move event stream (events.bin / drain_event_log).

Each EVENT_DTYPE record (clock, type, 6 labels) determines an exact tet-set
change; apply_event() maintains the running triangulation move-by-move at full
resolution. This lets us catch the birth / hop / death of a single deg-4 edge in
an FK neighbourhood -- the realized minimal disclination move -- which is too
transient for 150-sweep snapshots.

Move semantics (from move_geometry.replay_role_counts + SUPPORT_SPLIT):
  0 1->4 : A=tet(4), B=[new v]         rem tet(A); add tri(A)+v  (x4)
  1 2->3 : A=triangle(3), B=2 poles    rem A+pole (x2); add edge(A)+2poles (x3)
  2 3->2 : A=2 poles(edge), B=tri(3)   rem edge(B)+2poles (x3); add B+pole (x2)
  3 4->1 : A=[dead v], B=tet(4)         rem v+tri(B) (x4); add tet(B)
  4 4-4  : A=2 poles, B=4-cycle         rem A0A1+cycle-edge (x4); add B0B2+other (x4)
"""
import os
import sys
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "python"))
from discrete_differential_geometry.move_geometry import EVENT_DTYPE, SUPPORT_SPLIT


def _support(ev):
    nA, nB = SUPPORT_SPLIT[int(ev["type"])]
    labs = [int(x) for x in ev["labels"]]
    return labs[:nA], labs[nA:nA + nB]


def _tet(*vs):
    return tuple(sorted(vs))


def event_changes(ev):
    """(remove, add) tet lists for one accepted move."""
    t = int(ev["type"])
    A, B = _support(ev)
    if t == 0:                                   # 1->4
        return [_tet(*A)], [_tet(*tri, B[0]) for tri in combinations(A, 3)]
    if t == 1:                                   # 2->3
        return ([_tet(*A, B[0]), _tet(*A, B[1])],
                [_tet(*e, B[0], B[1]) for e in combinations(A, 2)])
    if t == 2:                                   # 3->2
        return ([_tet(*e, A[0], A[1]) for e in combinations(B, 2)],
                [_tet(*B, A[0]), _tet(*B, A[1])])
    if t == 3:                                   # 4->1
        return ([_tet(A[0], *tri) for tri in combinations(B, 3)], [_tet(*B)])
    # 4-4: removed edge A0A1 (deg 4), link 4-cycle B; new diagonal B0-B2
    rem = [_tet(A[0], A[1], B[i], B[(i + 1) % 4]) for i in range(4)]
    add = [_tet(B[0], B[2], A[0], B[1]), _tet(B[0], B[2], B[1], A[1]),
           _tet(B[0], B[2], A[1], B[3]), _tet(B[0], B[2], B[3], A[0])]
    return rem, add


def apply_event(tets, ev, strict=True):
    rem, add = event_changes(ev)
    for r in rem:
        if strict and r not in tets:
            raise KeyError(f"remove miss type={int(ev['type'])} {r}")
        tets.discard(r)
    tets.update(add)


def _verify():
    import discrete_differential_geometry as ddg
    sys.path.insert(0, _HERE + "/../")            # scripts/ for fk_skeleton
    from fk_skeleton import edges_from_facets
    R = os.path.join(_HERE, "..", "..")
    cell = os.path.join(R, "data/tcp_reference/T3_R_m2_N7248.mfd")
    ddg.set_random_seed(7)
    ref = ddg.Manifold.load(cell, 3)
    native = float(edges_from_facets(ref.facets())[1].mean())
    m = ddg.Manifold.load(cell, 3)
    p = ddg.SamplerParams(num_facets_target=ref.num_facets, num_facets_coef=0.1,
                          hinge_degree_target=native, num_hinges_coef=2.0,
                          hinge_degree_variance_coef=0.0,
                          codim3_degree_variance_coef=0.0,
                          hinge_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, 0.0)                  # free -> lots of moves incl 4-4
    v = s.manifold
    s.enable_event_log(64.0)
    s.drain_event_log()
    tets = {tuple(sorted(int(x) for x in t)) for t in v.facets()}
    ok = True
    for c in range(12):
        s.run(sweeps=20)
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            print("  overflow!"); ok = False; break
        try:
            for e in ev:
                apply_event(tets, e)
        except KeyError as err:
            print(f"  chunk {c}: {err}"); ok = False; break
        actual = {tuple(sorted(int(x) for x in t)) for t in v.facets()}
        match = tets == actual
        ok &= match
        print(f"  chunk {c}: {len(ev):4d} events, tets={len(tets)} "
              f"{'MATCH' if match else f'MISMATCH symdiff={len(tets ^ actual)}'}")
        if not match:
            break
    print(f"\nreplay verification {'PASSED' if ok else 'FAILED'}")


if __name__ == "__main__":
    _verify()
