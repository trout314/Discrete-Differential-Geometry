#!/usr/bin/env python3
"""Exact test battery for the parallel-transport / holonomy toolkit
(development.Quat, TransportContext, chain transport). Everything is
an exact equality check -- no tolerances.

  T1  quaternion algebra: lift/to_matrix roundtrip, composition
      homomorphism, inverse.
  T2  BC chain step: spinor lift has norm 6, cos = -2/3, axis == the
      screw axis from chain_axis; k-step transport reproduces the
      Chebyshev values cos(k phi) = T_k(-2/3) exactly.
  T3  edge-loop holonomy on the pristine m2 crystal: the loop around
      an interior edge of degree k is a rotation about that edge by
      the deficit angle, cos = T_k(1/3) exactly, axis parallel to the
      developed edge. Identifies the dock loop-holonomy spectrum:
      deg-5 -> 241/243 (7.36 deg), deg-6 -> 329/729 (63.17 deg).
  T4  Wilson-line composability: transport A->B->C equals transport
      of the concatenated path (global canonical frames).
  T5  spin lifts: the edge-loop's SU(2) lift and the 2x(edge loop)
      composition -- the double cover in action.
"""
import os
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import development as dev
from discrete_differential_geometry.development import Quat

REF = os.path.join(_HERE, "..", "..",
                   "data/tcp_reference/T3_R_m2_N7248.mfd")


def cheb(k, c):
    """T_k(c) exactly (Fraction)."""
    t0, t1 = Fraction(1), Fraction(c)
    for _ in range(k - 1):
        t0, t1 = t1, 2 * c * t1 - t0
    return t1 if k >= 1 else t0


def main():
    ok = True

    # ---- T1: quaternion algebra ----
    q1 = Quat(1, 2, 1, 0)
    q2 = Quat(3, -1, 2, 5)
    for q in (q1, q2, q1 * q2):
        assert Quat.lift(q.to_matrix()).same_rotation(q)
    M12 = (q1 * q2).to_matrix()
    R1, R2 = q1.to_matrix(), q2.to_matrix()
    prod = [tuple(sum(R1[i][k] * R2[k][j] for k in range(3))
                  for j in range(3)) for i in range(3)]
    assert M12 == prod, "Hamilton product is not a homomorphism?!"
    assert (q1 * q1.conj()).same_rotation(Quat(1, 0, 0, 0))
    print("T1 quaternion algebra: PASS")

    # ---- T2: chain step + Chebyshev ----
    seq = list(range(20))
    pos = dev.develop_chain(seq)
    qs = dev.chain_step_quat(seq, pos)
    assert qs.norm == 6, f"chain step norm {qs.norm} != 6"
    assert qs.cos_phi() == Fraction(-2, 3)
    ax_screw = dev.chain_axis(seq, pos)
    ax_q = qs.axis()
    assert dev.cos2_between(ax_q, ax_screw) == 1, \
        "quaternion axis != screw axis"
    for k in (2, 3, 5, 8):
        assert (qs ** k).cos_phi() == cheb(k, Fraction(-2, 3)), \
            f"k={k} transport fails Chebyshev"
    print(f"T2 chain step: PASS (q_step = {qs}, norm 6; "
          f"Chebyshev k=2,3,5,8 exact)")

    # ---- T3: edge-loop holonomy on the crystal ----
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    ctx = dev.TransportContext([tuple(int(x) for x in f) for f in F])
    # star map: edge -> tets
    star = defaultdict(list)
    for f in F:
        t = frozenset(int(x) for x in f)
        vs = sorted(t)
        for i in range(4):
            for j in range(i + 1, 4):
                star[(vs[i], vs[j])].append(t)

    def loop_of(edge):
        tets = star[edge]
        loop = [tets[0]]
        used = {tets[0]}
        while len(loop) < len(tets):
            cur = loop[-1]
            nxt = next(t for t in tets
                       if t not in used and len(t & cur) == 3)
            loop.append(nxt)
            used.add(nxt)
        return loop + [loop[0]]

    results = {}
    for edge, tets in star.items():
        k = len(tets)
        if k in results or k not in (5, 6):
            continue
        loop = loop_of(edge)
        q, cphi, ax = ctx.holonomy(loop)
        assert cphi == cheb(k, Fraction(1, 3)), \
            f"deg-{k} loop cos {cphi} != T_{k}(1/3) {cheb(k, Fraction(1,3))}"
        # axis parallel to the developed edge
        base = ctx.canon_placement(loop[0])
        evec = dev._sub(base[edge[1]], base[edge[0]])
        assert dev.cos2_between(ax, evec) == 1, "axis not along the edge"
        results[k] = (q, cphi)
        if len(results) == 2:
            break
    print("T3 edge loops: PASS --")
    for k, (q, cphi) in sorted(results.items()):
        deg = float(np.degrees(np.arccos(float(cphi))))
        print(f"    deg-{k} edge: cos phi = {cphi} ({deg:.2f} deg), "
              f"spin lift {q}")

    # ---- T4: Wilson-line composability ----
    t0 = frozenset(int(x) for x in F[0])
    pathAB = [t0]
    for _ in range(4):
        pathAB.append(ctx.dual[pathAB[-1]][0])
    pathBC = [pathAB[-1]]
    for _ in range(4):
        pathBC.append(ctx.dual[pathBC[-1]][1])
    qAB = ctx.wilson_line(pathAB)
    qBC = ctx.wilson_line(pathBC)
    qAC = ctx.wilson_line(pathAB + pathBC[1:])
    assert (qBC * qAB).same_rotation(qAC), "composition failed"
    print("T4 Wilson-line composition: PASS")

    # ---- T5: spin lifts / double cover ----
    k5 = 5
    edge5 = next(e for e, ts in star.items() if len(ts) == 5)
    loop5 = loop_of(edge5)
    q5 = ctx.wilson_line(loop5)
    q10 = ctx.wilson_line(loop5 + loop5[1:])
    assert q10.same_rotation(q5 * q5)
    # the double loop's rotation angle: cos(2 delta) exact
    assert q10.cos_phi() == 2 * q5.cos_phi() ** 2 - 1
    print(f"T5 spin lifts: PASS -- deg-5 loop lift {q5} "
          f"(a {'>' if q5.a > 0 else '<='} 0), double loop {q10}")

    print("\nALL TRANSPORT TESTS PASS (exact)")


if __name__ == "__main__":
    main()
