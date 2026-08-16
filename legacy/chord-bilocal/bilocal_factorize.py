"""Bilocal step-1/2: does the cost of two head regions FACTORIZE, and at
what separation?

The action splits as S = g(f1,f3) + (local terms), with

    g(f1,f3) = c_v (f3 - F3T)^2 + c_h (f1 - 6 f3/e*)^2

the only terms NONLINEAR in the f-vector. For two regions with f-offsets
dA, dB the global part NEVER factorizes -- it carries an exact,
s-INDEPENDENT cross term

    X(dA,dB) = 2 c_v dA3 dB3 + 2 c_h dwA dwB,   dw = df1 - 6 df3/e*

(the same quadratic identity wf0Compile exploits). The knot-collider work
sidestepped this by excluding the pins from dS_between; here we keep them
and subtract X analytically. Two tests:

  T1 non-conserving pair (knot at A + knot at B, dA = dB = (+1,+1)):
     residual R(s) = dS_AB - dS_A - dS_B - X  must -> 0 beyond contact.
     Its LOCAL twin (all pin content stripped) must do the same. This
     both measures the local interaction range s_min and validates X.

  T2 conserving pair (create knot at A, ANNIHILATE a knot at B;
     dA + dB = 0 exactly): the joint pin content g(f+dA+dB) - g(f) is
     identically zero, so the pair should cost its pure local part at
     ANY separation while either half alone pays the full pin. This is
     the design claim -- a conserving bilocal move is pin-free.

Separation is BC-chain index (knot_collider convention: 4 chain steps =
1 slide; chord-sharing contact at s = 4). Every measurement is applied
through the sampler and exactly unwound.
"""
import json
import os
import sys
import time

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, os.path.join(_ROOT, "scripts", "defect_dynamics"))
sys.argv = ["bilocal", "unused", "unused"]

import f0_worm as F  # noqa: E402  (action + Live wiring, identical to the worm)
import discrete_differential_geometry as ddg  # noqa: E402

C_V, C_H = 0.1, 1.0            # num_facets_coef, num_hinges_coef
MFD = os.environ.get("MFD", os.path.join(SCR,
                                         "quench_down5q_wOFF.final.mfd"))
SEPS = [int(x) for x in os.environ.get(
    "SEPS", "1,2,3,4,5,6,7,8,10,12,16,20,24,32,40,48").split(",")]
OUT = os.environ.get("OUT", os.path.join(SCR, "bilocal_factorize.json"))


def g(f1, f3):
    return C_V * (f3 - F.F3T) ** 2 + C_H * (f1 - 6.0 * f3 / F.ETARGET) ** 2


def cross(dA, dB):
    """X(dA,dB): the exact, s-independent global cross term."""
    (a1, a3), (b1, b3) = dA, dB
    wa, wb = a1 - 6.0 * a3 / F.ETARGET, b1 - 6.0 * b3 / F.ETARGET
    return 2.0 * C_V * a3 * b3 + 2.0 * C_H * wa * wb


def chain_face(chain, j):
    L = len(chain)
    return [chain[(j + 1) % L], chain[(j + 2) % L], chain[(j + 3) % L]]


def apexes_of(L_, face):
    fs = set(face)
    ts = [t for t in L_.v2t[face[0]] if fs <= set(t)]
    return sorted({x for t in ts for x in t} - fs), ts


def next_window(L_, w):
    """manifold.d nextChainWindow: drop w[0], exit through (w1,w2,w3),
    adopt the apex opposite w[0]. Deterministic and reversible, so the
    directed-window graph is a permutation -- chains never branch."""
    ap, ts = apexes_of(L_, (w[1], w[2], w[3]))
    if len(ts) != 2:
        return None
    other = [a for a in ap if a != w[0]]
    if len(other) != 1:
        return None
    return [w[1], w[2], w[3], other[0]]


def bc_chain(L_, w0, maxlen=40000):
    """Cyclic vertex sequence of the BC chain through window w0. Window j
    is chain[j..j+3]; a knot at j has chord (chain[j], chain[j+4])."""
    w = list(w0)
    chain = list(w)
    while len(chain) < maxlen:
        nxt = next_window(L_, w)
        if nxt is None:
            return None
        w = nxt
        k = len(chain) - 3
        if w == list(w0):
            return chain[:k]
        chain.append(w[3])
    return None


def knot_op(L_, chain, j):
    """(center, cocenter) of the 2->3 that creates the window-j knot, or
    None if the chain face is not flippable here."""
    L = len(chain)
    face = chain_face(chain, j)
    ap, ts = apexes_of(L_, face)
    want = sorted({chain[j % L], chain[(j + 4) % L]})
    if len(ts) != 2 or ap != want:
        return None
    if tuple(want) in L_.edeg:          # apex edge already present
        return None
    return (tuple(sorted(face)), tuple(want))


def bfs_from(L_, srcs, cap=14):
    """Vertex-graph distances from a source set (Live.v2t adjacency)."""
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


class Meas:
    """Apply an op list through the sampler, capture (dS, df), unwind."""

    def __init__(self, s, L_):
        self.s, self.L = s, L_

    def run(self, ops):
        s, L_ = self.s, self.L
        S0 = s.current_objective
        f0 = [int(x) for x in s.manifold.f_vector]
        done = []
        try:
            for cen, coc in ops:
                L_.do(cen, coc)
                done.append((cen, coc))
        except Exception as e:
            for cen, coc in reversed(done):
                L_.do(coc, cen)
            assert abs(s.current_objective - S0) < 1e-9
            return None, None, str(e)
        f1 = [int(x) for x in s.manifold.f_vector]
        dS = s.current_objective - S0
        df = (f1[1] - f0[1], f1[3] - f0[3])
        for cen, coc in reversed(done):
            L_.do(coc, cen)
        assert abs(s.current_objective - S0) < 1e-9, "undo not exact"
        return dS, df, None

    def local(self, dS, df):
        """Strip ALL pin content: local = dS - [g(f+df) - g(f)]."""
        f = [int(x) for x in self.s.manifold.f_vector]
        return dS - (g(f[1] + df[0], f[3] + df[1]) - g(f[1], f[3]))


def main():
    t0 = time.time()
    m = ddg.Manifold.load(MFD, 3)
    s, L = F.fresh(m)
    fv = [int(x) for x in s.manifold.f_vector]
    print(f"manifold {os.path.basename(MFD)}  f={fv}  "
          f"gap {fv[1] - 6*fv[3]/F.ETARGET:+.3f}  S={s.current_objective:.3f}")

    # --- BC chain -------------------------------------------------------
    chain = None
    for v0 in sorted(L.v2t)[:20]:
        for seed in sorted(L.v2t[v0])[:4]:
            chain = bc_chain(L, list(seed))
            if chain and len(chain) > 200:
                break
        if chain and len(chain) > 200:
            break
    assert chain, "no closing BC chain found"
    Lc = len(chain)
    print(f"BC chain: length {Lc} (seed tet {seed})")

    M = Meas(s, L)
    rows = []
    jA = 0
    while knot_op(L, chain, jA) is None and jA < 40:
        jA += 1
    opA = knot_op(L, chain, jA)
    assert opA, "no flippable chain window found for A"
    dSA_solo, dfA, _ = M.run([opA])
    print(f"\nhead A at window {jA}: knot dS {dSA_solo:+.4f} df {dfA} "
          f"(local {M.local(dSA_solo, dfA):+.4f})")

    # graph distance from A's window, computed once
    winA = [chain[(jA + i) % Lc] for i in range(5)]
    distA = bfs_from(L, winA)

    def dgraph(jB):
        w = [chain[(jB + i) % Lc] for i in range(5)]
        return min((distA.get(v, 99) for v in w), default=99)

    # === T1: non-conserving pair (knot + knot) ==========================
    print("\n=== T1  knot(A) + knot(B): residual after subtracting the "
          "analytic cross term ===")
    print(f"{'s':>4} {'dgrf':>5} {'dS_A':>9} {'dS_B':>9} {'dS_AB':>10} "
          f"{'X':>8} {'resid':>11} {'resid_loc':>11}")
    for sep in SEPS:
        jB = (jA + sep) % Lc
        opB = knot_op(L, chain, jB)
        if opB is None:
            rows.append({"test": "T1", "s": sep, "status": "no-op-B"})
            print(f"{sep:>4}  (window not flippable)")
            continue
        dA, dfA_, eA = M.run([opA])
        dB, dfB_, eB = M.run([opB])
        dAB, dfAB, eAB = M.run([opA, opB])
        if eAB is not None or dAB is None:
            rows.append({"test": "T1", "s": sep, "status": "blocked",
                         "err": eAB})
            print(f"{sep:>4}  BLOCKED (joint move invalid): {eAB}")
            continue
        X = cross(dfA_, dfB_)
        resid = dAB - dA - dB - X
        rloc = (M.local(dAB, dfAB) - M.local(dA, dfA_)
                - M.local(dB, dfB_))
        rows.append({"test": "T1", "s": sep, "dgraph": dgraph(jB),
                     "dS_A": dA, "dS_B": dB,
                     "dS_AB": dAB, "X": X, "resid": resid,
                     "resid_local": rloc, "dfA": dfA_, "dfB": dfB_,
                     "status": "ok"})
        print(f"{sep:>4} {dgraph(jB):>5} {dA:>9.4f} {dB:>9.4f} "
              f"{dAB:>10.4f} {X:>8.4f} {resid:>11.3e} {rloc:>11.3e}")

    # === T2: conserving pair (create A, annihilate B) ===================
    print("\n=== T2  knot(A) + UNknot(B), dA+dB = 0: the pair should be "
          "PIN-FREE at every s ===")
    print(f"{'s':>4} {'dgrf':>5} {'dS_A':>9} {'dS_B':>9} {'dS_AB':>10} "
          f"{'pin_AB':>8} {'X':>8} {'resid':>11}")
    for sep in SEPS:
        jB = (jA + sep) % Lc
        opB = knot_op(L, chain, jB)
        if opB is None:
            continue
        # base state: knot already at B
        cenB, cocB = opB
        Spre = s.current_objective
        L.do(cenB, cocB)
        unB = (cocB, cenB)                      # the 3->2 that removes it
        dA, dfA_, eA = M.run([opA])
        dB, dfB_, eB = M.run([unB])
        dAB, dfAB, eAB = M.run([opA, unB])
        ok = None not in (dA, dB, dAB)
        if ok:
            X = cross(dfA_, dfB_)
            resid = dAB - dA - dB - X
            pin_AB = dAB - M.local(dAB, dfAB)   # pin content of the pair
            rows.append({"test": "T2", "s": sep, "dgraph": dgraph(jB),
                         "dS_A": dA, "dS_B": dB,
                         "dS_AB": dAB, "pin_AB": pin_AB, "X": X,
                         "resid": resid, "dfAB": dfAB, "status": "ok"})
            print(f"{sep:>4} {dgraph(jB):>5} {dA:>9.4f} {dB:>9.4f} "
                  f"{dAB:>10.4f} {pin_AB:>8.1e} {X:>8.4f} {resid:>11.3e}")
        else:
            rows.append({"test": "T2", "s": sep, "status": "blocked"})
            print(f"{sep:>4}  BLOCKED")
        L.do(cocB, cenB)                        # restore
        assert abs(s.current_objective - Spre) < 1e-9
    json.dump({"mfd": MFD, "jA": jA, "chain_len": Lc, "rows": rows},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}  ({time.time() - t0:.0f}s)")


main()
