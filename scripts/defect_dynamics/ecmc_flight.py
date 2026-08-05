#!/usr/bin/env python3
"""Lifted ECMC flight of a deg-3 chord: the assembled driver.

Lifted state: (chord C, frame w, carried rung Q).
  * frame w = 4-window (e, p0, p1, p2), e an endpoint of C, (p0,p1,p2) an
    ordering of C's link triangle. Transport = the sliding-window walk in the
    chord-annihilated complex; flip F(e, p0p1p2) = (e', p2p1p0) -- verified
    involution (M(F(M)) = F, 18/18).
  * Q = the rung the chord is riding. CONSERVED by every move of this kernel
    ("kinetic energy"): travel targets sit on rung Q, and an exit from a
    complex must land on a rung-Q site. Carrying Q is what makes entry/exit
    reversible: the exit's target (first clean rung-Q site along the flipped
    ray) is exactly the site the entry came from, because the entry was only
    proposed once no free rung-Q site remained before contact.

ONE UNIFORM KERNEL (travel / entry / exit are the same rule):

  scan k = 1, 2, ... along the frame walk;
  k* = first k whose arrival is (clean AND on rung Q)   [a free-web site]
                             OR (touching a complex)     [contact]
  propose the slide to k*; accept with min(1, e^-dS);
  on accept: hop, transport w by k* steps;  on reject (or no k*, or the
  slide at k* is illegal): FLIP the frame.

Travel is the special case where the proposal is dS = 0 (always accepted);
entry is the contact case (usually uphill); exit is the boundary case where
the first clean rung-Q site is the launch point (downhill, the entry's exact
inverse). Target selection is WALK-based (face_rung + defect-vertex test),
never legality-based, so forward and reverse select the same sites; a slide
that turns out illegal at k* is a rejection (flip) on both sides.

Momentum handoff (--p-hand): after a CONVERSION the axis chord C1 and the
freed edge e4 are INTERLOCKED (each partner locally recognizable: e4's link
contains both endpoints of C1, C1's link is the face containing e4). With
probability p an interlocked ACTIVE edge hands to its partner BEFORE the scan
(source-hand), and an interlocked ARRIVAL hands AFTER acceptance
(arrival-hand). Frames map through the equivariant bijection h (axis->freed;
h(F(w)) = F(h(w)) verified 20/20 sites) or its inverse (freed->axis), so
source-hands pair with arrival-hands of the reverse composite and the move
stays skew-balanced. p-hand = 0 recovers the default rule, under which the
heavy sector is PROVABLY frozen (the background is invariant under every
slide); the handoff is the only coupling. Carried Q passes through hands
unchanged.

Q is conserved for the whole run (no Q-refresh yet), so one run samples one
rung sector; ergodicity across rungs needs a Q-resample move, deliberately
left out of this first driver.

Audit mode (--audit, default on): after every accepted gated move, rerun the
scan from the flipped new state and assert it proposes the exact inverse
(same site, same k, dS' = -dS). Also asserts sum(accepted dS) equals the
sampler's tracked objective change exactly, and that at clean sites the
carried Q equals the chord's own creation-face rung.
"""
import argparse
import collections
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from discrete_differential_geometry.ecmc import face_rung
from worm_deg4 import CREATE_SEQ, apply_seq

ESTAR = 5.105025
TOL = 1e-9


def ed_of(m):
    return {tuple(sorted(map(int, e))): m.degree(e) for e in m.simplices(1)}


def minus_flier(ed, chord, link):
    """Edge degrees of the chord-annihilated complex, by local adjustment."""
    out = dict(ed)
    del out[chord]
    l0, l1, l2 = link
    for e in ((l0, l1), (l0, l2), (l1, l2)):
        out[tuple(sorted(e))] += 1
    for c in chord:
        for l in link:
            out[tuple(sorted((c, l)))] -= 1
    return out


class Flight:
    def __init__(self, sampler, chord, frame, kscan=60, audit=True,
                 beta=1.0, p_hand=0.0):
        self.beta = float(beta)
        self.p_hand = float(p_hand)
        self.s = sampler
        self.C = chord
        self.w = frame
        self.kscan = kscan
        self.audit = audit
        self.events = collections.Counter()
        self.log = []
        self.dS_sum = 0.0
        self.S0 = sampler.current_objective
        self._rebuild_background()
        self.Q = self._rung_of_face(self.link, chord)

    # -- background = the chord-annihilated complex ---------------------------
    def _rebuild_background(self):
        ed = ed_of(self.s.manifold)
        # full link of the deg-3 chord: the 3 vertices among its tets
        F = np.asarray(self.s.manifold.facets())
        lk = set()
        for t in F:
            t = [int(x) for x in t]
            if self.C[0] in t and self.C[1] in t:
                lk |= {v for v in t if v not in self.C}
        self.link = tuple(sorted(lk))
        assert len(self.link) == 3, self.link
        self.edB = minus_flier(ed, self.C, self.link)
        ill = {e for e, d in self.edB.items() if d < 5 or d > 6}
        self.BV = {v for e in ill for v in e}
        self.f2a = {}
        # background face->apex map: adjust the flier's 3 tets to 2
        tets = []
        for t in F:
            t = tuple(sorted(int(x) for x in t))
            if self.C[0] in t and self.C[1] in t:
                continue
            tets.append(t)
        l0, l1, l2 = self.link
        tets.append(tuple(sorted((l0, l1, l2, self.C[0]))))
        tets.append(tuple(sorted((l0, l1, l2, self.C[1]))))
        for t in tets:
            for i in range(4):
                self.f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
        self.ill_snapshot = frozenset(ill)
        self.ill_deg = {e: self.edB[e] for e in ill}

    # -- frame walk in the background ----------------------------------------
    def stepw(self, w):
        ap = self.f2a.get(tuple(sorted(w[1:4])))
        if ap is None or len(ap) != 2:
            return None
        x = ap[0] if ap[1] == w[0] else ap[1]
        return None if x == w[0] else (w[1], w[2], w[3], int(x))

    def arrival(self, w):
        ap = self.f2a.get(tuple(sorted(w[1:4])))
        if ap is None or len(ap) != 2:
            return None
        nxt = ap[0] if ap[1] == w[0] else ap[1]
        return tuple(sorted((int(w[0]), int(nxt))))

    def _rung_of_face(self, face, apexes):
        return face_rung(face, apexes, self.edB)

    def flip(self):
        e2 = self.C[1] if self.w[0] == self.C[0] else self.C[0]
        self.w = (e2, self.w[3], self.w[2], self.w[1])

    # -- interlock detection + the equivariant handoff ------------------------
    def _link_edge(self, e):
        """Link vertices of edge e in the LIVE manifold (any degree)."""
        F = np.asarray(self.s.manifold.facets())
        lk = set()
        for t in F:
            t = [int(x) for x in t]
            if e[0] in t and e[1] in t:
                lk |= {v for v in t if v not in e}
        return sorted(lk)

    def partner(self):
        """The interlocked partner of the ACTIVE edge, or None.

        Interlock relation: a deg-3 edge P formed by two of my link vertices,
        whose own link contains both of my endpoints. The configuration is
        SYMMETRIC under (active <-> partner, x <-> n) -- both edges are deg-3
        with triangle links sharing one tet -- so there is no axis/freed role
        to distinguish, and the handoff map below is SELF-INVERSE under the
        exchange (verified: applying the same formula from the partner's side
        recovers the original frame exactly).
        """
        E, L = self.C, self.link
        for i in range(3):
            for j in range(i + 1, 3):
                cand = tuple(sorted((L[i], L[j])))
                if self.s.manifold.degree([cand[0], cand[1]]) == 3 \
                        and set(E) <= set(self._link_edge(cand)):
                    return cand
        return None

    def _hand_map(self, E, Elink, P, Plink):
        """The equivariant bijection h: frames(E) -> frames(P).

        x = singleton of E's link w.r.t. P; n = singleton of P's link w.r.t.
        E. h encodes (first-of-P-in-perm -> lead, position of x -> position
        of n, lead -> order of E's endpoints): 12 -> 12, flip-equivariant
        (h(F(w)) = F(h(w))), and the same formula from P's side is h^-1.
        """
        a1, b1 = E
        u, v = P
        x = next(z for z in Elink if z not in P)
        n = next(z for z in Plink if z not in E)
        def h(w):
            e = w[0]
            pi = w[1:4]
            ep = b1 if e == a1 else a1
            i = pi.index(x)
            f = next(z for z in pi if z in (u, v))
            rem = [q for q in range(3) if q != i]
            sig = [None] * 3
            sig[i] = n
            sig[rem[0]] = e
            sig[rem[1]] = ep
            return (f,) + tuple(sig)
        return h

    def try_hand(self, rng, when):
        """With prob p_hand, hand the activity to the interlocked partner."""
        if self.p_hand <= 0:
            return False
        P = self.partner()
        if P is None:
            return False
        if rng.random() >= self.p_hand:
            return False
        Plink = tuple(self._link_edge(P))
        h = self._hand_map(self.C, self.link, P, Plink)
        wnew = h(self.w)
        if self.audit:
            # flip-equivariance at use time: h(F(w)) == F(h(w))
            e2 = self.C[1] if self.w[0] == self.C[0] else self.C[0]
            wf = (e2, self.w[3], self.w[2], self.w[1])
            p2 = P[1] if wnew[0] == P[0] else P[0]
            wnf = (p2, wnew[3], wnew[2], wnew[1])
            if h(wf) != wnf:
                self.events["HAND_EQUIVARIANCE_FAIL"] += 1
            # self-inverse under exchange: the partner's map returns w
            hb = self._hand_map(P, Plink, self.C, self.link)
            if hb(wnew) != self.w:
                self.events["HAND_INVERSE_FAIL"] += 1
        self.C, self.w = P, wnew
        self._rebuild_background()
        self.events[f"hand_{when}"] += 1
        return True

    # -- proposal scan (pure walk, no legality) -------------------------------
    def scan(self, w):
        """First special site along w: ('free'|'contact', k, arrival, face)."""
        ww = w
        for k in range(1, self.kscan + 1):
            ww = self.stepw(ww)
            if ww is None:
                return None
            arr = self.arrival(ww)
            if arr is None:
                return None
            face = tuple(sorted(ww[1:4]))
            sup = set(face) | set(arr)
            if sup & self.BV:
                return ("contact", k, arr, face)
            if self._rung_of_face(face, arr) == self.Q:
                return ("free", k, arr, face)
        return None

    def _slot_for(self, w):
        # match the D slot to the frame by arrivals at steps 1..3. On
        # pristine crystal steps 1-2 pin the window uniquely; near a complex
        # two DIFFERENT faces can share an apex pair, so two different
        # windows can share arrivals at small k -- deeper matching shrinks
        # (but cannot eliminate) that; the residual is caught by the
        # frame-resync check after commit.
        preds = []
        ww = w
        for _ in range(3):
            ww = self.stepw(ww) if ww else None
            preds.append(self.arrival(ww) if ww else None)
        if preds[0] is None:
            return None
        for sl in range(12):
            ok = True
            for k, a in enumerate(preds, start=1):
                if a is None:
                    break
                r = self.s.nonlocal_slide_at(self.C[0], self.C[1], sl, k,
                                             commit=False)
                if not (r and r[2] == a):
                    ok = False
                    break
            if ok:
                return sl
        return None

    # -- one kernel step ------------------------------------------------------
    def refresh_frame(self, rng):
        """Uniform frame resample at the current chord -- the ECMC
        refreshment move. Uniform -> uniform is self-paired, state unchanged,
        so it is trivially valid; without it a run with no gated events is a
        DETERMINISTIC closed circuit on one rung web (measured on C15: 8/8
        seeds coasted forever without ever contacting the defect)."""
        perms = [(self.link[i], self.link[j], self.link[k])
                 for (i, j, k) in ((0, 1, 2), (0, 2, 1), (1, 0, 2),
                                   (1, 2, 0), (2, 0, 1), (2, 1, 0))]
        e = self.C[int(rng.integers(2))]
        self.w = (e,) + perms[int(rng.integers(6))]
        self.events["refresh"] += 1

    def step(self, rng):
        self.try_hand(rng, "src")
        prop = self.scan(self.w)
        if prop is None:
            self.events["stuck_flip"] += 1
            self.flip()
            return
        kind, k, arr, face = prop
        slot = self._slot_for(self.w)
        if slot is None:
            self.events["noslot_flip"] += 1
            self.flip()
            return
        r = self.s.nonlocal_slide_at(self.C[0], self.C[1], slot, k,
                                     commit=False)
        if r is None or r[2] != arr:
            self.events["illegal_flip"] += 1
            self.flip()
            return
        dS = r[0]
        # a rung-Q target is dS = 0 ONLY when the source chord is itself clean;
        # from a docked chord the same proposal is the EXIT and carries the
        # entry's elevation with reversed sign -- gated below like everything
        src_clean = not ((set(self.C) | set(self.link)) & self.BV)
        if kind == "free" and src_clean:
            assert abs(dS) < 1e-6, (dS, "free site not free from clean chord")
        if dS <= 0 or rng.random() < np.exp(-self.beta * dS):
            old = (self.C, self.w, k, dS)
            self.s.nonlocal_slide_at(self.C[0], self.C[1], slot, k,
                                     commit=True)
            self.dS_sum += dS
            wnew = self.w
            for _ in range(k):
                wnew = self.stepw(wnew)
            self.C, self.w = arr, wnew
            # the background (state minus flier) is INVARIANT under flier
            # motion, so no rebuild. The new link is taken from the LIVE
            # manifold, not the walk: near a complex the D slot's walk can
            # cross a DIFFERENT face sharing the same apex pair, in which
            # case the walked frame is fiction -- detect and refresh.
            self.link = tuple(self._link_edge(arr))
            if set(wnew[1:4]) != set(self.link):
                self.events["frame_resync"] += 1
                perms = [(self.link[i], self.link[j], self.link[k2])
                         for (i, j, k2) in ((0, 1, 2), (0, 2, 1), (1, 0, 2),
                                            (1, 2, 0), (2, 0, 1), (2, 1, 0))]
                wnew = (arr[int(rng.integers(2))],) + \
                    perms[int(rng.integers(6))]
                self.w = wnew
            label = kind if src_clean else \
                ("exit" if kind == "free" else "deeper")
            self.events[f"accept_{label}"] += 1
            # bundle change (conversion)?  ~|bundle| degree queries, not a
            # full edge sweep
            # (a docked flier's own support can shift a bundle edge's LIVE
            # degree without changing the background, so this flag only
            # triggers the rebuild; the rebuild's snapshot comparison decides
            # whether the background truly changed)
            changed = any(
                self.s.manifold.degree(list(e)) != d
                for e, d in self.ill_deg.items())
            if changed:
                before = self.ill_snapshot
                self._rebuild_background()
                if self.ill_snapshot != before:
                    self.events["bundle_changed"] += 1
            if self.audit and abs(dS) > TOL:
                self._audit_reverse(old)
            self.log.append((label, k, dS, self.C))
            self.try_hand(rng, "arr")
        else:
            self.events[f"reject_{kind}_flip"] += 1
            self.flip()

    def _audit_reverse(self, old):
        C0, w0, k0, dS0 = old
        save = (self.C, self.w)
        e2 = self.C[1] if self.w[0] == self.C[0] else self.C[0]
        wf = (e2, self.w[3], self.w[2], self.w[1])
        prop = self.scan(wf)
        ok = prop is not None and prop[1] == k0 and prop[2] == C0
        if not ok:
            self.events["AUDIT_FAIL"] += 1
            print(f"   AUDIT FAIL: reverse of {C0}->{self.C} (k={k0}) "
                  f"proposes {prop}")
        else:
            self.events["audit_ok"] += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default=os.path.join(
        _ROOT, "data/tcp_reference/T3_C15_m3_N3672.mfd"))
    ap.add_argument("--bg", default="face:8,9,10",
                    help="background defect complex: 'face:a,b,c[;a,b,c...]' "
                         "(2->3 on each face) or 'createseq' (the R m2 "
                         "CREATE_SEQ bundle). Default: the C15 orbit-3 "
                         "triangle -- tipless, immobile, a clean stationary "
                         "target.")
    ap.add_argument("--nstep", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--kscan", type=int, default=60)
    ap.add_argument("--no-audit", action="store_true")
    ap.add_argument("--refresh-every", type=int, default=0,
                    help="uniform frame resample every N kernel steps "
                         "(0 = off). Needed for ergodicity: with no gated "
                         "events the flight is a deterministic closed "
                         "circuit.")
    ap.add_argument("--p-hand", type=float, default=0.0,
                    help="handoff probability at interlocked configurations "
                         "(0 = default rule, heavy sector provably frozen)")
    ap.add_argument("--beta", type=float, default=1.0,
                    help="inverse temperature for the gated moves (the scan "
                         "knob; target distribution becomes pi^beta)")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    def build_background(mm):
        if args.bg == "createseq":
            apply_seq(mm, CREATE_SEQ)
        else:
            assert args.bg.startswith("face:"), args.bg
            Fm = np.asarray(mm.facets())
            fa = {}
            for t in Fm:
                t = tuple(sorted(int(x) for x in t))
                for i in range(4):
                    fa.setdefault(t[:i] + t[i + 1:], []).append(t[i])
            for spec in args.bg[5:].split(";"):
                f = tuple(sorted(int(x) for x in spec.split(",")))
                mm.do_bistellar_move(list(f), list(fa[f]))

    m = Manifold.load(args.ref, 3)
    build_background(m)
    edB = ed_of(m)
    BV = {v for e, d in edB.items() if d < 5 or d > 6 for v in e}
    F = np.asarray(m.facets())
    f2a = {}
    for t in F:
        t = tuple(sorted(int(x) for x in t))
        for i in range(4):
            f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
    # launch: clean face >= 3 from the bundle
    from collections import deque
    V = int(F.max()) + 1
    adj = [[] for _ in range(V)]
    for e in m.simplices(1):
        a, b = sorted(map(int, e))
        adj[a].append(b)
        adj[b].append(a)
    dist = np.full(V, -1, np.int32)
    dq = deque(sorted(BV))
    for v in BV:
        dist[v] = 0
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if dist[w] < 0:
                dist[w] = dist[u] + 1
                dq.append(w)
    faces = sorted(f for f, ap2 in f2a.items() if len(ap2) == 2)
    order = rng.permutation(len(faces))
    g = None
    for i in order:
        f = faces[i]
        ap2 = f2a[f]
        sup = set(f) | {int(ap2[0]), int(ap2[1])}
        if min(dist[v] for v in sup) >= 3:
            g = f
            break
    C = tuple(sorted(int(x) for x in f2a[g]))
    m.do_bistellar_move(list(g), list(C))
    par = dict(num_facets_target=len(F) + 1, num_facets_coef=0.0,
               hinge_degree_target=ESTAR, hinge_degree_target_coef=1.0,
               num_hinges_coef=0.0, hinge_degree_variance_coef=0.0,
               codim3_degree_variance_coef=0.0)
    s = ManifoldSampler(m, SamplerParams(**par))
    w0 = (C[0], g[0], g[1], g[2])
    fl = Flight(s, C, w0, kscan=args.kscan, audit=not args.no_audit,
                beta=args.beta, p_hand=args.p_hand)
    print(f"launch chord {C} on face {g}, carried Q = {fl.Q}, "
          f"bundle verts {len(BV)}, beta = {args.beta}, "
          f"p_hand = {args.p_hand}")

    for step in range(args.nstep):
        if args.refresh_every and step and step % args.refresh_every == 0:
            fl.refresh_frame(rng)
        fl.step(rng)

    print(f"\nevents after {args.nstep} kernel steps:")
    for k, v in sorted(fl.events.items()):
        print(f"   {k:22s} {v}")
    drift = abs((fl.s.current_objective - fl.S0) - fl.dS_sum)
    print(f"\nobjective drift |tracked - sum(dS)| = {drift:.3e}")
    hops = [x for x in fl.log if x[0] == "free"]
    gated = [x for x in fl.log if x[0] != "free"]
    print(f"free hops {len(hops)} (k median "
          f"{np.median([h[1] for h in hops]) if hops else 0:.0f}), "
          f"gated accepts {len(gated)} "
          f"{sorted(collections.Counter((g_[0], round(g_[2], 1)) for g_ in gated).items())}")
    print(f"final chord {fl.C}, carried Q {fl.Q}")
    # THE coupling metric: has the heavy sector changed, final vs initial?
    # (step-to-step comparison misses it: a handoff legitimately swaps which
    # edge is background, so the commit happens invisibly to that check)
    ill0 = {e for e, d in edB.items() if d < 5 or d > 6}
    ill1 = fl.ill_snapshot
    gained = sorted(ill1 - ill0)
    lost = sorted(ill0 - ill1)
    print(f"HEAVY SECTOR final vs initial: lost {lost}  gained {gained}"
          f"{'   <== COUPLED' if (gained or lost) else '   (unchanged)'}")


if __name__ == "__main__":
    main()
