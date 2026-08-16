#!/usr/bin/env python3
"""The deg-4 WORM move (deg-3-catalyzed transport step) -- reference
implementation / oracle, built along the same lines as worm_slide.py.

MOVE CLASS (one proposal = one 2->3 + one 3->2, f-vector neutral):
  * ANCHOR: a deg-4 edge e whose degree becomes legal (4->5/4->6) or whose
    edge is removed (4->0).  Proposal draws e uniformly from the deg-4 edge
    set; the move preserves that set's size, so the 1/n4 factors cancel.
  * LANDING: exactly ONE new deg-4 edge f, outside e's hinge-flip
    octahedron in the start state AND with e outside f's octahedron in the
    end state (escape symmetry, so the reverse is in the class).
  * CATALYST: the deg-3 edge->degree multiset is unchanged, or exactly one
    deg-3 removed and one created (the catalyst relocating).
  * NOTHING ELSE: the full illegal edge->degree MAP is unchanged outside
    {e, f, catalyst pair} (map, not set -- degree changes among retained
    illegal edges disqualify).

The anchor is determined by the transition (the unique healed deg-4), so a
given end state is reachable only through candidates at one anchor; the
Hastings ratio is

    alpha = min(1, exp(-lam*dS) * (k_r * n_f) / (k_f * n_r)),

with n = candidate count at the anchor, k = multiplicity of the end state
among them (by exact net facet diff).  dS is the FULL shape action (EDQ +
zleg*U(n6) + cimp*m^2) in lam=1 units.

Modes:
  --closure-test SNAP   exhaustive inverse-closure + k_f/k_r audit over a
                        subset of anchors (detailed-balance foundation)
  --walk SNAP           directed walk regression (worm_walk behaviour
                        through the class machinery)
FRAME: none -- graph identity only.
"""
import argparse
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import discrete_differential_geometry as ddg
from . import tip_retract_search as T

R_BALL, R_CORE = T.R_BALL, T.R_CORE


def _legal(d):
    return d in (5, 6)


class Live:
    """Incrementally-maintained global state alongside the Manifold: tet
    set, vertex->tets index, edge degrees.  Lets candidate enumeration
    build local patches in O(ball) instead of O(N)."""

    def __init__(self, m, driver=None):
        self.m = m
        self.driver = driver          # callable(cen, coc); default: manifold
        self.v2t = defaultdict(set)
        self.edeg = Counter()
        for t in (tuple(sorted(int(x) for x in f)) for f in
                  np.asarray(m.facets())):
            self._add(t)

    def _add(self, t):
        for v in t:
            self.v2t[v].add(t)
        for e in combinations(t, 2):
            self.edeg[e] += 1

    def _rm(self, t):
        for v in t:
            self.v2t[v].discard(t)
        for e in combinations(t, 2):
            self.edeg[e] -= 1
            if self.edeg[e] == 0:
                del self.edeg[e]

    def do(self, cen, coc):
        """Apply a bistellar move through the manifold AND the index."""
        if self.driver is not None:
            self.driver(list(cen), list(coc))
        else:
            self.m.do_bistellar_move(list(cen), list(coc))
        if len(cen) == 3:                       # 2->3
            a, b, c = cen
            x, y = coc
            for t in ((a, b, c, x), (a, b, c, y)):
                self._rm(tuple(sorted(t)))
            for t in ((x, y, a, b), (x, y, b, c), (x, y, a, c)):
                self._add(tuple(sorted(t)))
        elif len(cen) == 2:                      # 3->2
            u, w = cen
            a, b, c = coc
            for t in ((u, w, a, b), (u, w, b, c), (u, w, a, c)):
                self._rm(tuple(sorted(t)))
            for t in ((a, b, c, u), (a, b, c, w)):
                self._add(tuple(sorted(t)))
        elif len(cen) == 4:                      # 1->4 (coc = [new label])
            (v,) = coc
            self._rm(tuple(sorted(cen)))
            for f in combinations(cen, 3):
                self._add(tuple(sorted(f + (v,))))
        else:                                    # 4->1 (cen = [vertex])
            (v,) = cen
            for f in combinations(coc, 3):
                self._rm(tuple(sorted(f + (v,))))
            self._add(tuple(sorted(coc)))

    def ball(self, seeds, R):
        seen = set(seeds)
        frontier = set(seeds)
        for _ in range(R):
            nxt = set()
            for v in frontier:
                for t in self.v2t[v]:
                    nxt |= set(t)
            nxt -= seen
            seen |= nxt
            frontier = nxt
        return seen

    def region(self, Vball):
        out = set()
        for v in Vball:
            for t in self.v2t[v]:
                if set(t) <= Vball:
                    out.add(frozenset(t))
        return out

    def deg4(self):
        return [e for e, d in self.edeg.items() if d == 4]


def as_bistellar(mv):
    if mv[0] == "23":
        _, tri, x, y = mv
        return tuple(tri), (x, y)
    _, (u, w), a, b, c = mv
    return (u, w), (a, b, c)


def net_key(mv1, mv2):
    """Canonical net facet diff of the composite (cancellation-aware)."""
    net = Counter()

    def acc(mv, sgn=1):
        if mv[0] == "23":
            _, (a, b, c), x, y = mv[0], mv[1], mv[2], mv[3]
            net[tuple(sorted((a, b, c, x)))] -= sgn
            net[tuple(sorted((a, b, c, y)))] -= sgn
            for p, q in ((a, b), (b, c), (a, c)):
                net[tuple(sorted((x, y, p, q)))] += sgn
        else:
            _, (u, w), a, b, c = mv
            for p, q in ((a, b), (b, c), (a, c)):
                net[tuple(sorted((u, w, p, q)))] -= sgn
            net[tuple(sorted((a, b, c, u)))] += sgn
            net[tuple(sorted((a, b, c, w)))] += sgn
    acc(mv1)
    acc(mv2)
    return frozenset((t, c) for t, c in net.items() if c != 0)


def candidates(L, anchor, max_ds=None):
    """All class candidates at `anchor` in the live state.

    Returns [(mv1, mv2, dS, landing, key), ...].  Enumeration mirrors the
    pair census: first move touching the anchor, second the disturbed
    region; class filters applied exactly as documented above."""
    Vball = L.ball(anchor, R_BALL)
    Vcore = frozenset(L.ball(anchor, R_CORE))
    P = T.Patch(L.region(Vball))
    init_edeg = dict(P.edeg)
    if init_edeg.get(anchor, 0) != 4:
        return []
    init_ill = {e: d for e, d in init_edeg.items() if not _legal(d)}
    octa_a = set(anchor)
    for t in P.tets:
        if set(anchor) <= t:
            octa_a |= t

    def disturbed():
        f = set(anchor)
        for ed in P.touched:
            if P.edeg.get(ed, 0) != init_edeg.get(ed, 0):
                f |= set(ed)
        return f

    def apply_mv(mv):
        if mv[0] == "23":
            _, tri, x, y = mv
            P.apply_23(tri, x, y)
        else:
            _, (u, w), a, b, c = mv
            P.apply_32(u, w, a, b, c)

    def undo_mv(mv):
        if mv[0] == "23":
            _, tri, x, y = mv
            P.undo_23(tri, x, y)
        else:
            _, (u, w), a, b, c = mv
            P.undo_32(u, w, a, b, c)

    out = []
    for m1 in P.moves(Vcore, set(anchor)):
        apply_mv(m1)
        for m2 in P.moves(Vcore, disturbed()):
            if m1[0] == m2[0]:
                undo_mv_needed = False          # need one 23 + one 32
            else:
                undo_mv_needed = True
            if not undo_mv_needed:
                continue
            apply_mv(m2)
            ok = False
            ad = P.edeg.get(anchor, 0)
            if ad in (0, 5, 6):                 # anchor healed or removed
                cur_ill = {e: d for e, d in P.edeg.items()
                           if not _legal(d)}
                gone4 = [e for e, d in init_ill.items()
                         if d == 4 and cur_ill.get(e) != 4]
                new4 = [e for e, d in cur_ill.items()
                        if d == 4 and init_ill.get(e) != 4]
                gone3 = [e for e, d in init_ill.items()
                         if d == 3 and cur_ill.get(e) != 3]
                new3 = [e for e, d in cur_ill.items()
                        if d == 3 and init_ill.get(e) != 3]
                others_ok = all(
                    cur_ill.get(e) == d for e, d in init_ill.items()
                    if e not in gone4 + gone3) and all(
                    init_ill.get(e) == d for e, d in cur_ill.items()
                    if e not in new4 + new3)
                if (anchor in gone4
                        and 1 <= len(gone4) == len(new4) <= 2
                        and len(gone3) == len(new3) <= 1 and others_ok):
                    # symmetric escape: some landing f outside the anchor's
                    # octahedron AND anchor outside f's octahedron in x'
                    land = None
                    for f in new4:
                        if set(f) <= octa_a:
                            continue
                        octa_f = set(f)
                        for t in P.tets:
                            if set(f) <= t:
                                octa_f |= t
                        if not set(anchor) <= octa_f:
                            land = f
                            break
                    if land is not None:
                        ds = T.LAM * T.dS_full(
                            init_edeg, P.edeg, frozenset(P.touched))
                        if max_ds is None or ds <= max_ds:
                            out.append((m1, m2, ds, land,
                                        net_key(m1, m2), tuple(gone4),
                                        tuple(new4)))
            undo_mv(m2)
        undo_mv(m1)
    return out


def apply_pair(L, m1, m2):
    for mv in (m1, m2):
        cen, coc = as_bistellar(mv)
        L.do(cen, coc)


def undo_pair(L, m1, m2):
    for mv in (m2, m1):
        cen, coc = as_bistellar(mv)
        L.do(coc, cen)


def closure_test(args):
    """For a subset of anchors: every candidate transition must have its
    inverse in the class at the landing (k_r >= 1), with sane k_f/k_r."""
    m = ddg.Manifold.load(args.snap, 3)
    L = Live(m)
    d4 = sorted(L.deg4())
    rng = np.random.default_rng(args.seed)
    anchors = [d4[i] for i in
               rng.choice(len(d4), size=min(args.n_anchors, len(d4)),
                          replace=False)]
    print(f"{os.path.basename(args.snap)}: {len(d4)} deg-4 edges; "
          f"closure test on {len(anchors)} anchors")
    n_trans = n_closed = 0
    kf_h = Counter()
    kr_h = Counter()
    fails = []
    for ai, e in enumerate(anchors):
        cands = candidates(L, e)
        bykey = defaultdict(list)
        for c in cands:
            bykey[c[4]].append(c)
        for key, group in bykey.items():
            m1, m2, ds, f, _, gone4, new4 = group[0]
            n_trans += 1
            # forward proposal weight: sum over ALL anchors in gone4
            qf = 0.0
            k_f_tot = 0
            for a in gone4:
                ca = candidates(L, a)
                k = sum(1 for c in ca if c[4] == key)
                k_f_tot += k
                if ca:
                    qf += k / len(ca)
            apply_pair(L, m1, m2)
            inv = frozenset((t, -c) for t, c in key)
            qr = 0.0
            k_r_tot = 0
            ds_ok = True
            for b in new4:
                cb = candidates(L, b)
                k = sum(1 for c in cb if c[4] == inv)
                k_r_tot += k
                if cb:
                    qr += k / len(cb)
                ds_ok = ds_ok and all(abs(c[2] + ds) < 1e-9
                                      for c in cb if c[4] == inv)
            undo_pair(L, m1, m2)
            kf_h[k_f_tot] += 1
            kr_h[k_r_tot] += 1
            if qr > 0 and ds_ok:
                n_closed += 1
            else:
                fails.append((e, f, k_f_tot, k_r_tot, ds_ok))
        print(f"  anchor {ai + 1}/{len(anchors)} {e}: "
              f"{len(cands)} candidates, {len(bykey)} transitions",
              flush=True)
    print(f"\ntransitions: {n_trans}   closed (k_r>=1, dS antisymmetric): "
          f"{n_closed}")
    print(f"k_f histogram: {dict(sorted(kf_h.items()))}")
    print(f"k_r histogram: {dict(sorted(kr_h.items()))}")
    for x in fails[:10]:
        print(f"  FAIL: {x}")
    if n_trans and n_closed == n_trans:
        print("PASS: the class is inverse-closed with antisymmetric dS -- "
              "Hastings ratio well-defined")
    else:
        print("FAIL")
        sys.exit(1)


def walk(args):
    """Directed-walk regression through the class machinery."""
    m = ddg.Manifold.load(args.snap, 3)
    L = Live(m)
    head = tuple(args.head)
    comp0 = Counter(d for d in L.edeg.values() if not _legal(d))
    total_ds = 0.0
    applied = []
    for k in range(args.steps):
        cands = candidates(L, head, max_ds=args.max_ds)
        if not cands:
            print(f"step {k + 1}: no candidates at {head}; stop")
            break
        prev = applied[-1][2] if applied else None
        pool = [c for c in cands if c[3] != prev] or cands
        m1, m2, ds, f = min(pool, key=lambda c: c[2])[:4]
        apply_pair(L, m1, m2)
        applied.append((m1, m2, head))
        total_ds += ds
        print(f"step {k + 1}: {head} -> {f}  dS={ds:+.3f} "
              f"({len(cands)} candidates)")
        head = f
    comp1 = Counter(d for d in L.edeg.values() if not _legal(d))
    print(f"\ncomposition preserved: {comp0 == comp1}   "
          f"total dS {total_ds:+.3f}")
    for m1, m2, _ in reversed(applied):
        undo_pair(L, m1, m2)
    comp2 = Counter(d for d in L.edeg.values() if not _legal(d))
    print(f"undo: composition restored: {comp2 == comp0}")


def mh_step(L, rng):
    """One Metropolis-Hastings worm proposal on live state L (whose driver
    routes moves through a sampler when attached).

    Proposal: anchor uniform over deg-4 edges, candidate uniform at the
    anchor.  q(x->x') = (1/n4) sum_{a in gone4} k_f(a)/n_f(a); the class
    preserves n4 so that factor cancels.  Returns a diagnostics dict."""
    d4 = L.deg4()
    e = d4[int(rng.integers(len(d4)))]
    cands = candidates(L, e)
    if not cands:
        return {"status": "no-candidates", "anchor": e}
    ci = int(rng.integers(len(cands)))
    m1, m2, ds, f, key, gone4, new4 = cands[ci]
    qf = 0.0
    for a in gone4:
        ca = cands if a == e else candidates(L, a)
        k = sum(1 for c in ca if c[4] == key)
        if ca:
            qf += k / len(ca)
    apply_pair(L, m1, m2)
    inv = frozenset((t, -c) for t, c in key)
    qr = 0.0
    for b in new4:
        cb = candidates(L, b)
        k = sum(1 for c in cb if c[4] == inv)
        if cb:
            qr += k / len(cb)
    out = {"status": "rejected", "anchor": e, "landing": f, "dS": ds,
           "qf": qf, "qr": qr}
    if qr == 0:
        undo_pair(L, m1, m2)
        out["warn"] = "no reverse candidate"      # closure says impossible
        return out
    alpha = float(np.exp(-ds)) * qr / qf
    out["alpha"] = min(1.0, alpha)
    if rng.random() < alpha:
        out["status"] = "accepted"
    else:
        undo_pair(L, m1, m2)
    return out


def mh_test(args):
    """Integration: worm proposals through a LIVE sampler -- the tracked
    objective must move by exactly dS on accepts, the cocycle must stay
    attached, and the composite kernel (worms + thermal sweeps) must leave
    no objective drift vs a from-scratch rebuild."""
    from discrete_differential_geometry import cocycle as coc
    lam, et = T.LAM, T.ESTAR
    ddg.set_random_seed(args.seed)
    m = ddg.Manifold.load(args.snap, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, hinge_degree_target=et,
        num_facets_coef=0.1, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * et / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(T.ZLEG * lam, T.CIMP * lam, tilt=[0.0] * 5)
    cocpath = args.snap.replace(".mfd", ".cocycle.npz")
    have_coc = os.path.exists(cocpath)
    if have_coc:
        e0, w0, _ = coc.load_cocycle(cocpath)
        s.enable_cocycle(np.asarray(e0), np.asarray(w0))
        s.check_cocycle()
    v = s.manifold
    L = Live(v, driver=s.do_bistellar_move)
    rng = np.random.default_rng(args.seed)
    stats = Counter()
    worst = 0.0
    print(f"state: {os.path.basename(args.snap)}  N3={v.num_facets}  "
          f"n4={len(L.deg4())}  cocycle={'on' if have_coc else 'off'}")
    for k in range(args.props):
        obj0 = s.current_objective
        r = mh_step(L, rng)
        stats[r["status"]] += 1
        if r.get("warn"):
            stats["warn"] += 1
            print(f"  prop {k:3d}: WARN {r['warn']} at {r['anchor']}")
        if r["status"] != "accepted":
            continue
        err = abs((s.current_objective - obj0) - r["dS"])
        worst = max(worst, err)
        assert err < 1e-9, (f"objective moved by {s.current_objective-obj0} "
                            f"!= dS {r['dS']}")
        if have_coc:
            s.check_cocycle()
        print(f"  prop {k:3d}: ACCEPT {r['anchor']} -> {r['landing']}  "
              f"dS={r['dS']:+.3f} alpha={r['alpha']:.3f} "
              f"(qf={r['qf']:.3f} qr={r['qr']:.3f})  obj-err {err:.1e}"
              + ("  cocycle OK" if have_coc else ""))
    print(f"\n{args.props} proposals: {dict(stats)}  "
          f"worst |dObj - dS| = {worst:.2e}")
    # composite kernel: thermal sweeps after worms, then rebuild check
    s.run(sweeps=2)
    if have_coc:
        s.check_cocycle()
    m2 = v.dup()
    s2 = ddg.ManifoldSampler(m2, p)
    s2.set_n6_potential(T.ZLEG * lam, T.CIMP * lam, tilt=[0.0] * 5)
    drift = abs(s.current_objective - s2.current_objective)
    print(f"composite kernel (worms + 2 sweeps): objective drift {drift:.2e}"
          + (", cocycle OK" if have_coc else ""))
    assert drift < 1e-6
    print("PASS")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snap", default=os.path.join(
        _ROOT, "data/mgas/lam35r_snap15000.mfd"))
    ap.add_argument("--closure-test", action="store_true")
    ap.add_argument("--mh-test", action="store_true")
    ap.add_argument("--props", type=int, default=40)
    ap.add_argument("--walk", action="store_true")
    ap.add_argument("--head", type=int, nargs=2, default=[3282, 3318])
    ap.add_argument("--steps", type=int, default=6)
    ap.add_argument("--max-ds", type=float, default=1.0)
    ap.add_argument("--n-anchors", type=int, default=10)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    print("[frame] none -- graph identity only")
    if args.closure_test:
        closure_test(args)
    elif args.mh_test:
        mh_test(args)
    elif args.walk:
        walk(args)
    else:
        main.__doc__ and print(__doc__)


if __name__ == "__main__":
    main()
