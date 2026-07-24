#!/usr/bin/env python3
"""Worm program stage 2: minimal worm cycles and path-independence, as an
APPLY/UNDO depth-first search on a single D-core manifold.

v2 (2026-07-24): no copies, no python recounts. Per explored state:
  * illegal edges from the D core's incremental degree map (~1.6 ms);
  * link cycles via the O(degree) ridge walk (~11 us; hint tets maintained
    from each move's known support, slow-path fallback if stale);
  * moves applied AND undone through the validated targeted-move API
    (undo of a bistellar move = center/cocenter swapped; undo of a 4-4 =
    the 4-4 on the added diagonal with the interleaved link cycle);
  * a python Overlay tracks the net edge-degree diff vs the crystal
    incrementally (worm_moves exact deltas), so endpoint classification
    is O(changed edges).

Grammar per branch: [2-3 create] -> up to --walks 4-4 steps -> [3-2 close].
Legal endpoints (D-verified empty illegal set) are classified by their net
transformation (degree changes, junction transitions, dS_shape); every
state is content-hashed to find distinct sequences reaching identical
states (the groupoid relations).

Usage: worm_cycles.py [REF.mfd] [--walks N] [--json OUT] [--profile]
"""
import argparse
import json
import os
import sys
import time
from collections import defaultdict
from contextlib import contextmanager

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs
from worm_catalog import canon_sig

ESTAR = 5.105025

TIMES = defaultdict(float)
COUNTS = defaultdict(int)


@contextmanager
def tick(name):
    t0 = time.perf_counter()
    try:
        yield
    finally:
        TIMES[name] += time.perf_counter() - t0
        COUNTS[name] += 1


def profile_report():
    tot = sum(TIMES.values())
    print(f"\n── profile ({tot:.1f}s in instrumented phases) ──")
    print(f"{'phase':12s} {'total_s':>8s} {'calls':>8s} {'ms/call':>8s} {'share':>6s}")
    for k in sorted(TIMES, key=lambda k: -TIMES[k]):
        print(f"{k:12s} {TIMES[k]:8.2f} {COUNTS[k]:8d} "
              f"{1000 * TIMES[k] / COUNTS[k]:8.3f} {100 * TIMES[k] / tot:5.1f}%")


# ---------------------------------------------------------------------------
# incremental net-diff overlay vs the crystal
# ---------------------------------------------------------------------------


class Overlay:
    def __init__(self, base):
        self.base = base                    # edge(sorted tuple) -> crystal deg
        self.cur = {}                       # only where != base (None=absent)

    def deg(self, key):
        if key in self.cur:
            return self.cur[key]
        return self.base.get(key)

    def apply(self, deltas):
        rec = []
        for e, (old, new) in deltas.items():
            k = tuple(sorted(e))
            rec.append((k, self.cur.get(k, "-")))
            if new == self.base.get(k):
                self.cur.pop(k, None)
            else:
                self.cur[k] = new
        return rec

    def revert(self, rec):
        for k, prev in reversed(rec):
            if prev == "-":
                self.cur.pop(k, None)
            else:
                self.cur[k] = prev

    def net(self):
        return {k: (self.base.get(k), v) for k, v in self.cur.items()}


class DV:
    """frozenset-keyed degree view over an Overlay (for worm_moves deltas)."""

    def __init__(self, ov):
        self.ov = ov

    def __getitem__(self, fs):
        return self.ov.deg(tuple(sorted(fs)))


# ---------------------------------------------------------------------------


def link_cycle(m, a, b, hints):
    key = (a, b) if a < b else (b, a)
    hint = hints.get(key)
    if hint is not None:
        try:
            with tick("d_link"):
                return m.edge_link_cycle(a, b, hint).tolist()
        except RuntimeError:
            pass
    with tick("d_link_slow"):
        pairs = m.edge_link(a, b).tolist()
    adj = defaultdict(list)
    for x, y in pairs:
        adj[x].append(y)
        adj[y].append(x)
    verts = sorted(adj)
    v0 = verts[0]
    cyc = [v0, min(adj[v0])]
    while len(cyc) < len(verts):
        nxt = [x for x in adj[cyc[-1]] if x != cyc[-2]]
        if not nxt:
            return None
        cyc.append(nxt[0])
    return cyc


def set_hints(hints, tets, touched=None):
    # cover EVERY edge of every new tet (a move can invalidate the hint of an
    # edge whose degree it does not change -- the removed tets still carried it)
    rec = []
    seen = set()
    for t in tets:
        for i in range(4):
            for j in range(i + 1, 4):
                k = (t[i], t[j]) if t[i] < t[j] else (t[j], t[i])
                if k in seen:
                    continue
                seen.add(k)
                rec.append((k, hints.get(k, "-")))
                hints[k] = tuple(t)
    return rec


def revert_hints(hints, rec):
    for k, prev in reversed(rec):
        if prev == "-":
            hints.pop(k, None)
        else:
            hints[k] = prev


# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--walks", type=int, default=2,
                    help="max 4-4 walk steps per branch (default 2)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--profile", action="store_true")
    args = ap.parse_args()
    path = args.ref or best_refs(REF_GLOB)["r"]
    m = ddg.Manifold.load(path, 3)
    F0 = np.asarray(m.facets())
    faces0, edeg0, vedges0 = wm.build_tables(F0)
    em0 = {tuple(sorted(e)): d for e, d in edeg0.items()}
    n6m0 = {v: wm.vertex_counters(v, edeg0, vedges0) for v in vedges0}
    print(f"reference: {path}  N3={len(F0)}  walk budget {args.walks}")

    buckets = defaultdict(list)
    for face, d, e, valid in wm.two_three_sites(F0, faces0, edeg0):
        if valid:
            buckets[canon_sig(face, d, e, edeg0, vedges0)].append((face, d, e))
    reps = [s[0] for _, s in sorted(buckets.items())]
    print(f"{len(reps)} creation-class representatives\n")

    endpoints = defaultdict(list)
    legal_cycles = []
    ov = Overlay(em0)
    hints = {}
    fdiff = {}                    # tet -> +1 (added) / -1 (removed) vs crystal
    stats = dict(n_seq=0, n_legal=0)

    def fdiff_apply(remove, add):
        rec = []
        for t, s in [(t, -1) for t in remove] + [(t, +1) for t in add]:
            rec.append((t, fdiff.get(t, 0)))
            c = fdiff.get(t, 0) + s
            if c == 0:
                fdiff.pop(t, None)
            else:
                fdiff[t] = c
        return rec

    def fdiff_revert(rec):
        for t, prev in reversed(rec):
            if prev == 0:
                fdiff.pop(t, None)
            else:
                fdiff[t] = prev

    def state_hash():
        with tick("hash"):
            return hash(frozenset(fdiff.items()))

    def vertex_trans_from_net(net):
        dv = defaultdict(lambda: [0, 0])
        for (a, b), (d0, d1) in net.items():
            for v in (a, b):
                if d0 is not None:
                    dv[v][0] -= d0 >= 6
                    dv[v][1] -= d0 not in (5, 6)
                if d1 is not None:
                    dv[v][0] += d1 >= 6
                    dv[v][1] += d1 not in (5, 6)
        out = {}
        for v, (dn6, dm) in dv.items():
            b = n6m0.get(v, (0, 0))
            a_ = (b[0] + dn6, b[1] + dm)
            if b != a_:
                out[v] = (b, a_)
        return out

    def record_legal(pathdesc):
        with tick("classify"):
            net = ov.net()
            if not net:
                legal_cycles.append(dict(path=list(pathdesc), n_changed=0,
                                         dS=0.0, net6=0, sig=("identity",)))
                return
            x = ESTAR - int(ESTAR)
            dS = 0.0
            for k, (d0, d1) in net.items():
                if d0 is not None:
                    dS -= (ESTAR / 6.0) * ((d0 - ESTAR) ** 2 - x * (1 - x))
                if d1 is not None:
                    dS += (ESTAR / 6.0) * ((d1 - ESTAR) ** 2 - x * (1 - x))
            vtr = vertex_trans_from_net(net)
            for v, (b, a_) in vtr.items():
                dS += 0.6 * (wm.U_zleg(a_[0]) - wm.U_zleg(b[0]))
            net6 = sum(1 for d0, d1 in net.values() if d1 == 6) \
                - sum(1 for d0, d1 in net.values() if d0 == 6)
            sig = (tuple(sorted((d0 or 0, d1 or 0)
                                for d0, d1 in net.values())),
                   tuple(sorted(f"{wm.class_name(*b)}->{wm.class_name(*a_)}"
                                for b, a_ in vtr.values())))
            legal_cycles.append(dict(path=list(pathdesc), n_changed=len(net),
                                     dS=round(dS, 4), net6=net6, sig=sig))

    def explore(pathdesc, budget):
        endpoints[state_hash()].append(tuple(pathdesc))
        with tick("ill_overlay"):
            ill = [(k, v) for k, v in ov.cur.items()
                   if v is not None and v not in (5, 6)]
        if not ill:
            with tick("d_illegal"):        # authoritative check, rare
                pairs, _ = m.illegal_edges()
            assert len(pairs) == 0, "overlay disagrees with D core"
            stats["n_legal"] += 1
            record_legal(pathdesc)
            return
        e3 = [k for k, v in ill if v == 3]
        e4 = [k for k, v in ill if v == 4]
        for ea, eb in e3:
            cyc = link_cycle(m, ea, eb, hints)
            if cyc is None or len(cyc) != 3:
                continue
            link = sorted(cyc)
            with tick("d_hasmove"):
                ok = m.has_bistellar_move([ea, eb], link)
            if not ok:
                continue
            deltas = wm.delta_three_two(frozenset((ea, eb)), frozenset(link),
                                        DV(ov))
            with tick("d_move"):
                m.do_bistellar_move([ea, eb], link)
            rec = ov.apply(deltas)
            newt = [tuple(sorted(link + [ea])), tuple(sorted(link + [eb]))]
            oldt = [tuple(sorted([link[0], link[1], ea, eb])),
                    tuple(sorted([link[1], link[2], ea, eb])),
                    tuple(sorted([link[2], link[0], ea, eb]))]
            frec = fdiff_apply(oldt, newt)
            hrec = set_hints(hints, newt)
            stats["n_seq"] += 1
            pathdesc.append(("3-2", (ea, eb)))
            explore(pathdesc, budget)
            pathdesc.pop()
            with tick("d_move"):
                m.do_bistellar_move(link, [ea, eb])
            revert_hints(hints, hrec)
            fdiff_revert(frec)
            ov.revert(rec)
        if budget <= 0:
            return
        for ea, eb in e4:
            cyc = link_cycle(m, ea, eb, hints)
            if cyc is None or len(cyc) != 4:
                continue
            for dg in (0, 1):
                diag = (cyc[0], cyc[2]) if dg == 0 else (cyc[1], cyc[3])
                with tick("d_hasmove"):
                    ok = m.has_hinge_move([ea, eb], cyc, dg)
                if not ok:
                    continue
                deltas = wm.delta_four_four(frozenset((ea, eb)), list(cyc),
                                            frozenset(diag), DV(ov))
                with tick("d_move"):
                    m.do_hinge_move([ea, eb], cyc, dg)
                rec = ov.apply(deltas)
                b_, d_ = ((cyc[1], cyc[3]) if dg == 0 else (cyc[0], cyc[2]))
                tets = [tuple(sorted((diag[0], diag[1], ea, b_))),
                        tuple(sorted((diag[0], diag[1], b_, eb))),
                        tuple(sorted((diag[0], diag[1], eb, d_))),
                        tuple(sorted((diag[0], diag[1], d_, ea)))]
                oldt = [tuple(sorted((ea, eb, cyc[i], cyc[(i + 1) % 4])))
                        for i in range(4)]
                frec = fdiff_apply(oldt, tets)
                hrec = set_hints(hints, tets)
                stats["n_seq"] += 1
                pathdesc.append(("4-4", (ea, eb), tuple(sorted(diag))))
                explore(pathdesc, budget - 1)
                pathdesc.pop()
                inv_cyc = [ea, b_, eb, d_]
                with tick("d_move"):
                    m.do_hinge_move(sorted(diag), inv_cyc, 0)
                revert_hints(hints, hrec)
                fdiff_revert(frec)
                ov.revert(rec)

    t0 = time.time()
    for face, d, e in reps:
        fl = sorted(face)
        deltas = wm.delta_two_three(frozenset(face), d, e, DV(ov))
        with tick("d_move"):
            m.do_bistellar_move(fl, [d, e])
        rec = ov.apply(deltas)
        tets = [tuple(sorted((fl[0], fl[1], d, e))),
                tuple(sorted((fl[1], fl[2], d, e))),
                tuple(sorted((fl[2], fl[0], d, e)))]
        oldt = [tuple(sorted(fl + [d])), tuple(sorted(fl + [e]))]
        frec = fdiff_apply(oldt, tets)
        hrec = set_hints(hints, tets)
        stats["n_seq"] += 1
        pd = [("2-3", tuple(fl), (d, e))]
        explore(pd, args.walks)
        with tick("d_move"):
            m.do_bistellar_move([d, e], fl)
        revert_hints(hints, hrec)
        fdiff_revert(frec)
        ov.revert(rec)
        assert not ov.cur and not fdiff, "state not clean after undo"

    print(f"sequences applied (D core, apply/undo): {stats['n_seq']} "
          f"in {time.time() - t0:.1f}s; legal endpoints: {stats['n_legal']}")
    ident = sum(1 for c in legal_cycles if c["n_changed"] == 0)
    nontriv = [c for c in legal_cycles if c["n_changed"] > 0]
    print(f"  identity closures: {ident}")
    print(f"  NONTRIVIAL legal->legal cycles: {len(nontriv)}")
    bysig = defaultdict(list)
    for c in nontriv:
        bysig[c["sig"]].append(c)
    print(f"  distinct net-transformation classes: {len(bysig)}\n")
    for sig, cs in sorted(bysig.items(), key=lambda kv: -len(kv[1]))[:12]:
        c = cs[0]
        print(f"x{len(cs):4d}  edges changed {c['n_changed']:2d}  "
              f"net6 {c['net6']:+d}  dS_shape {c['dS']:+.4f}")
        print(f"      degree changes {sig[0]}")
        print(f"      junctions: {', '.join(sig[1]) or '(none)'}")
        print(f"      example: {' ; '.join(str(s) for s in c['path'])}")

    dup = {h: ps for h, ps in endpoints.items() if len(set(ps)) > 1}
    print(f"\npath-independence: {len(dup)} states reached by >1 distinct "
          f"sequence")
    for h, ps in sorted(dup.items(), key=lambda kv: -len(set(kv[1])))[:4]:
        u = sorted(set(ps))
        print(f"  state reached {len(u)} ways, e.g.:")
        for p in u[:3]:
            print(f"    {' ; '.join(str(s) for s in p)}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(dict(
                n_seq=stats["n_seq"], n_legal=stats["n_legal"], identity=ident,
                classes=[dict(count=len(cs),
                              **{k: v for k, v in cs[0].items() if k != "sig"},
                              degsig=[list(x) for x in cs[0]["sig"][0]],
                              juncsig=list(cs[0]["sig"][1]))
                         for cs in bysig.values()],
                n_multi_path_states=len(dup)), f, indent=1, default=str)
        print(f"\nwrote {os.path.abspath(args.json)}")
    if args.profile:
        profile_report()


if __name__ == "__main__":
    main()
