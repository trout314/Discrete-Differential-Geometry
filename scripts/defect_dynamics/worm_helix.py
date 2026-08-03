#!/usr/bin/env python3
"""Worm-along-BC-helix: find an advancing move motif, propagate it around a
torus-wrapping Boerdijk-Coxeter chain, and attempt a sector-changing closure.

Background (2026-07-24, CORRECTED 2026-08-02): the search behind this was
NOT exhaustive. worm_cycles.py seeded one DFS per canon_sig bucket, and that
signature fuses orbits -- 47 buckets for 102 face orbits on R, with follow-up
move sets differing inside 28 of the 34 fused ones -- so 55 starting orbits
were never explored. Re-verified with all 102 seeds: no nontrivial
legal->legal worm cycle within SIX moves (any grammar). The eight-move bound
is NOT established; walks-5 and walks-6 have yet to be re-run.

That rigidity is what motivates this script -- the conjectured spin-ice
pattern in which contractible worms are trivial and only torus-WRAPPING worms
act (changing the web winding sector W). The conjecture is unaffected in
spirit but now rests on a six-move bound, not eight. The deterministic BC chain
(sliding window: drop oldest vertex, exit the opposite face, adopt the apex --
after trout314/quantum-random-walks) closes into wrapping orbits in the
crystal, providing the track. Those orbits are now ENUMERATED rather than
stumbled on: R has exactly 14 chain classes, of which 4 wind along a pure
axis (L=315 winding (0,0,4) and L=2286 winding (0,0,10), each as a chiral
pair). Pick one with --chain-class; it is recorded in the output.

Chain-aligned worms are CHAIN-INTERNAL: the face between chain tets k, k+1 is
(v_{k+1}, v_{k+2}, v_{k+3}) with apexes (v_k, v_{k+4}), so 2-3/3-2 moves on
the chain touch only chain vertices, and a worm state is coded EXACTLY by its
edge-degree overlay in chain-relative indices. Motif = path segment between
two states with identical relative code at different chain offsets.

Stages:
  1. select a wrapping chain CLASS (exact lookup, default the shortest
     pure-axis one; see tools/chain_select.py);
  2. DFS in the chain tube (2-3 on faces of a sliding vertex window + 3-2 on
     worm deg-3 edges) from a chain-creation, detecting code repetition along
     the path -> MOTIF (relative move list, period p);
  3. re-instantiate the motif mechanically at successive offsets (D core),
     verifying the code repeats each period, around one full wrap;
  4. attempt closure to legality; report the net transformation.

Usage: worm_helix.py [REF.mfd] [--depth D] [--budget B] [--win W] [--json OUT]

FRAME NOTE: registry coordinates are used for TOPOLOGICAL
WINDING only (integer-valued, frame-robust; CONVENTIONS.md
sec 6).
"""
import argparse
import json
import os
import sys
import time
from collections import defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts", "../../tools"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from cocycle_check import reference_frac_positions
from crystal_grains import REF_GLOB, best_refs
from chain_select import add_chain_args

ESTAR = 5.105025


# ---------------------------------------------------------------------------
# BC chain
# ---------------------------------------------------------------------------


def bc_orbit(m, window):
    """Follow the sliding-window walk until the window repeats; return the
    cyclic vertex sequence v[0..L-1] (window k = v[k..k+3], indices mod L)."""
    seen = {}
    wins = []
    w = list(window)
    while True:
        key = tuple(w)
        if key in seen:
            assert seen[key] == 0, "entered cycle off the start (impossible?)"
            break
        seen[key] = len(wins)
        wins.append(tuple(w))
        face = w[1:]
        a1, a2 = m.face_apexes(*face)
        w = face + [a2 if a1 == w[0] else a1]
    return [win[0] for win in wins]


def orbit_winding(verts, rp, period):
    """Net winding (box periods) of the cyclic vertex sequence."""
    wind = np.zeros(3)
    L = len(verts)
    for i in range(L):
        d = rp[verts[(i + 1) % L]] - rp[verts[i]]
        d -= np.round(d / period) * period
        wind += d
    return np.round(wind / period).astype(int)


def find_axis_orbit(path, selector="axis"):
    """The BC chain class winding along a single lattice axis -- LOOKED UP,
    not searched for. Returns (ChainClasses, class index, vertex sequence,
    winding of that exact chain).

    This used to take up to 40 random seed frames and return whichever
    pure-axis chain it stumbled on. Two problems, both now gone: which axis
    chain you got depended on an RNG seed and was never recorded, and if none
    of the 40 tries hit, it silently returned the SHORTEST chain it had seen
    -- not an axis chain at all -- which the caller then reported as if it
    were. The chain classes are enumerated exactly now
    (``symmetry.CrystalSymmetry``), so the axis classes are a lookup, the
    choice is recorded, and their absence raises instead of degrading.
    """
    from chain_select import ChainClasses
    cc = ChainClasses(path)
    k = cc.select(selector)
    return cc, k, cc.vertices(k), cc.representative_winding(k)


# ---------------------------------------------------------------------------
# chain-relative worm search
# ---------------------------------------------------------------------------


class ChainWorm:
    """DFS for an advancing motif in the chain tube."""

    def __init__(self, m, verts, budget=8, win=7, depth=14):
        self.m = m
        self.v = verts
        self.L = len(verts)
        self.pos_of = {}                    # vertex -> chain index (local map)
        self.budget = budget
        self.win = win
        self.depth = depth
        self.overlay = {}                   # edge(sorted verts) -> (old, new)

    def vidx(self, x):
        return self.pos_of.get(x)

    def code(self):
        """Chain-relative code of the overlay, or None if not chain-internal.
        Returns (base_index, canonical tuple)."""
        items = []
        idxs = []
        for (a, b), (d0, d1) in self.overlay.items():
            ia, ib = self.vidx(a), self.vidx(b)
            if ia is None or ib is None:
                return None
            idxs += [ia, ib]
            items.append((min(ia, ib), max(ia, ib), d0, d1))
        if not items:
            return (0, ())
        base = min(idxs)
        rel = tuple(sorted((i - base, j - base, d0, d1)
                           for i, j, d0, d1 in items))
        return (base, rel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--depth", type=int, default=14,
                    help="max moves in the motif-search DFS")
    ap.add_argument("--budget", type=int, default=8,
                    help="max concurrent illegal edges")
    ap.add_argument("--win", type=int, default=7,
                    help="face-candidate window width (chain vertices)")
    ap.add_argument("--mcell", type=int, default=3)
    ap.add_argument("--json", default=None)
    ap.add_argument("--propagate", action="store_true",
                    help="drive the worm around the full orbit and close")
    add_chain_args(ap)
    args = ap.parse_args()
    path = args.ref or best_refs(REF_GLOB)["r"]
    m = ddg.Manifold.load(path, 3)
    F = np.asarray(m.facets())
    rp = np.asarray(reference_frac_positions("r", args.mcell), float)
    period = float(args.mcell)

    print(f"reference: {path}  N3={len(F)}")
    cc, kcls, verts, wind = find_axis_orbit(path, args.chain_class)
    print(f"  {cc.summary_line(kcls)}")
    L = len(verts)
    wrap = L / max(1, int(abs(wind).max()))
    print(f"orbit: length {L}, winding {wind.tolist()} box periods "
          f"(~{wrap:.0f} tets/wrap)")

    # base edge-degree map (python tables once)
    faces0, edeg0, vedges0 = wm.build_tables(F)
    em0 = {tuple(sorted(e)): d for e, d in edeg0.items()}

    # local chain-index map around the search region (start of orbit)
    pos_of = {}
    for i in range(min(L, 4 * args.depth + 16)):
        pos_of.setdefault(verts[i], i)      # first occurrence wins

    overlay = {}

    def deg(a, b):
        k = (a, b) if a < b else (b, a)
        if k in overlay:
            return overlay[k][1]
        return em0.get(k)

    def ov_apply(deltas):
        rec = []
        for e, (old, new) in deltas.items():
            k = tuple(sorted(e))
            rec.append((k, overlay.get(k, "-")))
            base = em0.get(k)
            if new == base:
                overlay.pop(k, None)
            else:
                overlay[k] = (base, new)
        return rec

    def ov_revert(rec):
        for k, prev in reversed(rec):
            if prev == "-":
                overlay.pop(k, None)
            else:
                overlay[k] = prev

    def code():
        items = []
        idxs = []
        for (a, b), (d0, d1) in overlay.items():
            ia, ib = pos_of.get(a), pos_of.get(b)
            if ia is None or ib is None:
                return None
            idxs += [ia, ib]
            items.append((min(ia, ib), max(ia, ib), d0, d1))
        if not items:
            return (0, ())
        base = min(idxs)
        return (base, tuple(sorted((i - base, j - base, d0, d1)
                                   for i, j, d0, d1 in items)))

    # DV adapter for worm_moves deltas
    class DVo:
        def __getitem__(self, fs):
            a, b = sorted(fs)
            return deg(a, b)

    dv = DVo()

    # ---- stage 2: DFS for the motif ----
    found = []

    def n_ill():
        return sum(1 for _, (d0, d1) in overlay.items()
                   if d1 is not None and d1 not in (5, 6))

    def dfs(pathrel, pathcodes, moves_left):
        if found:
            return
        c = code()
        if c is not None:
            for pi, (pb, prel) in enumerate(pathcodes[:-1] if pathcodes else []):
                if prel == c[1] and prel != () and c[0] > pb:
                    found.append(dict(period=c[0] - pb,
                                      start_at=pi,
                                      motif=pathrel[pi:],
                                      base_from=pb,
                                      code=prel))
                    return
            pathcodes = pathcodes + [c]
        if moves_left <= 0:
            return
        # candidate faces: triples of chain vertices in a window near the worm
        base = c[0] if c and c[1] else 0
        lo = max(0, base - 1)
        cand_idx = list(range(lo, min(lo + args.win, len(verts))))
        cand_v = [verts[i] for i in cand_idx]
        # 3-2 closes on worm deg-3 edges (chain-internal)
        for (a, b), (d0, d1) in list(overlay.items()):
            if d1 != 3:
                continue
            lk = m.edge_link(a, b).tolist()
            link = sorted({x for pr in lk for x in pr})
            if len(link) != 3 or not m.has_bistellar_move([a, b], link):
                continue
            deltas = wm.delta_three_two(frozenset((a, b)), frozenset(link), dv)
            m.do_bistellar_move([a, b], link)
            rec = ov_apply(deltas)
            dfs(pathrel + [("3-2", pos_of.get(a), pos_of.get(b))],
                pathcodes, moves_left - 1)
            if found:
                return                      # keep state committed
            m.do_bistellar_move(link, [a, b])
            ov_revert(rec)
        for tri in combinations(cand_v, 3):
            try:
                ap1, ap2 = m.face_apexes(*tri)
            except RuntimeError:
                continue
            if pos_of.get(ap1) is None or pos_of.get(ap2) is None:
                continue                    # chain-internal moves only
            if not m.has_bistellar_move(list(tri), [ap1, ap2]):
                continue
            deltas = wm.delta_two_three(frozenset(tri), ap1, ap2, dv)
            bad = sum((nw is not None and nw not in (5, 6))
                      - (o is not None and o not in (5, 6))
                      for o, nw in deltas.values())
            if n_ill() + bad > args.budget:
                continue
            fl = sorted(tri)
            m.do_bistellar_move(fl, [ap1, ap2])
            rec = ov_apply(deltas)
            dfs(pathrel + [("2-3", tuple(pos_of.get(x) for x in fl),
                            (pos_of.get(ap1), pos_of.get(ap2)))],
                pathcodes, moves_left - 1)
            if found:
                return                      # keep state committed
            m.do_bistellar_move([ap1, ap2], fl)
            ov_revert(rec)

    t0 = time.time()
    dfs([], [], args.depth)
    if not found:
        print(f"\nno motif found within depth {args.depth} "
              f"({time.time() - t0:.1f}s) -- raise --depth/--win/--budget")
        return
    mo = found[0]
    print(f"\nMOTIF FOUND ({time.time() - t0:.1f}s): period {mo['period']} "
          f"chain steps, {len(mo['motif'])} moves/period")
    print(f"  worm code: {mo['code']}")
    for mv in mo["motif"]:
        print(f"    {mv}")
    if not args.propagate:
        if args.json:
            with open(args.json, "w") as f:
                json.dump(dict(orbit_len=L, winding=wind.tolist(),
                               period=mo["period"], motif=mo["motif"],
                               code=[list(x) for x in mo["code"]]),
                          f, indent=1, default=str)
            print(f"wrote {os.path.abspath(args.json)}")
        return

    # ================= stage 3-4: propagation + closure =================
    steady_rel = tuple(mo["code"])
    fdiff = {}

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

    nmv = [0]

    def apply_23(tri, a1, a2):
        fl = sorted(tri)
        deltas = wm.delta_two_three(frozenset(tri), a1, a2, dv)
        m.do_bistellar_move(fl, [a1, a2])
        ovr = ov_apply(deltas)
        fr = fdiff_apply(
            [tuple(sorted(fl + [a1])), tuple(sorted(fl + [a2]))],
            [tuple(sorted((fl[0], fl[1], a1, a2))),
             tuple(sorted((fl[1], fl[2], a1, a2))),
             tuple(sorted((fl[2], fl[0], a1, a2)))])
        nmv[0] += 1
        return ("23", fl, a1, a2, ovr, fr)

    def apply_32(a, b, link):
        lk = sorted(link)
        deltas = wm.delta_three_two(frozenset((a, b)), frozenset(lk), dv)
        m.do_bistellar_move([a, b], lk)
        ovr = ov_apply(deltas)
        fr = fdiff_apply(
            [tuple(sorted([lk[0], lk[1], a, b])),
             tuple(sorted([lk[1], lk[2], a, b])),
             tuple(sorted([lk[2], lk[0], a, b]))],
            [tuple(sorted(lk + [a])), tuple(sorted(lk + [b]))])
        nmv[0] += 1
        return ("32", a, b, lk, ovr, fr)

    def undo_move(r):
        if r[0] == "23":
            m.do_bistellar_move([r[2], r[3]], r[1])
        else:
            m.do_bistellar_move(r[3], [r[1], r[2]])
        ov_revert(r[4])
        fdiff_revert(r[5])
        nmv[0] -= 1

    def remap(base):
        pos_of.clear()
        for i in range(max(0, base - 3), base + args.win + 10):
            pos_of[verts[i % L]] = i

    def n_ill():
        return sum(1 for _, (d0, d1) in overlay.items()
                   if d1 is not None and d1 not in (5, 6))

    nodes = [0]

    def seq_dfs(depth, accept, base_lo):
        if accept():
            return []
        if depth <= 0:
            return None
        nodes[0] += 1
        if nodes[0] > 100000:
            return None
        # 3-2 on worm deg-3 edges
        for (a, b), (d0, d1) in list(overlay.items()):
            if d1 != 3:
                continue
            lk = m.edge_link(a, b).tolist()
            link = sorted({x for pr in lk for x in pr})
            if len(link) != 3 or not m.has_bistellar_move([a, b], link):
                continue
            r = apply_32(a, b, link)
            sub = seq_dfs(depth - 1, accept, base_lo)
            if sub is not None:
                return [r] + sub
            undo_move(r)
        # 2-3 on chain-window triples
        cand_v = [verts[i % L] for i in range(base_lo, base_lo + args.win + 2)]
        for tri in combinations(cand_v, 3):
            try:
                ap1, ap2 = m.face_apexes(*tri)
            except RuntimeError:
                continue
            if not m.has_bistellar_move(list(tri), [ap1, ap2]):
                continue
            deltas = wm.delta_two_three(frozenset(tri), ap1, ap2, dv)
            grow = sum((nw is not None and nw not in (5, 6))
                       - (o is not None and o not in (5, 6))
                       for o, nw in deltas.values())
            if n_ill() + grow > args.budget:
                continue
            r = apply_23(tri, ap1, ap2)
            sub = seq_dfs(depth - 1, accept, base_lo)
            if sub is not None:
                return [r] + sub
            undo_move(r)
        return None

    # ---- start from the motif search's committed state ----
    def at_steady_past(bmin):
        def acc():
            c = code()
            return c is not None and c[1] == steady_rel and c[0] >= bmin
        return acc

    c0 = code()
    assert c0 is not None and c0[1] == steady_rel, "state not at steady code"
    base = c0[0]
    # reconstruct the facet diff once, exactly
    fs0 = set(map(tuple, np.sort(F, axis=1).tolist()))
    fs1 = set(map(tuple, np.sort(np.asarray(m.facets()), axis=1).tolist()))
    for t in fs1 - fs0:
        fdiff[t] = 1
    for t in fs0 - fs1:
        fdiff[t] = -1
    print(f"\n[propagate] starting from committed state: base {base}, "
          f"facet diff {len(fdiff)} tets")

    # ---- drive around the full orbit; fast-path = replay last period ----
    import gc as _gc
    from discrete_differential_geometry import _dlang as _dl
    bf = mo["base_from"]
    last_rel = []
    for mv in mo["motif"]:
        if mv[0] == "2-3":
            if None in mv[1] or mv[2][0] is None or mv[2][1] is None:
                last_rel = None
                break
            last_rel.append(("23", tuple(i - bf for i in mv[1]),
                             mv[2][0] - bf, mv[2][1] - bf))
        else:
            if mv[1] is None or mv[2] is None:
                last_rel = None
                break
            last_rel.append(("32", mv[1] - bf, mv[2] - bf))
    if last_rel is not None:
        print(f"[propagate] fast-path seeded from motif "
              f"({len(last_rel)} moves, rebased -{bf})")
    period_hist = []
    target_total = L          # one full orbit returns worm to launch site
    t2 = time.time()
    fails = 0
    while base < target_total:
        prev_base = base
        remap(base)
        advanced = None
        if last_rel is not None:
            recs = []
            ok = True
            for mv in last_rel:
                if mv[0] == "23":
                    tri = [verts[(base + i) % L] for i in mv[1]]
                    a1 = verts[(base + mv[2]) % L]
                    a2 = verts[(base + mv[3]) % L]
                    if not m.has_bistellar_move(sorted(tri), [a1, a2]):
                        ok = False
                        break
                    recs.append(apply_23(tuple(tri), a1, a2))
                else:
                    a = verts[(base + mv[1]) % L]
                    b = verts[(base + mv[2]) % L]
                    lk = m.edge_link(a, b).tolist() if True else None
                    link = sorted({x for pr in lk for x in pr})
                    if len(link) != 3 or \
                            not m.has_bistellar_move([a, b], link):
                        ok = False
                        break
                    recs.append(apply_32(a, b, link))
            if ok:
                c = code()
                if c is not None and c[1] == steady_rel and c[0] > base:
                    advanced = c[0]
            if advanced is None:
                for r in reversed(recs):
                    undo_move(r)
        if advanced is None:
            nodes[0] = 0
            seq = seq_dfs(len(mo["motif"]) + 2, at_steady_past(base + 1),
                          base)
            if seq is None:
                fails += 1
                print(f"[propagate] STUCK at base {base} "
                      f"({nmv[0]} moves so far)")
                break
            advanced = code()[0]
            # record the successful sequence relative to prev_base
            last_rel = []
            for r in seq:
                if r[0] == "23":
                    idx = [pos_of.get(x) for x in r[1]]
                    i1, i2 = pos_of.get(r[2]), pos_of.get(r[3])
                    if None in idx or i1 is None or i2 is None:
                        last_rel = None
                        break
                    last_rel.append(("23",
                                     tuple(i - prev_base for i in idx),
                                     i1 - prev_base, i2 - prev_base))
                else:
                    i1, i2 = pos_of.get(r[1]), pos_of.get(r[2])
                    if i1 is None or i2 is None:
                        last_rel = None
                        break
                    last_rel.append(("32", i1 - prev_base, i2 - prev_base))
        period_hist.append(advanced - prev_base)
        base = advanced
        if len(period_hist) % 50 == 0:
            _gc.collect()
            _dl._lib.ddg_gc_collect()
        if len(period_hist) % 100 == 0:
            print(f"[propagate] base {base}/{target_total} "
                  f"({nmv[0]} moves, {time.time() - t2:.0f}s)")

    if base >= target_total:
        print(f"[propagate] WRAPPED the orbit: base {base} >= {L}, "
              f"{nmv[0]} moves, {time.time() - t2:.0f}s; closing...")
        remap(base)

        def legal_now():
            return n_ill() == 0

        nodes[0] = 0
        seq = seq_dfs(args.depth, lambda: legal_now(), base)
        if seq is None:
            print("[close] could not reach legality within depth "
                  f"{args.depth}")
        else:
            pairs, _ = m.illegal_edges()
            assert len(pairs) == 0, "D core disagrees on legality"
            deg_net = {k: v for k, v in overlay.items()}
            print(f"\n════ CLOSED after {nmv[0]} moves ════")
            print(f"  degree net (overlay): {len(deg_net)} edges changed")
            for k, (d0, d1) in sorted(deg_net.items())[:12]:
                print(f"    {k}: {d0} -> {d1}")
            print(f"  FACET TRAIL: {len(fdiff)} tets differ from crystal")
            if not deg_net and not fdiff:
                print("  => IDENTITY (wrap acted trivially)")
            else:
                print("  => NONTRIVIAL LEGAL CYCLE"
                      + (" (dS = 0: pure degree-preserving "
                         "retriangulation)" if not deg_net else ""))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(dict(orbit_len=L, winding=wind.tolist(),
                           period=mo["period"], motif_moves=len(mo["motif"]),
                           total_moves=nmv[0], final_base=base,
                           n_deg_net=len(overlay), n_facet_trail=len(fdiff),
                           period_hist=period_hist[:50]),
                      f, indent=1, default=str)
        print(f"wrote {os.path.abspath(args.json)}")


if __name__ == "__main__":
    main()
