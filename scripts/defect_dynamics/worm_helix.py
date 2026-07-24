#!/usr/bin/env python3
"""Worm-along-BC-helix: find an advancing move motif, propagate it around a
torus-wrapping Boerdijk-Coxeter chain, and attempt a sector-changing closure.

Background (2026-07-24): exhaustive search proved the r-phase legal manifold
has NO nontrivial legal->legal worm cycle within 8 moves (any grammar) -- the
spin-ice pattern where contractible worms are trivial and only torus-WRAPPING
worms act (changing the web winding sector W). The deterministic BC chain
(sliding window: drop oldest vertex, exit the opposite face, adopt the apex --
after trout314/quantum-random-walks) closes into wrapping orbits in the
crystal (e.g. length 2286, winding (0,0,10)), providing the track.

Chain-aligned worms are CHAIN-INTERNAL: the face between chain tets k, k+1 is
(v_{k+1}, v_{k+2}, v_{k+3}) with apexes (v_k, v_{k+4}), so 2-3/3-2 moves on
the chain touch only chain vertices, and a worm state is coded EXACTLY by its
edge-degree overlay in chain-relative indices. Motif = path segment between
two states with identical relative code at different chain offsets.

Stages:
  1. find a wrapping orbit (pure-axis winding preferred);
  2. DFS in the chain tube (2-3 on faces of a sliding vertex window + 3-2 on
     worm deg-3 edges) from a chain-creation, detecting code repetition along
     the path -> MOTIF (relative move list, period p);
  3. re-instantiate the motif mechanically at successive offsets (D core),
     verifying the code repeats each period, around one full wrap;
  4. attempt closure to legality; report the net transformation.

Usage: worm_helix.py [REF.mfd] [--depth D] [--budget B] [--win W] [--json OUT]
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
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from cocycle_check import reference_frac_positions
from crystal_grains import REF_GLOB, best_refs

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


def find_axis_orbit(m, F, rp, period, tries=40, seed=1):
    """Random starts until an orbit with pure-axis winding is found."""
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(tries):
        t = F[rng.integers(len(F))]
        v = bc_orbit(m, [int(x) for x in rng.permutation(t)])
        w = orbit_winding(v, rp, period)
        nz = np.nonzero(w)[0]
        if len(nz) == 1:
            return v, w
        if best is None or len(v) < len(best[0]):
            best = (v, w)
    return best


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
    args = ap.parse_args()
    path = args.ref or best_refs(REF_GLOB)["r"]
    m = ddg.Manifold.load(path, 3)
    F = np.asarray(m.facets())
    rp = np.asarray(reference_frac_positions("r", args.mcell), float)
    period = float(args.mcell)

    print(f"reference: {path}  N3={len(F)}")
    verts, wind = find_axis_orbit(m, F, rp, period)
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
            m.do_bistellar_move(link, [a, b])
            ov_revert(rec)
            if found:
                return
        for tri in combinations(cand_v, 3):
            try:
                ap1, ap2 = m.face_apexes(*tri)
            except RuntimeError:
                continue
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
            m.do_bistellar_move([ap1, ap2], fl)
            ov_revert(rec)
            if found:
                return

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
    if args.json:
        with open(args.json, "w") as f:
            json.dump(dict(orbit_len=L, winding=wind.tolist(),
                           period=mo["period"], motif=mo["motif"],
                           code=[list(x) for x in mo["code"]]),
                      f, indent=1, default=str)
        print(f"wrote {os.path.abspath(args.json)}")


if __name__ == "__main__":
    main()
