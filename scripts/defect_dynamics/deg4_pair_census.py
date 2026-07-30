#!/usr/bin/env python3
"""Pair census: the COMPLETE 2-move reaction vocabulary at deg-4 fragment
tips.

For every deg-4 edge with a valence-1 endpoint in the illegal-edge graph
(fragment tips, including both ends of isolated edges), enumerate ALL valid
(2->3 | 3->2) x (2->3 | 3->2) composites -- first move touching the tip
edge's vertices, second touching the disturbed region -- with NO goal
filter.  Every composite that changes the illegal-edge content is recorded:

  * net reaction signature: the tip edge's transition tagged separately,
    plus the sorted (old,new) degree transitions of every other edge whose
    illegal status changed;
  * exact full-action dS (EDQ + zleg*U(n6) + cimp*m^2, x lam);
  * net tet count (f-vector neutrality);
  * for each NEW deg-4 edge: harmonic-chart displacement of its midpoint
    from the tip midpoint (cells), and whether it lies inside the tip's
    hinge-flip octahedron (landing there is flip-reachable, not transport).

This is the design table for the directional worm: the reactions that
exist, ranked by cost, with their landing geometry.
Writes <out>.jsonl (one record per composite) + a summary to stdout.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
import tip_retract_search as T
import defect_viewer as dv

CELL = 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snap", default=os.path.join(
        _ROOT, "data/mgas/lam35r_snap15000.mfd"))
    ap.add_argument("--out", default=os.path.join(
        _ROOT, "data/rxn_lam035_m4/pair_census"))
    ap.add_argument("--max-tips", type=int, default=0,
                    help="0 = all tips")
    args = ap.parse_args()

    m = ddg.Manifold.load(args.snap, 3)
    fac = np.asarray(m.facets())
    F = [tuple(sorted(int(x) for x in t)) for t in fac]
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    d4 = [e for e, d in ill.items() if d == 4]
    adjv = defaultdict(list)
    for a, b in d4:
        adjv[a].append(b)
        adjv[b].append(a)
    val = {v: len(nb) for v, nb in adjv.items()}
    tips = [e for e in d4 if val[e[0]] == 1 or val[e[1]] == 1]
    if args.max_tips:
        tips = tips[:args.max_tips]
    X, P = dv.positions(args.snap, fac)
    print(f"{os.path.basename(args.snap)}: {len(d4)} deg-4, "
          f"{len(tips)} fragment tips to census", flush=True)

    def legal(d):
        return d in (5, 6)

    fh = open(args.out + ".jsonl", "w")
    nrec = 0
    sig_stats = defaultdict(lambda: dict(n=0, ds=[], esc=0, fneu=0))

    for ti, e in enumerate(tips):
        Vball = T.ball(F, e, T.R_BALL)
        Vcore = frozenset(T.ball(F, e, T.R_CORE))
        region = [t for t in F if set(t) <= Vball]
        Pch = T.Patch([frozenset(t) for t in region])
        init_edeg = dict(Pch.edeg)
        init_tets = len(Pch.tets)
        octa = set(e)
        for t in Pch.tets:
            if set(e) <= t:
                octa |= t
        # tip midpoint in the chart (min-image local -- fine at this scale)
        exy = (X[e[0]] + np.round((X[e[1]] - X[e[0]]) / P) * -P + X[e[1]]) / 2
        emid = X[e[0]] + 0.5 * ((X[e[1]] - X[e[0]])
                                - np.round((X[e[1]] - X[e[0]]) / P) * P)

        def disturbed():
            f = set(e)
            for ed in Pch.touched:
                if Pch.edeg.get(ed, 0) != init_edeg.get(ed, 0):
                    f |= set(ed)
            return f

        def apply_mv(mv):
            if mv[0] == "23":
                _, tri, x, y = mv
                Pch.apply_23(tri, x, y)
            else:
                _, (u, w), a, b, c = mv
                Pch.apply_32(u, w, a, b, c)

        def undo_mv(mv):
            if mv[0] == "23":
                _, tri, x, y = mv
                Pch.undo_23(tri, x, y)
            else:
                _, (u, w), a, b, c = mv
                Pch.undo_32(u, w, a, b, c)

        m1s = Pch.moves(Vcore, set(e))
        for m1 in m1s:
            apply_mv(m1)
            for m2 in Pch.moves(Vcore, disturbed()):
                apply_mv(m2)
                # net transitions over touched edges
                trans = []
                for ed in Pch.touched:
                    a, b = init_edeg.get(ed, 0), Pch.edeg.get(ed, 0)
                    if a != b and (not legal(a) or not legal(b)
                                   or a == 0 or b == 0):
                        # keep only transitions with illegal/creation content
                        if not (legal(a) and legal(b)):
                            trans.append((ed, a, b))
                # does the composite change illegal content at all?
                ill_ch = [(ed, a, b) for ed, a, b in trans
                          if (a and not legal(a)) != (b and not legal(b))
                          or (a and not legal(a) and b and not legal(b))]
                if ill_ch:
                    tip_tr = next(((a, b) for ed, a, b in trans if ed == e),
                                  None)
                    rest = sorted((a, b) for ed, a, b in trans
                                  if ed != e and not
                                  (legal(a) and legal(b)) and a != b
                                  and ((a and not legal(a)) or
                                       (b and not legal(b))))
                    ds = T.LAM * T.dS_full(init_edeg, Pch.edeg,
                                           frozenset(Pch.touched))
                    dtets = len(Pch.tets) - init_tets
                    new4 = []
                    for ed, a, b in trans:
                        if b == 4 and a != 4:
                            d0 = X[ed[0]] - emid
                            d0 -= np.round(d0 / P) * P
                            d1 = X[ed[1]] - emid
                            d1 -= np.round(d1 / P) * P
                            mid = (d0 + d1) / 2 / CELL
                            new4.append(dict(
                                e=list(ed),
                                in_octa=set(ed) <= octa,
                                disp=[round(float(x), 4) for x in mid],
                                r=round(float(np.linalg.norm(mid)), 4)))
                    sig = (str(tip_tr), tuple(f"{a}>{b}" for a, b in rest))
                    esc = any(not n["in_octa"] for n in new4)
                    rec = dict(tip=list(e), m1=repr(m1), m2=repr(m2),
                               tip_tr=tip_tr, rest=[f"{a}>{b}"
                                                    for a, b in rest],
                               dS=round(float(ds), 4), dtets=dtets,
                               new4=new4)
                    fh.write(json.dumps(rec) + "\n")
                    nrec += 1
                    s = sig_stats[sig]
                    s["n"] += 1
                    s["ds"].append(ds)
                    s["esc"] += esc
                    s["fneu"] += (dtets == 0)
                undo_mv(m2)
            undo_mv(m1)
        if (ti + 1) % 10 == 0:
            print(f"  tip {ti + 1}/{len(tips)}: {nrec} reactions so far",
                  flush=True)
    fh.close()
    print(f"\nwrote {args.out}.jsonl  ({nrec} reaction records)", flush=True)

    # ---- summary ------------------------------------------------------
    print(f"\n=== reaction vocabulary (top 30 signatures) ===")
    print(f"{'n':>7s} {'tip transition':>15s} {'other illegal transitions':32s} "
          f"{'dS med':>8s} {'dS min':>8s} {'%fneu':>6s} {'%escape':>8s}")
    rows = sorted(sig_stats.items(), key=lambda kv: -kv[1]["n"])
    for (tip_tr, rest), s in rows[:30]:
        d = np.array(s["ds"])
        print(f"{s['n']:7d} {tip_tr:>15s} {' '.join(rest) or '-':32s} "
              f"{np.median(d):8.3f} {d.min():8.3f} "
              f"{100 * s['fneu'] / s['n']:6.0f} {100 * s['esc'] / s['n']:8.0f}")
    # highlight: cheapest reactions with an escaped new deg-4
    print(f"\n=== cheapest signatures with an ESCAPED new deg-4 (transport "
          f"candidates) ===")
    cand = [(k, s) for k, s in sig_stats.items() if s["esc"] > 0]
    cand.sort(key=lambda kv: np.median(kv[1]["ds"]))
    for (tip_tr, rest), s in cand[:15]:
        d = np.array(s["ds"])
        print(f"{s['n']:7d} {tip_tr:>15s} {' '.join(rest) or '-':32s} "
              f"med {np.median(d):7.3f}  min {d.min():7.3f}  "
              f"esc {100 * s['esc'] / s['n']:3.0f}%")


if __name__ == "__main__":
    main()
