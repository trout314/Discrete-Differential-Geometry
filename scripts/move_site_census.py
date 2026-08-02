#!/usr/bin/env python3
"""Classify the 2->3 move sites of a triangulation into local-isomorphism
classes, and report each class's curvature ledger and a representative site.

Every face of a 3-manifold triangulation whose two apexes are not already
joined is a legal 2->3; in a pristine TCP crystal that is EVERY face. The
resulting defect's boundary distance map is universal (see
``ball_boundary.py``) -- it depends only on the move's combinatorics -- so
what distinguishes one site from another is entirely (a) the curvature
LEDGER, i.e. which of the ten affected edges start at which degree, and
(b) the exterior embedding. This script enumerates both.

The ledger is what carries the physics. Independently of the site:

  * the new edge always has degree 3, and the three face edges each drop
    one degree, so the CORE gains exactly +2*pi of deficit;
  * the six spokes each gain one degree, so the SHELL loses exactly
    6*theta;

but *which* edges those are varies. A face edge at degree 5 goes to 4 (an
illegal, positive-disclination edge); one at degree 6 goes to 5 and stays
legal -- so the illegal signature is (3,4,4,4) or (3,4,4). A spoke at
degree 5 goes to 6, creating a NEW geodesic-permeable line (only degree
>= 6 edges can carry PL minimizers); a spoke already at 6 goes to 7,
widening an existing line's shadow wedge from 63.17 to 133.70 degrees.

Classification is by Weisfeiler-Leman refinement of the edge-degree
labelled 1-skeleton, canonicalised over the 12 symmetries of the abstract
configuration (S3 on the face x S2 on the apexes). WL colours are
automorphism invariants, so each class is a UNION of orbits and the count
reported is a rigorous LOWER bound on the number of inequivalent sites --
two sites in different classes are certainly inequivalent, two in the
same class are not certainly equivalent.

Usage:
  python scripts/move_site_census.py data/tcp_reference/T3_R_m3_N24462.mfd
  python scripts/move_site_census.py <mfd> --full --json sites.json
"""
import os, sys, json, argparse, itertools, collections

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "python"))

from discrete_differential_geometry import Manifold

#: the 12 symmetries of the abstract 2->3: S3 on (A,B,C) x S2 on (P,Q)
SYMS = [(fp, ap) for fp in itertools.permutations(range(3))
        for ap in ((0, 1), (1, 0))]


def load(path):
    m = Manifold.load(path, 3)
    facets = [tuple(sorted(map(int, f))) for f in m.facets()]
    deg = {tuple(sorted(map(int, e))): m.degree(e) for e in m.simplices(1)}
    nbr = collections.defaultdict(list)
    for (u, v) in deg:
        nbr[u].append(v); nbr[v].append(u)
    f2a = {}
    for t in facets:
        for i in range(4):
            f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
    return m, facets, deg, nbr, f2a


def sites(deg, f2a):
    return [(f, tuple(sorted(a))) for f, a in f2a.items()
            if len(a) == 2 and tuple(sorted(a)) not in deg]


def wl(deg, nbr, max_rounds=20):
    """WL refinement on the edge-degree-labelled 1-skeleton; returns
    (colour map, per-round class counts) stopping at the fixed point."""
    col = {v: len(ns) for v, ns in nbr.items()}
    hist = [len(set(col.values()))]
    for _ in range(max_rounds):
        new = {v: (col[v], tuple(sorted((col[u], deg[(min(v, u), max(v, u))])
                                        for u in nbr[v])))
               for v in col}
        pal = {c: i for i, c in enumerate(sorted(set(new.values())))}
        col = {v: pal[new[v]] for v in new}
        hist.append(len(pal))
        if hist[-1] == hist[-2]:
            break
    return col, hist


def invariant(face, apex, col, deg):
    """Canonical key: lexicographic min over the 12 relabellings."""
    e = lambda x, y: deg[(min(x, y), max(x, y))]
    best = None
    for fp, ap in SYMS:
        a, b, c = face[fp[0]], face[fp[1]], face[fp[2]]
        p, q = apex[ap[0]], apex[ap[1]]
        k = (col[a], col[b], col[c], col[p], col[q],
             e(a, b), e(a, c), e(b, c),
             e(a, p), e(b, p), e(c, p), e(a, q), e(b, q), e(c, q))
        if best is None or k < best:
            best = k
    return best


def ledger(face, apex, deg):
    """(sorted face-edge degrees, sorted spoke degrees) before the move."""
    e = lambda x, y: deg[(min(x, y), max(x, y))]
    return (tuple(sorted(e(*p) for p in itertools.combinations(face, 2))),
            tuple(sorted(e(x, y) for x in face for y in apex)))


def describe(fe, sp):
    return {"face_edges": list(fe), "spokes": list(sp),
            "illegal_after": [3] + [4] * sum(1 for d in fe if d == 5),
            "new_deg6": sum(1 for d in sp if d == 5),
            "new_deg7": sum(1 for d in sp if d == 6)}


def main(path, full, jsonout):
    m, facets, deg, nbr, f2a = load(path)
    S = sites(deg, f2a)
    dh = dict(sorted(collections.Counter(deg.values()).items()))
    Zh = dict(sorted(collections.Counter(len(nbr[v]) for v in nbr).items()))
    print(f"{os.path.basename(path)}: {len(nbr)} vertices, {len(facets)} tets, "
          f"{len(deg)} edges {dh}")
    print(f"  Z-census {Zh}")
    print(f"  faces {len(f2a)}, legal 2->3 sites {len(S)} "
          f"({'every face' if len(S) == len(f2a) else f'{100*len(S)/len(f2a):.1f}% of faces'})")

    col, hist = wl(deg, nbr)
    print(f"  WL vertex colours per round: {hist} -> {hist[-1]} stable classes")

    groups = collections.defaultdict(list)
    for f, a in S:
        groups[invariant(f, a, col, deg)].append((f, a))
    led = collections.defaultdict(list)
    for f, a in S:
        led[ledger(f, a, deg)].append((f, a))
    print(f"  move-site classes (WL lower bound): {len(groups)}"
          f"   distinct curvature ledgers: {len(led)}")

    print(f"\n  {'n':>8} {'%':>6}  {'face edges':<12} {'spokes':<22} "
          f"{'illegal after':<15} {'new6':>5} {'new7':>5}  representative")
    for (fe, sp), mem in sorted(led.items(), key=lambda kv: -len(kv[1])):
        d = describe(fe, sp)
        f, a = mem[0]
        print(f"  {len(mem):>8} {100*len(mem)/len(S):>5.1f}%  {str(fe):<12} "
              f"{str(sp):<22} {str(tuple(d['illegal_after'])):<15} "
              f"{d['new_deg6']:>5} {d['new_deg7']:>5}  face {f} apexes {a}")

    if full:
        print(f"\n  all {len(groups)} WL classes:")
        print(f"  {'#':>4} {'n':>7}  {'face edges':<12} {'spokes':<22} "
              f"representative")
        for i, (key, mem) in enumerate(
                sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))):
            f, a = mem[0]
            fe, sp = ledger(f, a, deg)
            print(f"  {i:>4} {len(mem):>7}  {str(fe):<12} {str(sp):<22} "
                  f"face {f} apexes {a}")

    if jsonout:
        recs = {"crystal": os.path.basename(path), "n_sites": len(S),
                "wl_rounds": hist, "n_classes": len(groups),
                "classes": []}
        for key, mem in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            f, a = mem[0]
            fe, sp = ledger(f, a, deg)
            recs["classes"].append(dict(count=len(mem), face=list(f),
                                        apexes=list(a), **describe(fe, sp)))
        with open(jsonout, "w") as fh:
            json.dump(recs, fh, indent=1)
        print(f"\n  wrote {jsonout}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("crystal", help="path to a .mfd triangulation")
    ap.add_argument("--full", action="store_true",
                    help="list every WL class, not just the ledger summary")
    ap.add_argument("--json", dest="jsonout", default=None,
                    help="write class representatives here")
    args = ap.parse_args()
    p = args.crystal if os.path.isabs(args.crystal) \
        else os.path.join(_ROOT, args.crystal)
    main(p, args.full, args.jsonout)
