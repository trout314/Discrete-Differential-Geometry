#!/usr/bin/env python3
"""Aut-orbit counts per crystal, refined by the simplex's own decoration.

Every decoration used here is a function of the edge-degree labelling, and Aut
preserves edge degrees, so the decoration is CONSTANT ON AN ORBIT -- refining
the orbit count by it is a genuine partition, not a bucketing (the distinction
notes/memory/crystal-symmetry-group.md insists on: WL colours give bounds,
orbits are exact). Row totals therefore reproduce the plain orbit counts
exactly, which is asserted.

    vertex   Z = coordination (the FK class: Z12 / Z14 / Z15 / Z16)
    edge     hinge degree (5 = minor site, 6 = on the disclination network)
    face     sorted triple of its three edges' degrees. The three edges of a
             triangle are interchangeable under its own symmetry, so the
             MULTISET is the complete invariant -- no canonicalisation needed.
    tet      the six edge degrees. Here the multiset is NOT complete: two
             degree-6 edges of a tet may share a vertex or be opposite, and
             those are inequivalent. So the key is canonicalised over the 24
             vertex permutations of the tet, which is the same thing as
             naming the ISOMORPHISM CLASS OF THE DEGREE-6 SUBGRAPH of K4 --
             i.e. how the disclination lines thread that cell.

Usage:
    python scripts/orbit_decoration_census.py [--m3] [--json out.json]
"""
import argparse
import collections
import itertools
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))

from discrete_differential_geometry import Manifold
from discrete_differential_geometry.symmetry import CrystalSymmetry

E6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

# Reference cell per crystal (the crystal_gas_scan.LIBRARY sizes), plus a
# smaller-supercell set used to check the counts are box-independent.
CELLS = [("a15", 5, "T3_A15_m5_N5750.mfd"), ("c15", 4, "T3_C15_m4_N8704.mfd"),
         ("c14", 5, "T3_C14_m5_N8500.mfd"), ("c36", 4, "T3_C36_m4_N8704.mfd"),
         ("sigma", 4, "T3_SIGMA_m4_N11008.mfd"), ("z", 6, "T3_Z_m6_N8640.mfd"),
         ("mu", 4, "T3_MU_m4_N14208.mfd"), ("r", 2, "T3_R_m2_N7248.mfd"),
         ("p", 3, "T3_P_m3_N8640.mfd"), ("delta", 3, "T3_DELTA_m3_N8640.mfd")]
CELLS_M3 = [("a15", 3, "T3_A15_m3_N1242.mfd"), ("c15", 3, "T3_C15_m3_N3672.mfd"),
            ("c14", 3, "T3_C14_m3_N1836.mfd"),
            ("sigma", 3, "T3_SIGMA_m3_N4644.mfd"), ("r", 3, "T3_R_m3_N24462.mfd")]

# Isomorphism classes of subgraphs of K4, keyed by (n_edges, sorted vertex
# degrees of the subgraph). These are the ways the six-edge (disclination)
# network can thread a single tetrahedron.
K4_SHAPES = {
    (0, (0, 0, 0, 0)): "empty",
    (1, (0, 0, 1, 1)): "1 edge",
    (2, (0, 1, 1, 2)): "2 adj",       # two 6-edges sharing a vertex
    (2, (1, 1, 1, 1)): "2 opp",       # perfect matching: opposite edges
    (3, (0, 2, 2, 2)): "triangle",
    (3, (1, 1, 1, 3)): "star",        # three 6-edges at one vertex
    (3, (1, 1, 2, 2)): "path",
    (4, (2, 2, 2, 2)): "4-cycle",
    (4, (1, 2, 2, 3)): "paw",         # triangle + pendant
    (5, (2, 2, 3, 3)): "K4-e",
    (6, (3, 3, 3, 3)): "K4",
}


def tet_key(degs):
    """Canonical form of a tet's six edge degrees, over the 24 relabellings."""
    best = None
    for perm in itertools.permutations(range(4)):
        pos = {}
        for k, (i, j) in enumerate(E6):
            a, b = perm[i], perm[j]
            pos[(a, b) if a < b else (b, a)] = degs[k]
        cand = tuple(pos[e] for e in E6)
        if best is None or cand < best:
            best = cand
    return best


def tet_shape(degs):
    """Name the degree-6 subgraph of K4 inside this tet."""
    six = [E6[k] for k, d in enumerate(degs) if d >= 6]
    vd = collections.Counter()
    for a, b in six:
        vd[a] += 1
        vd[b] += 1
    sig = (len(six), tuple(sorted(vd.get(v, 0) for v in range(4))))
    return K4_SHAPES.get(sig, f"?{sig}")


def decorate(path):
    """(sym, {kind: [decoration of each orbit]}, {kind: orbit sizes})."""
    mfd = Manifold.load(path, 3)
    tets = np.asarray(mfd.facets())
    sym = CrystalSymmetry.for_manifold_path(path)

    edeg = collections.Counter()
    for tv in tets:
        tv = tuple(int(x) for x in tv)
        for i, j in E6:
            a, b = tv[i], tv[j]
            edeg[(a, b) if a < b else (b, a)] += 1
    adj = collections.defaultdict(set)
    for u, w in edeg:
        adj[u].add(w)
        adj[w].add(u)

    def edge_of(a, b):
        return edeg[(a, b) if a < b else (b, a)]

    out, sizes = {}, {}
    for kind in ("vertex", "edge", "face", "tet"):
        reps = sym.orbit_representatives(kind)
        sizes[kind] = sym.orbit_sizes(kind)
        lab = []
        for r in reps:
            if kind == "vertex":
                lab.append(f"Z{len(adj[int(r)])}")
            elif kind == "edge":
                lab.append(str(edge_of(int(r[0]), int(r[1]))))
            elif kind == "face":
                t = tuple(sorted(edge_of(*p) for p in
                                 itertools.combinations(sorted(int(x) for x in r), 2)))
                lab.append("".join(map(str, t)))
            else:
                v = sorted(int(x) for x in r)
                degs = [edge_of(v[i], v[j]) for i, j in E6]
                lab.append((tet_key(degs), tet_shape(degs)))
        out[kind] = lab
    return sym, out, sizes, edeg, adj, tets


def tally(labels, sizes):
    """{decoration: (n_orbits, n_simplices)} preserving first-seen order."""
    n = collections.Counter()
    tot = collections.Counter()
    for lb, sz in zip(labels, sizes):
        n[lb] += 1
        tot[lb] += sz
    return n, tot


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--m3", action="store_true",
                    help="use the smaller supercells (box-independence check)")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    cells = CELLS_M3 if args.m3 else CELLS

    data = {}
    for name, m, fn in cells:
        path = os.path.join(_ROOT, "data", "tcp_reference", fn)
        sym, lab, sizes, edeg, adj, tets = decorate(path)
        data[name] = dict(m=m, order=sym.order, lab=lab, sizes=sizes,
                          fv=[int(x) for x in Manifold.load(path, 3).f_vector])
        # -- assertions: the refinement must be exact, not approximate
        for kind in lab:
            assert len(lab[kind]) == sym.n_orbits(kind), (name, kind)
            assert sum(sizes[kind]) == data[name]["fv"][
                {"vertex": 0, "edge": 1, "face": 2, "tet": 3}[kind]], (name, kind)
        # -- cross-check the Z census against the raw structure
        zc = collections.Counter(f"Z{len(adj[v])}" for v in adj)
        zo = collections.Counter()
        for lb, sz in zip(lab["vertex"], sizes["vertex"]):
            zo[lb] += sz
        assert zc == zo, (name, zc, zo)

    def show(kind, title, keyfn, colfmt):
        cols = sorted({k for d in data.values() for k in d["lab"][kind]},
                      key=keyfn)
        print(f"\n### {title}  (cells = # Aut-orbits)")
        w = max(len(colfmt(c)) for c in cols) + 1
        print(f"{'crystal':7s} {'tot':>4s} | " +
              " ".join(f"{colfmt(c):>{w}s}" for c in cols))
        print("-" * (14 + (w + 1) * len(cols)))
        for name, d in data.items():
            n, _ = tally(d["lab"][kind], d["sizes"][kind])
            print(f"{name:7s} {sum(n.values()):4d} | " +
                  " ".join(f"{n.get(c, 0) or '.':>{w}}" for c in cols))

    print(f"Aut-orbit counts refined by decoration "
          f"({'m3 check set' if args.m3 else 'library cells'})")
    show("vertex", "VERTICES by coordination Z",
         lambda c: int(c[1:]), lambda c: c)
    show("edge", "EDGES by hinge degree", int, lambda c: f"deg{c}")
    show("face", "TRIANGLES by sorted edge-degree triple",
         lambda c: c, lambda c: c)
    show("tet", "TETS by the six edge degrees (canonical; = shape of the "
                "degree-6 subgraph)",
         lambda c: (sum(1 for x in c[0] if x >= 6), c[0]),
         lambda c: f"{''.join(map(str, c[0]))}\n{c[1]}".split("\n")[0])

    print("\n  tet column legend (n6 = # degree-6 edges in the cell):")
    cols = sorted({k for d in data.values() for k in d["lab"]["tet"]},
                  key=lambda c: (sum(1 for x in c[0] if x >= 6), c[0]))
    for c in cols:
        n6 = sum(1 for x in c[0] if x >= 6)
        print(f"    {''.join(map(str, c[0]))}  n6={n6}  {c[1]}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({k: {"m": v["m"], "order": v["order"], "fv": v["fv"],
                           "lab": {kk: [list(x) if isinstance(x, tuple) else x
                                        for x in vv]
                                   for kk, vv in v["lab"].items()},
                           "sizes": v["sizes"]}
                       for k, v in data.items()}, fh, indent=1)
        print(f"\n-> {args.json}")


if __name__ == "__main__":
    main()
