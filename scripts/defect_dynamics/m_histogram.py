#!/usr/bin/env python3
"""m-histogram: what impurity valence does an ISOLATED defect carry, versus a
FUSED one?

This sets m* in the offset penalty V(m) = lambda * max(0, m - m*)^2
(sampler.VertexPot.impOffset). The design intent:

    isolated defect   pays NOTHING for its own internal structure
    two defects touching   pays little, so collisions are affordable
    deep burial       still forbidden

Note which way the overlap has to go. m* works if there is a gap between
"one or two defects locally" and "many piled up". Overlap between the
isolated and the touching distributions is GOOD -- it is what makes a contact
free. (An earlier note of mine had this backwards.) What must separate is
small-vs-large, not isolated-vs-contact.

So: m(v) = #incident edges of degree outside {5,6}, histogrammed by the SIZE
of the complex the vertex belongs to, over thermal gas snapshots. Complexes
are components of the subgraph induced on defect vertices by all edges -- the
convention that reproduces defect_catalog's complex counts and its
elementary sig=(3,4,4) population (the illegal-edge-graph convention
fragments them; see notes/memory/flight-contact-barrier.md).

Reported per class: the m distribution, its max, and -- the actionable part --
the fraction of that class's vertices that would be FREE (m <= m*) at each
candidate m*, plus the residual energy each class would still pay.

Usage:
    python scripts/defect_dynamics/m_histogram.py data/mgas/lam40*_snap*.mfd \
        data/mgas/l40s20*_snap*.mfd --out data/ecmc_ab/m_histogram.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_ROOT, "python"))

import discrete_differential_geometry as ddg
from discrete_differential_geometry.vertex_fields import edges_and_degrees

MSTARS = [0, 1, 2, 3, 4]
# classes by complex size in ILLEGAL EDGES; the elementary (3,4,4) flicker has 3
CLASSES = [(1, 3, "elementary (<=3 ill edges)"),
           (4, 7, "small fused (4-7)"),
           (8, 15, "medium (8-15)"),
           (16, 10 ** 9, "large / buried (>=16)")]


def components(pairs, verts):
    adj = collections.defaultdict(set)
    for u, v in pairs:
        adj[u].add(v)
        adj[v].add(u)
    seen, out = set(), []
    for s0 in verts:
        if s0 in seen:
            continue
        stack, cur = [s0], set()
        while stack:
            x = stack.pop()
            if x in cur:
                continue
            cur.add(x)
            stack.extend(adj[x] - cur)
        seen |= cur
        out.append(cur)
    return out


def one(path):
    F = np.asarray(ddg.Manifold.load(path, 3).facets())
    eu, ecnt, deg, V = edges_and_degrees(F)
    ill = (ecnt < 5) | (ecnt > 6)
    ille = [tuple(map(int, e)) for e in eu[ill]]
    if not ille:
        return []
    m = collections.Counter()
    for a, b in ille:
        m[a] += 1
        m[b] += 1
    defect = set(m)
    ind = [tuple(map(int, e)) for e in eu
           if int(e[0]) in defect and int(e[1]) in defect]
    comps = components(ind, defect)
    rows = []
    for c in comps:
        n_ill = sum(1 for a, b in ille if a in c and b in c)
        for v in c:
            rows.append({"m": m[v], "csize": n_ill})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = []
    for p in args.snapshots:
        try:
            rows += one(p)
        except Exception as ex:                        # noqa: BLE001
            print(f"  SKIP {os.path.basename(p)}: {ex}")
    if not rows:
        print("no defect vertices found")
        return
    print(f"{len(args.snapshots)} snapshots, {len(rows):,} defect vertices\n")

    mmax = max(r["m"] for r in rows)
    print(f"{'class':<28} {'n_vert':>7}  " +
          "  ".join(f"m={k}" for k in range(1, mmax + 1)) + "   max")
    summary = {}
    per_class = {}
    for lo, hi, lab in CLASSES:
        g = [r["m"] for r in rows if lo <= r["csize"] <= hi]
        if not g:
            continue
        per_class[lab] = g
        h = collections.Counter(g)
        cells = "  ".join(f"{100*h[k]/len(g):4.0f}%" for k in range(1, mmax + 1))
        print(f"{lab:<28} {len(g):>7,}  {cells}   {max(g)}")
        summary[lab] = {"n": len(g), "max": int(max(g)),
                        "mean": float(np.mean(g)),
                        "hist": {str(k): h[k] for k in sorted(h)}}

    print(f"\nfraction of each class FREE (m <= m*) under "
          f"V(m) = lambda*max(0, m-m*)^2:")
    print(f"{'class':<28} " + "  ".join(f"m*={s}" for s in MSTARS))
    for lab, g in per_class.items():
        cells = "  ".join(f"{100*np.mean([x <= s for x in g]):5.0f}%"
                          for s in MSTARS)
        print(f"{lab:<28} {cells}")

    print(f"\nresidual energy per vertex, lambda=1 "
          f"(mean of max(0, m-m*)^2) -- the SEPARATION m* buys:")
    print(f"{'class':<28} " + "  ".join(f"m*={s}" for s in MSTARS))
    for lab, g in per_class.items():
        cells = "  ".join(f"{np.mean([max(0, x - s) ** 2 for x in g]):6.2f}"
                          for s in MSTARS)
        print(f"{lab:<28} {cells}")
        summary[lab]["residual"] = {
            str(s): float(np.mean([max(0, x - s) ** 2 for x in g]))
            for s in MSTARS}

    # the ratio that matters: burial cost / elementary cost at each m*
    el = per_class.get("elementary (<=3 ill edges)")
    bl = per_class.get("large / buried (>=16)")
    if el and bl:
        print(f"\ncontrast (large/elementary residual) -- higher = m* separates "
              f"burial from an isolated defect better:")
        for s in MSTARS:
            e = np.mean([max(0, x - s) ** 2 for x in el])
            b = np.mean([max(0, x - s) ** 2 for x in bl])
            print(f"   m*={s}:  {b/e:6.1f}x" if e > 0
                  else f"   m*={s}:  inf (elementary is exactly free)")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"summary": summary, "n_snapshots": len(args.snapshots)},
                      f, indent=1, default=float)
        print(f"\nwrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
