#!/usr/bin/env python3
"""Dopant-dopant pair correlations in doped TCP crystals.

Companion to the doping experiments (tcp_melt.py --tilt): given a doped
crystal state, find the dopant vertices (a Z-class absent from the perfect
host, e.g. Z16 in A15) and measure whether they attract (precipitation /
second-phase nucleation), repel (solid solution with screened interactions),
or ignore each other (ideal dilute gas), via BFS graph distances.

Observables, each against a shuffled null (random vertex sets of the same
size, M draws):
  * pair-distance histogram -> g(r) = observed pairs at r / null mean
  * nearest-neighbour distance distribution (most sensitive at low counts)

Usage:
    python scripts/dopant_pairs.py data/dope_states/a15big_z16_mu*.mfd \
        --dopant-class Z16 --json out/dopant_pairs.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets

CLASS_N6 = {"Z12": 0, "Z14": 2, "Z15": 3, "Z16": 4}


def vertex_classes(facets):
    """Per-vertex (n6, m_impure) and adjacency list."""
    eu, edeg, V = edges_from_facets(facets)
    n6 = np.zeros(V, int)
    imp = np.zeros(V, int)
    adj = [[] for _ in range(V)]
    for (a, b), d in zip(eu, edeg):
        adj[a].append(b)
        adj[b].append(a)
        if d >= 6:
            n6[a] += 1
            n6[b] += 1
        if d < 5 or d > 6:
            imp[a] += 1
            imp[b] += 1
    return n6, imp, adj


def bfs_dists(src, adj, V):
    d = np.full(V, -1, int)
    d[src] = 0
    q = collections.deque([src])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if d[w] < 0:
                d[w] = d[u] + 1
                q.append(w)
    return d


def cluster_census(verts, adj):
    """Connected components of the dopant-adjacency graph: (size, edges) counts.
    edges == size => single cycle; size-1 => tree/chain; more => dense knot."""
    vs = set(verts)
    seen, comps = set(), collections.Counter()
    for s in verts:
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w in vs and w not in seen:
                    seen.add(w)
                    stack.append(w)
        cs = set(comp)
        ne = sum(1 for u in comp for w in adj[u] if w in cs) // 2
        comps[(len(comp), ne)] += 1
    return comps


def pair_stats(verts, adj, V, rmax):
    """(pair-distance histogram up to rmax, NN distances) for a vertex set."""
    verts = list(verts)
    hist = np.zeros(rmax + 1, int)
    nn = []
    for i, s in enumerate(verts):
        d = bfs_dists(s, adj, V)
        ds = d[verts[i + 1:]]
        for x in ds[(ds >= 0) & (ds <= rmax)]:
            hist[x] += 1
        dall = d[[v for v in verts if v != s]]
        if len(dall):
            nn.append(int(dall[dall >= 0].min()))
    return hist, np.array(nn)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--dopant-class", default="Z16", choices=list(CLASS_N6))
    ap.add_argument("--rmax", type=int, default=10)
    ap.add_argument("--shuffles", type=int, default=200)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(7)
    out = {}
    for path in args.samples:
        facets = np.asarray(ddg.Manifold.load(path, 3).facets())
        n6, imp, adj = vertex_classes(facets)
        V = len(n6)
        dop = np.where((n6 == CLASS_N6[args.dopant_class]) & (imp == 0))[0]
        nd = len(dop)
        name = os.path.basename(path)
        if nd < 4:
            print(f"[{name}] only {nd} dopants — too few for statistics, skipped")
            out[name] = {"n_dopants": nd}
            continue

        hist, nn = pair_stats(dop, adj, V, args.rmax)
        null_h = np.zeros((args.shuffles, args.rmax + 1), int)
        null_nn = []
        for k in range(args.shuffles):
            rand = rng.choice(V, size=nd, replace=False)
            h, n = pair_stats(rand, adj, V, args.rmax)
            null_h[k] = h
            null_nn.append(n.mean())
        nm, ns = null_h.mean(0), null_h.std(0)
        g = np.where(nm > 0, hist / np.maximum(nm, 1e-12), np.nan)
        z = np.where(ns > 0, (hist - nm) / np.maximum(ns, 1e-12), np.nan)
        nn_null_m, nn_null_s = np.mean(null_nn), np.std(null_nn)

        print(f"[{name}] V={V} dopants={nd} ({nd / V:.4f})")
        print(f"  r:      " + " ".join(f"{r:>6d}" for r in range(1, args.rmax + 1)))
        print(f"  pairs:  " + " ".join(f"{hist[r]:>6d}" for r in range(1, args.rmax + 1)))
        print(f"  g(r):   " + " ".join(f"{g[r]:>6.2f}" for r in range(1, args.rmax + 1)))
        print(f"  z(r):   " + " ".join(f"{z[r]:>6.1f}" for r in range(1, args.rmax + 1)))
        print(f"  <NN>={nn.mean():.2f} vs null {nn_null_m:.2f}±{nn_null_s:.2f} "
              f"(z={(nn.mean() - nn_null_m) / max(nn_null_s, 1e-12):+.1f}; "
              f"negative = clustering)")
        cc = cluster_census(dop, adj)
        n_multi = sum((sz * c) for (sz, ne), c in cc.items() if sz > 1)
        mean_sz = nd / sum(cc.values())
        print(f"  clusters: " + ", ".join(
            f"(sz{sz},e{ne})x{c}" for (sz, ne), c in sorted(cc.items()))
            + f"  | mean size {mean_sz:.2f}, {n_multi}/{nd} in multiplets")
        out[name] = {
            "V": V, "n_dopants": nd, "pairs": hist.tolist(),
            "null_mean": nm.tolist(), "null_std": ns.tolist(),
            "g": [None if np.isnan(x) else x for x in g],
            "nn_mean": float(nn.mean()), "nn_null_mean": float(nn_null_m),
            "nn_null_std": float(nn_null_s),
            "clusters": {f"{sz},{ne}": c for (sz, ne), c in sorted(cc.items())},
        }

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
