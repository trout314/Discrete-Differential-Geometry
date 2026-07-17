#!/usr/bin/env python3
"""Exact local-crystallography matcher: label vertices by reference TCP phase.

Answers "is this patch of the triangulation POSITIONALLY a piece of crystal X
(A15 / C15 / C14 / sigma / R / P / mu / ...)?" — a strictly stronger question
than the Z-class census in fk_skeleton.py, which only measures composition.

Method (purely combinatorial, no geometry):
  For every vertex v, take the radius-r ball in the 1-skeleton and compute a
  canonical CERTIFICATE of the rooted, edge-degree-decorated ball: iterated
  Weisfeiler-Leman color refinement where a vertex's initial color is its
  distance from v, and every edge carries its FULL-manifold edge degree (so
  the certificate sees one ring beyond the ball through the degrees). Two
  vertices get the same certificate iff (up to WL power, which with these
  decorations is effectively exact here) their balls are isomorphic as
  decorated complexes.

  A reference dictionary for phase X = the set of certificates of ALL vertices
  of a perfect T^3 crystal of X (from scripts/tcp_reference.py). A perfect
  crystal has as many distinct certificates as Wyckoff-orbit classes (printed
  as a sanity check). A sample vertex "matches X at radius r" iff its
  certificate is in X's dictionary. Patches = connected components of vertices
  matching the same phase.

  r=1 is composition-level (decorated link); r=2 pins the local lattice and
  distinguishes e.g. cubic vs hexagonal Laves stacking; r=3 is essentially
  unique to the crystal. Matching is EXACT (no tolerance): one flipped edge
  inside the ball breaks the match, so fractions are conservative lower
  bounds on "crystalline" content.

Usage:
    python scripts/crystal_match.py sample1.mfd sample2.mfd --radius 2
    python scripts/crystal_match.py state.mfd --radius 1 2 --json out/match.json
"""
import argparse
import collections
import glob
import hashlib
import json
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets

REF_GLOB = os.path.join(_ROOT, "data", "tcp_reference", "T3_*_m*.mfd")


def _h(obj):
    return hashlib.blake2b(repr(obj).encode(), digest_size=10).hexdigest()


def adjacency(facets):
    eu, edeg, V = edges_from_facets(facets)
    adj = [[] for _ in range(V)]
    for (a, b), d in zip(eu, edeg):
        adj[a].append((b, int(d)))
        adj[b].append((a, int(d)))
    return [sorted(a) for a in adj], V


def certificate(root, adj, r):
    """Canonical certificate of the decorated radius-r ball rooted at root."""
    dist = {root: 0}
    frontier = [root]
    for d in range(1, r + 1):
        nxt = []
        for u in frontier:
            for w, _ in adj[u]:
                if w not in dist:
                    dist[w] = d
                    nxt.append(w)
        frontier = nxt
    nodes = list(dist)
    inball = set(nodes)
    # Within-ball incidences carry the FULL-manifold edge degree; edges leaving
    # the ball are summarized per-vertex by their degree multiset (cheap look
    # one ring further out without growing the ball).
    inc = {u: [(w, dg) for w, dg in adj[u] if w in inball] for u in nodes}
    out = {u: tuple(sorted(dg for w, dg in adj[u] if w not in inball))
           for u in nodes}
    col = {u: _h((dist[u], out[u])) for u in nodes}
    for _ in range(r + 2):
        col = {u: _h((col[u], tuple(sorted((dg, col[w]) for w, dg in inc[u]))))
               for u in nodes}
    return _h((col[root], tuple(sorted(col.values()))))


def phase_dict(facets, r):
    """All distinct vertex certificates of a (perfect) reference crystal."""
    adj, V = adjacency(facets)
    return collections.Counter(certificate(v, adj, r) for v in range(V))


def load_facets(path):
    return np.asarray(ddg.Manifold.load(path, 3).facets())


def best_refs(ref_glob):
    """Largest-m reference file per structure name."""
    best = {}
    for p in glob.glob(ref_glob):
        mm = re.match(r"T3_([A-Z0-9]+)_m(\d+)_N\d+\.mfd", os.path.basename(p))
        if not mm:
            continue
        name, m = mm.group(1).lower(), int(mm.group(2))
        if name not in best or m > best[name][0]:
            best[name] = (m, p)
    return {n: p for n, (m, p) in sorted(best.items())}


def patch_sizes(matched_vertices, adj):
    """Connected-component sizes among matched vertices (1-skeleton induced)."""
    mset = set(matched_vertices)
    seen, sizes = set(), []
    for s in matched_vertices:
        if s in seen:
            continue
        stack, comp = [s], 0
        seen.add(s)
        while stack:
            u = stack.pop()
            comp += 1
            for w, _ in adj[u]:
                if w in mset and w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(comp)
    return sorted(sizes, reverse=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("samples", nargs="+", help=".mfd files to label")
    ap.add_argument("--radius", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--refs", default=REF_GLOB,
                    help="glob of reference crystal .mfd files")
    ap.add_argument("--json", default=None, help="write full results here")
    args = ap.parse_args()

    refs = best_refs(args.refs)
    if not refs:
        sys.exit(f"no reference crystals match {args.refs} — run tcp_reference.py")
    results = {"refs": refs, "radii": args.radius, "samples": {}}

    dicts = {}   # (phase, r) -> set of certs
    for r in args.radius:
        print(f"[refs r={r}]", end=" ")
        for name, path in refs.items():
            cnt = phase_dict(load_facets(path), r)
            dicts[(name, r)] = set(cnt)
            print(f"{name}:{len(cnt)} orbits", end="  ")
        print()
        # cross-phase certificate overlaps (shared local order, e.g. Laves pair)
        for a in refs:
            for b in refs:
                if a < b and dicts[(a, r)] & dicts[(b, r)]:
                    print(f"  note: {a}/{b} share {len(dicts[(a, r)] & dicts[(b, r)])} "
                          f"certificate(s) at r={r}")

    for spath in args.samples:
        facets = load_facets(spath)
        adj, V = adjacency(facets)
        sres = {}
        for r in args.radius:
            certs = [certificate(v, adj, r) for v in range(V)]
            labels = {}
            for v, c in enumerate(certs):
                hit = [n for n in refs if c in dicts[(n, r)]]
                if hit:
                    labels[v] = hit
            uniq = collections.Counter(tuple(h) for h in labels.values())
            frac = {n: sum(1 for h in labels.values() if n in h) / V for n in refs}
            row = {"V": V, "frac_matched_any": len(labels) / V,
                   "frac_by_phase": frac,
                   "label_combos": {"+".join(k): v for k, v in uniq.items()}}
            for n in refs:
                mv = [v for v, h in labels.items() if h == [n]]
                ps = patch_sizes(mv, adj)
                row[f"patches_{n}"] = {"n": len(ps), "largest": ps[0] if ps else 0}
            sres[f"r{r}"] = row
            fr = " ".join(f"{n}:{frac[n]:.4f}" for n in refs)
            print(f"[{os.path.basename(spath)} r={r}] any:{len(labels)/V:.4f}  {fr}  "
                  + " ".join(f"big_{n}:{row[f'patches_{n}']['largest']}"
                             for n in refs if row[f"patches_{n}"]["largest"]))
        results["samples"][spath] = sres

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
