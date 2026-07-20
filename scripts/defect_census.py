#!/usr/bin/env python3
"""General crystal-defect census by directional crystal-patch growth.

A defect is a vertex that does not fit onto the ideal crystal --- but detected
so it stays LOCALIZED to the misplaced vertex, not its perturbed neighbourhood.
Naive "vertex whose radius-r environment isn't crystalline" spreads: one point
defect contaminates every certificate within radius r, so at finite dopant
density the whole crystal flags. Instead we GROW the crystal patch and heal
inward (Aaron's algorithm, one-sided accretion):

  local environment  = (own Z-class, multiset of (neighbour Z-class, edge deg))
  1. SEED  crystal = vertices whose FULL local environment is a native one.
  2. GROW  absorb an unknown vertex into the crystal iff the multiset of its
     ALREADY-CRYSTAL neighbours is a sub-multiset of some native environment
     of its class (treating unknown/defect neighbours as don't-cares) -- so a
     good vertex next to a dopant still fits on its crystal-facing side.
  3. Repeat to a fixpoint; a fully-surrounded crystal vertex whose complete
     environment isn't native is demoted (catches early over-absorption).
  4. Leftover = defect CORES; connected components = dopant complexes.

Class-agnostic: a native Z14 in its Wyckoff site fits, a doped Z14 in the
wrong site does not, so it separates native from excess even when the dopant
IS a native class -- the case the Z-class census can't handle. r=1 local
descriptor => fast (no radius-2 WL certificates).

Usage:
    python scripts/defect_census.py STATE.mfd --ref c15
    python scripts/defect_census.py 'data/r_solubility/r_sol_Z14_mu4_r*_final.mfd' --ref r
"""
import argparse
import glob
import os
import sys
from collections import Counter, defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import discrete_differential_geometry as ddg
from tcp_melt import CRYSTALS
from dopant_pairs import vertex_classes
from fk_skeleton import edges_from_facets

N6_NAME = {-1: "illegal", 0: "Z12", 2: "Z14", 3: "Z15", 4: "Z16"}


def patch_sizes(vertices, adj):
    """Connected-component sizes of the subgraph induced by `vertices` on the
    plain neighbour-list adjacency `adj` (largest first)."""
    vset = set(vertices)
    seen = set()
    sizes = []
    for s in vertices:
        if s in seen:
            continue
        stack, n = [s], 0
        seen.add(s)
        while stack:
            u = stack.pop()
            n += 1
            for w in adj[u]:
                if w in vset and w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(n)
    return sorted(sizes, reverse=True)


def descriptors(facets):
    """Per-vertex helpers: vclass(v), adj, and nbr(v, mask) -> sorted multiset
    of (neighbour class, connecting-edge degree) over neighbours w with
    mask(w) True (mask=None => all)."""
    n6, imp, adj = vertex_classes(facets)
    eu, edeg, _ = edges_from_facets(facets)
    edge_d = {}
    for (a, b), d in zip(eu, edeg):
        edge_d[(int(a), int(b))] = int(d)

    def vclass(v):
        return -1 if imp[v] else int(n6[v])

    def ed(v, w):
        return edge_d[(v, w)] if v < w else edge_d[(w, v)]

    def nbr(v, mask=None):
        ws = adj[v] if mask is None else [w for w in adj[v] if mask(w)]
        return tuple(sorted((vclass(w), ed(v, w)) for w in ws))
    return vclass, adj, nbr, len(n6)


def defect_cores(facets, ref_facets):
    rv, radj, rnbr, RV = descriptors(ref_facets)
    native_envs = set()                       # (class, full neighbour tuple)
    native_by_class = defaultdict(list)       # class -> list of Counter
    for u in range(RV):
        env = rnbr(u)
        native_envs.add((rv(u), env))
        native_by_class[rv(u)].append(Counter(env))

    vclass, adj, nbr, V = descriptors(facets)
    CRYSTAL, UNKNOWN = 1, 0
    label = np.zeros(V, dtype=np.int8)
    for v in range(V):                        # 1. seed: exact native env
        if (vclass(v), nbr(v)) in native_envs:
            label[v] = CRYSTAL

    def fits(v):                              # crystal-facing sub-match
        cn = Counter(nbr(v, mask=lambda w: label[w] == CRYSTAL))
        return any(all(cn[k] <= M[k] for k in cn)
                   for M in native_by_class.get(vclass(v), []))

    changed = True                            # 2-3. grow + demote to fixpoint
    while changed:
        changed = False
        for v in range(V):
            if label[v] == UNKNOWN and any(label[w] == CRYSTAL for w in adj[v]) \
                    and fits(v):
                label[v] = CRYSTAL
                changed = True
        for v in range(V):                    # demote over-absorbed cores
            if label[v] == CRYSTAL and all(label[w] == CRYSTAL for w in adj[v]) \
                    and (vclass(v), nbr(v)) not in native_envs:
                label[v] = UNKNOWN
                changed = True

    dv = [v for v in range(V) if label[v] == UNKNOWN]
    return dv, vclass, adj, V


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("states", nargs="+")
    ap.add_argument("--ref", required=True, choices=list(CRYSTALS))
    args = ap.parse_args()
    files = []
    for s in args.states:
        files += sorted(glob.glob(s)) if any(c in s for c in "*?[") else [s]
    ref_facets = np.asarray(ddg.Manifold.load(CRYSTALS[args.ref], 3).facets())
    print(f"reference: {args.ref}")
    for f in files:
        facets = np.asarray(ddg.Manifold.load(f, 3).facets())
        dv, vclass, adj, V = defect_cores(facets, ref_facets)
        comp = Counter(N6_NAME.get(vclass(v), f"n6={vclass(v)}") for v in dv)
        sizes = patch_sizes(list(dv), adj)
        szh = Counter(sizes)
        print(f"\n{os.path.basename(f)}")
        print(f"  defect cores: {len(dv)}/{V} ({100*len(dv)/V:.2f}%)   "
              f"classes: " + " ".join(f"{k}:{v}" for k, v in comp.most_common()))
        print(f"  complexes: {len(sizes)} clusters, max {sizes[0] if sizes else 0}"
              f"   sizes[size×count]: "
              + " ".join(f"{k}×{v}" for k, v in sorted(szh.items())))


if __name__ == "__main__":
    main()
