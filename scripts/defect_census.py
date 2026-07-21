#!/usr/bin/env python3
"""Crystal-defect census: a reporting front-end over the crystal_grains identifier.

Defect cores are the vertices that do NOT sit in a genuine crystalline grain, as
determined by the single authoritative crystallinity detector, crystal_grains.py
(a single-valued covering map onto the reference TCP crystal). A vertex is
crystalline iff it develops into a registry-consistent grain of >= --min-size;
everything else is a defect core. Because the covering map preserves BOTH the
FK Z-class and all six edge degrees AND demands translational registry, it flags
EVERY departure from the ideal crystal -- a native Z-class atom in the WRONG site
(forced onto the reference vertex at its lattice position, the mismatch drops it
out), dislocations / coherent off-registry order where every atom is locally
native yet the lattice fails to close, and foreign phases -- while staying
localized (the FK-shell gate breaks only for the vertices a defect actually
touches, not their neighbours). See crystal_grains.py for the method.

This script adds reporting on top of that vertex set: the defect cores broken
down by Z-class (Z12/Z14/Z15/Z16/illegal) and clustered into connected complexes
(dopant cores / grain boundaries), plus the crystalline fraction and grain sizes.

Usage:
    python scripts/defect_census.py STATE.mfd --ref c15
    python scripts/defect_census.py STATE.mfd --ref r --min-size 30
    python scripts/defect_census.py 'data/r_solubility/r_sol_Z14_mu4_r*_final.mfd' --ref r
"""
import argparse
import glob
import os
import sys
from collections import Counter

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import discrete_differential_geometry as ddg
import crystal_grains as cg
from tcp_melt import CRYSTALS
from tcp_reference import STRUCTURES
from dopant_pairs import vertex_classes

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
    """Per-vertex helpers for defect-core REPORTING: vclass(v) (FK Z-class, -1 if
    illegal) and the plain neighbour adjacency. Vertex ids match crystal_grains:
    both relabel a facet array by the identical np.unique(F) sort, so a defect id
    from `defect_cores` indexes straight into these."""
    n6, imp, adj = vertex_classes(facets)

    def vclass(v):
        return -1 if imp[v] else int(n6[v])

    return vclass, adj, len(n6)


def defect_cores(facets, ref_facets, ref_name, min_size):
    """Defect cores = vertices NOT in a translational-crystalline grain of
    >= min_size, via the crystal_grains covering map. Returns
    (dv, vclass, adj, V, grain_sizes)."""
    rst = cg.build_struct(ref_facets)
    refs = {ref_name: rst}
    idx = cg.ref_index(refs)
    ns_of = {}                                            # translation-invariant site
    if ref_name in STRUCTURES:
        ns = len(STRUCTURES[ref_name][1])
        if ns and rst["V"] % ns == 0:
            ns_of[ref_name] = ns
    st = cg.build_struct(facets)
    grain_of_tet, sig_of_tet, phase_of_grain = cg.find_grains(st, refs, idx)
    interior = cg.interior_vertices(st, grain_of_tet, sig_of_tet,
                                    phase_of_grain, ns_of)
    crystalline = set()
    grain_sizes = []
    for _, verts in cg.grain_components(st, interior):
        if len(verts) >= min_size:
            crystalline.update(verts)
            grain_sizes.append(len(verts))
    vclass, adj, V = descriptors(facets)
    dv = [v for v in range(V) if v not in crystalline]
    return dv, vclass, adj, V, sorted(grain_sizes, reverse=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("states", nargs="+")
    ap.add_argument("--ref", required=True, choices=list(CRYSTALS))
    ap.add_argument("--min-size", type=int, default=30,
                    help="smallest crystalline grain kept; vertices in "
                         "sub-threshold grains count as defects. Default 30 "
                         "clears the coincidental-grain ceiling of a melted null "
                         "(~18 at V~1e3, and the ceiling grows only ~ln V); for "
                         "V >> 1e5 raise it, or calibrate against a matched melt.")
    args = ap.parse_args()
    files = []
    for s in args.states:
        files += sorted(glob.glob(s)) if any(c in s for c in "*?[") else [s]
    ref_facets = np.asarray(ddg.Manifold.load(CRYSTALS[args.ref], 3).facets())
    print(f"reference: {args.ref}   min-size: {args.min_size}")
    for f in files:
        facets = np.asarray(ddg.Manifold.load(f, 3).facets())
        dv, vclass, adj, V, grain_sizes = defect_cores(
            facets, ref_facets, args.ref, args.min_size)
        comp = Counter(N6_NAME.get(vclass(v), f"n6={vclass(v)}") for v in dv)
        sizes = patch_sizes(list(dv), adj)
        szh = Counter(sizes)
        ncry = V - len(dv)
        print(f"\n{os.path.basename(f)}")
        print(f"  crystalline: {ncry}/{V} ({100 * ncry / V:.2f}%)   "
              f"grains (>= {args.min_size}): {len(grain_sizes)}   "
              f"largest {grain_sizes[0] if grain_sizes else 0}")
        print(f"  defect cores: {len(dv)}/{V} ({100 * len(dv) / V:.2f}%)   "
              f"classes: " + " ".join(f"{k}:{v}" for k, v in comp.most_common()))
        print(f"  complexes: {len(sizes)} clusters, max {sizes[0] if sizes else 0}"
              f"   sizes[size×count]: "
              + " ".join(f"{k}×{v}" for k, v in sorted(szh.items())))


if __name__ == "__main__":
    main()
