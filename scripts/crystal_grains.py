#!/usr/bin/env python3
"""Crystalline-grain detection by covering-map development.

Identifies the connected components of the vertex graph that are *consistent
with a piece of a reference TCP crystal* -- i.e. genuine crystalline grains,
not merely regions of locally-crystal-like coordination.

Why not per-vertex local matching (crystal_match.py / defect_census.py):
  A single FK coordination polyhedron (Z12/Z14/Z15/Z16) is the universal local
  motif of *all* tetrahedrally-close-packed matter, crystalline AND amorphous,
  so any radius-r ball test is either too weak (an isolated liquid icosahedron
  matches -- defect_census r=1) or too brittle (one perturbed edge in a ~60-
  vertex r=2 ball breaks the match, so a 62%-ordered melt reads 0% crystal --
  crystal_match r=2, exact). Crystallinity is *translational registry* --
  vertices sitting in the crystal's actual periodic arrangement over an extended
  region -- which is inherently non-local.

Method (combinatorial covering map / "development"):
  A grain is a maximal region carrying a single-valued LOCAL ISOMORPHISM
  (covering map) onto the reference crystal. We work at the tetrahedron level:
    * A sample tet maps to a reference tet under a vertex correspondence that
      preserves decorations (vertex FK-class + all six hinge/edge degrees).
    * SEED: pick a sample tet, match it to some reference tet (fixing a
      registry + orientation).
    * DEVELOP (BFS across shared faces): crossing a face into a neighbour tet
      FORCES its image (three shared vertices keep their map; the apex maps to
      the reference apex across the corresponding reference face). Accept the
      neighbour iff its decorations match the forced reference tet; a mismatch
      or a single-valued conflict (holonomy) makes that face a grain BOUNDARY.
  Consistency propagation is what recovers r>=2 (and phase-exact)
  discrimination from an r=1 local test -- e.g. cubic C15 vs hexagonal C14
  share every link but their stackings make the map conflict within a step or
  two. Yet a covered vertex needs only its OWN link clean (not a pristine
  r=2 ball), so the map grows right up to defects: phase-exact AND tolerant.

  A vertex is INTERIOR-crystalline iff its entire star develops consistently
  into one grain (single-valued image) -- provably link-isomorphic to its
  reference vertex. Grains = connected components of interior vertices sharing
  a grain id. Everything else is defect/boundary. `--min-size` drops
  coincidental sub-threshold grains (calibrate against a fully-melted null,
  which yields ~0 grains).

  Registry is tracked only into the finite reference torus (single-valued; no
  continuous positions). Dislocation-type defects, where local structure is
  perfect but global registry fails to close, correctly appear as grain
  boundaries rather than being classified.

Usage:
    python scripts/crystal_grains.py STATE.mfd
    python scripts/crystal_grains.py STATE.mfd --ref c15 r --min-size 20
    python scripts/crystal_grains.py 'data/melt_test/*.mfd' --json out/grains.json
"""
import argparse
import collections
import glob
import os
import sys
from itertools import permutations

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes
from crystal_match import best_refs, REF_GLOB
from tcp_reference import STRUCTURES

# tet vertex-index pairs for the six edges
_E6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def build_struct(facets):
    """Combinatorial structure of a triangulation, in RELABELLED vertex ids
    (0..V-1, matching vertex_classes/edges_from_facets). Returns a dict with:
      tets   : (nT,4) int array of tet vertices
      typ    : (V,) int, vertex FK-class (-1 impure, else n6 code 0/2/3/4)
      edg    : dict (u<w) -> hinge degree
      face2  : dict sorted-triple -> list of (tet_id, apex_vertex)
      tetid  : dict sorted-4-tuple -> tet_id
      adj    : relabelled neighbour lists
      pure   : (V,) bool, not impure (used to seed interiors first)
    """
    F = np.asarray(facets, np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    tets = inv.reshape(F.shape)                       # relabelled 0..V-1
    n6, imp, adj = vertex_classes(facets)
    typ = np.where(imp > 0, -1, n6).astype(int)
    eu, edeg, V = edges_from_facets(facets)
    edg = {(int(a), int(b)): int(d) for (a, b), d in zip(eu, edeg)}
    face2 = collections.defaultdict(list)
    tetid = {}
    for t, tv in enumerate(tets):
        tv = tuple(int(x) for x in tv)
        tetid[tuple(sorted(tv))] = t
        for apex in tv:
            face = tuple(sorted(v for v in tv if v != apex))
            face2[face].append((t, apex))
    return dict(tets=tets, typ=typ, edg=edg, face2=face2, tetid=tetid,
                adj=adj, pure=(imp == 0), V=V, nT=len(tets))


def _edge(st, u, w):
    return st["edg"][(u, w)] if u < w else st["edg"][(w, u)]


def correspondences(st, sT, rt, rst):
    """All bijections {sample vert -> ref vert} taking sample tet `sT` onto ref
    tet `rt` (4-tuples) preserving vertex class + all six edge degrees."""
    styp, rtyp = st["typ"], rst["typ"]
    out = []
    for p in permutations(rt):
        if any(styp[sT[i]] != rtyp[p[i]] for i in range(4)):
            continue
        if all(_edge(st, sT[i], sT[j]) == _edge(rst, p[i], p[j]) for i, j in _E6):
            out.append({sT[i]: p[i] for i in range(4)})
    return out


def _tet_ok(st, rst, sT, sig):
    """Do sample tet `sT` (4 verts) and its forced ref image agree on all
    decorations under correspondence `sig`?"""
    if any(st["typ"][v] != rst["typ"][sig[v]] for v in sT):
        return False
    return all(_edge(st, sT[i], sT[j]) == _edge(rst, sig[sT[i]], sig[sT[j]])
               for i, j in _E6)


def canon(st, tv):
    """Canonical decorated-K4 form of a tet (for the reference seed index)."""
    best = None
    for p in permutations(tv):
        key = (tuple(int(st["typ"][v]) for v in p),
               tuple(_edge(st, p[i], p[j]) for i, j in _E6))
        if best is None or key < best:
            best = key
    return best


def ref_index(refs):
    """canon -> list of (phase, ref_tet_id): one representative per (canon,
    phase) is enough (symmetry-equivalent tets develop congruently)."""
    idx = collections.defaultdict(list)
    seen = set()
    for phase, rst in refs.items():
        for t, tv in enumerate(rst["tets"]):
            c = canon(rst, tuple(int(x) for x in tv))
            if (phase, c) in seen:
                continue
            seen.add((phase, c))
            idx[c].append((phase, t))
    return idx


def develop(st, rst, seed_tet, sig0):
    """Grow a single-valued covering map from `seed_tet` (sample tet id) with
    initial correspondence `sig0`. Returns {sample tet id -> sig} for every tet
    consistently reached (the grain, at tet granularity)."""
    seed_rt = tuple(sorted(sig0[int(v)] for v in st["tets"][seed_tet]))
    rtid0 = rst["tetid"][seed_rt]
    assign = {seed_tet: (rtid0, sig0)}
    q = collections.deque([seed_tet])
    while q:
        t = q.popleft()
        rtid, sig = assign[t]
        tv = [int(v) for v in st["tets"][t]]
        for apex in tv:
            face = [v for v in tv if v != apex]
            snb = [x for x in st["face2"][tuple(sorted(face))] if x[0] != t]
            if not snb:
                continue
            t2, apex2 = snb[0]
            if t2 in assign:
                continue                                   # already fixed
            rface = tuple(sorted(sig[v] for v in face))
            rnb = [x for x in rst["face2"][rface] if x[0] != rtid]
            if not rnb:
                continue
            rtid2, rapex2 = rnb[0]
            sig2 = {v: sig[v] for v in face}
            sig2[apex2] = rapex2
            if _tet_ok(st, rst, [int(v) for v in st["tets"][t2]], sig2):
                assign[t2] = (rtid2, sig2)
                q.append(t2)
    return assign


def find_grains(st, refs, idx):
    """Label sample tets by grain. Returns (grain_of_tet, sig_of_tet,
    phase_of_grain). Seeds interiors first (pure tets), so one seed grabs a
    whole grain; boundary/defect tets are reached and validated during a
    neighbour's development or left unassigned."""
    grain_of_tet = {}
    sig_of_tet = {}
    phase_of_grain = []
    order = sorted(range(st["nT"]),
                   key=lambda t: -int(st["pure"][st["tets"][t]].sum()))
    for t in order:
        if t in grain_of_tet:
            continue
        tv = tuple(int(v) for v in st["tets"][t])
        cands = idx.get(canon(st, tv), [])
        best = None
        for phase, rtid in cands:
            rst = refs[phase]
            rt = tuple(int(x) for x in rst["tets"][rtid])
            for sig0 in correspondences(st, tv, rt, rst):
                a = develop(st, rst, t, sig0)
                if best is None or len(a) > len(best[1]):
                    best = (phase, a)
        if best is None:
            continue                                       # unseedable -> defect
        gid = len(phase_of_grain)
        phase, assign = best
        phase_of_grain.append(phase)
        for tid, (_, sig) in assign.items():
            if tid not in grain_of_tet:
                grain_of_tet[tid] = gid
                sig_of_tet[tid] = sig
    return grain_of_tet, sig_of_tet, phase_of_grain


def interior_vertices(st, grain_of_tet, sig_of_tet, phase_of_grain, ns_of):
    """v -> grain id for vertices whose FULL star develops consistently into one
    grain with a single-valued reference image MODULO the unit-cell lattice
    (=> link-isomorphic to the reference vertex, in a consistent registry).

    Registry is compared by unit-cell SITE = ref_vertex % ns (ns = atoms per
    cell), which is translation-invariant: reaching v around the sample torus at
    reference images differing by a lattice translation (same site) is genuine
    crystal, not a conflict. Only a true stacking/phase inconsistency changes the
    site -- that's what excludes C14 tested against C15, and isolated motifs."""
    tets_of = collections.defaultdict(list)
    for t, tv in enumerate(st["tets"]):
        for v in tv:
            tets_of[int(v)].append(t)
    interior = {}
    for v, tl in tets_of.items():
        gids = {grain_of_tet.get(t, -1) for t in tl}
        if len(gids) != 1 or -1 in gids:
            continue                                       # boundary/defect
        gid = gids.pop()
        ns = ns_of.get(phase_of_grain[gid], 0)
        sites = {sig_of_tet[t][v] % ns if ns else sig_of_tet[t][v] for t in tl}
        if len(sites) == 1:                                # consistent registry
            interior[v] = gid
    return interior


def grain_components(st, interior):
    """Connected components (1-skeleton) of interior vertices sharing a grain
    id. Returns list of (grain_id, [vertices]), largest first."""
    seen, comps = set(), []
    for s in interior:
        if s in seen:
            continue
        gid = interior[s]
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in st["adj"][u]:
                if w not in seen and interior.get(w, -2) == gid:
                    seen.add(w)
                    stack.append(w)
        comps.append((gid, comp))
    comps.sort(key=lambda c: -len(c[1]))
    return comps


def analyze(facets, refs, idx, min_size, ns_of):
    st = build_struct(facets)
    grain_of_tet, sig_of_tet, phase_of_grain = find_grains(st, refs, idx)
    interior = interior_vertices(st, grain_of_tet, sig_of_tet,
                                 phase_of_grain, ns_of)
    comps = grain_components(st, interior)
    kept = [(g, c) for g, c in comps if len(c) >= min_size]
    by_phase = collections.Counter()
    for g, c in kept:
        by_phase[phase_of_grain[g]] += len(c)
    return dict(
        V=st["V"], n_interior=len(interior),
        n_grains=len(kept), grain_sizes=[len(c) for _, c in kept],
        by_phase=dict(by_phase),
        phase_of=[phase_of_grain[g] for g, _ in kept])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("states", nargs="+")
    ap.add_argument("--ref", nargs="+", default=None,
                    help="reference phases to try (default: all in "
                         "data/tcp_reference)")
    ap.add_argument("--refs-glob", default=REF_GLOB)
    ap.add_argument("--min-size", type=int, default=10,
                    help="drop crystalline components smaller than this")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    ref_paths = best_refs(args.refs_glob)
    if args.ref:
        ref_paths = {n: p for n, p in ref_paths.items() if n in args.ref}
    if not ref_paths:
        sys.exit("no reference crystals found (run tcp_reference.py)")
    print("references:", ", ".join(f"{n}" for n in ref_paths))
    refs = {n: build_struct(np.asarray(ddg.Manifold.load(p, 3).facets()))
            for n, p in ref_paths.items()}
    # atoms per unit cell per phase -> translation-invariant registry site.
    ns_of = {}
    for n in ref_paths:
        ns = len(STRUCTURES[n][1]) if n in STRUCTURES else 0
        if ns and refs[n]["V"] % ns == 0:
            ns_of[n] = ns
        else:
            print(f"  warn: {n} vertex count {refs[n]['V']} not a multiple of "
                  f"ns={ns}; registry check falls back to exact ref vertex")
    idx = ref_index(refs)

    files = []
    for s in args.states:
        files += sorted(glob.glob(s)) if any(c in s for c in "*?[") else [s]
    results = {}
    for f in files:
        facets = np.asarray(ddg.Manifold.load(f, 3).facets())
        r = analyze(facets, refs, idx, args.min_size, ns_of)
        results[f] = r
        sizes = collections.Counter(r["grain_sizes"])
        pha = " ".join(f"{k}:{v}" for k, v in sorted(r["by_phase"].items()))
        print(f"\n{os.path.basename(f)}")
        print(f"  interior-crystalline: {r['n_interior']}/{r['V']} "
              f"({100 * r['n_interior'] / r['V']:.2f}%)   by phase: {pha or '-'}")
        print(f"  grains (>= {args.min_size}): {r['n_grains']}   "
              f"largest {r['grain_sizes'][0] if r['grain_sizes'] else 0}   "
              f"sizes[size x count]: "
              + " ".join(f"{k}x{v}" for k, v in sorted(sizes.items())))

    if args.json:
        import json
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
