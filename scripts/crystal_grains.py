#!/usr/bin/env python3
"""Crystalline-grain detection by covering-map development.

Identifies the connected components of the vertex graph that are *consistent
with a piece of a reference TCP crystal* -- i.e. genuine crystalline grains,
not merely regions of locally-crystal-like coordination.

Why not per-vertex local matching (the superseded radius-r ball approach):
  A single FK coordination polyhedron (Z12/Z14/Z15/Z16) is the universal local
  motif of *all* tetrahedrally-close-packed matter, crystalline AND amorphous,
  so any radius-r ball test is either too weak (an isolated liquid icosahedron
  matches -- a radius-1 environment test) or too brittle (one perturbed edge in
  a ~60-vertex radius-2 ball breaks the match, so a 62%-ordered melt reads 0%
  crystal -- an exact radius-2 certificate). Crystallinity is *translational
  registry* -- vertices sitting in the crystal's actual periodic arrangement
  over an extended region -- which is inherently non-local. This covering-map
  detector is the single authoritative crystalline/defect identifier; defect
  census (defect_census.py) is a reporting front-end over it.

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

  A vertex is INTERIOR-crystalline iff (GATE) it sits in an intact Frank-Kasper
  coordination shell -- every spoke edge has hinge degree 5 or 6 with exactly the
  twelve mandatory fivefold disclinations -- AND (REGISTRY) the star tets the
  development ASSIGNED all share one grain and one translation-invariant unit-cell
  site. Star tets left unassigned because they straddle a defective neighbour are
  don't-cares, so a good vertex keeps its crystallinity on its crystal-facing side
  and the defect stays localized to the vertices a Pachner move actually touched
  (the shell gate breaks for exactly those) instead of spreading by a ring. Grains
  = connected components of interior vertices sharing a grain id. Everything else
  is defect/boundary. `--min-size` drops coincidental sub-threshold grains
  (calibrate against a fully-melted null, which yields ~0 grains).

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
import re
import sys
from itertools import permutations

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes
from tcp_reference import STRUCTURES

# tet vertex-index pairs for the six edges
_E6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

REF_GLOB = os.path.join(_ROOT, "data", "tcp_reference", "T3_*_m*.mfd")


def best_refs(ref_glob):
    """Largest-m reference .mfd per structure name -> {name: path}."""
    best = {}
    for p in glob.glob(ref_glob):
        mm = re.match(r"T3_([A-Z0-9]+)_m(\d+)_N\d+\.mfd", os.path.basename(p))
        if not mm:
            continue
        name, m = mm.group(1).lower(), int(mm.group(2))
        if name not in best or m > best[name][0]:
            best[name] = (m, p)
    return {n: p for n, (m, p) in sorted(best.items())}


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


def _drop_subsumed_grains(st, grain_of_tet, phase_of_grain):
    """Un-label redundant re-seed grains: a grain whose vertex set is contained in
    that of a single larger grain is the SAME crystal region developed a second
    time (the covering map re-seeds wherever a defect locally blocks the flood,
    picking an arbitrary symmetry-related frame). We drop its grain labels so its
    tets become unassigned (don't-cares); its vertices then keep only the frame of
    the grain that subsumes them and no longer straddle two frames of one crystal.
    This only edits the grain->tet bookkeeping, never the manifold itself.

    Only subset-of-a-SINGLE-larger-grain is dropped -- a genuine second domain
    (real grain boundary) has interior vertices of its own that no other grain
    covers, and holonomy-fragmented cross-phase grains (many similar-sized,
    mutually overlapping but non-nesting patches) are left intact, so this does
    not weaken phase discrimination."""
    verts = collections.defaultdict(set)
    tets = collections.defaultdict(list)
    for t, g in grain_of_tet.items():
        tets[g].append(t)
        for v in st["tets"][t]:
            verts[g].add(int(v))
    order = sorted(verts, key=lambda g: len(verts[g]))     # smallest first
    for i, g in enumerate(order):
        for h in order[i + 1:]:                            # only larger-or-equal
            if h in tets and g != h and verts[g] <= verts[h]:
                for t in tets[g]:                          # subsumed -> unassign
                    del grain_of_tet[t]
                del tets[g]
                break


def find_grains(st, refs, idx):
    """Label sample tets by grain. Returns (grain_of_tet, sig_of_tet,
    phase_of_grain). Seeds interiors first (pure tets), so one seed grabs a
    whole grain; boundary/defect tets are reached and validated during a
    neighbour's development or left unassigned. Redundant re-seed grains (subsumed
    by a larger grain) are un-labelled so their vertices carry a single frame."""
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
    _drop_subsumed_grains(st, grain_of_tet, phase_of_grain)
    return grain_of_tet, sig_of_tet, phase_of_grain


def _locally_fk(st, v):
    """Does v sit in a clean Frank-Kasper coordination shell? True iff every spoke
    edge (v,w) has hinge degree 5 or 6 and EXACTLY twelve of them are 5.

    The link of v is a triangulated 2-sphere, and the degree of a link node w
    equals hinge(v,w) (both count the tets around edge (v,w)); Euler forces any
    all-5/6 triangulated 2-sphere to carry exactly twelve fivefold vertices (the
    mandatory disclinations, Z-class = 12 + #sixes). A single Pachner move breaks
    this for *exactly* the vertices it touches -- a subdivided vertex gains a
    low-degree spoke, a flipped vertex trades a fivefold spoke for a sixfold so
    the count drops below 12 -- while every neighbour's shell is untouched. So it
    is a local, reference-free defect gate that does NOT spread by a ring. We
    deliberately test only spoke degrees here: neighbour FK-class and link-edge
    degrees both shift when a neighbour is defective and would leak the defect
    outward (that leak was the original one-ring over-flagging bug)."""
    deg = [_edge(st, v, w) for w in st["adj"][v]]
    return all(d in (5, 6) for d in deg) and deg.count(5) == 12


def interior_vertices(st, grain_of_tet, sig_of_tet, phase_of_grain, ns_of):
    """v -> grain id for interior-crystalline vertices. Two conditions, the
    strict-seed / loose-heal pair:

      GATE (strict, local): v has an intact FK shell (`_locally_fk`). This is what
        makes a genuine wrong-degree vertex a defect while sparing its good
        neighbours -- the defect stays localized to the vertices a move actually
        touched instead of spreading to their 1-ring.
      REGISTRY (loose heal): the star tets the covering map ASSIGNED all agree on
        one grain and one unit-cell SITE. Unassigned star tets -- those the strict
        development refused because they straddle a defective neighbour -- are
        DON'T-CARES: a good vertex on the crystal-facing side of a defect keeps
        its crystallinity from the tets that do develop. v needs at least one
        assigned tet (no registry at all => not interior).

    A single grain id suffices here because `find_grains` already un-labelled the
    redundant re-seed grains that would otherwise split one crystal across two
    frames (see `_drop_subsumed_grains`); what remains multi-grain at a vertex is
    a genuine boundary. SITE = ref_vertex % ns (ns = atoms per cell) is
    translation-invariant, so wrapping the sample torus to reference images
    differing by a lattice translation (same site) is genuine crystal, not a
    conflict; only a true stacking/phase inconsistency changes the site."""
    tets_of = collections.defaultdict(list)
    for t, tv in enumerate(st["tets"]):
        for v in tv:
            tets_of[int(v)].append(t)
    interior = {}
    for v, tl in tets_of.items():
        if not _locally_fk(st, v):
            continue                                       # broken shell -> defect
        assigned = [t for t in tl if t in grain_of_tet]
        if not assigned:
            continue                                       # no registry -> defect
        gids = {grain_of_tet[t] for t in assigned}
        if len(gids) != 1:
            continue                                       # straddles grains
        gid = gids.pop()
        ns = ns_of.get(phase_of_grain[gid], 0)
        sites = {sig_of_tet[t][v] % ns if ns else sig_of_tet[t][v]
                 for t in assigned}
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
