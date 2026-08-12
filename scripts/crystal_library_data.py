#!/usr/bin/env python3
"""Assemble everything known about the TCP crystal library into one JSON.

Feeds scripts/crystal_library_page.py (the HTML catalog). Everything here is
either read from the validated reference builds or computed exactly -- nothing
is transcribed by hand, so the page cannot drift from the code.

Per crystal:
  * cell metadata: prototype, atoms/cell, lattice matrix + derived
    a/b/c/alpha/beta/gamma, CN_mean, e_native = 6 - 12/CN, Z census
  * flat-pin arithmetic: e* - e_native, the extensive pin gap and the forced
    defect debt it implies
  * exact Aut(K): order, point+centering order |Aut|/m^3, orbit counts for
    k = 0..3, and the decoration refinement (Z class / hinge degree / face
    edge-degree triple / tet six-edge subgraph shape)
  * the defect-gas verdict at c_imp*
  * two 3D scenes -- one unit cell with its coordination shells and frame,
    one bulk block -- as flat typed arrays, with every edge tagged by hinge
    degree so the viewer can pick out the degree-6 (disclination) network

SCENE GEOMETRY. Degrees are computed on the T^3 quotient at the display
multiplier, which is exact: the crystal is periodic, so any m >= 2 gives the
true bulk hinge degree. For the chunk scene, edges that wrap the torus are
dropped (an edge wraps iff its fractional displacement exceeds half the box),
leaving a genuine finite piece whose interior degrees are exactly right. The
cell scene shows one unit cell's atoms with their coordination shells drawn
whole, each bond taken from an exact (site pair, lattice offset) table rather
than by min-imaging, plus the cell frame.

Usage:
    python scripts/crystal_library_data.py --out data/crystal_gas/library.json
"""
import argparse
import collections
import itertools
import json
import math
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, os.path.join(_ROOT, "scripts", "defect_dynamics"))

from discrete_differential_geometry import Manifold
from discrete_differential_geometry.symmetry import CrystalSymmetry
from tcp_reference import STRUCTURES, build_t3_triangulation
from cocycle_check import reference_frac_positions
from crystal_gas_scan import E_FLAT, LIBRARY, pin_gap
from orbit_decoration_census import E6, tet_key, tet_shape

# prototype / common name, and the display block size K (K^3 unit cells),
# chosen to land each chunk scene in ~300-1300 atoms.
META = {
    "a15":   ("Cr3Si (A15, beta-W)",          4),
    "c15":   ("MgCu2 (C15 Laves)",            3),
    "c14":   ("MgZn2 (C14 Laves)",            4),
    "c36":   ("MgNi2 (C36 Laves)",            3),
    "sigma": ("CrFe (sigma)",                 3),
    "z":     ("Zr4Al3 (Z)",                   5),
    "mu":    ("W6Fe7 (mu)",                   3),
    "r":     ("Co-Cr-Mo (R)",                 2),
    "p":     ("Cr-Ni-Mo (P)",                 2),
    "delta": ("MoNi (delta)",                 2),
}


def cell_geometry(L):
    """(a, b, c, alpha, beta, gamma) in lattice units / degrees."""
    a, b, c = (float(np.linalg.norm(L[i])) for i in range(3))
    ang = []
    for i, j in ((1, 2), (0, 2), (0, 1)):
        cos = float(L[i] @ L[j] / (np.linalg.norm(L[i]) * np.linalg.norm(L[j])))
        ang.append(math.degrees(math.acos(max(-1.0, min(1.0, cos)))))
    return a, b, c, ang[0], ang[1], ang[2]


def build_scene(name, K):
    """Positions (Cartesian), edges with hinge degree, Z per vertex, for a
    K^3 block of `name`."""
    L = STRUCTURES[name][0]
    ns = len(STRUCTURES[name][1])
    fac, nv = build_t3_triangulation(name, K)
    edeg = collections.Counter()
    for tv in fac:
        tv = tuple(int(x) for x in tv)
        for i, j in E6:
            a, b = tv[i], tv[j]
            edeg[(a, b) if a < b else (b, a)] += 1
    adj = collections.defaultdict(set)
    for u, w in edeg:
        adj[u].add(w)
        adj[w].add(u)
    frac = reference_frac_positions(name, K)          # period K, lattice units
    Z = np.array([len(adj[v]) for v in range(nv)], np.int16)
    return L, ns, frac, edeg, adj, Z, nv


def scene_chunk(L, frac, edeg, Z, K):
    """The full K^3 block; edges that wrap the torus are dropped."""
    xyz = frac @ L
    keep_e, deg = [], []
    for (u, w), d in edeg.items():
        dd = frac[w] - frac[u]
        if np.any(np.abs(dd) > K / 2.0):              # wraps the box
            continue
        keep_e.append((u, w))
        deg.append(d)
    return xyz, np.array(keep_e, np.int32), np.array(deg, np.int8), Z


def perturbed_sites(name):
    """The unit-cell site coordinates EXACTLY as tcp_reference perturbs them,
    but WITHOUT the final modulo.

    reference_frac_positions returns (sites + cell) % m. The deterministic
    ~1e-6 joggle pushes a coordinate sitting at 0 very slightly negative, and
    the modulo then wraps it to ~m -- harmless on the torus, fatal here: it
    makes a corner atom look like it lives m cells away, and every lattice
    offset computed against it comes out wrong by m. Rebuild the sites the
    way tcp_reference does and stop before the wrap."""
    sites = STRUCTURES[name][1]
    rng = np.random.default_rng(12345)               # same seed, same draw
    return sites + 1e-6 * rng.standard_normal(sites.shape)


def bond_table(ns, sites, frac, edeg, K):
    """{(site_u, site_w, lattice offset): hinge degree}.

    The crystal is periodic with period ONE unit cell, so every bond in the
    whole structure is one of these (site pair, offset) classes. Reading them
    off the K-torus once makes the cell drawing exact: an image pair is bonded
    iff its class is in this table, with no distance cutoff to tune."""
    tab = {}
    for (u, w), d in edeg.items():
        su, sw = u % ns, w % ns
        dd = frac[w] - frac[u]
        dd -= K * np.round(dd / K)                    # min image, period K
        off = np.round(dd - (sites[sw] - sites[su])).astype(int)
        tab[(su, sw, tuple(off))] = d
        tab[(sw, su, tuple(-off))] = d
    return tab


def scene_cell(L, ns, sites, frac, edeg, Z, K):
    """One unit cell's atoms with every coordination shell drawn whole.

    A cell-contents-only drawing is useless for a TCP phase: at CN ~ 13.4
    almost every bond leaves the cell, so the picture is a handful of stubs.
    What reads is the cell's atoms plus the complete shell around each."""
    tab = bond_table(ns, sites, frac, edeg, K)
    nodes = [(s, (0, 0, 0)) for s in range(ns)]       # the cell's own atoms
    idx = {k: i for i, k in enumerate(nodes)}

    def node(s, off):
        k = (s, tuple(int(x) for x in off))
        if k not in idx:
            idx[k] = len(nodes)
            nodes.append(k)
        return idx[k]

    # Every bond leaving a cell atom, drawn to the exact image it reaches.
    # Reading the image straight off the bond table means no min-imaging and
    # so no stray far-flung atoms -- the failure mode of placing each
    # neighbour once, relative to whichever core atom happened to reach it
    # first, which both dropped bonds and blew up the bounding box.
    edges, deg = [], []
    for (sa, sb, off), d in tab.items():
        if sa >= ns:
            continue
        edges.append((idx[(sa, (0, 0, 0))], node(sb, off)))
        deg.append(d)
    pos = np.array([sites[s] + np.array(off, float) for s, off in nodes])
    xyz = pos @ L
    zz = np.array([Z[s] for s, _ in nodes], np.int16)
    # cell frame: 12 edges of the parallelepiped, carried separately so the
    # viewer can draw it as a frame rather than as a bond
    corners = np.array(list(itertools.product([0, 1], repeat=3)), float) @ L
    box = [(i, j) for i in range(8) for j in range(i + 1, 8)
           if bin(i ^ j).count("1") == 1]
    return (xyz, np.array(edges, np.int32).reshape(-1, 2),
            np.array(deg, np.int8), zz, corners, np.array(box, np.int32), ns)


def pack(xyz, edges, deg, Z, box_xyz=None, box_e=None, ncore=0):
    """Compact JSON payload: positions as milli-units, edges flat.

    Atoms and any frame share ONE normalisation (centroid and scale taken over
    both), else the frame floats free of the structure it is supposed to
    enclose."""
    allp = xyz if box_xyz is None else np.vstack([xyz, box_xyz])
    ctr = allp.mean(0)
    sc = float(np.abs(allp - ctr).max()) or 1.0
    q = np.round((xyz - ctr) / sc * 1000).astype(int)
    out = {"p": q.ravel().tolist(), "n": int(len(xyz)),
           "e": edges.ravel().tolist(), "d": deg.tolist(),
           "z": Z.tolist(), "scale": round(sc, 4)}
    if box_xyz is not None:
        out["bp"] = np.round((box_xyz - ctr) / sc * 1000).astype(int).ravel().tolist()
        out["be"] = box_e.ravel().tolist()
    if ncore:
        out["ncore"] = int(ncore)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(
        _ROOT, "data", "crystal_gas", "library.json"))
    ap.add_argument("--verdicts", default=os.path.join(
        _ROOT, "data", "crystal_gas", "verdicts.json"))
    args = ap.parse_args()

    verd = {}
    if os.path.exists(args.verdicts):
        for r in json.load(open(args.verdicts)):
            verd.setdefault(r["crystal"], []).append(r)

    lib = {}
    for name in META:
        proto, K = META[name]
        fn, mcell = LIBRARY[name]
        path = os.path.join(_ROOT, "data", "tcp_reference", fn)
        L, sites, cn, census = STRUCTURES[name]
        ns = len(sites)
        e_nat, gap, nforced = pin_gap(name, mcell)
        a, b, c, al, be, ga = cell_geometry(np.asarray(L, float))

        # ---- exact symmetry + decorated orbit refinement
        sym = CrystalSymmetry.for_manifold_path(path)
        mfd = Manifold.load(path, 3)
        tets = np.asarray(mfd.facets())
        ed = collections.Counter()
        for tv in tets:
            tv = tuple(int(x) for x in tv)
            for i, j in E6:
                x, y = tv[i], tv[j]
                ed[(x, y) if x < y else (y, x)] += 1
        nbr = collections.defaultdict(set)
        for u, w in ed:
            nbr[u].add(w)
            nbr[w].add(u)

        def edge_of(x, y):
            return ed[(x, y) if x < y else (y, x)]

        orb = {}
        for kind in ("vertex", "edge", "face", "tet"):
            reps = sym.orbit_representatives(kind)
            sizes = sym.orbit_sizes(kind)
            tal = collections.Counter()
            for r, sz in zip(reps, sizes):
                if kind == "vertex":
                    key = f"Z{len(nbr[int(r)])}"
                elif kind == "edge":
                    key = f"deg{edge_of(int(r[0]), int(r[1]))}"
                elif kind == "face":
                    vs = sorted(int(x) for x in r)
                    key = "".join(str(edge_of(*p))
                                  for p in sorted(itertools.combinations(vs, 2)))
                    key = "".join(sorted(key))
                else:
                    vs = sorted(int(x) for x in r)
                    dg = [edge_of(vs[i], vs[j]) for i, j in E6]
                    key = f"{''.join(map(str, tet_key(dg)))} ({tet_shape(dg)})"
                tal[key] += 1
            orb[kind] = {"total": sym.n_orbits(kind), "by": dict(tal)}

        # ---- 3D scenes
        Lm, _, frac, edeg, adjm, Z, nv = build_scene(name, K)
        cx, ce, cd, cz, bxyz, bedg, ncore = scene_cell(
            Lm, ns, perturbed_sites(name), frac, edeg, Z, K)
        kx, ke, kd, kz = scene_chunk(Lm, frac, edeg, Z, K)

        gas = None
        rows = sorted(verd.get(name, []), key=lambda r: r["cimp"])
        best = next((r for r in rows if r["verdict"] == "GAS"), None)
        if best:
            gas = {k: best[k] for k in
                   ("cimp", "n_ill", "dfrac", "ncomp", "top1", "max_m",
                    "turnover", "drift_z", "e_mean", "paid", "pin_note")}
        lib[name] = {
            "proto": proto, "ns": ns, "cn": cn, "e_native": e_nat,
            "census": census, "cell": {"a": a, "b": b, "c": c,
                                       "alpha": al, "beta": be, "gamma": ga},
            "mcell": mcell, "file": fn,
            "fv": [int(x) for x in mfd.f_vector],
            "gap": gap, "nforced": nforced, "de": E_FLAT - e_nat,
            "aut": sym.order, "point": sym.order // mcell ** 3,
            "orb": orb, "gas": gas,
            "K": K,
            "scene_cell": pack(cx, ce, cd, cz, bxyz, bedg, ncore),
            "scene_chunk": pack(kx, ke, kd, kz),
            "series": [{k: r[k] for k in
                        ("cimp", "n_ill", "dfrac", "ncomp", "top1", "max_m",
                         "turnover", "drift_z", "e_mean", "paid", "verdict")}
                       for r in rows],
        }
        print(f"{name:6s} K={K} cell={ncore}+{len(cx)-ncore} atoms/"
              f"{len(ce)} bonds  chunk={len(kx)} atoms/{len(ke)} bonds  "
              f"|Aut|={sym.order}", flush=True)

    payload = {"e_flat": E_FLAT, "crystals": lib}
    with open(args.out, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"\n-> {args.out}  ({os.path.getsize(args.out)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
