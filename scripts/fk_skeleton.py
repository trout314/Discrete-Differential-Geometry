#!/usr/bin/env python3
"""Frank-Kasper / disclination-skeleton census of seed triangulations.

Physics context (see notes + memory: QC/phason/TCP program): the flat edge
degree 2*pi/arccos(1/3) ~ 5.1043 is irrational, so a near-flat triangulation at
minimal degree variance must be a degree-5 majority threaded by a network of
degree>=6 edges -- combinatorially a tetrahedrally-close-packed (Frank-Kasper)
structure, whose ordered representatives include FK crystals and quasicrystals.
This script measures, on saved .mfd seeds, how much FK/TCP local order the
equilibrium ensemble develops as the VDV coupling g approaches the glass wall.

Exact combinatorial backbone: for every vertex v the link is a triangulated
S^2, and the degree of edge {v,u} equals the degree of u in link(v); Euler for
the link gives the LOCAL SUM RULE  sum_{e ni v} (6 - deg e) = 12  exactly.
Hence a vertex whose incident edges are all degree 5 has exactly 12 of them
(icosahedral link, Z12), and pure-{5,6} vertices with two/three/four 6-edges
are the FK coordinations Z14/Z15/Z16. A pure-{5,6} vertex with exactly ONE
6-edge cannot exist (no S^2 triangulation with degree census 5^12 6^1), so
"disclination lines cannot end" is a theorem in the pure-{5,6} sector; line
endpoints require a compensating deg<=4 or deg>=7 edge at the same vertex.

Measured per seed:
  * edge-degree census (fractions by degree, EDV);
  * vertex Z-class census: Z12/Z14/Z15/Z16, other pure-{5,6}, rest;
  * the 6-skeleton (edges deg>=6): valence histogram of its vertices
    (valence-1 = broken lines, forbidden in ideal TCP), connected components,
    largest-component share, cyclomatic number (loops);
  * composition hyperuniformity: window variance of the per-vertex 6-edge
    charge q(v) = (1/2) sum_{e ni v} 1[deg e >= 6] over fixed-size BFS windows
    vs a value-shuffled null (same methodology as scripts/hyperuniformity.py);
    ratio -> 0 with window size = hyperuniform composition, ~1 = random,
    > 1 = clustered/phase-separating.

Output: printed per-family table + out/fk_skeleton.json.

Usage:
    python scripts/fk_skeleton.py                 # default N=1e4 g-ladders
    python scripts/fk_skeleton.py --reps 4 --centers 100   # quicker pass
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold

# (tag, family stem, g) -- k=2 ladders at N=1e4 across the three edge pins.
# g is the VDV coupling beta/N; None = no VDV term (base family).
ENSEMBLES = [
    ("ED5p0043 g=0",    "S3_N1e4_1e-1_ED5p0043_2",            0.0),
    ("ED5p0043 g=1e-3", "S3_N1e4_1e-1_ED5p0043_2_VDVs_1e-3",  1e-3),
    ("ED5p0043 g=2e-3", "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3",  2e-3),
    ("ED5p0043 g=4e-3", "S3_N1e4_1e-1_ED5p0043_2_VDVs_4e-3",  4e-3),
    ("ED5p0043 g=8e-3", "S3_N1e4_1e-1_ED5p0043_2_VDVs_8e-3",  8e-3),
    ("ED5p0043 g=1e-2", "S3_N1e4_1e-1_ED5p0043_2_VDVs_1e-2",  1e-2),
    ("ED5p1043 g=0",    "S3_N1e4_1e-1_ED5p1043_2",            0.0),
    ("ED5p1043 g=1e-3", "S3_N1e4_1e-1_ED5p1043_2_VDVs_1e-3",  1e-3),
    ("ED5p1043 g=2e-3", "S3_N1e4_1e-1_ED5p1043_2_VDVs_2e-3",  2e-3),
    ("ED5p1043 g=4e-3", "S3_N1e4_1e-1_ED5p1043_2_VDVs_4e-3",  4e-3),
    ("ED5p2043 g=0",    "S3_N1e4_1e-1_ED5p2043_2",            0.0),
    ("ED5p2043 g=1e-3", "S3_N1e4_1e-1_ED5p2043_2_VDVs_1e-3",  1e-3),
    ("ED5p2043 g=2e-3", "S3_N1e4_1e-1_ED5p2043_2_VDVs_2e-3",  2e-3),
    ("ED5p2043 g=4e-3", "S3_N1e4_1e-1_ED5p2043_2_VDVs_4e-3",  4e-3),
    # Hinge-degree-variance penalties: the coupling that directly targets the
    # TCP/FK composition (edge degrees), alone and on top of VDVs_2e-3.
    ("HDV 5e-2",        "S3_N1e4_1e-1_ED5p0043_2_HDVs_5e-2",  5e-2),
    ("HDV 1e-1",        "S3_N1e4_1e-1_ED5p0043_2_HDVs_1e-1",  1e-1),
    ("HDV 4e-1",        "S3_N1e4_1e-1_ED5p0043_2_HDVs_4e-1",  4e-1),
    ("VDV2e-3+HDV.032", "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p032", 0.032),
    ("VDV2e-3+HDV.101", "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p101", 0.101),
    ("VDV2e-3+HDV.228", "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p228", 0.228),
    ("VDV2e-3+HDV.405", "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p405", 0.405),
]

EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def load_edges(path):
    """Load a seed; return (eu, edeg, V): unique edges (relabo 0..V-1), degrees."""
    m = Manifold.load(path, 3)
    F = np.asarray(m.facets(), np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    T = inv.reshape(F.shape)
    V = len(lab)
    epairs = np.sort(np.vstack([T[:, [i, j]] for i, j in EDGE_PAIRS]), axis=1)
    eu, edeg = np.unique(epairs, axis=0, return_counts=True)
    return eu, edeg, V


def vertex_class_census(eu, edeg, V):
    """Per-vertex counts of incident edges by degree class + Z-classification."""
    def incident(mask):
        c = np.zeros(V, dtype=np.int64)
        np.add.at(c, eu[mask, 0], 1)
        np.add.at(c, eu[mask, 1], 1)
        return c

    n_le4 = incident(edeg <= 4)
    n5 = incident(edeg == 5)
    n6 = incident(edeg == 6)
    n_ge7 = incident(edeg >= 7)

    # Exact local sum rule (vertex link = S^2): sum(6 - deg) over incident = 12.
    charge = np.zeros(V, dtype=np.int64)
    np.add.at(charge, eu[:, 0], 6 - edeg)
    np.add.at(charge, eu[:, 1], 6 - edeg)
    assert np.all(charge == 12), "link sum rule violated -- not a 3-manifold?"

    pure56 = (n_le4 == 0) & (n_ge7 == 0)
    fz = {}
    for name, k in (("Z12", 0), ("Z14", 2), ("Z15", 3), ("Z16", 4)):
        fz[name] = float(np.mean(pure56 & (n6 == k)))
    n_broken = int(np.sum(pure56 & (n6 == 1)))  # forbidden by S^2 combinatorics
    fz["pure56_other"] = float(np.mean(pure56)) - sum(fz.values())
    fz["impure"] = float(np.mean(~pure56))
    return fz, n_broken


def skeleton_stats(eu, edeg, V):
    """Structure of the deg>=6 disclination network."""
    sk = eu[edeg >= 6]
    if len(sk) == 0:
        return dict(n_edges=0, frac_val1=np.nan, frac_val2=np.nan,
                    n_comp=0, largest_frac=np.nan, cyclomatic=0)
    val = np.zeros(V, dtype=np.int64)
    np.add.at(val, sk[:, 0], 1)
    np.add.at(val, sk[:, 1], 1)
    touched = val > 0
    vh = np.bincount(val[touched])
    ntouch = int(touched.sum())

    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in sk:
        a, b = int(a), int(b)
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = {}
    for a, b in sk:
        r = find(int(a))
        roots[r] = roots.get(r, 0) + 1
    comp_sizes = sorted(roots.values(), reverse=True)
    return dict(
        n_edges=int(len(sk)),
        frac_val1=float(vh[1] / ntouch) if len(vh) > 1 else 0.0,
        frac_val2=float(vh[2] / ntouch) if len(vh) > 2 else 0.0,
        n_comp=len(comp_sizes),
        largest_frac=float(comp_sizes[0] / len(sk)),
        cyclomatic=int(len(sk) - ntouch + len(comp_sizes)),
    )


def bfs_orders(eu, V, n_centers, rng):
    adj = [[] for _ in range(V)]
    for a, b in eu:
        adj[a].append(b)
        adj[b].append(a)
    mmax = V // 6
    orders = []
    for _ in range(n_centers):
        src = int(rng.integers(V))
        seen = np.zeros(V, bool)
        seen[src] = True
        order = [src]
        frontier = [src]
        while len(order) < mmax and frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if not seen[w]:
                        seen[w] = True
                        nxt.append(w)
            rng.shuffle(nxt)
            order.extend(nxt)
            frontier = nxt
        orders.append(np.array(order[:mmax]))
    return orders, mmax


def window_ratio(q, orders, mgrid, rng, n_shuf=4):
    """Var over centers of window charge, observed / value-shuffled null."""
    def var(field):
        sums = np.array([np.cumsum(field[o])[mgrid - 1] for o in orders])
        return sums.var(axis=0)

    vr = var(q)
    vs = np.mean([var(rng.permutation(q)) for _ in range(n_shuf)], axis=0)
    return vr / vs


def analyze_seed(path, n_centers, rng, n_shuf=4):
    eu, edeg, V = load_edges(path)
    f1 = len(eu)
    row = dict(V=V, f1=f1,
               edv=float(np.var(edeg)),
               mean_deg=float(np.mean(edeg)),
               p_le4=float(np.mean(edeg <= 4)),
               p5=float(np.mean(edeg == 5)),
               p6=float(np.mean(edeg == 6)),
               p_ge7=float(np.mean(edeg >= 7)),
               deg_hist=np.bincount(edeg).tolist())
    fz, n_broken = vertex_class_census(eu, edeg, V)
    row.update({f"f{k}": v for k, v in fz.items()})
    row["fFK"] = fz["Z12"] + fz["Z14"] + fz["Z15"] + fz["Z16"]
    row["n_broken56"] = n_broken
    row.update(skeleton_stats(eu, edeg, V))
    # Null for the skeleton geometry: permute degree labels over the SAME
    # 1-skeleton (keeps density and graph, kills degree-degree correlations).
    nulls = [skeleton_stats(eu, rng.permutation(edeg), V) for _ in range(n_shuf)]
    for k in ("frac_val1", "largest_frac", "cyclomatic"):
        row[k + "_null"] = float(np.mean([s[k] for s in nulls]))
    row["frac_val1_ratio"] = (row["frac_val1"] / row["frac_val1_null"]
                              if row["frac_val1_null"] else np.nan)

    q6 = np.zeros(V)
    sk6 = edeg >= 6
    np.add.at(q6, eu[sk6, 0], 0.5)
    np.add.at(q6, eu[sk6, 1], 0.5)
    orders, mmax = bfs_orders(eu, V, n_centers, rng)
    mgrid = np.unique(np.geomspace(8, mmax, 10).astype(int))
    row["mgrid"] = mgrid.tolist()
    row["win_ratio"] = window_ratio(q6, orders, mgrid, rng).tolist()
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=8, help="replicas per family")
    ap.add_argument("--centers", type=int, default=200, help="BFS window centers")
    ap.add_argument("--out", default=os.path.join(_ROOT, "out", "fk_skeleton.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    results = {}
    scalar_keys = ["edv", "p_le4", "p5", "p6", "p_ge7",
                   "fZ12", "fZ14", "fZ15", "fZ16", "fpure56_other", "fimpure",
                   "fFK", "frac_val1", "frac_val2", "largest_frac",
                   "n_comp", "n_edges", "cyclomatic", "n_broken56",
                   "frac_val1_null", "largest_frac_null", "cyclomatic_null",
                   "frac_val1_ratio"]
    print(f"{'family':>18} {'EDV':>6} {'p_le4':>6} {'p5':>6} {'p6':>6} "
          f"{'p_ge7':>6} {'pure56':>6} {'val1':>6} {'val1/nul':>8} "
          f"{'ratio_last':>10}")
    for tag, stem, g in ENSEMBLES:
        paths = sorted(glob.glob(os.path.join(_ROOT, "seeds", stem + "_s0*.mfd")))
        paths = paths[:args.reps]
        if not paths:
            print(f"{tag:>18}  -- no seeds found ({stem})", file=sys.stderr)
            continue
        t0 = time.time()
        rows = [analyze_seed(p, args.centers, rng) for p in paths]
        agg = {"g": g, "stem": stem, "n_seeds": len(rows),
               "mgrid": rows[0]["mgrid"], "seeds": rows}
        for k in scalar_keys:
            vals = np.array([r[k] for r in rows], dtype=float)
            agg[k] = float(np.nanmean(vals))
            agg[k + "_sem"] = (float(np.nanstd(vals, ddof=1) / np.sqrt(len(vals)))
                               if len(vals) > 1 else 0.0)
        W = np.array([r["win_ratio"] for r in rows])
        agg["win_ratio_mean"] = W.mean(0).tolist()
        agg["win_ratio_sem"] = (W.std(0, ddof=1) / np.sqrt(len(W))).tolist() \
            if len(W) > 1 else (W.std(0) * 0).tolist()
        results[tag] = agg
        pure56 = 1.0 - agg["fimpure"]
        print(f"{tag:>18} {agg['edv']:6.3f} {agg['p_le4']:6.3f} {agg['p5']:6.3f} "
              f"{agg['p6']:6.3f} {agg['p_ge7']:6.3f} {pure56:6.3f} "
              f"{agg['frac_val1']:6.3f} {agg['frac_val1_ratio']:8.3f} "
              f"{agg['win_ratio_mean'][-1]:10.3f}   ({time.time()-t0:.0f}s)",
              flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
