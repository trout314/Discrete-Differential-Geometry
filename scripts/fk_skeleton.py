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
  * vertex Z-class census: Z12/Z14/Z15/Z16, other pure-{5,6}, rest -- with Z16
    split into the FK Friauf polyhedron (T_d) and the D2 isomer, which has the
    same n_6 = 4 but is not a Frank-Kasper coordination (tools/link_classes.py);
    fFK_strict excludes it, fFK (historical, n_6-bucketed) does not;
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
# Census core promoted to the package (2026-08 cleanup, Phase 1a); re-exported
# here so legacy ``from fk_skeleton import ...`` keeps working.
from discrete_differential_geometry.fk_skeleton import (   # noqa: F401
    EDGE_PAIRS, edges_from_facets, load_edges, vertex_classes,
    vertex_class_census, skeleton_stats, bfs_orders, window_ratio)

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

def analyze_seed(path, n_centers, rng, n_shuf=4):
    facets = np.asarray(Manifold.load(path, 3).facets(), np.int64)
    eu, edeg, V = edges_from_facets(facets)
    f1 = len(eu)
    row = dict(V=V, f1=f1,
               edv=float(np.var(edeg)),
               mean_deg=float(np.mean(edeg)),
               p_le4=float(np.mean(edeg <= 4)),
               p5=float(np.mean(edeg == 5)),
               p6=float(np.mean(edeg == 6)),
               p_ge7=float(np.mean(edeg >= 7)),
               deg_hist=np.bincount(edeg).tolist())
    fz, n_broken = vertex_class_census(eu, edeg, V, facets=facets)
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
                   "fZ16_Td", "fZ16_D2", "fFK_strict",
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
