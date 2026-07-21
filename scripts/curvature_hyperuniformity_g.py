#!/usr/bin/env python3
"""Curvature hyperuniformity vs coupling: does the deficit field uniformize
approaching (and beyond) the glass transition?

Extends scripts/hyperuniformity.py (same estimators, importable + parameterized)
to a full g-ladder at N=1e4 plus the beyond-wall ANNEALED TCP states. Fields:

  q_R(v) = sum_{e ni v} delta_e / 2   (curvature charge, delta = 2pi - theta*deg)
  q_V(v) = D_v / 4                    (volume charge)

Estimator 1 (window): Var over 300 random centers of the total charge in BFS
windows grown to exactly m vertices, divided by the same variance with q values
shuffled over vertices. Hyperuniform <=> ratio -> 0 as m grows; ratio ~ 1 =
Poisson-like; > 1 = clustered/over-disperse.

Estimator 2 (spectral): graph-Laplacian mode power S(lambda) = |phi_n . dq|^2
vs shuffle, lambda ~ k^2; hyperuniform <=> ratio -> 0 as lambda -> 0.

Physics context (memory: phason/TCP program): PN defect binding measured in the
annealed states is a screening mechanism; screened Coulomb systems are
hyperuniform. So the signed curvature field may uniformize with g (and with
annealing) even while the defect DENSITY clusters into domain walls.
"""
import glob
import json
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import graph_hyperuniformity as gh
from discrete_differential_geometry.vertex_fields import FIELDS, edges_and_degrees

NCENT, NSHUF, NSPEC = 300, 4, 41

# (tag, glob pattern, n_seeds)
ENSEMBLES = [
    ("eq g=1e-3",   "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_1e-3_s0*.mfd", 8),
    ("eq g=2e-3",   "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_s0*.mfd", 8),
    ("eq g=4e-3",   "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_4e-3_s0*.mfd", 8),
    ("eq g=8e-3",   "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_8e-3_s0*.mfd", 8),
    ("eq g=1e-2",   "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_1e-2_s0*.mfd", 8),
    ("eq combo 0p405", "seeds/S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3_HDVs_0p405_s0*.mfd", 8),
    ("annealed d300",  "data/fk_anneal/N1e4_x100/rate300_rep*.final.mfd", 2),
    ("annealed d3000", "data/fk_anneal/N1e4_x100/rate3000_rep*.final.mfd", 2),
    ("annealed d30000","data/fk_anneal/N1e4_x100/rate30000_rep*.final.mfd", 2),
]


def main():
    rng = np.random.default_rng(0)
    results = {}
    print(f"{'ensemble':>18} {'winR(mid)':>15} {'winR(max)':>15} {'specR(low)':>15}")
    for tag, pat, reps in ENSEMBLES:
        paths = sorted(glob.glob(os.path.join(_ROOT, pat)))
        paths = paths[:(16 if "seeds/" in pat else reps)]
        if not paths:
            print(f"{tag:>18}  -- missing {pat}", file=sys.stderr)
            continue
        t0 = time.time()
        ratios, specs, lam0 = [], [], None
        for path in paths:
            fac = np.asarray(Manifold.load(path, 3).facets(), np.int64)
            qR = FIELDS["curvature_charge"](fac)
            eu, ecnt, deg, V = edges_and_degrees(fac)
            adj = gh.adjacency(eu, V)
            mmax = V // 6
            mgrid = np.unique(np.geomspace(8, mmax, 12).astype(int))
            orders = [gh.bfs_order(adj, V, int(rng.integers(V)), mmax, rng)
                      for _ in range(NCENT)]
            ratios.append(gh.window_ratio(qR, orders, mgrid, NSHUF, rng))
            lam, sr = gh.spectral_ratio(eu, V, qR - qR.mean(),
                                        k=NSPEC, nshuf=32, rng=rng)
            specs.append(sr)
            lam0 = lam if lam0 is None else lam0
        M = np.array(ratios)
        S = np.array(specs)
        mid = int(np.searchsorted(mgrid, mgrid[-1] // 4))
        out = dict(mgrid=mgrid.tolist(),
                   ratio_mean=M.mean(0).tolist(),
                   ratio_sem=(M.std(0, ddof=1) / np.sqrt(len(M))).tolist()
                   if len(M) > 1 else (M.std(0) * 0).tolist(),
                   spec_lam=lam0.tolist(),
                   spec_ratio_mean=S.mean(0).tolist(),
                   spec_ratio_sem=(S.std(0, ddof=1) / np.sqrt(len(S))).tolist()
                   if len(S) > 1 else (S.std(0) * 0).tolist(),
                   n_seeds=len(paths))
        results[tag] = out
        lo = S[:, :5].mean(1)
        print(f"{tag:>18} {out['ratio_mean'][mid]:9.3f}±{out['ratio_sem'][mid]:5.3f}"
              f" {out['ratio_mean'][-1]:9.3f}±{out['ratio_sem'][-1]:5.3f}"
              f" {lo.mean():9.3f}±{lo.std(ddof=1)/np.sqrt(len(lo)) if len(lo)>1 else 0:5.3f}"
              f"   ({time.time()-t0:.0f}s)", flush=True)

    outp = os.path.join(_ROOT, "out", "curvature_hyperuniformity_g.json")
    with open(outp, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {outp}")


if __name__ == "__main__":
    main()
