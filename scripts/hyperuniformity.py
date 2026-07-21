#!/usr/bin/env python3
"""Graph-proxy hyperuniformity test: is a vertex field more uniform than random?

Thin CLI over the package estimators (discrete_differential_geometry.
graph_hyperuniformity) and the shared field library (vertex_fields) -- no
coordinates, topology only. Two estimators, each vs a charge-permutation shuffle
null (same graph/windows, field values permuted over vertices -- isolates spatial
organization):

1. Window variance: BFS balls grown to EXACTLY m vertices; Var over random
   centers of the window charge Q = sum q(v), vs m. Hyperuniform: ratio decays
   with m; thermal/Poisson: ~ 1; clustered: > 1.
2. Spectral (low-k): matrix-free low-pass power ratio (`lowpass_ratio`, robust on
   crystals); --full-spectrum for the eigensolver spectrum instead.

Default field is curvature_charge (the Regge scalar-curvature density = the
Hamiltonian-constraint quantity). Pools inputs by --group regex.

Usage:
    python scripts/hyperuniformity.py 'seeds/S3_N1e4_*_s0*.mfd' \
        --field curvature_charge --out out/hu.json
"""
import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold
from discrete_differential_geometry.vertex_fields import FIELDS, edges_and_degrees
from discrete_differential_geometry import graph_hyperuniformity as gh


def hu_one(path, field, ncent, nshuf, nspec, rng, full_spectrum):
    fac = np.asarray(Manifold.load(path, 3).facets())
    q = FIELDS[field](fac)
    eu, ecnt, deg, V = edges_and_degrees(fac)
    adj = gh.adjacency(eu, V)
    mmax = V // 6
    mgrid = np.unique(np.geomspace(8, mmax, 12).astype(int))
    orders = [gh.bfs_order(adj, V, int(rng.integers(V)), mmax, rng)
              for _ in range(ncent)]
    wr = gh.window_ratio(q, orders, mgrid, nshuf, rng)
    dq = q - q.mean()
    out = {"mgrid": mgrid, "window_ratio": wr,
           "lowk_spectral": gh.lowpass_ratio(eu, V, dq, rng=rng)}
    if full_spectrum:
        out["spec_lambda"], out["spec_ratio"] = gh.spectral_ratio(
            eu, V, dq, k=nspec, rng=rng)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+", help=".mfd files (globs ok)")
    ap.add_argument("--field", choices=list(FIELDS), default="curvature_charge")
    ap.add_argument("--group", default=r"_s\d+.*$",
                    help="regex removed from basename to pool replicas")
    ap.add_argument("--ncent", type=int, default=300)
    ap.add_argument("--nshuf", type=int, default=4)
    ap.add_argument("--nspec", type=int, default=41)
    ap.add_argument("--full-spectrum", action="store_true",
                    help="also compute the eigensolver spectrum (slow; thrashes "
                         "on crystals -- lowpass low-k ratio is the robust default)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    files = []
    for s in args.snapshots:
        files += sorted(glob.glob(s)) if any(c in s for c in "*?[") else [s]
    rng = np.random.default_rng(0)
    pooled = defaultdict(list)
    for path in files:
        label = re.sub(args.group, "", os.path.splitext(os.path.basename(path))[0])
        r = hu_one(path, args.field, args.ncent, args.nshuf, args.nspec, rng,
                   args.full_spectrum)
        pooled[label].append(r)
        print(f"  {os.path.basename(path)}: window m={r['mgrid'][0]}:"
              f"{r['window_ratio'][0]:.3f} .. m={r['mgrid'][-1]}:"
              f"{r['window_ratio'][-1]:.3f}  low-k spectral={r['lowk_spectral']:.3f}",
              flush=True)

    results = {}
    for label, runs in sorted(pooled.items()):
        mgrid = runs[0]["mgrid"]
        wr = np.mean([r["window_ratio"] for r in runs], axis=0)
        lowk = float(np.mean([r["lowk_spectral"] for r in runs]))
        results[label] = {"mgrid": mgrid.tolist(), "window_ratio": wr.tolist(),
                          "lowk_spectral": lowk, "n": len(runs)}
        if "spec_ratio" in runs[0]:
            results[label]["spec_lambda"] = runs[0]["spec_lambda"].tolist()
            results[label]["spec_ratio"] = runs[0]["spec_ratio"].tolist()
        print(f"[{label:>16}] window ratio m={mgrid[0]}:{wr[0]:.3f} .. "
              f"m={mgrid[-1]}:{wr[-1]:.3f}   low-k spectral={lowk:.3f}   "
              f"({len(runs)} states)")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        json.dump(results, open(args.out, "w"), indent=1)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
