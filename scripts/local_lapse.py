#!/usr/bin/env python3
"""Measured combinatorial lapse N(v) from exact per-move counters (Stage 2).

Uses the sampler's per-vertex move-attribution counters (Option A: each event
distributes total weight 1 uniformly over its support — the 5 bistellar-ball
vertices = one glued 4-simplex, or the 6 support vertices of a 4-4 hinge move;
see sampler.MoveCounters and notes/combinatorial-lapse-plan.md). Supersedes the
per-sweep facet-diff proxy in local_lapse_proto.py, which saturates at high
churn and cannot see rejected/invalid attempts.

Per stable vertex (present at both window endpoints), with D_v the vertex
degree (endpoint average) and D_v/4 its share of the 3-volume:

  N(v)   = deposited 4-volume rate / local 3-volume
         = (acc_bistellar + 2*acc_hinge) / window / (D_v/4)
           (a 4-4 hinge move = 2 stacked 4-simplices)
  K(v)   = valid/proposed   -- kinematic availability (a legal move exists)
  A(v)   = accepted/valid   -- energetic (Metropolis) acceptance
  R(v)   = proposed rate per unit volume -- proposal pressure

Reports conditional means vs (even) vertex degree, and the spatial
autocorrelation of the degree-residual dN on the 1-skeleton (dynamical
heterogeneity). Aggregates over replicas.

Examples
--------
    python scripts/local_lapse.py                      # three-phase N=1e3 default
    python scripts/local_lapse.py --seeds 'seeds/S3_N1e4_1e-1_ED5p0043_4_VDVs_8e-3_s00*.mfd' \
        --window 100 --out out/lapse_8e-3_N1e4
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from seed_utils import load_seed_metadata
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

DEFAULT_PHASES = [
    ("extended_2e-3", "seeds/S3_N1e3_1e-1_ED5p0043_2_VDVs_2e-3_s0*.mfd"),
    ("hub_beta0_flatpin", "seeds/S3_N1e3_1e-1_ED5p1043_2_s0*.mfd"),
    ("glassadj_8e-3_k4", "seeds/S3_N1e3_1e-1_ED5p0043_4_VDVs_8e-3_s0*.mfd"),
]
HINGE_VOL = 2.0    # 4-4 hinge move = 2 stacked 4-simplices (pentachoron convention)


def params_from_metadata(path):
    md = load_seed_metadata(path)
    return SamplerParams(
        num_facets_target=int(md["num_facets_target"]),
        hinge_degree_target=float(md["hinge_degree_target"]),
        num_facets_coef=float(md["num_facets_coef"]),
        num_hinges_coef=float(md["num_hinges_coef"]),
        hinge_degree_variance_coef=float(md["hinge_degree_variance_coef"]),
        codim3_degree_variance_coef=float(md["codim3_degree_variance_coef"]),
    )


def degrees(F):
    lab, cnt = np.unique(F, return_counts=True)
    return dict(zip(lab.tolist(), cnt.tolist()))


def measure_seed(path, burn, window):
    """Run one replica window; return per-stable-vertex fields + the end facets."""
    s = ManifoldSampler(Manifold.load(path, 3), params_from_metadata(path))
    s.run(sweeps=burn)
    deg0 = degrees(np.asarray(s.manifold.facets(), np.int64))
    s.track_move_counts()
    s.reset_move_counts()
    s.run(sweeps=window)
    mc = s.move_counts()
    F1 = np.asarray(s.manifold.facets(), np.int64)
    deg1 = degrees(F1)

    idx = {v: i for i, v in enumerate(mc["vertex"].tolist())}
    stable = sorted(set(deg0) & set(deg1) & set(idx))
    sel = np.array([idx[v] for v in stable])
    D = np.array([(deg0[v] + deg1[v]) / 2.0 for v in stable])
    prop = mc["proposed"][sel]
    valid = mc["valid"][sel]
    acc_b = mc["accepted_bistellar"][sel]
    acc_h = mc["accepted_hinge"][sel]

    vol4 = D / 4.0
    fields = dict(
        vertex=np.array(stable), deg=D,
        N=(acc_b + HINGE_VOL * acc_h) / window / vol4,
        K=np.divide(valid, prop, out=np.zeros_like(valid), where=prop > 0),
        A=np.divide(acc_b + acc_h, valid, out=np.zeros_like(valid), where=valid > 0),
        R=prop / window / vol4,
    )
    return fields, F1


def cond_by_degree(deg, x, min_count=3):
    out = {}
    for d, v in zip(deg, x):
        out.setdefault(int(2 * round(d / 2)), []).append(v)
    return {b: (float(np.mean(v)), float(np.std(v) / np.sqrt(len(v))), len(v))
            for b, v in sorted(out.items()) if len(v) >= min_count}


def spatial_corr(fields, F1, rng, rmax=6, nsrc=150):
    """Autocorrelation of the degree-residual N field along the 1-skeleton."""
    deg, N = fields["deg"], fields["N"]
    mu = {b: c[0] for b, c in cond_by_degree(deg, N, min_count=1).items()}
    dN = np.array([n - mu[int(2 * round(d / 2))] for n, d in zip(N, deg)])
    idx = {v: i for i, v in enumerate(fields["vertex"].tolist())}
    pr = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    ee = np.unique(np.sort(np.vstack([F1[:, [i, j]] for i, j in pr]), axis=1), axis=0)
    ee = np.array([[idx[a], idx[b]] for a, b in ee if a in idx and b in idx])
    n = len(dN)
    G = coo_matrix((np.ones(2 * len(ee)),
                    (np.r_[ee[:, 0], ee[:, 1]], np.r_[ee[:, 1], ee[:, 0]])),
                   shape=(n, n)).tocsr()
    src = rng.choice(n, min(nsrc, n), replace=False)
    dist = shortest_path(G, method="D", unweighted=True, directed=False, indices=src)
    var = dN.var()
    corr = {}
    for r in range(1, rmax + 1):
        pairs = [(i, j) for si, i in enumerate(src) for j in np.where(dist[si] == r)[0]]
        if len(pairs) < 30:
            break
        a = np.array([dN[i] for i, _ in pairs])
        b = np.array([dN[j] for _, j in pairs])
        corr[r] = float(np.mean(a * b) / var) if var > 0 else 0.0
    return corr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default=None,
                    help="Seed glob for a single family (default: three-phase N=1e3 preset).")
    ap.add_argument("--label", default=None, help="Label for --seeds mode.")
    ap.add_argument("--burn", type=int, default=20, help="Burn-in sweeps (default 20).")
    ap.add_argument("--window", type=int, default=200, help="Measurement sweeps (default 200).")
    ap.add_argument("--rng-seed", type=int, default=0)
    ap.add_argument("--max-reps", type=int, default=None,
                    help="Cap replicas per phase (default: all matching).")
    ap.add_argument("--out", default="out/local_lapse",
                    help="Output prefix for <out>.json (default out/local_lapse).")
    args = ap.parse_args()

    if args.seeds:
        phases = [(args.label or os.path.basename(args.seeds), args.seeds)]
    else:
        phases = DEFAULT_PHASES
    rng = np.random.default_rng(args.rng_seed)

    results = {}
    for name, pat in phases:
        paths = sorted(glob.glob(os.path.join(_ROOT, pat))) or sorted(glob.glob(pat))
        if not paths:
            print(f"[{name}] no seeds match {pat}", file=sys.stderr)
            continue
        if args.max_reps is not None:
            paths = paths[:args.max_reps]
        deg_all, fld_all = [], {k: [] for k in ("N", "K", "A", "R")}
        corr_reps = []
        for p in paths:
            fields, F1 = measure_seed(p, args.burn, args.window)
            deg_all.append(fields["deg"])
            for k in fld_all:
                fld_all[k].append(fields[k])
            corr_reps.append(spatial_corr(fields, F1, rng))
        deg = np.concatenate(deg_all)
        cond = {k: cond_by_degree(deg, np.concatenate(v)) for k, v in fld_all.items()}
        rs = sorted(set().union(*[set(c) for c in corr_reps]))
        corr = {r: [c[r] for c in corr_reps if r in c] for r in rs}
        N = np.concatenate(fld_all["N"])
        results[name] = dict(
            replicas=len(paths), window=args.window,
            N_mean=float(N.mean()), N_cv=float(N.std() / N.mean()),
            cond={k: {str(b): v for b, v in c.items()} for k, c in cond.items()},
            spatial_corr={str(r): dict(
                              mean=float(np.mean(v)),
                              sem=float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
                              lo=float(np.min(v)), hi=float(np.max(v)), n=len(v))
                          for r, v in corr.items()},
        )
        c1 = corr.get(1, [0.0])
        sem1 = np.std(c1, ddof=1) / np.sqrt(len(c1)) if len(c1) > 1 else 0.0
        print(f"[{name}] reps={len(paths)}  ⟨N⟩={N.mean():.3f} CV={N.std()/N.mean():.2f}  "
              f"corr(r=1)={np.mean(c1):+.3f}±{sem1:.3f}  "
              f"K̄={np.concatenate(fld_all['K']).mean():.3f} "
              f"Ā={np.concatenate(fld_all['A']).mean():.3f}", flush=True)

    out = args.out if os.path.isabs(args.out) else os.path.join(_ROOT, args.out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out + ".json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {out}.json")


if __name__ == "__main__":
    main()
