#!/usr/bin/env python3
"""Stage-1 prototype: measured combinatorial lapse N(v) from per-sweep facet diffs.

For each phase's certified N=1e3 seed, run the sampler AT the family objective and,
each sweep, diff the facet set against the previous sweep. Every changed tet
(appeared or disappeared) attributes one event to each of its 4 vertices. Then

    N(v) = events(v) / sweeps / (D_v / 4),

the accepted-activity density per unit local 3-volume (D_v = vertex degree,
averaged over window endpoints). Analyses: conditional mean N vs D_v, activity
share vs degree, and spatial autocorrelation of dN = N - <N|D> (residual after
removing the degree dependence) along the 1-skeleton -> heterogeneity length.
"""
import os, sys, json
import numpy as np
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python")); sys.path.insert(0, os.path.join(_ROOT, "tools"))
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from seed_utils import load_seed_metadata
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

PHASES = [
    ("extended g=2e-3",      "seeds/S3_N1e3_1e-1_ED5p0043_2_VDVs_2e-3_s000.mfd"),
    ("hub beta=0 flat-pin",  "seeds/S3_N1e3_1e-1_ED5p1043_2_s000.mfd"),
    ("glass-adj g=8e-3 k=4", "seeds/S3_N1e3_1e-1_ED5p0043_4_VDVs_8e-3_s000.mfd"),
]
BURN, WINDOW = 20, 200          # sweeps; window >> per-vertex event resolution, ~ few tau
rng = np.random.default_rng(0)

def facet_key(F):
    """Set of facet tuples (already sorted rows from the D core)."""
    return set(map(tuple, np.sort(F, axis=1).tolist()))

def degrees(F):
    """Per-vertex-label degree dict from a facet array."""
    lab, cnt = np.unique(F, return_counts=True)
    return dict(zip(lab.tolist(), cnt.tolist()))

results = {}
for name, seed in PHASES:
    path = os.path.join(_ROOT, seed)
    md = load_seed_metadata(path)
    params = SamplerParams(
        num_facets_target=int(md["num_facets_target"]),
        hinge_degree_target=float(md["hinge_degree_target"]),
        num_facets_coef=float(md["num_facets_coef"]),
        num_hinges_coef=float(md["num_hinges_coef"]),
        hinge_degree_variance_coef=float(md["hinge_degree_variance_coef"]),
        codim3_degree_variance_coef=float(md["codim3_degree_variance_coef"]),
    )
    s = ManifoldSampler(Manifold.load(path, 3), params)
    s.run(sweeps=BURN)
    F0 = np.asarray(s.manifold.facets(), np.int64)
    prev = facet_key(F0); deg0 = degrees(F0)
    events = {}
    changed_per_sweep = []
    for _ in range(WINDOW):
        s.run(sweeps=1)
        F = np.asarray(s.manifold.facets(), np.int64)
        cur = facet_key(F)
        diff = prev ^ cur                      # tets created or destroyed this sweep
        changed_per_sweep.append(len(diff))
        for tet in diff:
            for v in tet:
                events[v] = events.get(v, 0) + 1
        prev = cur
    F1 = F; deg1 = degrees(F1)
    # vertices present at both window endpoints (identity-stable subset)
    stable = sorted(set(deg0) & set(deg1))
    frac_stable = len(stable) / len(set(deg0) | set(deg1))
    D = np.array([(deg0[v] + deg1[v]) / 2.0 for v in stable])
    E = np.array([events.get(v, 0) for v in stable], float)
    Nlapse = E / WINDOW / (D / 4.0)

    # conditional mean N | degree (even-degree bins)
    cond = {}
    for d, n in zip(D, Nlapse):
        b = 2 * round(d / 2)
        cond.setdefault(b, []).append(n)
    cond = {b: (np.mean(v), np.std(v) / np.sqrt(len(v)), len(v)) for b, v in sorted(cond.items()) if len(v) >= 3}

    # spatial autocorrelation of residual dN on the 1-skeleton (window-end graph)
    idx = {v: i for i, v in enumerate(stable)}
    pr = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    ee = np.unique(np.sort(np.vstack([F1[:, [i, j]] for i, j in pr]), axis=1), axis=0)
    m = np.array([(a in idx) and (b in idx) for a, b in ee])
    ee = np.array([[idx[a], idx[b]] for a, b in ee[m]])
    n = len(stable)
    G = coo_matrix((np.ones(2 * len(ee)), (np.r_[ee[:, 0], ee[:, 1]], np.r_[ee[:, 1], ee[:, 0]])), shape=(n, n)).tocsr()
    # residual after removing degree dependence (else deg-deg correlation leaks in)
    mu = {b: c[0] for b, c in cond.items()}
    dN = np.array([nl - mu.get(2 * round(d / 2), nl) for nl, d in zip(Nlapse, D)])
    src = rng.choice(n, min(120, n), replace=False)
    dist = shortest_path(G, method="D", unweighted=True, directed=False, indices=src)
    corr = {}
    var = dN.var()
    for r in range(1, 8):
        pairs = [(i, j) for si, i in enumerate(src) for j in np.where(dist[si] == r)[0]]
        if len(pairs) < 30: break
        a = np.array([dN[i] for i, _ in pairs]); b = np.array([dN[j] for _, j in pairs])
        corr[r] = float(np.mean(a * b) / var) if var > 0 else 0.0

    results[name] = dict(
        n_stable=len(stable), frac_stable=round(frac_stable, 3),
        moves_per_sweep=round(np.mean(changed_per_sweep) / 2, 1),   # ~2 tets changed per move avg-ish
        N_mean=round(float(Nlapse.mean()), 4), N_cv=round(float(Nlapse.std() / Nlapse.mean()), 3),
        frozen_frac=round(float((E == 0).mean()), 3),
        cond_N_by_deg={int(b): (round(c[0], 4), round(c[1], 4), c[2]) for b, c in cond.items()},
        spatial_corr={r: round(c, 3) for r, c in corr.items()},
    )
    print(f"[{name}] stable={len(stable)}/{frac_stable:.0%}  ⟨N⟩={Nlapse.mean():.3f} "
          f"CV={Nlapse.std()/Nlapse.mean():.2f}  frozen={(E==0).mean():.1%}  "
          f"tet-changes/sweep={np.mean(changed_per_sweep):.0f}", flush=True)

OUT = os.path.join(_ROOT, "out")
os.makedirs(OUT,exist_ok=True); json.dump(results, open(os.path.join(OUT,"local_lapse_proto.json"),"w"), indent=1)
print("\nwrote out/local_lapse_proto.json")
