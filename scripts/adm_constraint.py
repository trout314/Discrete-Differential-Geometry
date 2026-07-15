#!/usr/bin/env python3
"""ADM constraint test (pass 1, N=1e3).

Identification under test: volume pin = discrete maximal slicing (trK ~ 0);
move activity = extrinsic-curvature fluctuation. Hamiltonian-constraint form:
does local intrinsic curvature R(v) relate affinely to local K^2(v),
  R(v) = a + lambda * K2(v) + residual,
with a ~ 2*Lambda_eff and lambda calibrating MC time to geometric units?

Per-vertex fields over short windows (W sweeps, quasi-static geometry):
  R  = sum_{e ni v} delta_e / (2 * D_v/4),  delta_e = 2*pi - theta * d_e
       (Regge deficit density per unit 3-volume; averaged over window endpoints)
  trK = volume_flux / W / (D_v/4)          (signed; global mean ~ 0 check)
  K2a = (lapse rate per volume)^2          (total-activity proxy)
  K2d = sum_{e ni v} (deficit rate)^2/(D_v/4)  (curvature-velocity proxy)
"""
import os, sys, glob, json, time
import numpy as np
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python")); sys.path.insert(0, os.path.join(_ROOT, "tools"))
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from discrete_differential_geometry import move_geometry as mg
from seed_utils import load_seed_metadata

THETA = mg.THETA_TET
ENSEMBLES = [
    ("pos_2e-3",  "S3_N1e3_1e-1_ED5p0043_2_VDVs_2e-3"),
    ("flat_2e-3", "S3_N1e3_1e-1_ED5p1043_2_VDVs_2e-3"),
    ("neg_2e-3",  "S3_N1e3_1e-1_ED5p2043_2_VDVs_2e-3"),
    ("pos_8e-3_k4", "S3_N1e3_1e-1_ED5p0043_4_VDVs_8e-3"),
    ("hub_b0",    "S3_N1e3_1e-1_ED5p1043_2"),
]
REPS, BURN, NWIN, W = 6, 20, 30, 5

def params_of(path):
    md = load_seed_metadata(path)
    return SamplerParams(num_facets_target=int(md["num_facets_target"]),
        hinge_degree_target=float(md["hinge_degree_target"]),
        num_facets_coef=float(md["num_facets_coef"]),
        num_hinges_coef=float(md["num_hinges_coef"]),
        hinge_degree_variance_coef=float(md["hinge_degree_variance_coef"]),
        codim3_degree_variance_coef=float(md["codim3_degree_variance_coef"]))

def fields_from_facets(F):
    """Per-vertex: degree D_v and deficit sum Sum_{e ni v} delta_e."""
    lab, cnt = np.unique(F, return_counts=True)
    D = dict(zip(lab.tolist(), cnt.tolist()))
    E = {}
    pr = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for i, j in pr:
        for a, b in np.sort(F[:, [i, j]], axis=1):
            k = (int(a), int(b)); E[k] = E.get(k, 0) + 1
    dsum = {v: 0.0 for v in D}
    for (a, b), d in E.items():
        delta = 2*np.pi - THETA*d
        dsum[a] += delta; dsum[b] += delta
    return D, dsum

results = {}
for tag, stem in ENSEMBLES:
    t0 = time.time()
    pool = {k: [] for k in ("R","trK","K2a","K2d")}
    trK_glob = []
    for path in sorted(glob.glob(f"{_ROOT}/seeds/{stem}_s0*.mfd"))[:REPS]:
        s = ManifoldSampler(Manifold.load(path, 3), params_of(path))
        s.run(sweeps=BURN)
        s.track_geometry()
        Fprev = np.asarray(s.manifold.facets(), np.int64)
        Dp, DSp = fields_from_facets(Fprev)
        for _ in range(NWIN):
            s.reset_geometry()
            s.run(sweeps=W)
            F = np.asarray(s.manifold.facets(), np.int64)
            Dn, DSn = fields_from_facets(F)
            vr = s.vertex_role_counts(); er = s.edge_role_counts()
            vlab = vr["vertex"]; vc = vr["counts"]
            # per-vertex deficit-velocity^2 from the edge ledger
            edot = er["counts"] @ mg.E_DELTA            # Δdeg per edge over window
            k2d_sum = {}
            for (a, b), dd in zip(er["edge"].tolist(), edot):
                val = (THETA * dd / W) ** 2
                k2d_sum[a] = k2d_sum.get(a, 0.0) + val
                k2d_sum[b] = k2d_sum.get(b, 0.0) + val
            vflux = mg.volume_flux(vc); act = mg.lapse(vc)
            bd = vc @ (mg.V_BIRTH + mg.V_DEATH)
            for i, v in enumerate(vlab.tolist()):
                if bd[i] > 0 or v not in Dp or v not in Dn:
                    continue
                vol = 0.5 * (Dp[v] + Dn[v]) / 4.0
                R = 0.5 * (DSp[v] + DSn[v]) / 2.0 / vol
                trk = vflux[i] / W / vol
                pool["R"].append(R); pool["trK"].append(trk)
                pool["K2a"].append((act[i] / W / vol) ** 2)
                pool["K2d"].append(k2d_sum.get(v, 0.0) / vol)
            trK_glob.append(vflux.sum() / W)          # net volume rate (tets/sweep)
            Fprev, Dp, DSp = F, Dn, DSn
    R = np.array(pool["R"]); trK = np.array(pool["trK"])
    K2a = np.array(pool["K2a"]); K2d = np.array(pool["K2d"])
    out = dict(n=len(R), R_mean=float(R.mean()), R_std=float(R.std()),
               trK_mean=float(trK.mean()), trK_std=float(trK.std()),
               trK_glob_mean=float(np.mean(trK_glob)), trK_glob_std=float(np.std(trK_glob)),
               K2a_mean=float(K2a.mean()), K2d_mean=float(K2d.mean()))
    for name, K2 in (("K2a", K2a), ("K2d", K2d)):
        lam, a = np.polyfit(K2, R, 1)
        res = R - (a + lam * K2)
        out[f"fit_{name}"] = dict(lam=float(lam), a=float(a),
                                  r=float(np.corrcoef(K2, R)[0, 1]),
                                  res_frac=float(res.var() / R.var()))
        # binned conditional for the figure
        q = np.quantile(K2, np.linspace(0, 1, 13))
        cb = []
        for lo, hi in zip(q[:-1], q[1:]):
            m = (K2 >= lo) & (K2 < hi)
            if m.sum() > 20:
                cb.append([float(K2[m].mean()), float(R[m].mean()),
                           float(R[m].std() / np.sqrt(m.sum()))])
        out[f"cond_{name}"] = cb
    # trK-R correlation (should be ~0 if slicing is maximal pointwise-on-average)
    out["corr_trK_R"] = float(np.corrcoef(trK, R)[0, 1])
    results[tag] = out
    fa = out["fit_K2a"]; fd = out["fit_K2d"]
    print(f"[{tag:>11}] n={len(R):>6}  <R>={R.mean():+.3f}  <trK>={trK.mean():+.2e} "
          f"(glob {np.mean(trK_glob):+.3f}±{np.std(trK_glob):.2f} tets/sweep)\n"
          f"    R~K2a: lam={fa['lam']:+.3f} a={fa['a']:+.3f} r={fa['r']:+.3f} resfrac={fa['res_frac']:.3f} | "
          f"R~K2d: lam={fd['lam']:+.3f} a={fd['a']:+.3f} r={fd['r']:+.3f} resfrac={fd['res_frac']:.3f} "
          f"({time.time()-t0:.0f}s)", flush=True)

OUT = os.path.join(_ROOT, "out"); os.makedirs(OUT, exist_ok=True)
json.dump(results, open(os.path.join(OUT,"adm_constraint.json"), "w"), indent=1)
print("wrote out/adm_constraint.json")
