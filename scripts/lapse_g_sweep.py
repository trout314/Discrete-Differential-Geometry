#!/usr/bin/env python3
"""xi_N(g): dynamical-heterogeneity length of the measured lapse vs VDV coupling.

Fixed edge (5.0043), fixed k=2, N=1e3; g = beta/N ladder {0, 1e-3, 2e-3, 4e-3,
8e-3}. For each g, measure the per-vertex lapse field over all replicas (exact
Option-A counters), compute the spatial autocorrelation of the degree-residual
dN on the 1-skeleton per replica, then aggregate corr(r) mean +/- SEM and fit a
correlation length xi from the positive decaying part. Errors on corr via SEM
over replicas; error on xi via jackknife over replicas.
"""
import os, sys, glob, json
import numpy as np
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
from local_lapse import measure_seed, spatial_corr

LADDER = [("0", "S3_N1e3_1e-1_ED5p0043_2_s0*.mfd", 0.0),
          ("1e-3", "S3_N1e3_1e-1_ED5p0043_2_VDVs_1e-3_s0*.mfd", 1e-3),
          ("2e-3", "S3_N1e3_1e-1_ED5p0043_2_VDVs_2e-3_s0*.mfd", 2e-3),
          ("4e-3", "S3_N1e3_1e-1_ED5p0043_2_VDVs_4e-3_s0*.mfd", 4e-3),
          ("8e-3", "S3_N1e3_1e-1_ED5p0043_2_VDVs_8e-3_s0*.mfd", 8e-3)]
BURN, WINDOW = 20, 200
rng = np.random.default_rng(0)

def xi_from_corr(cmean):
    """Corr-length from exp fit A*exp(-r/xi) over the leading positive run r>=1."""
    rs = sorted(cmean)
    pos = []
    for r in rs:                         # take the contiguous positive run from r=1
        if cmean[r] > 0: pos.append(r)
        else: break
    if len(pos) < 2:
        return np.nan
    x = np.array(pos, float); y = np.log(np.array([cmean[r] for r in pos]))
    slope = np.polyfit(x, y, 1)[0]
    return -1.0 / slope if slope < 0 else np.nan

results = {}
for tag, pat, g in LADDER:
    paths = sorted(glob.glob(os.path.join(_ROOT, "seeds", pat)))
    corr_reps, Nmeans, cvs = [], [], []
    for p in paths:
        fields, F1 = measure_seed(p, BURN, WINDOW)
        corr_reps.append(spatial_corr(fields, F1, rng))
        Nmeans.append(float(fields["N"].mean()))
        cvs.append(float(fields["N"].std() / fields["N"].mean()))
    rs = sorted(set().union(*[set(c) for c in corr_reps]))
    # per-r mean +/- SEM across replicas
    cmean = {r: float(np.mean([c[r] for c in corr_reps if r in c])) for r in rs}
    csem = {r: float(np.std([c[r] for c in corr_reps if r in c], ddof=1)
                     / np.sqrt(sum(r in c for c in corr_reps))) for r in rs}
    # xi jackknife over replicas
    xi_full = xi_from_corr(cmean)
    n = len(corr_reps); xis = []
    for i in range(n):
        sub = [corr_reps[j] for j in range(n) if j != i]
        cm = {r: np.mean([c[r] for c in sub if r in c]) for r in rs}
        xis.append(xi_from_corr(cm))
    xis = np.array([x for x in xis if np.isfinite(x)])
    xi_sem = (np.sqrt((len(xis) - 1) / len(xis) * np.sum((xis - xis.mean()) ** 2))
              if len(xis) > 1 else np.nan)
    results[tag] = dict(g=g, reps=n, N_mean=float(np.mean(Nmeans)), N_cv=float(np.mean(cvs)),
                        corr={str(r): [cmean[r], csem[r]] for r in rs},
                        xi=float(xi_full), xi_sem=float(xi_sem))
    print(f"g={tag:>4} (={g:.0e}): reps={n}  corr(1)={cmean[1]:+.3f}±{csem[1]:.3f}  "
          f"corr(2)={cmean.get(2,0):+.3f}  xi={xi_full:.2f}±{xi_sem:.2f}  "
          f"CV={np.mean(cvs):.3f}", flush=True)

OUT = os.path.join(_ROOT, "out"); os.makedirs(OUT, exist_ok=True)
json.dump(results, open(os.path.join(OUT,"lapse_g_sweep.json"), "w"), indent=1)
print("wrote out/lapse_g_sweep.json")
