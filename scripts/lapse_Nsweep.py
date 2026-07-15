#!/usr/bin/env python3
"""Cross-N heterogeneity: does the lapse correlation LENGTH grow with N?

N=1e4 full g-ladder (room for a length: diameter ~50) + N=1e5 g=2e-3 point.
N=1e3 comes from the earlier out/lapse_g_sweep.json. Fixed edge 5.0043, k=2.
"""
import os, sys, glob, json, time
import numpy as np
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
from local_lapse import measure_seed, spatial_corr

BURN, WINDOW = 20, 200

# cells: (N_tag, N, g_tag, g, filename-stem, reps)
CELLS = []
for gt, g in [("0", 0.0), ("1e-3", 1e-3), ("2e-3", 2e-3), ("4e-3", 4e-3),
              ("8e-3", 8e-3), ("1e-2", 1e-2)]:
    stem = "S3_N1e4_1e-1_ED5p0043_2" + ("" if g == 0 else f"_VDVs_{gt}")
    CELLS.append(("1e4", 10000, gt, g, stem, 12))
CELLS.append(("1e5", 100000, "2e-3", 2e-3, "S3_N1e5_1e-1_ED5p0043_2_VDVs_2e-3", 4))

def xi_from_corr(cmean):
    rs = sorted(cmean); pos = []
    for r in rs:
        if cmean[r] > 0: pos.append(r)
        else: break
    if len(pos) < 2: return np.nan
    x = np.array(pos, float); y = np.log([cmean[r] for r in pos])
    s = np.polyfit(x, y, 1)[0]
    return -1.0 / s if s < 0 else np.nan

rng = np.random.default_rng(0)
results = {}
for Ntag, N, gt, g, stem, reps in CELLS:
    paths = sorted(glob.glob(os.path.join(_ROOT, "seeds", stem + "_s0*.mfd")))[:reps]
    if not paths:
        print(f"[{Ntag} g={gt}] no seeds", flush=True); continue
    t0 = time.time()
    corr_reps, Nm, cvs = [], [], []
    for p in paths:
        f, F1 = measure_seed(p, BURN, WINDOW)
        corr_reps.append(spatial_corr(f, F1, rng, rmax=10, nsrc=250))
        Nm.append(float(f["N"].mean())); cvs.append(float(f["N"].std() / f["N"].mean()))
    rs = sorted(set().union(*[set(c) for c in corr_reps]))
    cmean = {r: float(np.mean([c[r] for c in corr_reps if r in c])) for r in rs}
    csem = {r: float(np.std([c[r] for c in corr_reps if r in c], ddof=1)
                     / np.sqrt(sum(r in c for c in corr_reps))) for r in rs}
    xi_full = xi_from_corr(cmean)
    n = len(corr_reps); xis = []
    for i in range(n):
        sub = [corr_reps[j] for j in range(n) if j != i]
        cm = {r: np.mean([c[r] for c in sub if r in c]) for r in rs}
        xis.append(xi_from_corr(cm))
    xis = np.array([x for x in xis if np.isfinite(x)])
    xi_sem = float(np.sqrt((len(xis) - 1) / len(xis) * np.sum((xis - xis.mean()) ** 2))) \
             if len(xis) > 1 else float("nan")
    results[f"{Ntag}_{gt}"] = dict(N=N, g=g, reps=n, N_mean=float(np.mean(Nm)),
                                   N_cv=float(np.mean(cvs)),
                                   corr={str(r): [cmean[r], csem[r]] for r in rs},
                                   xi=float(xi_full), xi_sem=xi_sem)
    print(f"[N={Ntag} g={gt:>4}] reps={n}  corr(1)={cmean[1]:+.3f}±{csem[1]:.3f} "
          f"corr(2)={cmean.get(2,0):+.3f} corr(3)={cmean.get(3,0):+.3f}  "
          f"xi={xi_full:.2f}±{xi_sem:.2f}  CV={np.mean(cvs):.3f}  ({time.time()-t0:.0f}s)", flush=True)

OUT = os.path.join(_ROOT, "out"); os.makedirs(OUT, exist_ok=True)
json.dump(results, open(os.path.join(OUT,"lapse_Nsweep.json"), "w"), indent=1)
print("wrote out/lapse_Nsweep.json")
