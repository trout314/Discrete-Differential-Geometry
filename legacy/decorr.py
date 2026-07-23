#!/usr/bin/env python3
"""Decorrelation diagnostic from the run5h ledger: is the defect configuration
mobile (long chains buy samples) or glassy (need many independent chains)?

(1) IAT of the scalar n_illegal(t) per chain.
(2) Spatial persistence: Jaccard overlap of the illegal-vertex SET vs lag; where
    it decays to the random baseline is the pattern decorrelation time.
(3) Survival of the initial defect cores; cross-chain basin separation.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

SP = sys.argv[1]
V = 10176


def load(f):
    recs = [json.loads(l) for l in open(f)]
    seen, out = set(), []
    for r in sorted(recs, key=lambda r: r["sweep"]):
        if r["sweep"] in seen:
            continue
        seen.add(r["sweep"])
        ill = set()
        for c in r.get("comps", []):
            ill |= set(c["members"])
        out.append((r["sweep"], r["n_illegal"], r.get("n_5knots", 0), ill))
    return out


def iat_samples(x):
    x = np.asarray(x, float)
    if x.std() == 0:
        return np.inf
    x = x - x.mean()
    n = len(x)
    acf = np.correlate(x, x, "full")[n - 1:] / (x @ x)
    tau = 1.0
    for k in range(1, n):
        if acf[k] < 0:
            break
        tau += 2 * acf[k]
        if k > 5 * tau:
            break
    return tau


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


files = sorted(glob.glob(f"{SP}/run5h/*.ts.jsonl"))
lagJ = defaultdict(list)          # lag(sweeps) -> jaccard list (pooled)
survival = defaultdict(list)      # lag -> frac of initial cores still illegal
iats, means = [], []
print(f"{'chain':10s} {'<n_ill>':>8s} {'IAT(sw)':>8s}  {'note'}")
for f in files:
    rec = load(f)
    sw = np.array([r[0] for r in rec])
    nill = np.array([r[1] for r in rec])
    sets = [r[3] for r in rec]
    tau = iat_samples(nill)
    dsw = np.median(np.diff(sw))
    tau_sw = tau * dsw if np.isfinite(tau) else np.inf
    iats.append(tau_sw); means.append(nill.mean())
    print(f"{os.path.basename(f)[:-9]:10s} {nill.mean():8.1f} {tau_sw:8.0f}"
          f"  {'FROZEN' if not np.isfinite(tau_sw) else ''}")
    # pooled Jaccard vs lag + core survival from t0
    for i in range(len(rec)):
        for j in range(i, len(rec)):
            lagJ[sw[j] - sw[i]].append(jaccard(sets[i], sets[j]))
    s0 = sets[0]
    for j in range(len(rec)):
        survival[sw[j] - sw[0]].append(len(s0 & sets[j]) / max(len(s0), 1))

base = np.mean([len(sets[0])])  # typical size
print(f"\nrandom-set Jaccard baseline ~ {base/(2*V):.4f}  (illegal-set size ~{base:.0f}, V={V})")
print(f"\n{'lag(sw)':>8s} {'Jaccard':>8s} {'core-surv':>9s}")
for lag in sorted(lagJ):
    if lag % 3000 == 0 or lag in (150, 750, 1500):
        j = np.mean(lagJ[lag])
        sv = np.mean(survival.get(lag, [np.nan]))
        print(f"{lag:8d} {j:8.3f} {sv:9.3f}")

# cross-chain: mean n_illegal spread + inter-chain final-set overlap
finals = [load(f)[-1][3] for f in files]
cross = [jaccard(finals[i], finals[j]) for i in range(len(finals))
         for j in range(i + 1, len(finals))]
print(f"\ncross-chain <n_illegal> = {np.mean(means):.1f} +/- {np.std(means):.1f}"
      f"   (spread/mean = {np.std(means)/np.mean(means):.2%})")
print(f"inter-chain illegal-set Jaccard = {np.mean(cross):.4f} (vs baseline "
      f"{base/(2*V):.4f}) -> {'SHARED cores' if np.mean(cross) > 0.05 else 'independent basins'}")
print(f"\nmedian within-chain IAT = {np.median([t for t in iats if np.isfinite(t)]):.0f} sweeps")
