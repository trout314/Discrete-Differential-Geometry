#!/usr/bin/env python3
"""Certify the strain-gas densities: multi-chain R-hat over the volume cells.

For each strain cell (crystal x vol-scale, channel on), pools the original
chain with the dedicated certification seeds and reports, on the late window
(post-burn second half) of n_ill and f0:
  * quantized_split_rhat (quantum = 1: both observables are integer counts),
    two-sided certification per CONVENTIONS;
  * pooled mean +- pooled SD;
  * per-chain late-window slope +- moving-block-bootstrap sigma (stationarity
    is claimed only if every chain's slope is consistent with zero).

Usage:
    python scripts/defect_dynamics/certify_strain_gas.py
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry.convergence import quantized_split_rhat

CELLS = {
    ("c15", "0.95"): ["c15_v0.95_cs"] + [f"c15_v0.95_cs_s{s}"
                                         for s in (195, 295, 395)],
    ("c15", "1.05"): ["c15_v1.05_cs"] + [f"c15_v1.05_cs_s{s}"
                                         for s in (196, 296, 396)],
    ("a15", "0.95"): ["a15_v0.95_cs"] + [f"a15_v0.95_cs_s{s}"
                                         for s in (197, 297, 397)],
    ("a15", "1.05"): ["a15_v1.05_cs"] + [f"a15_v1.05_cs_s{s}"
                                         for s in (198, 298, 398)],
}
RHAT_GATE = 1.01


def late_series(stem):
    path = os.path.join(_ROOT, "data", "crystal_gas", stem + ".obs.jsonl")
    if not os.path.exists(path):
        return None
    rows = [json.loads(l) for l in open(path)]
    data = [r for r in rows if "n_ill" in r and r.get("burned")]
    half = data[len(data) // 2:]
    return (np.array([r["n_ill"] for r in half], float),
            np.array([float(r["f"][0]) for r in half]))


def slope_bootstrap(x, nblock=8, nboot=400, seed=0):
    x = np.asarray(x, float)
    t = np.arange(len(x))
    slope = np.polyfit(t, x, 1)[0]
    blk = len(x) // nblock
    if blk < 2:
        return slope, float("nan")
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nboot):
        idx = np.concatenate([np.arange(s, s + blk) for s in
                              rng.integers(0, len(x) - blk, nblock)])
        bs.append(np.polyfit(np.arange(len(idx)), x[idx], 1)[0])
    return slope, float(np.std(bs))


def main():
    all_ok = True
    for (crystal, vol), stems in CELLS.items():
        chains = [late_series(s) for s in stems]
        missing = [s for s, c in zip(stems, chains) if c is None]
        chains = [c for c in chains if c is not None]
        print(f"\n=== {crystal} vol x{vol}: {len(chains)} chains"
              + (f" (missing: {missing})" if missing else ""))
        if len(chains) < 2:
            print("  not enough chains yet")
            all_ok = False
            continue
        n = min(len(c[0]) for c in chains)
        for name, k in (("n_ill", 0), ("f0", 1)):
            series = [c[k][-n:] for c in chains]
            rhat = quantized_split_rhat(series, 1.0)
            pooled = np.concatenate(series)
            verdict = "PASS" if rhat < RHAT_GATE else "FAIL"
            print(f"  {name:5s}: Rhat = {rhat:.4f} [{verdict} @ {RHAT_GATE}] "
                  f" pooled {pooled.mean():.1f} +- {pooled.std():.1f}")
            worst_z = 0.0
            for stem, s in zip(stems, series):
                sl, err = slope_bootstrap(s)
                z = sl / err if err > 0 else float("nan")
                worst_z = max(worst_z, abs(z))
                if abs(z) > 3:
                    print(f"    NONSTATIONARY {stem}: slope {sl:+.3f} "
                          f"+- {err:.3f} (z={z:+.1f})")
            print(f"    worst |slope z| across chains: {worst_z:.1f}")
            if rhat >= RHAT_GATE or worst_z > 3:
                all_ok = False
    print(f"\nOVERALL: {'CERTIFIED' if all_ok else 'NOT CERTIFIED'}")


if __name__ == "__main__":
    main()
