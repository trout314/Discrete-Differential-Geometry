"""Test importance weights for the default sampler (no Hastings correction).

The sampler uses pure Metropolis acceptance — no Hastings correction for
proposal asymmetry. Its stationary distribution is exp(-obj)*V, so 1/V
importance weights correct it back to exp(-obj).

Compares default sampler (with 1/V weights) against exact sampler.

Usage: python3 tests/run_naive_test.py [samples] [target] [spacing] [batch_size]
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import stats


def weighted_ks_2samp(values1, weights1, values2):
    order1 = np.argsort(values1)
    sv1, cw1 = values1[order1], np.cumsum(weights1[order1])
    cw1 /= cw1[-1]
    order2 = np.argsort(values2)
    sv2, cw2 = values2[order2], np.cumsum(np.ones(len(values2)))
    cw2 /= cw2[-1]
    all_vals = np.sort(np.concatenate([sv1, sv2]))
    ecdf1 = np.searchsorted(sv1, all_vals, side="right")
    ecdf1 = np.where(ecdf1 > 0, cw1[ecdf1 - 1], 0.0)
    ecdf2 = np.searchsorted(sv2, all_vals, side="right")
    ecdf2 = np.where(ecdf2 > 0, cw2[ecdf2 - 1], 0.0)
    return np.max(np.abs(ecdf1 - ecdf2))


WORKER_SCRIPT = '''
import json, sys
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams

target, num_samples, spacing = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
mode = sys.argv[4]  # "exact" or "naive"
start_mfd_path = sys.argv[5]
outpath = sys.argv[6]

params = SamplerParams(
    num_facets_target=target, hinge_degree_target=4.5,
    num_facets_coef=0.5, num_hinges_coef=0.1,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
)

start = Manifold.load(start_mfd_path, 3)
sampler = ManifoldSampler(start, params)

exact = mode == "exact"
sampler.run(moves=1000, exact=exact)  # burn-in

num_facets, num_edges, weights = [], [], []
for i in range(num_samples):
    sampler.run(moves=spacing, exact=exact)
    m = sampler.manifold
    fv = m.f_vector
    num_facets.append(int(fv[3]))
    num_edges.append(int(fv[1]))
    weights.append(float(m.importance_weight()))

json.dump({"num_facets": num_facets, "num_edges": num_edges, "weights": weights},
          open(outpath, "w"))
'''


def run_batch(target, num_samples, spacing, mode, outpath, start_mfd_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = "python:" + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", WORKER_SCRIPT,
         str(target), str(num_samples), str(spacing),
         mode, start_mfd_path, outpath],
        env=env, capture_output=True, text=True, timeout=600,
    )
    if result.returncode != 0:
        print(f"  Worker failed: {result.stderr[-500:]}")
        return None
    with open(outpath) as f:
        return json.load(f)


def main():
    num_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    target = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    spacing = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    batch_size = int(sys.argv[4]) if len(sys.argv) > 4 else 500

    num_batches = (num_samples + batch_size - 1) // batch_size

    print(f"Naive sampler test: target={target}, spacing={spacing}")
    print(f"Collecting {num_samples} samples in {num_batches} batches of {batch_size}")

    tmpdir = tempfile.mkdtemp()

    # Create shared starting manifold
    print("Growing starting manifold...", end=" ", flush=True)
    start_mfd_path = os.path.join(tmpdir, "start.mfd")
    env = os.environ.copy()
    env["PYTHONPATH"] = "python:" + env.get("PYTHONPATH", "")
    subprocess.run(
        [sys.executable, "-c", f"""
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
params = SamplerParams(num_facets_target={target}, hinge_degree_target=4.5,
    num_facets_coef=0.5, num_hinges_coef=0.1,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0)
seed = Manifold.standard_sphere(3)
grower = ManifoldSampler(seed, params)
grower.run(moves=5000, exact=False)
grower.manifold.save("{start_mfd_path}")
"""], env=env, check=True, capture_output=True, timeout=120)
    print("done\n")

    for mode in ["exact", "naive"]:
        all_facets, all_edges, all_weights = [], [], []
        label = "Exact" if mode == "exact" else "Naive"

        for b in range(num_batches):
            n = min(batch_size, num_samples - b * batch_size)
            outpath = os.path.join(tmpdir, f"{mode}_{b}.json")
            print(f"  {label} batch {b+1}/{num_batches} ({n} samples)...", end=" ", flush=True)
            data = run_batch(target, n, spacing, mode, outpath, start_mfd_path)
            if data is None:
                return 1
            all_facets.extend(data["num_facets"])
            all_edges.extend(data["num_edges"])
            all_weights.extend(data["weights"])
            print("done")

        np.savez(os.path.join(tmpdir, f"{mode}.npz"),
                 num_facets=all_facets, num_edges=all_edges, weights=all_weights)

    exact = np.load(os.path.join(tmpdir, "exact.npz"))
    naive = np.load(os.path.join(tmpdir, "naive.npz"))

    alpha = 0.01
    n = len(exact["num_facets"])
    ks_critical = np.sqrt(-0.5 * np.log(alpha / 2)) * np.sqrt(2.0 / n)

    print(f"\nResults (n={n}, KS critical={ks_critical:.4f}):")
    for key in ["num_facets", "num_edges"]:
        ks_weighted = weighted_ks_2samp(naive[key], naive["weights"], exact[key])
        ks_unweighted, _ = stats.ks_2samp(naive[key], exact[key])

        weighted_ok = ks_weighted < ks_critical
        unweighted_differs = ks_unweighted > ks_critical

        if weighted_ok and unweighted_differs:
            status = "PASS (weights necessary and correct)"
        elif weighted_ok and not unweighted_differs:
            status = "WEAK PASS (weights correct but test not discriminating)"
        else:
            status = "FAIL (weighted distribution rejected)"
        print(f"  {key}: KS_weighted={ks_weighted:.4f}, KS_unweighted={ks_unweighted:.4f} [{status}]")

    w = naive["weights"]
    print(f"\n  Weight stats (1/V): mean={np.mean(w):.6f}, std={np.std(w):.6f}, "
          f"min={np.min(w):.6f}, max={np.max(w):.6f}")
    print(f"  Weight CV (std/mean): {np.std(w)/np.mean(w):.2f}")
    return 0


if __name__ == "__main__":
    exit(main())
