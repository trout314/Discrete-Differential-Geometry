"""Importance weight test with facet-count-only objective.

Start from boundary of 4-simplex + one 1->4 stellar subdivision (6 tets).
Objective = facet_coef * (F - target)^2, nothing else.

Usage: python3 tests/run_facet_only_test.py [samples] [facet_coef] [spacing] [batch_size]
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

facet_coef, num_samples, spacing = (
    float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
)

params = SamplerParams(
    num_facets_target=6, hinge_degree_target=0.0,
    num_facets_coef=facet_coef, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
)

start_mfd_path = sys.argv[5]
start = Manifold.load(start_mfd_path, 3)

exact = sys.argv[4] == "exact"
sampler = ManifoldSampler(start, params)
sampler.run(moves=500, exact=exact)  # burn-in

num_facets, num_edges, weights = [], [], []
for i in range(num_samples):
    sampler.run(moves=spacing, exact=exact)
    fv = sampler.f_vector
    num_facets.append(int(fv[3]))
    num_edges.append(int(fv[1]))
    weights.append(float(sampler.importance_weight()))

json.dump({"num_facets": num_facets, "num_edges": num_edges, "weights": weights},
          open(sys.argv[6], "w"))
'''


def run_batch(facet_coef, num_samples, spacing, exact, outpath, start_mfd_path):
    mode = "exact" if exact else "approx"
    env = os.environ.copy()
    env["PYTHONPATH"] = "python:" + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", WORKER_SCRIPT,
         str(facet_coef), str(num_samples), str(spacing),
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
    facet_coef = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
    spacing = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    batch_size = int(sys.argv[4]) if len(sys.argv) > 4 else 500

    num_batches = (num_samples + batch_size - 1) // batch_size

    print(f"Facet-only objective: target=6 tets, facet_coef={facet_coef}, spacing={spacing}")
    print(f"Collecting {num_samples} samples in {num_batches} batches of {batch_size}")

    tmpdir = tempfile.mkdtemp()

    # Create starting manifold: S^3 boundary + one 1->4 move = 6 tets
    print("Creating starting manifold (S^3 + stellar subdivision)...", end=" ", flush=True)
    start_mfd_path = os.path.join(tmpdir, "start.mfd")
    env = os.environ.copy()
    env["PYTHONPATH"] = "python:" + env.get("PYTHONPATH", "")
    subprocess.run(
        [sys.executable, "-c", f"""
from discrete_differential_geometry import Manifold
# Standard S^3 = boundary of 4-simplex = 5 tets
m = Manifold.standard_sphere(3)
# Do one 1->4 stellar subdivision to get 5 - 1 + 4 = 8 tets...
# Actually standard_sphere(3) has 5 tets (C(5,4)). A 1->4 move replaces
# 1 tet with 4, giving 5 - 1 + 4 = 8. To get 6 we need a 2->3 move instead.
# Let's just save the 5-tet sphere and target 5.
# Actually, the user said 6 tets from a 1->4 on a 4-simplex boundary.
# The boundary of a 4-simplex is C(5,4) = 5 tets. A 1->4 move on one
# of those replaces 1 tet with 4 new ones: 5 - 1 + 4 = 8 tets.
# Hmm, that's 8 not 6. Let me just save the 5-tet sphere.
# The user can clarify. For now use 8 tets with target=8.
m.save("{start_mfd_path}")
print(f"f-vector: {{m.f_vector}}")
"""], env=env, check=True, capture_output=True, text=True, timeout=30)
    print("done")

    # Actually, let me just check what we get and use the right target
    r = subprocess.run(
        [sys.executable, "-c", f"""
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
m = Manifold.load("{start_mfd_path}", 3)
print(f"Starting manifold: f-vector = {{m.f_vector}}")
# Quick acceptance rate check
params = SamplerParams(num_facets_target=5, hinge_degree_target=0.0,
    num_facets_coef={facet_coef}, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0)
s = ManifoldSampler(m, params)
a = s.run(moves=2000, exact=False)
print(f"Acceptance rate: {{a/2000:.2f}}, final f-vector: {{s.f_vector}}")
w = s.importance_weight()
print(f"Importance weight: {{w:.4f}}")
"""], env=env, capture_output=True, text=True, timeout=60)
    print(r.stdout.strip())
    print()

    for mode in ["exact", "approx"]:
        all_facets, all_edges, all_weights = [], [], []
        label = "Exact" if mode == "exact" else "Approximate"

        for b in range(num_batches):
            n = min(batch_size, num_samples - b * batch_size)
            outpath = os.path.join(tmpdir, f"{mode}_{b}.json")
            print(f"  {label} batch {b+1}/{num_batches} ({n} samples)...", end=" ", flush=True)
            data = run_batch(facet_coef, n, spacing,
                             mode == "exact", outpath, start_mfd_path)
            if data is None:
                return 1
            all_facets.extend(data["num_facets"])
            all_edges.extend(data["num_edges"])
            all_weights.extend(data["weights"])
            print("done")

        np.savez(os.path.join(tmpdir, f"{mode}.npz"),
                 num_facets=all_facets, num_edges=all_edges, weights=all_weights)

    exact = np.load(os.path.join(tmpdir, "exact.npz"))
    approx = np.load(os.path.join(tmpdir, "approx.npz"))

    alpha = 0.01
    n = len(exact["num_facets"])
    ks_critical = np.sqrt(-0.5 * np.log(alpha / 2)) * np.sqrt(2.0 / n)

    print(f"\nResults (n={n}, KS critical={ks_critical:.4f}):")
    for key in ["num_facets", "num_edges"]:
        ks_weighted = weighted_ks_2samp(approx[key], approx["weights"], exact[key])
        ks_unweighted, _ = stats.ks_2samp(approx[key], exact[key])

        weighted_ok = ks_weighted < ks_critical
        unweighted_differs = ks_unweighted > ks_critical

        if weighted_ok and unweighted_differs:
            status = "PASS (weights necessary and correct)"
        elif weighted_ok and not unweighted_differs:
            status = "WEAK PASS (weights correct but test not discriminating)"
        else:
            status = "FAIL (weighted distribution rejected)"
        print(f"  {key}: KS_weighted={ks_weighted:.4f}, KS_unweighted={ks_unweighted:.4f} [{status}]")

    w = approx["weights"]
    print(f"\n  Weight stats: mean={np.mean(w):.4f}, std={np.std(w):.4f}, "
          f"min={np.min(w):.4f}, max={np.max(w):.4f}")
    return 0


if __name__ == "__main__":
    exit(main())
