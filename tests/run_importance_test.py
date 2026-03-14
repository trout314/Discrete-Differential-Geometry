"""Run importance weight validation with configurable sample count."""

import sys
import resource
import numpy as np
from scipy import stats

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams, gc_collect


def get_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def collect_samples(manifold, params, *, num_samples, burn_in_moves, spacing_moves, exact):
    sampler = ManifoldSampler(manifold, params)
    sampler.run(moves=burn_in_moves, exact=exact)

    num_tets = []
    num_edges = []
    weights = []

    for i in range(num_samples):
        sampler.run(moves=spacing_moves, exact=exact)
        fv = sampler.f_vector
        num_tets.append(fv[3])
        num_edges.append(fv[1])
        weights.append(sampler.importance_weight())

        if (i + 1) % 200 == 0:
            print(f"    {i+1}/{num_samples} samples, RSS={get_rss_mb():.0f} MB")

    return {
        "num_tets": np.array(num_tets, dtype=float),
        "num_edges": np.array(num_edges, dtype=float),
        "weights": np.array(weights),
    }


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


def main():
    num_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

    params = SamplerParams(
        num_facets_target=30, hinge_degree_target=4.5,
        num_facets_coef=0.5, num_hinges_coef=0.1,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    )

    seed = Manifold.standard_sphere(3)
    grower = ManifoldSampler(seed, params)
    grower.run(moves=2000, exact=False)
    start = grower.manifold

    burn_in = 1000
    spacing = 50

    print(f"Collecting {num_samples} exact samples...")
    exact_data = collect_samples(start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=True)

    print(f"Collecting {num_samples} approximate samples...")
    approx_data = collect_samples(start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=False)

    # Check results
    alpha = 0.01
    n = num_samples
    ks_critical = np.sqrt(-0.5 * np.log(alpha / 2)) * np.sqrt(2.0 / n)

    print(f"\nResults (n={n}, KS critical={ks_critical:.4f}):")
    for key in ["num_tets", "num_edges"]:
        ks_weighted = weighted_ks_2samp(approx_data[key], approx_data["weights"], exact_data[key])
        ks_unweighted, _ = stats.ks_2samp(approx_data[key], exact_data[key])
        status = "PASS" if ks_weighted < ks_critical else "FAIL"
        print(f"  {key}: KS_weighted={ks_weighted:.4f}, KS_unweighted={ks_unweighted:.4f} [{status}]")

    w = approx_data["weights"]
    print(f"\n  Weight stats: mean={np.mean(w):.4f}, std={np.std(w):.4f}, "
          f"min={np.min(w):.4f}, max={np.max(w):.4f}")
    print(f"  Final RSS: {get_rss_mb():.0f} MB")


if __name__ == "__main__":
    main()
