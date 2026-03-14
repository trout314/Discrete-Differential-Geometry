"""Explore how importance weights behave at different manifold sizes.

Runs exact and approximate samplers at num_facets_target = 20 and 50,
prints KS statistics and weight summary stats for comparison.
"""

import numpy as np
from scipy import stats

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams


def collect_samples(manifold, params, *, num_samples, burn_in_moves, spacing_moves, exact):
    """Run a sampler and collect observable snapshots."""
    sampler = ManifoldSampler(manifold, params)
    sampler.run(moves=burn_in_moves, exact=exact)

    num_tets = []
    num_edges = []
    vertex_deg_var = []
    weights = []

    for i in range(num_samples):
        sampler.run(moves=spacing_moves, exact=exact)
        m = sampler.manifold
        fv = m.f_vector
        num_tets.append(fv[3])
        num_edges.append(fv[1])

        verts = m.simplices(0)
        degs = np.array([m.degree(v) for v in verts])
        vertex_deg_var.append(np.var(degs))

        weights.append(m.importance_weight())

    return {
        "num_tets": np.array(num_tets, dtype=float),
        "num_edges": np.array(num_edges, dtype=float),
        "vertex_deg_var": np.array(vertex_deg_var),
        "weights": np.array(weights),
    }


def weighted_ecdf(values, weights):
    order = np.argsort(values)
    sorted_vals = values[order]
    sorted_weights = weights[order]
    cumw = np.cumsum(sorted_weights)
    cumw /= cumw[-1]
    return sorted_vals, cumw


def weighted_ks_2samp(values1, weights1, values2):
    sv1, cw1 = weighted_ecdf(values1, weights1)
    sv2, cw2 = weighted_ecdf(values2, np.ones(len(values2)))

    all_vals = np.sort(np.concatenate([sv1, sv2]))
    ecdf1 = np.searchsorted(sv1, all_vals, side="right")
    ecdf1 = np.where(ecdf1 > 0, cw1[ecdf1 - 1], 0.0)
    ecdf2 = np.searchsorted(sv2, all_vals, side="right")
    ecdf2 = np.where(ecdf2 > 0, cw2[ecdf2 - 1], 0.0)

    return np.max(np.abs(ecdf1 - ecdf2))


def run_experiment(num_facets_target):
    print(f"\n{'='*70}")
    print(f"  num_facets_target = {num_facets_target}")
    print(f"{'='*70}")

    params = SamplerParams(
        num_facets_target=num_facets_target,
        hinge_degree_target=4.5,
        num_facets_coef=0.5,
        num_hinges_coef=0.1,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.0,
    )

    seed = Manifold.standard_sphere(3)
    grower = ManifoldSampler(seed, params)
    grower.run(moves=3000, exact=False)
    start = grower.manifold
    print(f"  Starting manifold: {start.f_vector[3]} tets")

    num_samples = 1000
    burn_in = 1500
    spacing = 50

    print(f"  Collecting {num_samples} exact samples...")
    exact_data = collect_samples(
        start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=True,
    )

    print(f"  Collecting {num_samples} approximate samples...")
    approx_data = collect_samples(
        start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=False,
    )

    # Weight statistics
    w = approx_data["weights"]
    print(f"\n  Weight statistics:")
    print(f"    mean   = {np.mean(w):.4f}")
    print(f"    std    = {np.std(w):.4f}")
    print(f"    min    = {np.min(w):.4f}")
    print(f"    max    = {np.max(w):.4f}")
    print(f"    cv     = {np.std(w)/np.mean(w):.4f}")

    # KS results
    n = num_samples
    alpha = 0.01
    ks_critical = np.sqrt(-0.5 * np.log(alpha / 2)) * np.sqrt(2.0 / n)

    print(f"\n  {'Observable':<20} {'KS_weighted':>12} {'KS_unweighted':>14} {'critical':>10}")
    print(f"  {'-'*56}")

    for key in ["num_tets", "num_edges", "vertex_deg_var"]:
        exact_vals = exact_data[key]
        approx_vals = approx_data[key]

        ks_w = weighted_ks_2samp(approx_vals, w, exact_vals)
        ks_u, _ = stats.ks_2samp(approx_vals, exact_vals)

        better = "  <--" if ks_w < ks_u else ""
        print(f"  {key:<20} {ks_w:>12.4f} {ks_u:>14.4f} {ks_critical:>10.4f}{better}")


if __name__ == "__main__":
    run_experiment(20)
    run_experiment(50)
    print()
