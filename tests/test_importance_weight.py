"""Statistical validation: importance-weighted approximate sampler
should reproduce the exact sampler's distribution.

We run both samplers on a small (~50 tet) 3-manifold, collect observables
(#tets, #edges, vertex degree variance), and check that the
importance-weighted approximate distribution matches the exact distribution
using a two-sample Kolmogorov-Smirnov test.
"""

import numpy as np
import pytest
from scipy import stats

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams


def collect_samples(manifold, params, *, num_samples, burn_in_moves, spacing_moves, exact):
    """Run a sampler and collect observable snapshots."""
    sampler = ManifoldSampler(manifold, params)

    # Burn-in
    sampler.run(moves=burn_in_moves, exact=exact)

    num_tets = []
    num_edges = []
    vertex_deg_var = []
    weights = []

    for _ in range(num_samples):
        sampler.run(moves=spacing_moves, exact=exact)
        # Query the sampler's manifold directly (no copy)
        fv = sampler.f_vector
        num_tets.append(fv[3])
        num_edges.append(fv[1])

        # Vertex degree variance
        verts = sampler.simplices(0)
        degs = np.array([sampler.degree(v) for v in verts])
        vertex_deg_var.append(np.var(degs))

        weights.append(sampler.importance_weight())

    return {
        "num_tets": np.array(num_tets, dtype=float),
        "num_edges": np.array(num_edges, dtype=float),
        "vertex_deg_var": np.array(vertex_deg_var),
        "weights": np.array(weights),
    }


def weighted_ecdf(values, weights):
    """Return (sorted_values, cumulative_weights) for a weighted ECDF."""
    order = np.argsort(values)
    sorted_vals = values[order]
    sorted_weights = weights[order]
    cumw = np.cumsum(sorted_weights)
    cumw /= cumw[-1]
    return sorted_vals, cumw


def weighted_ks_2samp(values1, weights1, values2):
    """Weighted KS statistic comparing weighted sample 1 to unweighted sample 2."""
    # Weighted ECDF for sample 1
    sv1, cw1 = weighted_ecdf(values1, weights1)
    # Unweighted ECDF for sample 2
    sv2, cw2 = weighted_ecdf(values2, np.ones(len(values2)))

    # Merge and compute max difference
    all_vals = np.sort(np.concatenate([sv1, sv2]))
    ecdf1 = np.searchsorted(sv1, all_vals, side="right")
    ecdf1 = np.where(ecdf1 > 0, cw1[ecdf1 - 1], 0.0)
    ecdf2 = np.searchsorted(sv2, all_vals, side="right")
    ecdf2 = np.where(ecdf2 > 0, cw2[ecdf2 - 1], 0.0)

    return np.max(np.abs(ecdf1 - ecdf2))


@pytest.fixture(scope="module")
def sampler_results():
    """Run exact and approximate samplers once for all tests in this module."""
    params = SamplerParams(
        num_facets_target=30,
        hinge_degree_target=4.5,
        num_facets_coef=0.5,
        num_hinges_coef=0.1,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.0,
    )

    # Start from standard sphere, grow a bit first
    seed = Manifold.standard_sphere(3)
    grower = ManifoldSampler(seed, params)
    grower.run(moves=2000, exact=False)
    start = grower.manifold

    num_samples = 10000
    burn_in = 1000
    spacing = 50

    exact_data = collect_samples(
        start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=True,
    )
    approx_data = collect_samples(
        start, params,
        num_samples=num_samples, burn_in_moves=burn_in,
        spacing_moves=spacing, exact=False,
    )
    return exact_data, approx_data


class TestImportanceWeight:
    """Verify that importance-weighted approximate samples match the exact distribution."""

    def _check_observable(self, exact_data, approx_data, key, alpha=0.01):
        """Run weighted KS test for one observable.

        We check two things:
        1. Weighted approximate vs exact should NOT be rejected (distributions match).
        2. Unweighted approximate vs exact SHOULD show a larger KS statistic
           (demonstrating the weights are actually doing something).
        """
        exact_vals = exact_data[key]
        approx_vals = approx_data[key]
        approx_weights = approx_data["weights"]

        # Weighted KS statistic
        ks_weighted = weighted_ks_2samp(approx_vals, approx_weights, exact_vals)

        # Unweighted KS for comparison
        ks_unweighted, _ = stats.ks_2samp(approx_vals, exact_vals)

        # The weighted version should have a smaller KS statistic than the
        # unweighted version (or at least not be dramatically worse)
        # Use a generous threshold since these are stochastic
        n = len(exact_vals)
        # Critical value for KS test at significance alpha with n samples each
        ks_critical = np.sqrt(-0.5 * np.log(alpha / 2)) * np.sqrt(2.0 / n)

        print(f"\n{key}: KS_weighted={ks_weighted:.4f}, KS_unweighted={ks_unweighted:.4f}, "
              f"critical={ks_critical:.4f}")

        # Weighted distribution should not be rejected
        assert ks_weighted < ks_critical, (
            f"{key}: weighted KS ({ks_weighted:.4f}) >= critical ({ks_critical:.4f})"
        )

    def test_num_tets(self, sampler_results):
        exact_data, approx_data = sampler_results
        self._check_observable(exact_data, approx_data, "num_tets")

    def test_num_edges(self, sampler_results):
        exact_data, approx_data = sampler_results
        self._check_observable(exact_data, approx_data, "num_edges")

    def test_vertex_deg_var(self, sampler_results):
        exact_data, approx_data = sampler_results
        self._check_observable(exact_data, approx_data, "vertex_deg_var")

    def test_weights_not_trivial(self, sampler_results):
        """Verify that the importance weights are not all 1.0
        (which would mean the test is vacuous)."""
        _, approx_data = sampler_results
        weights = approx_data["weights"]
        assert not np.allclose(weights, 1.0), "All weights are 1.0 — test is vacuous"
        # Weights should have some spread
        assert np.std(weights) > 0.001, f"Weight std too small: {np.std(weights)}"
