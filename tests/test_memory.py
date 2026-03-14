"""Quick test to verify D GC memory doesn't grow unboundedly.

Runs a sampling loop similar to test_importance_weight.py and checks
that RSS stays bounded.
"""

import os
import resource

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams, gc_collect, gc_minimize


def get_rss_mb():
    """Get current RSS in MB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def main():
    params = SamplerParams(
        num_facets_target=30,
        hinge_degree_target=4.5,
        num_facets_coef=0.5,
        num_hinges_coef=0.1,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.0,
    )

    seed = Manifold.standard_sphere(3)
    sampler = ManifoldSampler(seed, params)
    sampler.run(moves=2000, exact=False)

    rss_start = get_rss_mb()
    print(f"Start RSS: {rss_start:.1f} MB")

    num_samples = 2000
    for i in range(num_samples):
        sampler.run(moves=50, exact=False)

        # Query the sampler (same as test_importance_weight)
        fv = sampler.f_vector
        verts = sampler.simplices(0)
        for v in verts:
            sampler.degree(v)
        w = sampler.importance_weight()

        # Periodic GC + status
        if (i + 1) % 200 == 0:
            gc_collect()
            gc_minimize()
            rss_now = get_rss_mb()
            print(f"  Sample {i+1}/{num_samples}: RSS={rss_now:.1f} MB (delta={rss_now - rss_start:+.1f} MB)")

    rss_end = get_rss_mb()
    print(f"End RSS: {rss_end:.1f} MB (delta={rss_end - rss_start:+.1f} MB)")

    # RSS should not have grown by more than 50 MB
    growth = rss_end - rss_start
    if growth > 50:
        print(f"FAIL: RSS grew by {growth:.1f} MB — likely a memory leak")
        return 1
    else:
        print(f"OK: RSS growth {growth:.1f} MB is within acceptable bounds")
        return 0


if __name__ == "__main__":
    exit(main())
