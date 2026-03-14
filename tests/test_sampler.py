"""Tests for ManifoldSampler Python bindings."""

import pytest

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams


class TestSampler:
    def test_create_and_run(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(
            num_facets_target=20,
            hinge_degree_target=4.5,
            num_facets_coef=0.1,
            num_hinges_coef=0.05,
            hinge_degree_variance_coef=0.2,
            codim3_degree_variance_coef=0.1,
        )
        sampler = ManifoldSampler(m, params)
        accepted = sampler.run(moves=100)
        assert accepted >= 0

        # Check manifold changed
        result = sampler.manifold
        assert result.dimension == 3

    def test_run_sweeps(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=10)
        sampler = ManifoldSampler(m, params)
        accepted = sampler.run(sweeps=1.0)
        assert accepted >= 0

    def test_callback_early_stop(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=10)
        sampler = ManifoldSampler(m, params)

        calls = []

        def cb(done, total):
            calls.append(done)
            return done >= 2000  # stop after 2000 moves

        accepted = sampler.run(moves=100000, callback=cb)
        # Should have stopped early
        assert len(calls) > 0

    def test_requires_sweeps_or_moves(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams()
        sampler = ManifoldSampler(m, params)
        with pytest.raises(ValueError):
            sampler.run()
        with pytest.raises(ValueError):
            sampler.run(sweeps=1.0, moves=100)
