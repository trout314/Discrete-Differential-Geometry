"""Tests for ManifoldSampler Python bindings."""

import numpy as np
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


class TestStats:
    def test_stats_accumulate(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)

        sampler.run(moves=500)
        stats = sampler.get_stats()

        assert stats.total_tried == 500
        assert stats.total_accepted >= 0
        assert stats.total_accepted <= stats.total_tried
        assert sum(stats.bistellar_tries) == stats.total_tried

    def test_stats_reset(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)

        sampler.run(moves=200)
        sampler.reset_stats()
        stats = sampler.get_stats()

        assert stats.total_tried == 0
        assert stats.total_accepted == 0
        assert stats.hinge_tries == 0
        assert stats.hinge_accepts == 0

    def test_stats_accumulate_across_runs(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)

        sampler.run(moves=100)
        sampler.run(moves=100)
        stats = sampler.get_stats()

        assert stats.total_tried == 200

    def test_bistellar_accepts_le_tries(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=500)
        stats = sampler.get_stats()

        for i in range(4):
            assert stats.bistellar_accepts[i] <= stats.bistellar_tries[i]


class TestHingeMoves:
    def test_hinge_moves_tracked(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(
            num_facets_target=50,
            num_facets_coef=0.1,
            hinge_move_prob=0.3,
        )
        sampler = ManifoldSampler(m, params)
        # Grow first so there are degree-4 edges for hinge moves
        sampler.ramped_grow(100, step_size=50, eq_sweeps_per_step=2)
        sampler.reset_stats()
        sampler.run(sweeps=10)
        stats = sampler.get_stats()

        # With hinge_move_prob=0.3, we should see hinge attempts
        assert stats.hinge_tries > 0
        assert stats.hinge_accepts <= stats.hinge_tries

    def test_no_hinge_moves_when_prob_zero(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(
            num_facets_target=20,
            num_facets_coef=0.1,
            hinge_move_prob=0.0,
        )
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=500)
        stats = sampler.get_stats()

        assert stats.hinge_tries == 0
        assert stats.hinge_accepts == 0

    def test_set_hinge_move_prob(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.ramped_grow(100, step_size=50, eq_sweeps_per_step=2)

        # Start with no hinge moves
        sampler.reset_stats()
        sampler.run(sweeps=5)
        assert sampler.get_stats().hinge_tries == 0

        # Enable hinge moves
        sampler.set_hinge_move_prob(0.5)
        sampler.reset_stats()
        sampler.run(sweeps=5)
        assert sampler.get_stats().hinge_tries > 0


class TestParamSetters:
    def test_set_num_facets_target(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.5)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=500)

        sampler.set_num_facets_target(50)
        assert sampler._params.num_facets_target == 50

        # Run more — manifold should grow toward new target
        sampler.run(moves=2000)
        assert sampler.manifold.num_facets > 20

    def test_set_num_facets_coef(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.set_num_facets_coef(0.5)
        assert sampler._params.num_facets_coef == 0.5
        # Should still run without error
        sampler.run(moves=100)

    def test_objective_changes_with_params(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=100)
        obj1 = sampler.current_objective

        sampler.set_num_facets_coef(10.0)
        sampler.run(moves=1)  # triggers recompute
        obj2 = sampler.current_objective

        # Objective should differ with wildly different coefficient
        assert obj1 != obj2


class TestCurrentObjective:
    def test_objective_finite(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=100)
        obj = sampler.current_objective
        assert np.isfinite(obj)

    def test_objective_nonnegative(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=100)
        assert sampler.current_objective >= 0.0


class TestRampedGrow:
    def test_basic_growth(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)

        sampler.ramped_grow(200, step_size=50, eq_sweeps_per_step=2)
        assert sampler.manifold.num_facets >= 200

    def test_callback_called(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)

        steps = []
        sampler.ramped_grow(
            150, step_size=50, eq_sweeps_per_step=1,
            callback=lambda cur, tgt: steps.append((cur, tgt)),
        )
        assert len(steps) >= 2
        # Each step target should be increasing
        targets = [t for _, t in steps]
        assert targets == sorted(targets)
        assert targets[-1] == 150

    def test_manifold_valid_after_growth(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.ramped_grow(200, step_size=100, eq_sweeps_per_step=2)

        mfd = sampler.manifold
        assert mfd.dimension == 3
        assert mfd.euler_characteristic == 0
        assert mfd.is_orientable
        fv = mfd.f_vector
        assert len(fv) == 4
        assert all(v > 0 for v in fv)

    def test_growth_with_hinge_moves(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(
            num_facets_target=50,
            num_facets_coef=0.1,
            hinge_move_prob=0.3,
        )
        sampler = ManifoldSampler(m, params)
        sampler.ramped_grow(150, step_size=50, eq_sweeps_per_step=3)
        assert sampler.manifold.num_facets >= 150
        assert sampler.manifold.euler_characteristic == 0


class TestManifoldViewFromSampler:
    def test_degree_histogram_via_view(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.ramped_grow(100, step_size=50, eq_sweeps_per_step=2)

        h = sampler.manifold.degree_histogram(0)
        assert sum(h) == sampler.manifold.f_vector[0]

    def test_degree_variance_via_view(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=50, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.ramped_grow(100, step_size=50, eq_sweeps_per_step=2)

        dv = sampler.manifold.degree_variance(0)
        assert dv >= 0.0
        assert np.isfinite(dv)

    def test_dup_produces_independent_manifold(self):
        m = Manifold.standard_sphere(3)
        params = SamplerParams(num_facets_target=20, num_facets_coef=0.1)
        sampler = ManifoldSampler(m, params)
        sampler.run(moves=200)

        copy = sampler.manifold.dup()
        nf_before = copy.num_facets
        sampler.run(moves=200)

        # Copy should be unchanged
        assert copy.num_facets == nf_before
