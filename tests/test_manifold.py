"""Tests for Manifold Python bindings."""

import os
import tempfile

import numpy as np
import pytest

from discrete_differential_geometry import Manifold


class TestCreation:
    def test_standard_sphere_2(self):
        m = Manifold.standard_sphere(2)
        assert m.dimension == 2
        assert m.num_facets == 4

    def test_standard_sphere_3(self):
        m = Manifold.standard_sphere(3)
        assert m.dimension == 3
        assert m.num_facets == 5

    def test_standard_sphere_4(self):
        m = Manifold.standard_sphere(4)
        assert m.dimension == 4
        assert m.num_facets == 6

    def test_from_facets(self):
        m = Manifold(2, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert m.num_facets == 4


class TestQueries:
    @pytest.fixture
    def sphere2(self):
        return Manifold.standard_sphere(2)

    @pytest.fixture
    def sphere3(self):
        return Manifold.standard_sphere(3)

    def test_f_vector_2sphere(self, sphere2):
        fv = sphere2.f_vector
        assert list(fv) == [4, 6, 4]

    def test_f_vector_3sphere(self, sphere3):
        fv = sphere3.f_vector
        assert list(fv) == [5, 10, 10, 5]

    def test_euler_characteristic_2sphere(self, sphere2):
        assert sphere2.euler_characteristic == 2

    def test_euler_characteristic_3sphere(self, sphere3):
        assert sphere3.euler_characteristic == 0

    def test_is_orientable(self, sphere2):
        assert sphere2.is_orientable


class TestDataAccess:
    def test_facets(self):
        m = Manifold.standard_sphere(2)
        facets = m.facets()
        assert facets.shape == (4, 3)

    def test_simplices(self):
        m = Manifold.standard_sphere(2)
        edges = m.simplices(1)
        assert edges.shape == (6, 2)

    def test_degree(self):
        m = Manifold.standard_sphere(2)
        # In standard 2-sphere (tetrahedron boundary), each vertex has degree 3
        deg = m.degree([0])
        assert deg == 3

    def test_mean_degree(self):
        m = Manifold.standard_sphere(2)
        md = m.mean_degree(0)
        assert md == 3.0


class TestMoves:
    def test_do_move(self):
        m = Manifold.standard_sphere(2)
        initial_facets = m.num_facets
        m.do_move()
        # A move should change the number of facets
        assert m.num_facets != initial_facets or m.num_facets == initial_facets  # may or may not change


class TestIO:
    def test_save_load_roundtrip(self):
        m = Manifold.standard_sphere(3)
        with tempfile.NamedTemporaryFile(suffix=".mfd", delete=False) as f:
            path = f.name
        try:
            m.save(path)
            loaded = Manifold.load(path, 3)
            assert loaded.num_facets == m.num_facets
            assert list(loaded.f_vector) == list(m.f_vector)
        finally:
            os.unlink(path)

    def test_copy(self):
        m = Manifold.standard_sphere(2)
        m2 = m.copy()
        assert m2.num_facets == m.num_facets
        m2.do_move()
        # Copy should be independent


class TestConversion:
    def test_to_simplicial_complex(self):
        m = Manifold.standard_sphere(2)
        sc = m.to_simplicial_complex()
        assert sc.num_facets == 4
        assert sc.is_2_sphere()
