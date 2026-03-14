"""Tests for SimplicialComplex Python bindings."""

import os
import tempfile

import numpy as np
import pytest

from discrete_differential_geometry import SimplicialComplex, join


class TestCreation:
    def test_empty(self):
        sc = SimplicialComplex()
        assert sc.num_facets == 0

    def test_from_triangles(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert sc.num_facets == 4

    def test_from_edges(self):
        sc = SimplicialComplex([[0, 1], [1, 2], [0, 2]])
        assert sc.num_facets == 3


class TestQueries:
    @pytest.fixture
    def tetrahedron(self):
        return SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])

    def test_f_vector(self, tetrahedron):
        fv = tetrahedron.f_vector
        assert list(fv) == [4, 6, 4]

    def test_euler_characteristic(self, tetrahedron):
        assert tetrahedron.euler_characteristic == 2

    def test_is_connected(self, tetrahedron):
        assert tetrahedron.is_connected

    def test_is_pure(self, tetrahedron):
        assert tetrahedron.is_pure

    def test_contains(self, tetrahedron):
        assert tetrahedron.contains([0, 1])
        assert tetrahedron.contains([0, 1, 2])
        assert not tetrahedron.contains([0, 1, 2, 3])

    def test_contains_facet(self, tetrahedron):
        assert tetrahedron.contains_facet([0, 1, 2])
        assert not tetrahedron.contains_facet([0, 1])


class TestMutation:
    def test_insert_remove(self):
        sc = SimplicialComplex()
        sc.insert_facet([1, 2])
        assert sc.num_facets == 1
        assert sc.contains_facet([1, 2])
        sc.remove_facet([1, 2])
        assert sc.num_facets == 0


class TestDataAccess:
    def test_simplices(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        verts = sc.simplices(0)
        assert verts.shape == (4, 1)
        edges = sc.simplices(1)
        assert edges.shape == (6, 2)

    def test_star(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        star = sc.star([0])
        assert len(star) == 3  # 3 triangles contain vertex 0

    def test_link(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        lk = sc.link([0])
        assert len(lk) == 3  # link of vertex 0 has 3 edges


class TestTopology:
    def test_connected_components(self):
        sc = SimplicialComplex([[0, 1], [1, 2], [3, 4]])
        comps = sc.connected_components()
        assert len(comps) == 2

    def test_is_circle(self):
        sc = SimplicialComplex([[0, 1], [1, 2], [2, 3], [0, 3]])
        assert sc.is_circle()

    def test_is_2_sphere(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert sc.is_2_sphere()

    def test_is_pure_of_dim(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert sc.is_pure_of_dim(2)
        assert not sc.is_pure_of_dim(1)

    def test_join(self):
        sc1 = SimplicialComplex([[0], [1]])
        sc2 = SimplicialComplex([[2], [3]])
        result = join(sc1, sc2)
        assert result.num_facets == 4


class TestIO:
    def test_save_load_roundtrip(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        with tempfile.NamedTemporaryFile(suffix=".sc", delete=False) as f:
            path = f.name
        try:
            sc.save(path)
            loaded = SimplicialComplex.load(path)
            assert loaded.num_facets == sc.num_facets
            assert list(loaded.f_vector) == list(sc.f_vector)
        finally:
            os.unlink(path)

    def test_copy(self):
        sc = SimplicialComplex([[0, 1, 2], [0, 1, 3]])
        sc2 = sc.copy()
        assert sc2.num_facets == sc.num_facets
        sc2.insert_facet([4, 5, 6])
        assert sc2.num_facets != sc.num_facets
