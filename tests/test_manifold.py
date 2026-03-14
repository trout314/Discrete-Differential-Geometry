"""Tests for Manifold Python bindings."""

import os
import tempfile

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

    def test_from_facets_dim2(self):
        m = Manifold(2, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert m.num_facets == 4

    def test_from_facets_dim3(self):
        # Standard 3-sphere = boundary of 4-simplex
        m = Manifold(3, [
            [0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4],
            [0, 2, 3, 4], [1, 2, 3, 4],
        ])
        assert m.num_facets == 5
        assert list(m.f_vector) == [5, 10, 10, 5]

    def test_from_facets_dim4(self):
        # Standard 4-sphere = boundary of 5-simplex
        m = Manifold(4, [
            [0, 1, 2, 3, 4], [0, 1, 2, 3, 5], [0, 1, 2, 4, 5],
            [0, 1, 3, 4, 5], [0, 2, 3, 4, 5], [1, 2, 3, 4, 5],
        ])
        assert m.num_facets == 6
        assert list(m.f_vector) == [6, 15, 20, 15, 6]


class TestQueries:
    @pytest.fixture
    def sphere2(self):
        return Manifold.standard_sphere(2)

    @pytest.fixture
    def sphere3(self):
        return Manifold.standard_sphere(3)

    @pytest.fixture
    def sphere4(self):
        return Manifold.standard_sphere(4)

    # -- f-vector --

    def test_f_vector_2sphere(self, sphere2):
        assert list(sphere2.f_vector) == [4, 6, 4]

    def test_f_vector_3sphere(self, sphere3):
        assert list(sphere3.f_vector) == [5, 10, 10, 5]

    def test_f_vector_4sphere(self, sphere4):
        # Boundary of 5-simplex: binom(6, k+1) for k=0..4
        assert list(sphere4.f_vector) == [6, 15, 20, 15, 6]

    # -- Euler characteristic --

    def test_euler_characteristic_2sphere(self, sphere2):
        assert sphere2.euler_characteristic == 2

    def test_euler_characteristic_3sphere(self, sphere3):
        assert sphere3.euler_characteristic == 0  # odd-dimensional sphere

    def test_euler_characteristic_4sphere(self, sphere4):
        assert sphere4.euler_characteristic == 2  # even-dimensional sphere

    # -- Orientability --

    def test_is_orientable_2sphere(self, sphere2):
        assert sphere2.is_orientable

    def test_is_orientable_3sphere(self, sphere3):
        assert sphere3.is_orientable

    def test_is_orientable_4sphere(self, sphere4):
        assert sphere4.is_orientable


class TestDataAccess:
    # -- Dim 2 --

    def test_facets_dim2(self):
        m = Manifold.standard_sphere(2)
        facets = m.facets()
        assert facets.shape == (4, 3)

    def test_simplices_dim2(self):
        m = Manifold.standard_sphere(2)
        edges = m.simplices(1)
        assert edges.shape == (6, 2)

    def test_degree_dim2(self):
        m = Manifold.standard_sphere(2)
        # In standard 2-sphere (tetrahedron boundary), each vertex has degree 3
        for v in range(4):
            assert m.degree([v]) == 3

    def test_mean_degree_dim2(self):
        m = Manifold.standard_sphere(2)
        assert m.mean_degree(0) == 3.0

    # -- Dim 3 --

    def test_facets_dim3(self):
        m = Manifold.standard_sphere(3)
        facets = m.facets()
        assert facets.shape == (5, 4)

    def test_simplices_dim3_vertices(self):
        m = Manifold.standard_sphere(3)
        verts = m.simplices(0)
        assert verts.shape == (5, 1)

    def test_simplices_dim3_edges(self):
        m = Manifold.standard_sphere(3)
        edges = m.simplices(1)
        assert edges.shape == (10, 2)

    def test_simplices_dim3_triangles(self):
        m = Manifold.standard_sphere(3)
        tris = m.simplices(2)
        assert tris.shape == (10, 3)

    def test_degree_dim3(self):
        m = Manifold.standard_sphere(3)
        # In standard 3-sphere (boundary of 4-simplex), each vertex has degree 4
        for v in range(5):
            assert m.degree([v]) == 4
        # Each edge has degree 3
        edges = m.simplices(1)
        for edge in edges:
            assert m.degree(edge) == 3

    def test_mean_degree_dim3(self):
        m = Manifold.standard_sphere(3)
        # Each vertex is in 4 tetrahedra: mean_degree(0) = 4
        assert m.mean_degree(0) == 4.0
        # Each edge is in 3 tetrahedra: mean_degree(1) = 3
        assert m.mean_degree(1) == 3.0
        # Each triangle is in 2 tetrahedra: mean_degree(2) = 2
        assert m.mean_degree(2) == 2.0

    # -- Dim 4 --

    def test_facets_dim4(self):
        m = Manifold.standard_sphere(4)
        facets = m.facets()
        assert facets.shape == (6, 5)

    def test_simplices_dim4_edges(self):
        m = Manifold.standard_sphere(4)
        edges = m.simplices(1)
        assert edges.shape == (15, 2)

    def test_simplices_dim4_triangles(self):
        m = Manifold.standard_sphere(4)
        tris = m.simplices(2)
        assert tris.shape == (20, 3)

    def test_simplices_dim4_tetrahedra(self):
        m = Manifold.standard_sphere(4)
        tets = m.simplices(3)
        assert tets.shape == (15, 4)

    def test_degree_dim4(self):
        m = Manifold.standard_sphere(4)
        # In standard 4-sphere (boundary of 5-simplex), each vertex has degree 5
        for v in range(6):
            assert m.degree([v]) == 5
        # Each edge has degree 4
        edges = m.simplices(1)
        for edge in edges:
            assert m.degree(edge) == 4

    def test_mean_degree_dim4(self):
        m = Manifold.standard_sphere(4)
        # In boundary of 5-simplex: mean_degree(d) = dim - d + 1
        assert m.mean_degree(0) == 5.0
        assert m.mean_degree(1) == 4.0
        assert m.mean_degree(2) == 3.0
        assert m.mean_degree(3) == 2.0


class TestMoves:
    def test_do_move_dim2(self):
        m = Manifold.standard_sphere(2)
        initial_facets = m.num_facets
        m.do_move()
        # On the standard 2-sphere there are no valid bistellar moves,
        # so the only option is a 1->3 stellar subdivision
        assert m.num_facets == initial_facets + 2

    def test_do_move_dim3(self):
        m = Manifold.standard_sphere(3)
        initial_facets = m.num_facets
        m.do_move()
        # Same logic: standard 3-sphere has no bistellar moves,
        # so a 1->4 stellar subdivision adds 3 facets
        assert m.num_facets == initial_facets + 3

    def test_do_move_dim4(self):
        m = Manifold.standard_sphere(4)
        initial_facets = m.num_facets
        m.do_move()
        # 1->5 stellar subdivision adds 4 facets
        assert m.num_facets == initial_facets + 4

    def test_multiple_moves_dim3(self):
        m = Manifold.standard_sphere(3)
        for _ in range(5):
            m.do_move()
        # After several moves, manifold should still have valid f-vector
        fv = m.f_vector
        assert len(fv) == 4
        assert all(v > 0 for v in fv)
        # Euler characteristic of a 3-manifold is always 0
        assert m.euler_characteristic == 0


class TestIO:
    def test_save_load_roundtrip_dim2(self):
        m = Manifold.standard_sphere(2)
        with tempfile.NamedTemporaryFile(suffix=".mfd", delete=False) as f:
            path = f.name
        try:
            m.save(path)
            loaded = Manifold.load(path, 2)
            assert loaded.num_facets == m.num_facets
            assert list(loaded.f_vector) == list(m.f_vector)
        finally:
            os.unlink(path)

    def test_save_load_roundtrip_dim3(self):
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

    def test_save_load_roundtrip_dim4(self):
        m = Manifold.standard_sphere(4)
        with tempfile.NamedTemporaryFile(suffix=".mfd", delete=False) as f:
            path = f.name
        try:
            m.save(path)
            loaded = Manifold.load(path, 4)
            assert loaded.num_facets == m.num_facets
            assert list(loaded.f_vector) == list(m.f_vector)
        finally:
            os.unlink(path)

    def test_save_load_after_moves_dim3(self):
        m = Manifold.standard_sphere(3)
        m.do_move()
        m.do_move()
        with tempfile.NamedTemporaryFile(suffix=".mfd", delete=False) as f:
            path = f.name
        try:
            m.save(path)
            loaded = Manifold.load(path, 3)
            assert loaded.num_facets == m.num_facets
            assert list(loaded.f_vector) == list(m.f_vector)
            assert loaded.euler_characteristic == m.euler_characteristic
        finally:
            os.unlink(path)

    def test_copy_dim3(self):
        m = Manifold.standard_sphere(3)
        m2 = m.copy()
        assert m2.num_facets == m.num_facets
        m2.do_move()
        # Copy should be independent — original unchanged
        assert m.num_facets == 5
        assert m2.num_facets != 5

    def test_copy_dim4(self):
        m = Manifold.standard_sphere(4)
        m2 = m.copy()
        assert m2.num_facets == m.num_facets
        m2.do_move()
        assert m.num_facets == 6
        assert m2.num_facets != 6

    def test_save_edge_graph_dim3(self):
        m = Manifold.standard_sphere(3)
        with tempfile.NamedTemporaryFile(suffix=".edge_graph", delete=False) as f:
            path = f.name
        try:
            m.save_edge_graph(path)
            with open(path) as f:
                lines = f.readlines()
            # 10 edges in standard 3-sphere
            assert len(lines) == 10
        finally:
            os.unlink(path)

    def test_save_dual_graph_dim3(self):
        m = Manifold.standard_sphere(3)
        with tempfile.NamedTemporaryFile(suffix=".dual_graph", delete=False) as f:
            path = f.name
        try:
            m.save_dual_graph(path)
            with open(path) as f:
                lines = f.readlines()
            # 5 tetrahedra, each shares 4 ridges, but each edge counted once: 10 edges
            assert len(lines) == 10
        finally:
            os.unlink(path)


class TestConversion:
    def test_to_simplicial_complex_dim2(self):
        m = Manifold.standard_sphere(2)
        sc = m.to_simplicial_complex()
        assert sc.num_facets == 4
        assert sc.is_2_sphere()

    def test_to_simplicial_complex_dim3(self):
        m = Manifold.standard_sphere(3)
        sc = m.to_simplicial_complex()
        assert sc.num_facets == 5
        assert sc.euler_characteristic == 0

    def test_to_simplicial_complex_dim4(self):
        m = Manifold.standard_sphere(4)
        sc = m.to_simplicial_complex()
        assert sc.num_facets == 6
        assert sc.euler_characteristic == 2
