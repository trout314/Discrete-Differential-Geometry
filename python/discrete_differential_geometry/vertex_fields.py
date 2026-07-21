"""Per-vertex scalar fields on a 3-manifold triangulation, for structure-factor /
hyperuniformity analysis.

Every field is a pure function ``facets -> ndarray(V,)`` in the ``np.unique(F)``
DENSE vertex order (0..V-1). That single labeling convention is shared with the
coordinate builder (:func:`cocycle.torus_positions`) and both estimator modules
(:mod:`structure_factor`, :mod:`graph_hyperuniformity`), so a field array, a
position array, and an adjacency all index the same vertices with no relabeling.

Fields (see :data:`FIELDS` for the CLI registry):

* ``n6``               -- #incident edges of degree >= 6 (cheap disclination-density proxy)
* ``curvature_charge`` -- Regge scalar-curvature charge q_R = 1/2 sum_e delta_e
* ``curvature_density``-- q_R / (D_v/4)  (per unit dual 3-volume)
* ``volume_charge``    -- D_v / 4
* ``defect_indicator`` -- 1 at illegal (impurity) vertices, else 0

The Regge deficit is delta_e = 2*pi - theta*deg(e) with theta = arccos(1/3) the
regular-tet dihedral angle; the flat edge degree is 2*pi/theta = 5.1043... (the
native TCP value), NOT 6. The naive combinatorial deficit sum_e (6 - deg e) is
identically 12 by the link sum rule (Euler on the S^2 link) and so is a dead,
constant field -- it is deliberately not offered here.
"""
import numpy as np

from .move_geometry import THETA_TET

_EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def edges_and_degrees(facets):
    """Combinatorics in dense labels: (eu (E,2) sorted u<w, ecnt (E,) edge
    degrees, deg (V,) vertex degrees, V). Shared preamble for every field."""
    F = np.asarray(facets, np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    T = inv.reshape(F.shape)
    V = len(lab)
    deg = np.bincount(T.ravel(), minlength=V)
    epairs = np.sort(np.vstack([T[:, [i, j]] for i, j in _EDGE_PAIRS]), axis=1)
    eu, ecnt = np.unique(epairs, axis=0, return_counts=True)
    return eu, ecnt, deg, V


def _accumulate(eu, w, V):
    """Sum a per-edge weight w onto both endpoints -> per-vertex field."""
    q = np.zeros(V)
    np.add.at(q, eu[:, 0], w)
    np.add.at(q, eu[:, 1], w)
    return q


def n6(facets):
    """#incident edges with degree >= 6 (disclination-line-density proxy)."""
    eu, ecnt, deg, V = edges_and_degrees(facets)
    return _accumulate(eu, (ecnt >= 6).astype(float), V)


def curvature_charge(facets):
    """Regge scalar-curvature CHARGE q_R(v) = 1/2 sum_{e ni v} delta_e,
    delta_e = 2*pi - theta*deg(e). The (integrated) Hamiltonian-constraint
    quantity: sum_v q_R = sum_e delta_e."""
    eu, ecnt, deg, V = edges_and_degrees(facets)
    return _accumulate(eu, (2 * np.pi - THETA_TET * ecnt) / 2, V)


def curvature_density(facets):
    """Regge scalar-curvature DENSITY R(v) = q_R(v) / (D_v/4) (deficit per unit
    dual 3-volume)."""
    eu, ecnt, deg, V = edges_and_degrees(facets)
    qR = _accumulate(eu, (2 * np.pi - THETA_TET * ecnt) / 2, V)
    return qR / (deg / 4.0)


def volume_charge(facets):
    """Volume charge D_v/4 (vertex degree / 4; global total pinned by the facet
    pin)."""
    eu, ecnt, deg, V = edges_and_degrees(facets)
    return deg / 4.0


def defect_indicator(facets):
    """1 at illegal (impurity) vertices -- those with an incident edge of degree
    outside {5, 6} -- else 0. Isolates the defect arrangement (the crystalline
    bulk is invisible), unlike n6/curvature which the periodic crystal dominates."""
    eu, ecnt, deg, V = edges_and_degrees(facets)
    bad = (ecnt < 5) | (ecnt > 6)
    q = np.zeros(V)
    np.add.at(q, eu[bad, 0], 1.0)
    np.add.at(q, eu[bad, 1], 1.0)
    return (q > 0).astype(float)


FIELDS = {
    "n6": n6,
    "curvature_charge": curvature_charge,
    "curvature_density": curvature_density,
    "volume_charge": volume_charge,
    "defect_indicator": defect_indicator,
}
