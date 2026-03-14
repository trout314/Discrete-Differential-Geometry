"""ctypes loader and thin wrappers for the D shared library."""

import ctypes
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Library discovery
# ---------------------------------------------------------------------------

def _find_library() -> ctypes.CDLL:
    """Locate and load ddg_dlang.so (or .dylib / .dll)."""
    suffixes = {
        "linux": ".so",
        "darwin": ".dylib",
        "win32": ".dll",
    }
    ext = suffixes.get(sys.platform, ".so")
    name = f"ddg_dlang{ext}"

    search_dirs = [
        # Same directory as this Python file (for installed packages)
        Path(__file__).parent,
        # builddir/ relative to project root
        Path(__file__).parent.parent.parent / "builddir",
    ]

    for d in search_dirs:
        candidate = d / name
        if candidate.exists():
            return ctypes.CDLL(str(candidate))

    raise OSError(
        f"Cannot find {name}. Searched: {[str(d) for d in search_dirs]}. "
        f"Build with: meson compile -C builddir ddg_dlang"
    )


_lib = _find_library()

# ---------------------------------------------------------------------------
# Error helper
# ---------------------------------------------------------------------------

_lib.ddg_last_error.argtypes = []
_lib.ddg_last_error.restype = ctypes.c_char_p


def _check_null(result, func, args):
    """ctypes errcheck: raise on NULL return."""
    if result is None:
        err = _lib.ddg_last_error()
        msg = err.decode() if err else "unknown error"
        raise RuntimeError(f"{func.__name__}: {msg}")
    return result


def _check_int(result, func, args):
    """ctypes errcheck: raise on negative return."""
    if result < 0:
        err = _lib.ddg_last_error()
        msg = err.decode() if err else "unknown error"
        raise RuntimeError(f"{func.__name__}: {msg}")
    return result


# ---------------------------------------------------------------------------
# Manifold lifecycle
# ---------------------------------------------------------------------------

_lib.ddg_manifold_standard_sphere.argtypes = [ctypes.c_int]
_lib.ddg_manifold_standard_sphere.restype = ctypes.c_void_p
_lib.ddg_manifold_standard_sphere.errcheck = _check_null

_lib.ddg_manifold_from_facets.argtypes = [
    ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_manifold_from_facets.restype = ctypes.c_void_p
_lib.ddg_manifold_from_facets.errcheck = _check_null

_lib.ddg_manifold_load.argtypes = [ctypes.c_char_p, ctypes.c_int]
_lib.ddg_manifold_load.restype = ctypes.c_void_p
_lib.ddg_manifold_load.errcheck = _check_null

_lib.ddg_manifold_copy.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_copy.restype = ctypes.c_void_p
_lib.ddg_manifold_copy.errcheck = _check_null

_lib.ddg_manifold_free.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_free.restype = None

# ---------------------------------------------------------------------------
# Manifold queries
# ---------------------------------------------------------------------------

_lib.ddg_manifold_dimension.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_dimension.restype = ctypes.c_int

_lib.ddg_manifold_num_facets.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_num_facets.restype = ctypes.c_long
_lib.ddg_manifold_num_facets.errcheck = _check_int

_lib.ddg_manifold_euler_characteristic.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_euler_characteristic.restype = ctypes.c_int

_lib.ddg_manifold_f_vector.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.c_int,
]
_lib.ddg_manifold_f_vector.restype = ctypes.c_int
_lib.ddg_manifold_f_vector.errcheck = _check_int

_lib.ddg_manifold_is_orientable.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_is_orientable.restype = ctypes.c_int
_lib.ddg_manifold_is_orientable.errcheck = _check_int

_lib.ddg_manifold_num_connected_components.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_num_connected_components.restype = ctypes.c_int
_lib.ddg_manifold_num_connected_components.errcheck = _check_int

# ---------------------------------------------------------------------------
# Manifold data
# ---------------------------------------------------------------------------

_lib.ddg_manifold_facets.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
_lib.ddg_manifold_facets.restype = ctypes.c_long
_lib.ddg_manifold_facets.errcheck = _check_int

_lib.ddg_manifold_simplices.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_manifold_simplices.restype = ctypes.c_long
_lib.ddg_manifold_simplices.errcheck = _check_int

_lib.ddg_manifold_degree.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_manifold_degree.restype = ctypes.c_long
_lib.ddg_manifold_degree.errcheck = _check_int

_lib.ddg_manifold_count_valid_moves.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_count_valid_moves.restype = ctypes.c_long
_lib.ddg_manifold_count_valid_moves.errcheck = _check_int

_lib.ddg_manifold_importance_weight.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_importance_weight.restype = ctypes.c_double

_lib.ddg_manifold_mean_degree.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_manifold_mean_degree.restype = ctypes.c_double

# ---------------------------------------------------------------------------
# Manifold I/O
# ---------------------------------------------------------------------------

_lib.ddg_manifold_save.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.ddg_manifold_save.restype = ctypes.c_int
_lib.ddg_manifold_save.errcheck = _check_int

_lib.ddg_manifold_save_edge_graph.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.ddg_manifold_save_edge_graph.restype = ctypes.c_int
_lib.ddg_manifold_save_edge_graph.errcheck = _check_int

_lib.ddg_manifold_save_dual_graph.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.ddg_manifold_save_dual_graph.restype = ctypes.c_int
_lib.ddg_manifold_save_dual_graph.errcheck = _check_int

# ---------------------------------------------------------------------------
# Manifold moves
# ---------------------------------------------------------------------------

_lib.ddg_manifold_do_pachner_move.argtypes = [ctypes.c_void_p]
_lib.ddg_manifold_do_pachner_move.restype = ctypes.c_int
_lib.ddg_manifold_do_pachner_move.errcheck = _check_int

# ---------------------------------------------------------------------------
# SimplicialComplex lifecycle
# ---------------------------------------------------------------------------

_lib.ddg_sc_create.argtypes = []
_lib.ddg_sc_create.restype = ctypes.c_void_p
_lib.ddg_sc_create.errcheck = _check_null

_lib.ddg_sc_from_facets.argtypes = [
    ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_int,
]
_lib.ddg_sc_from_facets.restype = ctypes.c_void_p
_lib.ddg_sc_from_facets.errcheck = _check_null

_lib.ddg_sc_load.argtypes = [ctypes.c_char_p]
_lib.ddg_sc_load.restype = ctypes.c_void_p
_lib.ddg_sc_load.errcheck = _check_null

_lib.ddg_sc_copy.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_copy.restype = ctypes.c_void_p
_lib.ddg_sc_copy.errcheck = _check_null

_lib.ddg_sc_free.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_free.restype = None

# ---------------------------------------------------------------------------
# SimplicialComplex queries
# ---------------------------------------------------------------------------

_lib.ddg_sc_num_facets.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_num_facets.restype = ctypes.c_long
_lib.ddg_sc_num_facets.errcheck = _check_int

_lib.ddg_sc_f_vector.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.c_int,
]
_lib.ddg_sc_f_vector.restype = ctypes.c_int
_lib.ddg_sc_f_vector.errcheck = _check_int

_lib.ddg_sc_contains.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_sc_contains.restype = ctypes.c_int
_lib.ddg_sc_contains.errcheck = _check_int

_lib.ddg_sc_contains_facet.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_sc_contains_facet.restype = ctypes.c_int
_lib.ddg_sc_contains_facet.errcheck = _check_int

# ---------------------------------------------------------------------------
# SimplicialComplex data
# ---------------------------------------------------------------------------

_lib.ddg_sc_facets.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sc_facets.restype = ctypes.c_long
_lib.ddg_sc_facets.errcheck = _check_int

_lib.ddg_sc_facets_of_dim.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sc_facets_of_dim.restype = ctypes.c_long
_lib.ddg_sc_facets_of_dim.errcheck = _check_int

_lib.ddg_sc_simplices.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sc_simplices.restype = ctypes.c_long
_lib.ddg_sc_simplices.errcheck = _check_int

_lib.ddg_sc_star.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sc_star.restype = ctypes.c_long
_lib.ddg_sc_star.errcheck = _check_int

_lib.ddg_sc_link.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sc_link.restype = ctypes.c_long
_lib.ddg_sc_link.errcheck = _check_int

# ---------------------------------------------------------------------------
# SimplicialComplex mutation
# ---------------------------------------------------------------------------

_lib.ddg_sc_insert_facet.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_sc_insert_facet.restype = ctypes.c_int
_lib.ddg_sc_insert_facet.errcheck = _check_int

_lib.ddg_sc_remove_facet.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_sc_remove_facet.restype = ctypes.c_int
_lib.ddg_sc_remove_facet.errcheck = _check_int

# ---------------------------------------------------------------------------
# SimplicialComplex I/O
# ---------------------------------------------------------------------------

_lib.ddg_sc_save.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.ddg_sc_save.restype = ctypes.c_int
_lib.ddg_sc_save.errcheck = _check_int

_lib.ddg_sc_save_edge_graph.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.ddg_sc_save_edge_graph.restype = ctypes.c_int
_lib.ddg_sc_save_edge_graph.errcheck = _check_int

# ---------------------------------------------------------------------------
# SimplicialComplex algorithms
# ---------------------------------------------------------------------------

_lib.ddg_sc_euler_characteristic.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_euler_characteristic.restype = ctypes.c_int

_lib.ddg_sc_is_connected.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_is_connected.restype = ctypes.c_int
_lib.ddg_sc_is_connected.errcheck = _check_int

_lib.ddg_sc_connected_components.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_void_p),
]
_lib.ddg_sc_connected_components.restype = ctypes.c_int
_lib.ddg_sc_connected_components.errcheck = _check_int

_lib.ddg_sc_is_pure.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_is_pure.restype = ctypes.c_int
_lib.ddg_sc_is_pure.errcheck = _check_int

_lib.ddg_sc_is_pure_of_dim.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sc_is_pure_of_dim.restype = ctypes.c_int
_lib.ddg_sc_is_pure_of_dim.errcheck = _check_int

_lib.ddg_sc_is_circle.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_is_circle.restype = ctypes.c_int
_lib.ddg_sc_is_circle.errcheck = _check_int

_lib.ddg_sc_is_2_sphere.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_is_2_sphere.restype = ctypes.c_int
_lib.ddg_sc_is_2_sphere.errcheck = _check_int

_lib.ddg_sc_is_2_torus.argtypes = [ctypes.c_void_p]
_lib.ddg_sc_is_2_torus.restype = ctypes.c_int
_lib.ddg_sc_is_2_torus.errcheck = _check_int

_lib.ddg_sc_is_orientable_surface_of_genus.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sc_is_orientable_surface_of_genus.restype = ctypes.c_int
_lib.ddg_sc_is_orientable_surface_of_genus.errcheck = _check_int

_lib.ddg_sc_join.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib.ddg_sc_join.restype = ctypes.c_void_p
_lib.ddg_sc_join.errcheck = _check_null

# ---------------------------------------------------------------------------
# Sampler
# ---------------------------------------------------------------------------

CALLBACK_TYPE = ctypes.CFUNCTYPE(
    ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_void_p,
)

_lib.ddg_sampler_create.argtypes = [
    ctypes.c_void_p,
    ctypes.c_int, ctypes.c_double,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double,
]
_lib.ddg_sampler_create.restype = ctypes.c_void_p
_lib.ddg_sampler_create.errcheck = _check_null

_lib.ddg_sampler_run.argtypes = [
    ctypes.c_void_p, ctypes.c_long, CALLBACK_TYPE, ctypes.c_void_p,
]
_lib.ddg_sampler_run.restype = ctypes.c_long
_lib.ddg_sampler_run.errcheck = _check_int

_lib.ddg_sampler_run_exact.argtypes = [
    ctypes.c_void_p, ctypes.c_long, CALLBACK_TYPE, ctypes.c_void_p,
]
_lib.ddg_sampler_run_exact.restype = ctypes.c_long
_lib.ddg_sampler_run_exact.errcheck = _check_int

_lib.ddg_sampler_run_naive.argtypes = [
    ctypes.c_void_p, ctypes.c_long, CALLBACK_TYPE, ctypes.c_void_p,
]
_lib.ddg_sampler_run_naive.restype = ctypes.c_long
_lib.ddg_sampler_run_naive.errcheck = _check_int

_lib.ddg_sampler_naive_importance_weight.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_naive_importance_weight.restype = ctypes.c_double

_lib.ddg_sampler_get_manifold.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_get_manifold.restype = ctypes.c_void_p
_lib.ddg_sampler_get_manifold.errcheck = _check_null

_lib.ddg_sampler_free.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_free.restype = None

# Sampler direct queries (avoid manifold copy)

_lib.ddg_sampler_f_vector.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.c_int,
]
_lib.ddg_sampler_f_vector.restype = ctypes.c_int
_lib.ddg_sampler_f_vector.errcheck = _check_int

_lib.ddg_sampler_importance_weight.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_importance_weight.restype = ctypes.c_double

_lib.ddg_sampler_simplices.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
]
_lib.ddg_sampler_simplices.restype = ctypes.c_long
_lib.ddg_sampler_simplices.errcheck = _check_int

_lib.ddg_sampler_degree.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
]
_lib.ddg_sampler_degree.restype = ctypes.c_long
_lib.ddg_sampler_degree.errcheck = _check_int

# ---------------------------------------------------------------------------
# GC control
# ---------------------------------------------------------------------------

_lib.ddg_gc_collect.argtypes = []
_lib.ddg_gc_collect.restype = None

_lib.ddg_gc_minimize.argtypes = []
_lib.ddg_gc_minimize.restype = None
