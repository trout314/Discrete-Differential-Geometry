"""ctypes loader and thin wrappers for the D shared library.

The library is rebuilt from the current source tree on import (a no-op when
already up to date), so scripts always run the latest code. Controls:

  DDG_BUILD=debug        Prefer builddir/ over builddir-release/ (default: release).
  DDG_LIBRARY=<path>     Explicit override. A directory is treated as a build
                         dir (kept fresh + loaded); a file is loaded as-is with
                         NO rebuild (use this to pin a specific binary).
  DDG_NO_AUTOBUILD=1     Skip the rebuild step (load whatever exists). Use when
                         the D toolchain isn't installed, or to pin the binary.
"""

import ctypes
import os
import subprocess
import sys
from pathlib import Path
from shutil import which

# ---------------------------------------------------------------------------
# Library discovery + build freshness
# ---------------------------------------------------------------------------

_EXT = {"linux": ".so", "darwin": ".dylib", "win32": ".dll"}.get(sys.platform, ".so")
_LIB_STEM = "ddg_dlang"          # meson shared_library() target name
_LIB_NAME = f"{_LIB_STEM}{_EXT}"  # output filename / ninja target


def _autobuild_disabled() -> bool:
    return os.environ.get("DDG_NO_AUTOBUILD", "").strip().lower() in (
        "1", "true", "yes", "on")


def _build_command(build_dir: Path):
    """Preferred build command for `build_dir`, or None if no tool is found."""
    if which("ninja"):
        # Ninja target is the output filename; near-instant when already fresh.
        return ["ninja", "-C", str(build_dir), _LIB_NAME]
    if which("meson"):
        return ["meson", "compile", "-C", str(build_dir), _LIB_STEM]
    return None


def _ensure_fresh(build_dir: Path) -> None:
    """Rebuild the shared library from current source if stale.

    A directory lock serializes the build so the equilibrium driver's fan-out
    of worker subprocesses can't race on it (the parent builds once; workers
    then see a no-op). A failed *build* raises rather than silently loading a
    stale library; a missing build *tool* only warns.
    """
    if _autobuild_disabled():
        return
    if not (build_dir / "build.ninja").exists():
        return  # not a source-tree build dir (e.g. installed package)

    lock_fd = None
    try:
        import fcntl
        lock_fd = os.open(str(build_dir / ".ddg_autobuild.lock"),
                          os.O_CREAT | os.O_RDWR, 0o644)
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
    except Exception:
        lock_fd = None  # locking unavailable (non-posix) — proceed unlocked

    try:
        cmd = _build_command(build_dir)
        if cmd is None:
            print(f"[ddg] WARNING: neither ninja nor meson found; loading "
                  f"existing {_LIB_NAME} from {build_dir} without a freshness "
                  f"check. Set DDG_NO_AUTOBUILD=1 to silence.", file=sys.stderr)
            return
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"[ddg] Rebuild of {_LIB_NAME} failed; refusing to load a "
                f"possibly-stale library.\n  cmd: {' '.join(cmd)}\n"
                f"{proc.stdout}\n{proc.stderr}")
        if "no work to do" not in proc.stdout:
            print(f"[ddg] rebuilt {_LIB_NAME} in {build_dir.name} "
                  f"(source changed)", file=sys.stderr)
    finally:
        if lock_fd is not None:
            os.close(lock_fd)


def _select_build_dir(root: Path):
    """Pick the source-tree build dir to use (release preferred), or None."""
    prefer_debug = os.environ.get("DDG_BUILD", "release").strip().lower() == "debug"
    release, debug = root / "builddir-release", root / "builddir"
    order = [debug, release] if prefer_debug else [release, debug]
    for d in order:
        if (d / "build.ninja").exists() or (d / _LIB_NAME).exists():
            return d
    return None


def _find_library() -> ctypes.CDLL:
    """Rebuild (if needed) and load ddg_dlang.so (or .dylib / .dll)."""
    # 1. Explicit override.
    env_override = os.environ.get("DDG_LIBRARY")
    if env_override:
        p = Path(env_override)
        if p.is_file():
            return ctypes.CDLL(str(p))            # pinned file: load as-is
        if p.is_dir():
            _ensure_fresh(p)
            candidate = p / _LIB_NAME
            if candidate.exists():
                return ctypes.CDLL(str(candidate))
        raise OSError(f"DDG_LIBRARY={env_override!r} does not point to {_LIB_NAME}")

    # 2. Source checkout: auto-build the preferred build dir, then load.
    root = Path(__file__).parent.parent.parent
    build_dir = _select_build_dir(root)
    if build_dir is not None:
        _ensure_fresh(build_dir)
        candidate = build_dir / _LIB_NAME
        if candidate.exists():
            return ctypes.CDLL(str(candidate))

    # 3. Installed package: .so shipped next to this module.
    pkg_local = Path(__file__).parent / _LIB_NAME
    if pkg_local.exists():
        return ctypes.CDLL(str(pkg_local))

    raise OSError(
        f"Cannot find {_LIB_NAME}. Build with "
        f"`meson compile -C builddir-release {_LIB_STEM}`, or set DDG_LIBRARY.")


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

_lib.ddg_manifold_save_with_comments.argtypes = [
    ctypes.c_void_p, ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_char_p), ctypes.c_int,
]
_lib.ddg_manifold_save_with_comments.restype = ctypes.c_int
_lib.ddg_manifold_save_with_comments.errcheck = _check_int

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

_lib.ddg_sampler_create_ext.argtypes = [
    ctypes.c_void_p,
    ctypes.c_int, ctypes.c_double,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_double,
]
_lib.ddg_sampler_create_ext.restype = ctypes.c_void_p
_lib.ddg_sampler_create_ext.errcheck = _check_null

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

_lib.ddg_sampler_set_callback_interval.argtypes = [ctypes.c_void_p, ctypes.c_long]
_lib.ddg_sampler_set_callback_interval.restype = ctypes.c_int
_lib.ddg_sampler_set_callback_interval.errcheck = _check_int

_lib.ddg_sampler_set_num_facets_target.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sampler_set_num_facets_target.restype = ctypes.c_int
_lib.ddg_sampler_set_num_facets_target.errcheck = _check_int

_lib.ddg_sampler_set_hinge_move_prob.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_hinge_move_prob.restype = ctypes.c_int
_lib.ddg_sampler_set_hinge_move_prob.errcheck = _check_int

_lib.ddg_sampler_set_num_facets_coef.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_num_facets_coef.restype = ctypes.c_int
_lib.ddg_sampler_set_num_facets_coef.errcheck = _check_int

_lib.ddg_sampler_set_num_hinges_coef.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_num_hinges_coef.restype = ctypes.c_int
_lib.ddg_sampler_set_num_hinges_coef.errcheck = _check_int

_lib.ddg_sampler_set_hinge_degree_variance_coef.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_hinge_degree_variance_coef.restype = ctypes.c_int
_lib.ddg_sampler_set_hinge_degree_variance_coef.errcheck = _check_int

_lib.ddg_sampler_set_codim3_degree_variance_coef.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_codim3_degree_variance_coef.restype = ctypes.c_int
_lib.ddg_sampler_set_codim3_degree_variance_coef.errcheck = _check_int

_lib.ddg_sampler_set_hinge_degree_target.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.ddg_sampler_set_hinge_degree_target.restype = ctypes.c_int
_lib.ddg_sampler_set_hinge_degree_target.errcheck = _check_int

_lib.ddg_sampler_reset_stats.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_reset_stats.restype = ctypes.c_int
_lib.ddg_sampler_reset_stats.errcheck = _check_int

_lib.ddg_sampler_track_move_counts.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sampler_track_move_counts.restype = ctypes.c_int
_lib.ddg_sampler_track_move_counts.errcheck = _check_int

_lib.ddg_sampler_reset_move_counts.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_reset_move_counts.restype = ctypes.c_int
_lib.ddg_sampler_reset_move_counts.errcheck = _check_int

_lib.ddg_sampler_move_counts.argtypes = [
    ctypes.c_void_p, ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
]
_lib.ddg_sampler_move_counts.restype = ctypes.c_long
_lib.ddg_sampler_move_counts.errcheck = _check_int

_lib.ddg_sampler_get_stats.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
    ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
    ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
    ctypes.c_int,
]
_lib.ddg_sampler_get_stats.restype = ctypes.c_int
_lib.ddg_sampler_get_stats.errcheck = _check_int

_lib.ddg_sampler_current_objective.argtypes = [ctypes.c_void_p]
_lib.ddg_sampler_current_objective.restype = ctypes.c_double

# ---------------------------------------------------------------------------
# Degree variance
# ---------------------------------------------------------------------------

_lib.ddg_manifold_degree_variance.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_manifold_degree_variance.restype = ctypes.c_double

_lib.ddg_sampler_degree_variance.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sampler_degree_variance.restype = ctypes.c_double

# ---------------------------------------------------------------------------
# Degree histogram
# ---------------------------------------------------------------------------

_lib.ddg_manifold_degree_histogram.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_long),
]
_lib.ddg_manifold_degree_histogram.restype = ctypes.c_long
_lib.ddg_manifold_degree_histogram.errcheck = _check_int

_lib.ddg_sampler_degree_histogram.argtypes = [
    ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_long),
]
_lib.ddg_sampler_degree_histogram.restype = ctypes.c_long
_lib.ddg_sampler_degree_histogram.errcheck = _check_int

# ---------------------------------------------------------------------------
# GC control
# ---------------------------------------------------------------------------

_lib.ddg_gc_collect.argtypes = []
_lib.ddg_gc_collect.restype = None

_lib.ddg_gc_minimize.argtypes = []
_lib.ddg_gc_minimize.restype = None

_lib.ddg_gc_stats.argtypes = [ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
_lib.ddg_gc_stats.restype = None
