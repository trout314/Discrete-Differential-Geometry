"""Shared utilities for seed triangulation management."""

import math
import re
import subprocess


def encode_float(x: float) -> str:
    """Encode a float for use in filenames.

    Examples: 1000 -> '1e3', 0.1 -> '1e-1', 5.1 -> '5p1', 0 -> '0'
    """
    if x == 0:
        return "0"
    # Check if it's an exact power-of-10 times a small integer
    # Try scientific notation for clean cases
    if x == int(x) and x != 0:
        x_int = int(x)
        # Check for power of 10
        if x_int > 0:
            log = math.log10(x_int)
            if log == int(log) and log >= 1:
                return f"1e{int(log)}"
        # Try n * 10^k form
        if x_int > 0:
            s = str(x_int)
            trailing = len(s) - len(s.rstrip("0"))
            if trailing >= 1:
                mantissa = x_int // (10 ** trailing)
                if mantissa * (10 ** trailing) == x_int and mantissa < 10:
                    return f"{mantissa}e{trailing}"
        return str(x_int)
    # Negative exponents: 0.1 -> 1e-1, 0.01 -> 1e-2, etc.
    if 0 < x < 1:
        log = math.log10(x)
        if log == int(log):
            return f"1e{int(log)}"
        # Try n * 10^-k: 0.5 -> 5e-1, 0.05 -> 5e-2
        for exp in range(1, 10):
            mantissa = x * (10 ** exp)
            if mantissa == int(mantissa) and int(mantissa) < 10:
                return f"{int(mantissa)}e-{exp}"
    # Fallback: use 'p' for decimal point
    s = f"{x:g}"
    return s.replace(".", "p")


def decode_float(s: str) -> float:
    """Decode a filename-encoded float back to a number."""
    if s == "0":
        return 0.0
    s = s.replace("p", ".")
    if "e" in s:
        return float(s)
    return float(s)


def build_seed_filename(topology: str, params, seed_index: int | None = None) -> str:
    """Construct a seed filename from topology and SamplerParams-like object.

    Parameters
    ----------
    topology : str
        e.g. 'S3'
    params : object
        Must have attributes: num_facets_target, num_facets_coef,
        hinge_degree_target, num_hinges_coef, hinge_degree_variance_coef,
        codim3_degree_variance_coef.
    seed_index : int, optional
        Replica index. When set, appends ``_s{index:03d}`` before ``.mfd``.
    """
    parts = [topology]

    # N term (always present)
    parts.append(f"N{encode_float(params.num_facets_target)}_{encode_float(params.num_facets_coef)}")

    # ED (edge/hinge degree) term
    if params.num_hinges_coef != 0:
        parts.append(f"ED{encode_float(params.hinge_degree_target)}_{encode_float(params.num_hinges_coef)}")

    # VDV (vertex degree variance = codim3 degree variance) term
    if params.codim3_degree_variance_coef != 0:
        parts.append(f"VDV_{encode_float(params.codim3_degree_variance_coef)}")

    # HDV (hinge degree variance) term — last when present
    if params.hinge_degree_variance_coef != 0:
        parts.append(f"HDV_{encode_float(params.hinge_degree_variance_coef)}")

    name = "_".join(parts)
    if seed_index is not None:
        name += f"_s{seed_index:03d}"
    return name + ".mfd"


def load_seed_metadata(path: str) -> dict:
    """Parse # key = value header lines from an .mfd file."""
    metadata = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                break
            # Strip leading '# '
            content = line[1:].strip()
            if "=" in content:
                key, _, value = content.partition("=")
                metadata[key.strip()] = value.strip()
    return metadata


def get_git_info() -> tuple[str, bool]:
    """Return (commit_hash, is_dirty) for the current repo."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
        status = subprocess.check_output(
            ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL
        ).decode().strip()
        return commit, len(status) > 0
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown", True


def build_metadata_comments(
    *,
    topology: str,
    dimension: int,
    initial_triangulation: str,
    num_facets_target: int,
    num_facets_coef: float,
    hinge_degree_target: float,
    num_hinges_coef: float,
    hinge_degree_variance_coef: float,
    codim3_degree_variance_coef: float,
    growth_step_size: int,
    eq_sweeps_per_step: int,
    equilibration_sweeps: int,
    manifold_view,
    objective: float,
) -> list[str]:
    """Build the list of comment strings for .save()."""
    git_commit, git_dirty = get_git_info()
    comments = []
    comments.append(f"git_commit = {git_commit}")
    comments.append(f"git_dirty = {'true' if git_dirty else 'false'}")
    comments.append(f"initial_triangulation = {initial_triangulation}")
    comments.append(f"topology = {topology}")
    comments.append(f"dimension = {dimension}")
    comments.append(f"num_facets_target = {num_facets_target}")
    comments.append(f"num_facets_coef = {num_facets_coef}")
    comments.append(f"hinge_degree_target = {hinge_degree_target}")
    comments.append(f"num_hinges_coef = {num_hinges_coef}")
    comments.append(f"hinge_degree_variance_coef = {hinge_degree_variance_coef}")
    comments.append(f"codim3_degree_variance_coef = {codim3_degree_variance_coef}")
    comments.append(f"growth_step_size = {growth_step_size}")
    comments.append(f"eq_sweeps_per_step = {eq_sweeps_per_step}")
    comments.append(f"equilibration_sweeps = {equilibration_sweeps}")

    mfd = manifold_view
    fv = list(mfd.f_vector)
    dim = mfd.dimension
    comments.append(f"actual_num_facets = {mfd.num_facets}")
    comments.append(f"actual_f_vector = {fv}")
    comments.append(f"actual_mean_hinge_degree = {mfd.mean_degree(dim - 2):.3f}")
    comments.append(f"actual_hinge_degree_variance = {mfd.degree_variance(dim - 2):.2f}")
    comments.append(f"actual_vertex_degree_variance = {mfd.degree_variance(0):.1f}")
    comments.append(f"final_objective = {objective:.1f}")

    return comments
