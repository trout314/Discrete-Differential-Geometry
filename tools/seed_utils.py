"""Shared utilities for seed triangulation management."""

import json
import math
import os
import re
import subprocess


def set_header_field(path: str, key: str, value: str) -> None:
    """Set ``# key = value`` in an .mfd header, replacing an existing line for
    that key or inserting one before the facet data. Atomic; preserves all other
    header lines and the facet body."""
    with open(path) as f:
        lines = f.readlines()
    end = 0
    while end < len(lines) and lines[end].lstrip().startswith("#"):
        end += 1
    header, body, out, replaced = lines[:end], lines[end:], [], False
    for ln in header:
        content = ln.lstrip()[1:].strip()
        if "=" in content and content.split("=", 1)[0].strip() == key:
            out.append(f"# {key} = {value}\n"); replaced = True
        else:
            out.append(ln)
    if not replaced:
        out.append(f"# {key} = {value}\n")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.writelines(out + body)
    os.replace(tmp, path)


def get_free_memory_gb() -> float:
    """Return available system memory in GB by reading /proc/meminfo.

    Uses MemAvailable, which accounts for reclaimable caches/buffers.
    Returns float('inf') on non-Linux systems where /proc/meminfo is absent.
    """
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    # Value is in kB
                    kb = int(line.split()[1])
                    return kb / (1024 * 1024)
    except (FileNotFoundError, OSError):
        return float("inf")
    return float("inf")


def get_total_memory_gb() -> float:
    """Return total physical system memory in GB (MemTotal from /proc/meminfo).

    Returns float('inf') on non-Linux systems where /proc/meminfo is absent.
    """
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])  # value is in kB
                    return kb / (1024 * 1024)
    except (FileNotFoundError, OSError):
        return float("inf")
    return float("inf")


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


def round_sig(x: float, sig: int = 3) -> float:
    """Round x to `sig` significant figures. Keeps clean grid values clean
    (0.1 -> 0.1) while taming messy ratios (0.101266 -> 0.101) to a short token."""
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sig - 1))


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

    # ED (edge/hinge degree) term. REQUIRED whenever a fixed-target quadratic
    # (VDQ/EDQ) is present -- the target d-bar_t defines both the quadratic
    # targets and the per-tet scaling of their tokens -- with a ZERO pin
    # coefficient allowed (e.g. ED5p1043_0 = target defined, pin off).
    vdq_coef = getattr(params, "codim3_degree_target_coef", 0) or 0
    edq_coef = getattr(params, "hinge_degree_target_coef", 0) or 0
    if params.num_hinges_coef != 0 or vdq_coef or edq_coef:
        parts.append(f"ED{encode_float(params.hinge_degree_target)}_{encode_float(params.num_hinges_coef)}")

    # VDV (vertex degree variance = codim3 degree variance) term.
    # Encoded as the SCALED coefficient beta/N ("VDVs_"): one Pachner move shifts
    # global VDV by ~O(1/N), so beta/N -- not raw beta -- is the coupling on which
    # the equilibrium VDV collapses across volumes. Raw beta is recoverable as
    # (beta/N) * num_facets_target from the metadata.
    if params.codim3_degree_variance_coef != 0:
        beta_over_n = params.codim3_degree_variance_coef / params.num_facets_target
        parts.append(f"VDVs_{encode_float(beta_over_n)}")

    # HDV (hinge/edge degree variance) term — SCALED coef/N ("HDVs_"), like VDVs_:
    # coef/N is the natural HDV coupling too (equipartition; a move shifts global
    # HDV by ~O(1/N_edges), N_edges ∝ N), so the same coupling gets the same label
    # across N. Rounded to 3 sig figs; raw coef recoverable from the metadata.
    if params.hinge_degree_variance_coef != 0:
        hdv_over_n = params.hinge_degree_variance_coef / params.num_facets_target
        parts.append(f"HDVs_{encode_float(round_sig(hdv_over_n, 3))}")

    # VDQ/EDQ (fixed-target "degree quadratic", strictly local) — tokens carry
    # the PER-TETRAHEDRON scaled coupling lambda = c * (elements per tet at the
    # TARGET edge degree): lambda_v = c_v*(6/dbar_t - 1), lambda_e = c_e*6/dbar_t.
    # Under the exact old->new coupling map these equal the old scaled labels
    # (VDVs_g <-> VDQ_g, HDVs_h <-> EDQ_h), so matched physics carries the same
    # number across modes, edges, and N. Raw per-element couplings c_v/c_e are
    # exact in metadata (codim3/hinge_degree_target_coef); conversion uses the
    # target dbar_t from the ED token (always present when these terms are on).
    if vdq_coef:
        lam_v = vdq_coef * (6.0 / params.hinge_degree_target - 1.0)
        parts.append(f"VDQ_{encode_float(round_sig(lam_v, 3))}")
    if edq_coef:
        lam_e = edq_coef * 6.0 / params.hinge_degree_target
        parts.append(f"EDQ_{encode_float(round_sig(lam_e, 3))}")

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


# ---------------------------------------------------------------------------
# Flattened provenance history
#
# Each seed's header carries a self-contained `history` field: a JSON list of
# "legs", one per contiguous run under a single build + objective, oldest first.
# A leg answers "how many sweeps at what objective", so a seed's full lineage is
# reconstructable without any ancestor file on disk. All seed writers funnel
# through build_metadata_comments(), which requires the source's prior history
# (or root="sphere") plus this run's legs -- so lineage can never be dropped
# silently. See read_history() to obtain the prior history from a source seed.
# ---------------------------------------------------------------------------

def obj_of(params) -> dict:
    """Compact objective dict for a provenance leg, from a SamplerParams-like."""
    obj = {
        "nf": params.num_facets_target, "nf_c": params.num_facets_coef,
        "ht": params.hinge_degree_target, "nh_c": params.num_hinges_coef,
        "hdv_c": params.hinge_degree_variance_coef,
        "vdv_c": params.codim3_degree_variance_coef,
    }
    # Fixed-target quadratics (VDQ/EDQ): only recorded when active, so legs of
    # old-mode runs are byte-identical to before.
    if getattr(params, "codim3_degree_target_coef", 0):
        obj["vdq_c"] = params.codim3_degree_target_coef
        obj["vdq_t"] = getattr(params, "codim3_degree_target", 0)
    if getattr(params, "hinge_degree_target_coef", 0):
        obj["edq_c"] = params.hinge_degree_target_coef
    return obj


def make_leg(op: str, obj: dict, sweeps, *, from_: str = "prev",
             commit: str = None, dirty: bool = None,
             tried=None, accepted=None, reconstructed: bool = False) -> dict:
    """Construct one provenance leg. `sweeps` is nominal; `tried`/`accepted` are
    the sampler move counts for this leg (None when not cleanly attributable)."""
    if commit is None:
        commit, dirty = get_git_info()
    leg = {
        "op": op, "from": from_,
        "commit": commit[:7] if commit != "unknown" else commit,
        "dirty": bool(dirty),
        "obj": obj, "sweeps": sweeps,
        "tried": tried, "accepted": accepted,
    }
    if reconstructed:
        leg["reconstructed"] = True
    return leg


def read_history(path: str) -> list:
    """Return the flattened provenance history (list of legs) recorded in an
    .mfd header, or [] if the file predates history tracking. Writers pass the
    result as build_metadata_comments(prior_history=...)."""
    h = load_seed_metadata(path).get("history")
    if not h:
        return []
    try:
        legs = json.loads(h)
        return legs if isinstance(legs, list) else []
    except (ValueError, TypeError):
        return []


def verify_history(path: str, manifold_view=None) -> tuple[bool, list]:
    """Check a seed's flattened history is well-formed. Returns (ok, issues)."""
    md = load_seed_metadata(path)
    raw = md.get("history")
    if not raw:
        return False, ["no history field"]
    try:
        legs = json.loads(raw)
    except (ValueError, TypeError):
        return False, ["history is not valid JSON"]
    issues = []
    if not legs:
        issues.append("empty history")
    else:
        if legs[0].get("from") != "sphere":
            issues.append(f"not rooted at sphere (first leg from={legs[0].get('from')!r})")
        for i, leg in enumerate(legs[1:], start=1):
            if leg.get("from") != "prev":
                issues.append(f"leg {i} from={leg.get('from')!r} (expected 'prev')")
            for k in ("op", "obj", "sweeps"):
                if k not in leg:
                    issues.append(f"leg {i} missing '{k}'")
    if manifold_view is not None:
        rec, act = md.get("actual_f_vector"), str(list(manifold_view.f_vector))
        if rec is not None and rec != act:
            issues.append(f"final f_vector mismatch: header {rec} vs manifold {act}")
    return (not issues), issues


def history_fields(history: list, note: str = None) -> list:
    """Comment lines encoding a flattened history: root/totals/note + JSON.

    Shared by build_metadata_comments (live writes) and the back-fill migration
    so both emit an identical schema."""
    rooted = bool(history) and history[0].get("from") == "sphere"
    tot_s = sum(l["sweeps"] for l in history if isinstance(l.get("sweeps"), (int, float)))
    tot_t = sum(l["tried"] for l in history if isinstance(l.get("tried"), (int, float)))
    out = [f"history_root = {'sphere' if rooted else 'incomplete'}",
           f"history_total_sweeps = {tot_s}",
           f"history_total_tried = {tot_t}"]
    if note:
        out.append(f"history_note = {note}")
    out.append("history = " + json.dumps(history, separators=(",", ":")))
    return out


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
    hinge_degree_target_coef: float = 0.0,
    codim3_degree_target_coef: float = 0.0,
    codim3_degree_target: float = 0.0,
    growth_step_size: int,
    eq_sweeps_per_step: int,
    equilibration_sweeps: int,
    manifold_view,
    objective: float,
    sampler_stats=None,
    legs: list = None,
    prior_history: list = None,
    root: str = None,
    history_note: str = None,
) -> list[str]:
    """Build the list of comment strings for .save().

    Lineage is mandatory: pass this run's `legs` (a non-empty list from
    make_leg) AND exactly one of `prior_history` (the source seed's history,
    from read_history) or `root="sphere"` (a from-scratch seed). The emitted
    `history` field is prior_history + legs, so it can never be silently
    dropped. `history_note` is free text for caveats (e.g. reconstructed
    ancestry).
    """
    if legs is None or len(legs) == 0:
        raise ValueError(
            "build_metadata_comments requires legs=[...]: declare this run's "
            "provenance legs (see make_leg). Lineage is not optional.")
    if (prior_history is None) == (root is None):
        raise ValueError(
            "build_metadata_comments requires exactly one of prior_history="
            "<list from read_history> or root='sphere'.")
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
    comments.append(f"hinge_degree_target_coef = {hinge_degree_target_coef}")
    comments.append(f"codim3_degree_target_coef = {codim3_degree_target_coef}")
    comments.append(f"codim3_degree_target = {codim3_degree_target}")
    comments.append(f"growth_step_size = {growth_step_size}")
    comments.append(f"eq_sweeps_per_step = {eq_sweeps_per_step}")
    comments.append(f"equilibration_sweeps = {equilibration_sweeps}")

    mfd = manifold_view
    fv = list(mfd.f_vector)
    dim = mfd.dimension
    # Record only observables a downstream reader of the .mfd file cannot cheaply
    # recompute for themselves:
    #   * the f-vector -- the fundamental simplex counts;
    #   * the degree *variances* -- these need the full degree distribution (the
    #     means are fixed by the f-vector, but the spread is not); the sampler
    #     tracks them for free, whereas a file reader would have to iterate;
    #   * valid_moves -- needs move enumeration.
    # Deliberately omitted: num_facets and the mean degrees (exact functions of
    # the f-vector), and euler_characteristic / is_orientable (topological
    # invariants -- unchanged by Pachner moves, so a change would be a bug, not
    # data worth storing; is_orientable is also O(n^2) and OOMs on large
    # manifolds, which is what froze the reburn).
    comments.append(f"actual_f_vector = {fv}")
    comments.append(f"actual_hinge_degree_variance = {mfd.degree_variance(dim - 2):.6f}")
    comments.append(f"actual_vertex_degree_variance = {mfd.degree_variance(0):.6f}")
    comments.append(f"valid_moves = {mfd.count_valid_moves()}")
    comments.append(f"final_objective = {objective:.6f}")

    if sampler_stats is not None:
        stats = sampler_stats
        accept_rate = stats.total_accepted / stats.total_tried if stats.total_tried > 0 else 0.0
        comments.append(f"eq_total_tried = {stats.total_tried}")
        comments.append(f"eq_total_accepted = {stats.total_accepted}")
        comments.append(f"eq_acceptance_rate = {accept_rate:.6f}")
        bt = list(stats.bistellar_tries)
        ba = list(stats.bistellar_accepts)
        comments.append(f"eq_bistellar_tries = {bt}")
        comments.append(f"eq_bistellar_accepts = {ba}")

    # Flattened provenance: prior lineage + this run's legs (see module header).
    history = (list(prior_history) if prior_history is not None else []) + list(legs)
    comments.extend(history_fields(history, history_note))
    return comments
