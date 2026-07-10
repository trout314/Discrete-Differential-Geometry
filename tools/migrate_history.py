#!/usr/bin/env python3
"""Back-fill the flattened `history` field into existing seed .mfd headers.

Every seed written before provenance tracking records only a single "leg"
(its own pass) plus a parent-filename pointer. This tool reconstructs each
seed's full sphere-rooted history and (with --apply) rewrites the header to be
self-contained, so no seed depends on an ancestor file existing on disk.

Three cases:
  * self-contained (initial_triangulation = standard_sphere): [grow, equilibrate]
    legs synthesized from the file's own flat metadata.
  * chained w/ surviving parent: recurse into the parent, append this leg.
  * orphaned (parent file gone): reconstruct the parent's grow-from-sphere legs
    from the child's objective (reburn preserves it) + the generation grid
    (growth_step_size = max(50, N//20), eq_sweeps_per_step = 5) + a per-N burn
    anchor read from any surviving grow-from-sphere seed. Reconstructed legs are
    tagged "reconstructed":true and a history_note records the situation.

Default is a DRY RUN (counts only, writes nothing). Pass --apply to rewrite.
"""

import argparse
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_utils import (history_fields, load_seed_metadata, make_leg,
                        read_history, verify_history)

SPHERE_FACETS = {2: 4, 3: 5, 4: 6}      # standardSphere(dim) facet count
MARCH_BUILD = "a63f4d6"                  # generation-grid batch build (reconstructed legs)


def _obj_from_md(md):
    g = lambda k, d=0.0: float(md.get(k, d))
    return {"nf": int(g("num_facets_target")), "nf_c": g("num_facets_coef"),
            "ht": g("hinge_degree_target"), "nh_c": g("num_hinges_coef"),
            "hdv_c": g("hinge_degree_variance_coef"),
            "vdv_c": g("codim3_degree_variance_coef")}


def _stats_from_md(md, sweeps):
    """(tried, accepted), but only when tried is consistent with sweeps*N.

    Produce chains reset_stats() per measured sample, so their flat
    eq_total_tried is just the last sample -- attaching it to a multi-hundred-
    sweep leg would be wrong. A move count within [0.25, 4]x of sweeps*N is a
    genuine whole-leg total (e.g. reburn); anything else is dropped to None.
    """
    t, a = md.get("eq_total_tried"), md.get("eq_total_accepted")
    if t is None:
        return None, None
    t, a = int(t), (int(a) if a is not None else None)
    n = int(float(md.get("num_facets_target", 0)))
    if sweeps and n and 0.25 <= t / (sweeps * n) <= 4.0:
        return t, a
    return None, None


def build_index(dirs):
    """basename -> [paths]; and the per-N burn anchor from grow-from-sphere seeds."""
    index, anchor = {}, {}
    for d in dirs:
        for p in glob.glob(os.path.join(d, "*.mfd")):
            index.setdefault(os.path.basename(p), []).append(p)
            md = load_seed_metadata(p)
            if md.get("initial_triangulation", "").startswith("standard_sphere") \
                    and "equilibration_sweeps" in md and "num_facets_target" in md:
                anchor.setdefault(int(float(md["num_facets_target"])),
                                  int(float(md["equilibration_sweeps"])))
    return index, anchor


def _anchor_burn(n, anchor):
    """Burn length for facet count n: exact anchor, else nearest in log-N."""
    if n in anchor:
        return anchor[n], True
    if not anchor:
        return 500, False
    nearest = min(anchor, key=lambda a: abs(math.log(a) - math.log(n)))
    return anchor[nearest], False


def _grow_legs(obj, dim, gs, eqps, burn, commit, dirty, *, reconstructed):
    """A [grow, equilibrate] pair rooting a triangulation at the sphere."""
    start = SPHERE_FACETS.get(dim, 5)
    steps = max(0, math.ceil((obj["nf"] - start) / gs)) if gs else 0
    legs = [make_leg("grow", obj, steps * eqps, from_="sphere",
                     commit=commit, dirty=dirty, reconstructed=reconstructed)]
    if burn:
        legs.append(make_leg("equilibrate", obj, burn, commit=commit, dirty=dirty,
                             reconstructed=reconstructed))
    return legs


def _own_legs(md, path, anchor):
    """Legs contributed by THIS file's own pass, from its flat metadata."""
    obj = _obj_from_md(md)
    dim = int(float(md.get("dimension", 3)))
    commit = md.get("git_commit", "unknown")[:7]
    dirty = md.get("git_dirty") == "true"
    gs = int(float(md.get("growth_step_size", 0)))
    eq = int(float(md.get("equilibration_sweeps", 0)))
    tried, accepted = _stats_from_md(md, eq)
    if gs > 0:                                     # grew from sphere this pass
        eqps = int(float(md.get("eq_sweeps_per_step", 5)))
        legs = _grow_legs(obj, dim, gs, eqps, eq, commit, dirty, reconstructed=True)
        if legs and eq and tried is not None:
            legs[-1]["tried"], legs[-1]["accepted"] = tried, accepted
        return legs
    op = "reburn" if os.sep + "seeds_reburned" + os.sep in path else "equilibrate"
    return [make_leg(op, obj, eq, from_="prev", commit=commit, dirty=dirty,
                     tried=tried, accepted=accepted, reconstructed=True)]


def _reconstruct_orphan(md, parent_name, anchor):
    """Reconstruct a gone parent's grow-from-sphere legs from the child's obj."""
    obj = _obj_from_md(md)
    dim = int(float(md.get("dimension", 3)))
    gs = max(50, obj["nf"] // 20)
    burn, exact = _anchor_burn(obj["nf"], anchor)
    legs = _grow_legs(obj, dim, gs, 5, burn, MARCH_BUILD, True, reconstructed=True)
    note = (f"ancestor '{parent_name}' not retained; earlier legs reconstructed "
            f"from generation grid; burn length "
            + ("from same-N anchor" if exact else "ESTIMATED (no same-N anchor)"))
    return legs, note


def resolve_history(path, index, anchor, memo, resolving):
    """Full sphere-rooted flattened history for `path` (list, note)."""
    if path in memo:
        return memo[path]
    existing = read_history(path)
    if existing:                                   # already migrated: idempotent
        memo[path] = (existing, None)
        return memo[path]
    md = load_seed_metadata(path)
    own = _own_legs(md, path, anchor)
    init = md.get("initial_triangulation", "")
    note = None
    if init.startswith("standard_sphere"):
        prior = []                                 # own legs already root at sphere
    else:
        cands = [c for c in index.get(init, []) if os.path.abspath(c) != os.path.abspath(path)]
        cands = [c for c in cands if os.path.abspath(c) not in resolving]
        if cands:
            prior, note = resolve_history(cands[0], index, anchor, memo,
                                          resolving | {os.path.abspath(path)})
        else:
            prior, note = _reconstruct_orphan(md, init, anchor)
    full = list(prior) + list(own)
    memo[path] = (full, note)
    return memo[path]


def apply_history(path, legs, note):
    """Append the flattened-history comment lines to an .mfd header, atomically,
    preserving all existing header lines and facet data. Idempotent: a file that
    already has a `history` field is left untouched."""
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and lines[i].lstrip().startswith("#"):
        if lines[i].lstrip()[1:].lstrip().startswith("history"):
            return False                        # already has history -> skip
        i += 1
    new = lines[:i] + [f"# {c}\n" for c in history_fields(legs, note)] + lines[i:]
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.writelines(new)
    os.replace(tmp, path)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dirs", nargs="+",
                    default=["seeds", "seeds_reburned", "seeds_uncertified"])
    ap.add_argument("--apply", action="store_true", help="Rewrite headers (default: dry run).")
    ap.add_argument("--sample", type=int, default=0, help="Print N resolved histories and exit.")
    args = ap.parse_args()

    index, anchor = build_index(args.dirs)
    all_paths = sorted(p for ps in index.values() for p in ps)
    print(f"indexed {len(all_paths)} seeds across {args.dirs}; "
          f"{len(anchor)} burn anchors", file=sys.stderr)

    memo = {}
    from collections import Counter
    stat = Counter()
    samples = []
    for p in all_paths:
        legs, note = resolve_history(p, index, anchor, memo, set())
        rooted = bool(legs) and legs[0].get("from") == "sphere"
        stat["rooted" if rooted else "unrooted"] += 1
        if note and "ESTIMATED" in note:
            stat["burn_estimated"] += 1
        elif note:
            stat["reconstructed_ancestry"] += 1
        if args.sample and len(samples) < args.sample:
            samples.append((p, legs, note, rooted))
        if args.apply and not args.sample:
            stat["written" if apply_history(p, legs, note) else "already_had_history"] += 1

    if args.sample:
        for p, legs, note, rooted in samples:
            print(f"\n{p}  (rooted={rooted})")
            for lg in legs:
                r = " [recon]" if lg.get("reconstructed") else ""
                print(f"   {lg['op']:<12} from={lg['from']:<6} sweeps={lg['sweeps']:<7} "
                      f"tried={lg['tried']} commit={lg['commit']}{r}")
            if note:
                print(f"   note: {note}")
        return

    print("\n=== dry-run summary ===")
    for k, v in stat.items():
        print(f"  {k:<22} {v}")
    print(f"  total                  {len(all_paths)}")
    if not args.apply:
        print("\n(dry run — nothing written; pass --apply to rewrite headers)")


if __name__ == "__main__":
    main()
