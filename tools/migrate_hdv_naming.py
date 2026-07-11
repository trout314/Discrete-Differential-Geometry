#!/usr/bin/env python3
"""Rename HDV seed families from raw-coef (_HDV_{coef}) to SCALED coef/N
(_HDVs_{coef/N}) tokens, mirroring the earlier VDV_->VDVs_ migration. coef/N is
the natural HDV coupling (equipartition: a move shifts global HDV by ~O(1/N_edges),
N_edges proportional to N), so scaled naming labels the coupling consistently
across N. The raw coef stays in the .mfd metadata (source of truth); only the
filename token changes. Also fixes initial_triangulation cross-references pointing
at any renamed file.

Default is a DRY RUN (prints the rename map, writes nothing). Pass --apply.
Idempotent: once renamed there are no _HDV_ files left to match.
"""

import argparse
import glob
import os
import re
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_utils import build_seed_filename, load_seed_metadata, set_header_field


def _params(md):
    g = lambda k, d=0.0: float(md.get(k, d))
    return SimpleNamespace(
        num_facets_target=int(g("num_facets_target")),
        num_facets_coef=g("num_facets_coef"),
        hinge_degree_target=g("hinge_degree_target"),
        num_hinges_coef=g("num_hinges_coef"),
        hinge_degree_variance_coef=g("hinge_degree_variance_coef"),
        codim3_degree_variance_coef=g("codim3_degree_variance_coef"))


def new_name(path):
    md = load_seed_metadata(path)
    base = os.path.basename(path)
    m = re.search(r"_s(\d+)\.mfd$", base)
    idx = int(m.group(1)) if m else None
    topo = base.split("_", 1)[0]
    return build_seed_filename(topo, _params(md), seed_index=idx)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds-dir", default="seeds")
    ap.add_argument("--apply", action="store_true", help="Rename (default: dry run).")
    args = ap.parse_args()

    hdv_files = sorted(glob.glob(os.path.join(args.seeds_dir, "*_HDV_*.mfd")))
    rename = {}  # old basename -> new basename
    for p in hdv_files:
        ob, nb = os.path.basename(p), new_name(p)
        if nb != ob:
            rename[ob] = nb

    stems = {}
    for ob, nb in rename.items():
        stems[re.sub(r"_s\d+\.mfd$", "", ob)] = re.sub(r"_s\d+\.mfd$", "", nb)
    print(f"{len(hdv_files)} _HDV_ files, {len(rename)} to rename "
          f"({len(stems)} families):")
    for o in sorted(stems):
        print(f"  {o}\n    -> {stems[o]}")

    if not args.apply:
        print("\n(dry run -- pass --apply to rename)")
        return

    # collisions guard
    news = list(rename.values())
    if len(set(news)) != len(news):
        sys.exit("ERROR: rename target collision; aborting.")

    for ob, nb in rename.items():
        os.rename(os.path.join(args.seeds_dir, ob), os.path.join(args.seeds_dir, nb))
    fixed = 0
    for p in glob.glob(os.path.join(args.seeds_dir, "*.mfd")):
        init = load_seed_metadata(p).get("initial_triangulation", "")
        if init in rename:
            set_header_field(p, "initial_triangulation", rename[init])
            fixed += 1
    print(f"\nrenamed {len(rename)} files; fixed {fixed} initial_triangulation refs")


if __name__ == "__main__":
    main()
