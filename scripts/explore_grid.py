#!/usr/bin/env python3
"""Explore phase of an explore-then-produce sweep: a SHORT, `--dry-run` map of
the beta/N x edge-target x N grid via the vdv x edge_deg x num_facets factorial
bracket. Runs the chains + convergence gate per cell and reports the verdict
(binding feature, ESS, pass/fail) but WRITES NOTHING to the library -- use it to
locate the certifiable frontier and confirm run lengths before committing to the
production `produce_grid.py` pass.

Founds each cell from the nearest existing family at that (N, edge). Writes a
per-cell summary CSV under --out-dir.

Examples
--------
    # full-grid short map
    python scripts/explore_grid.py
    # one N tier, longer, to disambiguate a marginal frontier
    python scripts/explore_grid.py --only-n 1e3 --burnin 800 --n-samples 300
"""

import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
import grid_sweep as G


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bracket", nargs="+", default=G.DEFAULT_BRACKET,
                   help="Factorial bracket coordinates (default: vdv edge_deg num_facets).")
    p.add_argument("--replicas", type=int, default=16)
    p.add_argument("--num-hinges-coef", nargs="+", type=float, default=[2.0],
                   help="Edge-pin stiffness k AXIS (one or more; default 2.0). "
                        "edge_deg mixing ~ 1/k -- raising k helps edge_deg but "
                        "eventually breaks facet mixing (k>=4).")
    p.add_argument("--hdv-coef", nargs="+", type=float, default=[0.0],
                   help="Hinge-degree-variance coef AXIS (one or more; default 0, "
                        "unpenalized). Distinct k/hdv make distinct families.")
    p.add_argument("--burnin", type=int, default=300)
    p.add_argument("--n-samples", type=int, default=150)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--only-n", nargs="+", default=None,
                   help="Restrict to these N tokens/values (e.g. 1e3 562).")
    p.add_argument("--only-edge", nargs="+", default=None,
                   help="Restrict to these edge tokens (e.g. 5p0043).")
    p.add_argument("--only-bon", nargs="+", type=float, default=None,
                   help="Restrict to these beta/N values (e.g. 0.001 0.002).")
    p.add_argument("--seeds-dir", default="seeds",
                   help="Library to found from (relative to repo root or absolute).")
    p.add_argument("--out-dir", default="data/grid_explore")
    args = p.parse_args()

    G.sweep(_ROOT, dry_run=True, bracket=args.bracket, replicas=args.replicas,
            burnin=args.burnin, nsamp=args.n_samples, thin=args.thin,
            k_values=args.num_hinges_coef, hdv_values=args.hdv_coef,
            seeds_dir=args.seeds_dir, out_root=args.out_dir,
            only_n=args.only_n, only_edge=args.only_edge, only_bon=args.only_bon)


if __name__ == "__main__":
    main()
