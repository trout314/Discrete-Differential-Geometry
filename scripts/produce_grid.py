#!/usr/bin/env python3
"""Produce phase of an explore-then-produce sweep: certify the beta/N x
edge-target x N grid at PRODUCTION length via the vdv x edge_deg x num_facets
factorial bracket, copying passing families into seeds/ (VDVs_ naming).

Founds each cell from the nearest existing family (progressive: later cells start
from just-produced closer-beta families). Idempotent -- `--produce` skips families
already in seeds/, so re-runs resume. By default prunes the known-dead corner
(beta/N>=0.008 at raised edge targets 5p1043/5p2043); --no-prune keeps them.

WRITES TO THE LIBRARY (seeds/). Writes a per-cell summary CSV under --out-dir.

Examples
--------
    # full production pass (55 cells after prune)
    python scripts/produce_grid.py
    # re-run a single marginal cell longer
    python scripts/produce_grid.py --only-n 562 --only-edge 5p2043 \
        --only-bon 0.001 --burnin 8000 --n-samples 2500
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
    p.add_argument("--replicas", type=int, default=32)
    p.add_argument("--num-hinges-coef", nargs="+", type=float, default=[2.0],
                   help="Edge-pin stiffness k AXIS (one or more; default 2.0). "
                        "Part of the objective/filename, so each k makes a distinct "
                        "family (no collision with the k=2 library).")
    p.add_argument("--hdv-coef", nargs="+", type=float, default=[0.0],
                   help="Hinge-degree-variance coef AXIS (one or more; default 0, "
                        "unpenalized). Each hdv value makes a distinct family "
                        "(adds an _HDV_ token to the filename).")
    p.add_argument("--burnin", type=int, default=5000)
    p.add_argument("--n-samples", type=int, default=1500)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--adaptive", action="store_true",
                   help="Retry a cell that FAILS on VDV-only under-burn once at "
                        "--retry-burnin/--retry-n-samples (genuine frontier breaks "
                        "are not retried). Lets --burnin be the shorter 'long' "
                        "setting while still clearing slow low-beta/N VDV cells.")
    p.add_argument("--retry-burnin", type=int, default=5000)
    p.add_argument("--retry-n-samples", type=int, default=1500)
    p.add_argument("--no-prune", action="store_true",
                   help="Run the known-dead corner (beta/N>=0.008 at raised edges) too.")
    p.add_argument("--only-n", nargs="+", default=None,
                   help="Restrict to these N tokens/values (e.g. 562).")
    p.add_argument("--only-edge", nargs="+", default=None,
                   help="Restrict to these edge tokens (e.g. 5p2043).")
    p.add_argument("--beta-over-n", nargs="+", type=float, default=None,
                   help="VDV coupling beta/N AXIS (one or more; default the standard "
                        "grid). Pass 0 for base (no-VDV) families.")
    p.add_argument("--only-bon", nargs="+", type=float, default=None,
                   help="Restrict the beta/N axis to these values (a filter).")
    p.add_argument("--seeds-dir", default="seeds",
                   help="Library to found from and copy into (repo-root-relative "
                        "or absolute).")
    p.add_argument("--out-dir", default="data/grid_produce")
    args = p.parse_args()

    kw = dict(beta_over_N=args.beta_over_n) if args.beta_over_n is not None else {}
    G.sweep(_ROOT, dry_run=False, bracket=args.bracket, replicas=args.replicas,
            burnin=args.burnin, nsamp=args.n_samples, thin=args.thin,
            k_values=args.num_hinges_coef, hdv_values=args.hdv_coef,
            adaptive=args.adaptive, retry_burnin=args.retry_burnin,
            retry_nsamp=args.retry_n_samples,
            prune=None if args.no_prune else G.prune_raised_frontier,
            seeds_dir=args.seeds_dir, out_root=args.out_dir,
            only_n=args.only_n, only_edge=args.only_edge, only_bon=args.only_bon, **kw)


if __name__ == "__main__":
    main()
