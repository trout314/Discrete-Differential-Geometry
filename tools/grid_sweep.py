"""Shared engine for beta/N x edge-target x N grid sweeps over
`equilibrium_vdv.py --produce`.

Used by `scripts/explore_grid.py` (short `--dry-run` map, writes nothing) and
`scripts/produce_grid.py` (production certify + copy into seeds/). Each grid cell
is one `--produce` invocation with the vdv x edge_deg x num_facets factorial
bracket, founded from the NEAREST existing family at that (N, edge). Because
founding re-globs seeds/ every cell, a produce sweep enjoys progressive founding
(later cells start from the just-created closer-beta family).
"""

import csv
import glob
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_utils import load_seed_metadata

# Standard grid (the certified N=1e4 grid, extended to small N).
N_TIERS = [("1e2", 100), ("178", 178), ("316", 316), ("562", 562), ("1e3", 1000)]
EDGES = [("5p0043", 5.0043), ("5p1043", 5.1043), ("5p2043", 5.2043)]
BETA_OVER_N = [0.001, 0.002, 0.004, 0.008, 0.01]
DEFAULT_BRACKET = ["vdv", "edge_deg", "num_facets"]


def founding_seed(root, ntok, edgetok, beta_target, seeds_dir="seeds"):
    """s000 of the existing family at (N, edge) whose beta is nearest beta_target
    (globs seeds_dir live, so newly produced families are candidates too)."""
    pat = os.path.join(root, seeds_dir, f"S3_N{ntok}_1e-1_ED{edgetok}_2*_s000.mfd")
    best, best_d = None, None
    for p in glob.glob(pat):
        try:
            b = float(load_seed_metadata(p).get("codim3_degree_variance_coef", 0.0))
        except Exception:
            continue
        d = abs(b - beta_target)
        if best_d is None or d < best_d:
            best, best_d = p, d
    return best


def _parse_verdict(out, dry_run):
    """(verdict, detail) from a --produce run's combined stdout+stderr."""
    if "not duplicating" in out:
        return "SKIP-EXISTS", ""
    if dry_run:
        line = next((l for l in out.splitlines() if "[dry-run]" in l), "")
        v = "PASS" if "[dry-run] PASS" in line else ("FAIL" if "[dry-run] FAIL" in line
                                                     else "ERR")
        m = re.search(r"binding qRhat=\S+.*", line)
        return v, (m.group(0)[:90] if m else line[:90])
    if "CONVERGED" in out:
        m = re.search(r"binding qRhat \S+, ESS \S+", out)
        return "PASS", (m.group(0) if m else "")
    m = re.search(r"(not converged in .*|ESS \d+.*run longer)", out)
    return "FAIL", (m.group(1)[:90] if m else "ERR")


def run_cell(root, seed, n, edgeval, beta, *, bracket, replicas, burnin, nsamp,
             thin, dry_run, seeds_dir, out_dir):
    """One grid cell = one `equilibrium_vdv.py --produce` invocation."""
    cmd = [sys.executable, os.path.join(root, "scripts", "equilibrium_vdv.py"),
           "--produce", "--seed-file", seed, "--n-target", str(n), "--beta", str(beta),
           "--bracket", *bracket, "--hinge-target", str(edgeval),
           "--num-hinges-coef", "2.0", "--num-facets-coef", "0.1",
           "--replicas", str(replicas), "--production-burnin", str(burnin),
           "--n-samples", str(nsamp), "--thin", str(thin),
           "--seeds-dir", os.path.join(root, seeds_dir), "--output-dir", out_dir]
    if dry_run:
        cmd.append("--dry-run")
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    sys.stdout.write(out); sys.stdout.flush()
    return _parse_verdict(out, dry_run)


def sweep(root, *, dry_run, bracket=DEFAULT_BRACKET, replicas, burnin, nsamp, thin,
          n_tiers=N_TIERS, edges=EDGES, beta_over_N=BETA_OVER_N, prune=None,
          seeds_dir="seeds", out_root="data/grid_sweep", only_n=None, only_edge=None,
          only_bon=None):
    """Run the grid. `prune(edgetok, bon) -> bool` skips cells; only_* (lists)
    restrict to a subset. Writes an incremental summary CSV under out_root and
    returns the collected rows."""
    out_root = out_root if os.path.isabs(out_root) else os.path.join(root, out_root)
    os.makedirs(out_root, exist_ok=True)
    rows = []
    for ntok, n in n_tiers:
        if only_n and ntok not in only_n and str(n) not in only_n:
            continue
        for edgetok, edgeval in edges:
            if only_edge and edgetok not in only_edge:
                continue
            for bon in beta_over_N:
                if only_bon and not any(abs(bon - float(b)) < 1e-9 for b in only_bon):
                    continue
                if prune and prune(edgetok, bon):
                    continue
                beta = bon * n
                seed = founding_seed(root, ntok, edgetok, beta, seeds_dir)
                if not seed:
                    print(f"SKIP N={ntok} ED={edgetok} b/N={bon}: no founding family",
                          flush=True)
                    continue
                cell = f"N{ntok}_ED{edgetok}_bN{bon:g}"
                print(f"\n===== {cell}  beta={beta:g}  "
                      f"found={os.path.basename(seed)} =====", flush=True)
                verdict, detail = run_cell(
                    root, seed, n, edgeval, beta, bracket=bracket, replicas=replicas,
                    burnin=burnin, nsamp=nsamp, thin=thin, dry_run=dry_run,
                    seeds_dir=seeds_dir, out_dir=os.path.join(out_root, cell))
                rows.append(dict(N=n, edge=edgeval, beta_over_N=bon, beta=beta,
                                 verdict=verdict, detail=detail,
                                 found=os.path.basename(seed)))
                with open(os.path.join(out_root, "summary.csv"), "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                    w.writeheader(); w.writerows(rows)

    print("\n\n================ SWEEP SUMMARY ================")
    for r in rows:
        print(f"{r['N']:>6} ED{r['edge']:>7} b/N={r['beta_over_N']:<6g} "
              f"{r['verdict']:>12}  {r['detail']}")
    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    nskip = sum(1 for r in rows if r["verdict"] == "SKIP-EXISTS")
    print(f"\n{npass} pass, {nskip} pre-existing, {len(rows)-npass-nskip} fail, "
          f"of {len(rows)} cells. CSV: {os.path.join(out_root, 'summary.csv')}")
    return rows


def prune_raised_frontier(edgetok, bon):
    """The known-dead corner: beta/N>=0.008 at raised edge targets."""
    return edgetok in ("5p1043", "5p2043") and bon >= 0.008
