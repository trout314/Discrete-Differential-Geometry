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

# Standard grid (the certified N=1e4 grid, extended to small N). Every axis below
# is a full grid dimension; k (edge-pin stiffness) and hdv (hinge-degree-variance
# coef) default to a single library value but accept a list, so the sweep covers
# the whole N x edge x beta/N x k x hdv product. Distinct k/hdv make distinct
# families (encoded in the .mfd filename), so no collision with the k=2/hdv=0 grid.
N_TIERS = [("1e2", 100), ("178", 178), ("316", 316), ("562", 562), ("1e3", 1000),
           ("1778", 1778), ("3162", 3162), ("5623", 5623)]
EDGES = [("5p0043", 5.0043), ("5p1043", 5.1043), ("5p2043", 5.2043)]
BETA_OVER_N = [0.001, 0.002, 0.004, 0.008, 0.01]
K_VALUES = [2.0]        # edge-pin stiffness num_hinges_coef (default: library k=2)
# HDV coupling is coef/N (equipartition; like beta/N for VDV), so the axis is
# specified as coef/N and raw hdv_coef = (coef/N) * N is set per cell.
HDV_OVER_N = [0.0]      # hinge-degree-variance coef/N (default: unpenalized)
DEFAULT_BRACKET = ["vdv", "edge_deg", "num_facets"]


def founding_seed(root, ntok, edgetok, beta_target, seeds_dir="seeds"):
    """s000 of the existing family at (N, edge) whose beta is nearest beta_target
    (globs seeds_dir live, so newly produced families are candidates too).
    k-agnostic: matches any edge-pin coef, since the founding config is only a
    starting triangulation -- the study objective re-pins at whatever k is set."""
    pat = os.path.join(root, seeds_dir, f"S3_N{ntok}_1e-1_ED{edgetok}_*_s000.mfd")
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


def retry_worthwhile(detail):
    """True if a produce FAIL is RUN-LENGTH-fixable rather than a genuine frontier
    break: an ESS 'run longer' (under-sampled), or VDV-only non-convergence (vdv
    in the bad list but none of the frontier coordinates edge_deg / num_facets /
    hdv). A frontier break (edge_deg/num_facets/hdv qRhat) is not retried."""
    if "run longer" in detail:
        return True
    m = re.search(r"not converged in \[(.*?)\]", detail)
    if not m:
        return False
    bad = m.group(1)
    return "vdv" in bad and not any(x in bad for x in
                                    ("edge_deg", "num_facets", "hdv"))


def run_cell(root, seed, n, edgeval, beta, *, bracket, replicas, burnin, nsamp,
             thin, dry_run, seeds_dir, out_dir, num_hinges_coef=2.0, hdv_coef=0.0,
             vdq_coef=0.0, edq_coef=0.0,
             max_workers=None, max_memory_gb=None, topology="S3"):
    """One grid cell = one `equilibrium_vdv.py --produce` invocation. Optional
    max_workers / max_memory_gb cap the concurrent-chain fan-out (cores / RAM).
    vdq_coef/edq_coef are the RAW per-element fixed-target couplings (the queue
    and grid axes use per-tet scaled values; convert before calling)."""
    cmd = [sys.executable, os.path.join(root, "scripts", "equilibrium_vdv.py"),
           "--produce", "--seed-file", seed, "--n-target", str(n), "--beta", str(beta),
           "--bracket", *bracket, "--hinge-target", str(edgeval),
           "--num-hinges-coef", str(num_hinges_coef), "--hdv-coef", str(hdv_coef),
           "--num-facets-coef", "0.1",
           "--replicas", str(replicas), "--production-burnin", str(burnin),
           "--n-samples", str(nsamp), "--thin", str(thin),
           "--seeds-dir", os.path.join(root, seeds_dir), "--output-dir", out_dir,
           "--topology", topology]
    if vdq_coef:
        cmd += ["--vdq-coef", str(vdq_coef)]
    if edq_coef:
        cmd += ["--edq-coef", str(edq_coef)]
    if max_workers is not None:
        cmd += ["--max-workers", str(max_workers)]
    if max_memory_gb is not None:
        cmd += ["--max-memory-gb", str(max_memory_gb)]
    if dry_run:
        cmd.append("--dry-run")
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    sys.stdout.write(out); sys.stdout.flush()
    return _parse_verdict(out, dry_run)


def sweep(root, *, dry_run, bracket=DEFAULT_BRACKET, replicas, burnin, nsamp, thin,
          n_tiers=N_TIERS, edges=EDGES, beta_over_N=BETA_OVER_N, k_values=K_VALUES,
          hdv_over_N=HDV_OVER_N, prune=None, seeds_dir="seeds",
          out_root="data/grid_sweep", only_n=None, only_edge=None, only_bon=None,
          adaptive=False, retry_burnin=5000, retry_nsamp=1500):
    """Run the N x edge x beta/N x k x hdv grid. `prune(edgetok, bon) -> bool`
    skips cells; only_* (lists) restrict the base axes; k_values is the edge-pin
    axis and hdv_over_N is the HDV-coupling axis in coef/N (raw hdv_coef = coef/N
    * N per cell). With adaptive=True, a cell that FAILS on a run-length-fixable
    cause is retried once at retry_burnin/retry_nsamp. Writes an incremental
    summary CSV under out_root and returns the rows."""
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
                for k in k_values:
                    for hdv_on in hdv_over_N:
                        hdv = hdv_on * n                      # raw coef = (coef/N)*N
                        seed = founding_seed(root, ntok, edgetok, beta, seeds_dir)
                        if not seed:
                            print(f"SKIP N={ntok} ED={edgetok} b/N={bon}: no founding",
                                  flush=True)
                            continue
                        cell = f"N{ntok}_ED{edgetok}_bN{bon:g}_k{k:g}_hdvN{hdv_on:g}"
                        print(f"\n===== {cell}  beta={beta:g}  "
                              f"found={os.path.basename(seed)} =====", flush=True)
                        verdict, detail = run_cell(
                            root, seed, n, edgeval, beta, bracket=bracket,
                            replicas=replicas, burnin=burnin, nsamp=nsamp, thin=thin,
                            dry_run=dry_run, seeds_dir=seeds_dir,
                            out_dir=os.path.join(out_root, cell),
                            num_hinges_coef=k, hdv_coef=hdv)
                        if adaptive and verdict == "FAIL" and retry_worthwhile(detail):
                            print(f"  [adaptive] {cell}: run-length-fixable -> retry "
                                  f"at {retry_burnin}/{retry_nsamp}", flush=True)
                            # Fresh out_dir: a chain whose --save-config already
                            # exists short-circuits (resume), so reusing the first
                            # attempt's staging would re-gate the SAME short run
                            # instead of running longer. _retry gets clean staging.
                            verdict, detail = run_cell(
                                root, seed, n, edgeval, beta, bracket=bracket,
                                replicas=replicas, burnin=retry_burnin,
                                nsamp=retry_nsamp, thin=thin, dry_run=dry_run,
                                seeds_dir=seeds_dir,
                                out_dir=os.path.join(out_root, cell + "_retry"),
                                num_hinges_coef=k, hdv_coef=hdv)
                            detail = "[retry] " + detail
                        rows.append(dict(N=n, edge=edgeval, beta_over_N=bon, beta=beta,
                                         k=k, hdv_over_n=hdv_on, hdv_raw=hdv,
                                         verdict=verdict, detail=detail,
                                         found=os.path.basename(seed)))
                        with open(os.path.join(out_root, "summary.csv"), "w",
                                  newline="") as f:
                            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                            w.writeheader(); w.writerows(rows)

    print("\n\n================ SWEEP SUMMARY ================")
    for r in rows:
        print(f"{r['N']:>6} ED{r['edge']:>7} b/N={r['beta_over_N']:<6g} "
              f"k={r['k']:<4g} hdvN={r['hdv_over_n']:<5g} {r['verdict']:>12}  {r['detail']}")
    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    nskip = sum(1 for r in rows if r["verdict"] == "SKIP-EXISTS")
    print(f"\n{npass} pass, {nskip} pre-existing, {len(rows)-npass-nskip} fail, "
          f"of {len(rows)} cells. CSV: {os.path.join(out_root, 'summary.csv')}")
    return rows


def prune_raised_frontier(edgetok, bon):
    """The known-dead corner: beta/N>=0.008 at raised edge targets."""
    return edgetok in ("5p1043", "5p2043") and bon >= 0.008
