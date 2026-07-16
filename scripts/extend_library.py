#!/usr/bin/env python3
"""Fill in the existing seed families at LARGER N (a long-running background job).

For every family already in seeds/, walk the geometric N ladder upward from that
family's current largest tier to --max-n, and at each new tier: GROW the family's
biggest seed up to the target facet count under the family's own objective
(grow_seed.py), then CERTIFY at that N via `equilibrium_vdv.py --produce` (the
multi-chain convergence gate), copying the passing family into seeds/. Founding is
progressive -- each tier grows from the just-produced next-smaller tier.

The family's objective is read from its largest existing seed's .mfd metadata and
SCALED to the new N: beta/N and hdv-coef/N are held constant (the VDVs_/HDVs_
couplings), so a produced tier lands on the SAME family stem. Edge target, edge-pin
k, and num-facets-coef are held fixed.

Resource-capped for a shared box: families run ONE AT A TIME (peak = one produce),
and each produce fans its chains out under --max-workers cores / --max-memory-gb
RAM (equilibrium_vdv auto-throttles concurrent chains to fit). Idempotent and
resumable: a tier whose s000 already exists is skipped, so re-launching resumes.
If a tier fails to certify, the family's higher tiers are skipped (they only mix
slower). Writes an incremental summary CSV under --out-dir.

Legacy raw `_VDV_` families (constant-beta diagonals, not constant-beta/N) are
skipped by default: the scaled naming can't give them a stable stem across N.
Pass --include-legacy to extend them anyway (each tier lands on a scaled stem).

Examples
--------
    # the whole job: fill every family up to N=1e4, 10 cores / 16 GB
    python scripts/extend_library.py --max-n 10000 --max-workers 10 --max-memory-gb 16
    # just one family, to test
    python scripts/extend_library.py --only-stem 1e-1_ED5p0043_2_VDVs_4e-3
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import grid_sweep as G
from seed_utils import encode_float, load_seed_metadata

# Geometric N ladder (10^(k/4)); tokens must match encode_float / existing files.
LADDER = [("1e2", 100), ("178", 178), ("316", 316), ("562", 562), ("1e3", 1000),
          ("1778", 1778), ("3162", 3162), ("5623", 5623), ("1e4", 10000),
          ("17783", 17783), ("31623", 31623), ("56234", 56234), ("1e5", 100000)]


def families(seeds_dir):
    """{stem -> {"tiers": {N: path}, "params": <largest-tier metadata dict>}}."""
    fam = {}
    for p in glob.glob(os.path.join(seeds_dir, "*_s000.mfd")):
        m = re.match(r"S3_N([0-9e]+)_(.*)_s000\.mfd", os.path.basename(p))
        if not m:
            continue
        ntok, stem = m.group(1), m.group(2)
        n = next((nn for tt, nn in LADDER if tt == ntok), None)
        if n is None:
            continue
        fam.setdefault(stem, {"tiers": {}})["tiers"][n] = p
    for stem, d in fam.items():
        big = d["tiers"][max(d["tiers"])]           # largest-tier seed = reference
        d["params"] = load_seed_metadata(big)
    return fam


def objective(md):
    """(edge, k, nfc, beta_over_n, hdv_over_n) from a reference seed's metadata."""
    g = lambda key, dflt=0.0: float(md.get(key, dflt))
    n_ref = g("num_facets_target")
    return dict(edge=g("hinge_degree_target"), k=g("num_hinges_coef"),
                nfc=g("num_facets_coef", 0.1),
                beta_over_n=g("codim3_degree_variance_coef") / n_ref,
                hdv_over_n=g("hinge_degree_variance_coef") / n_ref)


def grow(root, src, target_n, obj, beta, hdv, out, min_free_gb,
         vdq_coef=0.0, edq_coef=0.0):
    """Grow src up to target_n facets under the family objective. Returns out on
    success (or if it already exists), else None. vdq_coef/edq_coef are the RAW
    per-element fixed-target (VDQ/EDQ) couplings, 0 = off."""
    if os.path.exists(out):
        return out
    cmd = [sys.executable, os.path.join(root, "scripts", "grow_seed.py"),
           "--from", src, "--target-facets", str(target_n), "--out", out,
           "--hinge-target", str(obj["edge"]), "--num-facets-coef", str(obj["nfc"]),
           "--num-hinges-coef", str(obj["k"]), "--hdv-coef", str(hdv),
           "--vdv-coef", str(beta), "--min-free-gb", str(min_free_gb)]
    if vdq_coef:
        cmd += ["--vdq-coef", str(vdq_coef)]
    if edq_coef:
        cmd += ["--edq-coef", str(edq_coef)]
    r = subprocess.run(cmd)
    return out if (r.returncode == 0 and os.path.exists(out)) else None


def largest_source(root, seeds_dir, stem, below_n):
    """s000 of the family's largest tier with N < below_n (re-globbed live, so a
    just-produced tier is a candidate)."""
    best, best_n = None, -1
    for tok, n in LADDER:
        if n >= below_n:
            continue
        p = os.path.join(root, seeds_dir, f"S3_N{tok}_{stem}_s000.mfd")
        if os.path.exists(p) and n > best_n:
            best, best_n = p, n
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--max-n", type=int, default=10000, help="Fill tiers up to this N.")
    p.add_argument("--max-workers", type=int, default=10, help="Concurrent-chain core cap.")
    p.add_argument("--max-memory-gb", type=float, default=16.0, help="Concurrent-chain RAM cap.")
    p.add_argument("--replicas", type=int, default=32)
    p.add_argument("--burnin", type=int, default=5000)
    p.add_argument("--n-samples", type=int, default=1500)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--retry-burnin", type=int, default=8000,
                   help="Longer burn-in for the one adaptive retry of a run-length fail.")
    p.add_argument("--retry-n-samples", type=int, default=2500)
    p.add_argument("--grow-min-free-gb", type=float, default=6.0,
                   help="grow_seed aborts below this free RAM.")
    p.add_argument("--only-stem", nargs="+", default=None,
                   help="Restrict to these exact family stems.")
    p.add_argument("--include-legacy", action="store_true",
                   help="Also extend legacy raw _VDV_ families (scaled-renamed per N).")
    p.add_argument("--seeds-dir", default="seeds")
    p.add_argument("--out-dir", default="data/extend_library")
    p.add_argument("--scratch-dir",
                   default=os.path.join(_ROOT, "data", "extend_library", "_grown"))
    args = p.parse_args()

    seeds_abs = args.seeds_dir if os.path.isabs(args.seeds_dir) \
        else os.path.join(_ROOT, args.seeds_dir)
    out_root = args.out_dir if os.path.isabs(args.out_dir) \
        else os.path.join(_ROOT, args.out_dir)
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)

    fam = families(seeds_abs)
    # Build worklist: (stem, [target tiers strictly above current max, up to max-n]).
    work = []
    skipped_legacy, nothing = [], []
    for stem, d in sorted(fam.items()):
        if args.only_stem and stem not in args.only_stem:
            continue
        if "_VDV_" in stem and not args.include_legacy:
            skipped_legacy.append(stem)
            continue
        cur_max = max(d["tiers"])
        targets = [(tok, n) for tok, n in LADDER if n > cur_max and n <= args.max_n]
        if not targets:
            nothing.append(stem)
            continue
        work.append((stem, d, targets))
    # Front-load quick wins: families needing the fewest new tiers first.
    work.sort(key=lambda w: (len(w[2]), w[0]))

    print(f"extend_library: {len(work)} families to fill, "
          f"{len(nothing)} already complete, {len(skipped_legacy)} legacy skipped.")
    print(f"  cores<= {args.max_workers}, RAM<= {args.max_memory_gb} GB, "
          f"replicas={args.replicas}, burnin={args.burnin}, samples={args.n_samples}")
    if skipped_legacy:
        print(f"  legacy raw _VDV_ (use --include-legacy): {', '.join(skipped_legacy)}")
    total_tiers = sum(len(t) for _, _, t in work)
    print(f"  {total_tiers} (family x tier) cells queued.\n", flush=True)

    rows = []

    def record(**kw):
        rows.append(kw)
        with open(os.path.join(out_root, "summary.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    for stem, d, targets in work:
        obj = objective(d["params"])
        print(f"\n########## {stem}  edge={obj['edge']:g} k={obj['k']:g} "
              f"b/N={obj['beta_over_n']:g} hdv/N={obj['hdv_over_n']:g}  "
              f"targets={[t for t, _ in targets]} ##########", flush=True)
        for tok, n in targets:
            dest = os.path.join(seeds_abs, f"S3_N{tok}_{stem}_s000.mfd")
            if os.path.exists(dest):
                print(f"  [{tok}] already present -> skip", flush=True)
                record(stem=stem, N=n, verdict="EXISTS", detail="", grown="")
                continue
            src = largest_source(_ROOT, args.seeds_dir, stem, n)
            if not src:
                print(f"  [{tok}] no smaller source seed -> skip", flush=True)
                record(stem=stem, N=n, verdict="NO-SOURCE", detail="", grown="")
                break
            beta, hdv = obj["beta_over_n"] * n, obj["hdv_over_n"] * n
            grown = os.path.join(args.scratch_dir, f"grow_N{tok}_{stem}.mfd")
            print(f"  [{tok}] grow {os.path.basename(src)} -> {n} facets "
                  f"(beta={beta:g}, hdv={hdv:g})", flush=True)
            g = grow(_ROOT, src, n, obj, beta, hdv, grown, args.grow_min_free_gb)
            if not g:
                print(f"  [{tok}] GROW FAILED -> stop climbing this family", flush=True)
                record(stem=stem, N=n, verdict="GROW-FAIL", detail="", grown="")
                break
            cell = os.path.join(out_root, f"{stem}_N{tok}")
            verdict, detail = G.run_cell(
                _ROOT, g, n, obj["edge"], beta, bracket=G.DEFAULT_BRACKET,
                replicas=args.replicas, burnin=args.burnin, nsamp=args.n_samples,
                thin=args.thin, dry_run=False, seeds_dir=args.seeds_dir,
                out_dir=cell, num_hinges_coef=obj["k"], hdv_coef=hdv,
                max_workers=args.max_workers, max_memory_gb=args.max_memory_gb)
            if verdict == "FAIL" and G.retry_worthwhile(detail):
                print(f"  [{tok}] run-length-fixable -> retry at "
                      f"{args.retry_burnin}/{args.retry_n_samples}", flush=True)
                # Fresh out_dir: a chain whose --save-config already exists short-
                # circuits (resume), so the retry must use clean staging to
                # actually run longer rather than re-gate the same short run.
                verdict, detail = G.run_cell(
                    _ROOT, g, n, obj["edge"], beta, bracket=G.DEFAULT_BRACKET,
                    replicas=args.replicas, burnin=args.retry_burnin,
                    nsamp=args.retry_n_samples, thin=args.thin, dry_run=False,
                    seeds_dir=args.seeds_dir, out_dir=cell + "_retry",
                    num_hinges_coef=obj["k"], hdv_coef=hdv,
                    max_workers=args.max_workers, max_memory_gb=args.max_memory_gb)
                detail = "[retry] " + detail
            record(stem=stem, N=n, verdict=verdict, detail=detail,
                   grown=os.path.basename(g))
            if verdict not in ("PASS", "SKIP-EXISTS"):
                print(f"  [{tok}] {verdict} -> stop climbing this family", flush=True)
                break
            try:                                    # free the grown scratch config
                os.remove(g)
            except OSError:
                pass

    print("\n\n================ EXTEND SUMMARY ================")
    for r in rows:
        print(f"  N{r['N']:<6} {r['verdict']:>10}  {r['stem']}  {r['detail']}")
    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    nex = sum(1 for r in rows if r["verdict"] in ("EXISTS", "SKIP-EXISTS"))
    nfail = len(rows) - npass - nex
    print(f"\n{npass} newly certified, {nex} pre-existing, {nfail} fail/blocked, "
          f"of {len(rows)} cells. CSV: {os.path.join(out_root, 'summary.csv')}")


if __name__ == "__main__":
    main()
