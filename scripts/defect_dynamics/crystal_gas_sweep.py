#!/usr/bin/env python3
"""Drive crystal_gas_scan.py over the whole crystal library x c_imp grid.

One worker process per cell (the D loader takes a directory lock, so the
shared-library rebuild is serialised safely). Cells already on disk are
skipped unless --force.

Usage:
    python scripts/defect_dynamics/crystal_gas_sweep.py --jobs 6 \
        --out data/crystal_gas
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from crystal_gas_scan import LIBRARY

CIMPS = [0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0]
SCAN = os.path.join(_HERE, "crystal_gas_scan.py")


def run_cell(job):
    crystal, cimp, out, burn, span, chunk, seed, cocycle = job
    stem = os.path.join(out, f"{crystal}_c{cimp:.2f}")
    cmd = [sys.executable, SCAN, "--crystal", crystal, "--cimp", str(cimp),
           "--burn", str(burn), "--span", str(span), "--chunk", str(chunk),
           "--seed", str(seed), "--out", stem]
    if cocycle:
        cmd.append("--cocycle")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    tail = (p.stdout or p.stderr).strip().splitlines()
    return (crystal, cimp, p.returncode, round(time.time() - t0, 1),
            tail[-1] if tail else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "crystal_gas"))
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--burn", type=int, default=1000)
    ap.add_argument("--span", type=int, default=2000)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20260806)
    ap.add_argument("--cimps", type=float, nargs="*", default=CIMPS)
    ap.add_argument("--crystals", nargs="*", default=sorted(LIBRARY))
    ap.add_argument("--cocycle", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    jobs = []
    for crystal in args.crystals:
        for cimp in args.cimps:
            stem = os.path.join(args.out, f"{crystal}_c{cimp:.2f}")
            if not args.force and os.path.exists(stem + ".final.mfd"):
                continue
            jobs.append((crystal, cimp, args.out, args.burn, args.span,
                         args.chunk, args.seed, args.cocycle))
    print(f"{len(jobs)} cells, {args.jobs} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_cell, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            crystal, cimp, rc, dt, msg = f.result()
            flag = "ok " if rc == 0 else "FAIL"
            print(f"[{i}/{len(jobs)} {time.time()-t0:6.0f}s] {flag} "
                  f"{crystal:6s} c={cimp:<5} {dt:6.1f}s  {msg[:120]}",
                  flush=True)
    print(f"sweep done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
