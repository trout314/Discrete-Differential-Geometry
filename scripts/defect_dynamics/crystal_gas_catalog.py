#!/usr/bin/env python3
"""Production gas run + HTML defect catalog, per crystal, at its own c_imp*.

Phase 2 of the library survey. For every crystal that crystal_gas_report.py
certified as holding a gas, this re-runs the cell at production length WITH the
harmonic cocycle attached from sweep 0 -- mandatory, because post-hoc cocycle
regeneration FAILS on these gases (monodromy; see notes/memory/m2-only-gas.md)
and the 3D viewer has no positions without it -- and then renders
defect_catalog.py over the final snapshot.

Output per crystal: data/viz/<crystal>_gas_c<cimp>.final/index.html

Usage:
    python scripts/defect_dynamics/crystal_gas_catalog.py \
        --pick data/crystal_gas/verdicts.json --jobs 4
    python scripts/defect_dynamics/crystal_gas_catalog.py \
        --cell c15:0.5 --cell a15:0.45          # explicit override
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from crystal_gas_scan import E_FLAT, LIBRARY, pin_gap

SCAN = os.path.join(_HERE, "crystal_gas_scan.py")
CATALOG = os.path.join(_HERE, "defect_catalog.py")


def note_for(crystal, cimp):
    """Provenance stamp for the page header (CONVENTIONS sec 2)."""
    fname, mcell = LIBRARY[crystal]
    e_nat, gap, nforced = pin_gap(crystal, mcell)
    sign = "2->3" if gap > 0 else "3->2 / 1->4"
    return (f"host {crystal.upper()} m{mcell} ({fname}); "
            f"action 0.1(f3-f3ref)^2 + 1.0(f1-6f3/e*)^2 + {cimp}*sum_v m(v)^2 "
            f"(no EDQ term, no U(n6), no VDV/HDV); "
            f"e* = 2pi/arccos(1/3) = {E_FLAT:.7f} (FLAT), "
            f"e_native = {e_nat:.7f}, pin gap = {gap:+.2f} f1 "
            f"=> ~{nforced:.0f} forced {sign} moves; "
            f"single chain, cocycle attached at sweep 0; "
            f"PROVISIONAL (uncertified -- no two-sided R-hat)")


def run_one(job):
    crystal, cimp, out, burn, span, chunk, seed = job
    stem = os.path.join(out, f"{crystal}_gas_c{cimp:.2f}")
    t0 = time.time()
    p = subprocess.run(
        [sys.executable, SCAN, "--crystal", crystal, "--cimp", str(cimp),
         "--burn", str(burn), "--span", str(span), "--chunk", str(chunk),
         "--seed", str(seed), "--cocycle", "--out", stem],
        capture_output=True, text=True)
    if p.returncode:
        return crystal, cimp, "RUN-FAIL", (p.stderr or "")[-400:], None
    snap = stem + ".final.mfd"
    if not os.path.exists(snap[:-4] + ".cocycle.npz"):
        return crystal, cimp, "NO-COCYCLE", snap, None
    q = subprocess.run(
        [sys.executable, CATALOG, snap, "--note", note_for(crystal, cimp)],
        capture_output=True, text=True)
    if q.returncode:
        return crystal, cimp, "CATALOG-FAIL", (q.stderr or "")[-400:], None
    idx = os.path.join(_ROOT, "data", "viz",
                       os.path.basename(snap)[:-4], "index.html")
    return (crystal, cimp, f"ok {time.time()-t0:.0f}s",
            (p.stdout or "").strip().splitlines()[-1], idx)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pick", default=os.path.join(
        _ROOT, "data", "crystal_gas", "verdicts.json"),
        help="verdicts.json from crystal_gas_report.py --json")
    ap.add_argument("--cell", action="append", default=[],
                    help="explicit crystal:cimp (repeatable); overrides --pick")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "crystal_gas"))
    ap.add_argument("--burn", type=int, default=3000)
    ap.add_argument("--span", type=int, default=3000)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--seed", type=int, default=770806)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    if args.cell:
        picks = [(c.split(":")[0], float(c.split(":")[1])) for c in args.cell]
    else:
        rows = json.load(open(args.pick))
        best = {}
        for r in sorted(rows, key=lambda r: r["cimp"]):
            if r["verdict"] == "GAS" and r["crystal"] not in best:
                best[r["crystal"]] = r["cimp"]
        picks = sorted(best.items())
    if not picks:
        sys.exit("no crystal certified a gas -- nothing to catalog")

    print(f"{len(picks)} catalogs: " +
          ", ".join(f"{c}@{v:.2f}" for c, v in picks), flush=True)
    jobs = [(c, v, args.out, args.burn, args.span, args.chunk, args.seed)
            for c, v in picks]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        links = []
        for f in as_completed(futs):
            crystal, cimp, status, msg, idx = f.result()
            print(f"  {crystal:6s} c={cimp:<5} {status:14s} {msg[:110]}",
                  flush=True)
            if idx:
                links.append((crystal, cimp, idx))
    print("\n=== catalogs ===")
    for crystal, cimp, idx in sorted(links):
        print(f"  {crystal:6s} c_imp={cimp:.2f}  {idx}")


if __name__ == "__main__":
    main()
