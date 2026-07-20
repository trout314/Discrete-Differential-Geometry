#!/usr/bin/env python3
"""R dopant-solubility isotherms: sidecar-tracked dope_hold holds over a
species x chemical-potential grid, pinned at R's own q-bar.

For each dopant species (a minority legal Z-class R already contains) and each
chemical potential mu, hold perfect R at its own curvature under a mu tilt on
that class and equilibrate. The equilibrium EXCESS population n(mu) over the
native crystal count is the solubility isotherm; the analyzer
(r_isotherm.py) fits it and reads off the solubility limit, insertion gap,
and the curvature lever dq/dn.

Every run is a sidecar-tracked dope_hold segment (fresh R start), so the
grid is fully provenanced. Launches are memory-gated and staggered so the
grid coexists with other running fleets.

Usage:
    python scripts/r_solubility.py --species Z14 Z15 Z16 \
        --mu 0 0.5 1 2 3 4 6 8 --replicas 3 --sweeps 20000
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
from seed_utils import get_free_memory_gb, get_git_info

DOPE = os.path.join(_ROOT, "scripts", "dope_hold.py")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", default=["Z14", "Z15", "Z16"],
                    choices=["Z12", "Z14", "Z15", "Z16"])
    ap.add_argument("--mu", nargs="+", type=float,
                    default=[0, 0.5, 1, 2, 3, 4, 6, 8])
    ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--sweeps", type=int, default=20000)
    ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--struct", default="r")
    ap.add_argument("--edge-target", nargs="+", type=float, default=[None],
                    help="edge-degree pin target(s); default (None) = the "
                         "crystal's own q-bar. Set below native to test "
                         "whether a curvature-lowering dopant (Z12) becomes "
                         "soluble when the pin moves in its direction.")
    ap.add_argument("--out-dir", default="data/r_solubility")
    ap.add_argument("--base-seed", type=int, default=30000)
    ap.add_argument("--mem-floor-gb", type=float, default=10.0)
    ap.add_argument("--stagger", type=float, default=12.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    launched = []
    seed = args.base_seed
    for sp in args.species:
        for et in args.edge_target:
            for mu in args.mu:
                mtag = f"{mu:g}".replace(".", "p")
                ettag = "" if et is None else f"_et{et:g}".replace(".", "p")
                for k in range(args.replicas):
                    name = f"r_sol_{sp}{ettag}_mu{mtag}_r{k}"
                    base = f"{args.out_dir}/{name}"
                    # R's native classes are the host (0 2 3 4); dope_hold's
                    # n_dop counts the dopant class -> analyzer subtracts native.
                    cmd = ["python3", DOPE, "--structure", args.struct,
                           "--dopant-class", sp, "--mu", f"{mu:g}",
                           "--lam", "1.0", "--host-n6", "0", "2", "3", "4",
                           "--sweeps-total", str(args.sweeps), "--chunk",
                           str(args.chunk), "--seed", str(seed),
                           "--out-csv", f"{base}.csv",
                           "--save-final", f"{base}_final.mfd"]
                    if et is not None:
                        cmd += ["--edge-target", f"{et:g}"]
                    seed += 1
                    if args.dry_run:
                        if k == 0:
                            print(f"# {sp} et={et} mu={mu:g}: "
                                  + " ".join(cmd[2:]))
                        continue
                    while get_free_memory_gb() < args.mem_floor_gb:
                        print(f"  mem {get_free_memory_gb():.1f}GB<floor; wait",
                              flush=True)
                        time.sleep(30)
                    subprocess.Popen(
                        f"cd {_ROOT} && setsid nohup {' '.join(cmd)} "
                        f"> {base}.log 2>&1 &", shell=True)
                    launched.append(dict(species=sp, mu=mu, edge_target=et,
                                         replica=k, seed=seed-1, name=name))
                    print(f"launched {name} (seed {seed-1})", flush=True)
                    time.sleep(args.stagger)

    if args.dry_run:
        return
    with open(f"{args.out_dir}/manifest.json", "w") as f:
        json.dump(dict(campaign="r_solubility", struct=args.struct,
                       species=args.species, mu=args.mu,
                       replicas=args.replicas, sweeps=args.sweeps,
                       started=datetime.datetime.now().isoformat(
                           timespec="seconds"),
                       git=get_git_info(), launched=launched), f, indent=1)
    print(f"\n{len(launched)} holds launched; "
          f"manifest -> {args.out_dir}/manifest.json")


if __name__ == "__main__":
    main()
