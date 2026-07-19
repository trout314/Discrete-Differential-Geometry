#!/usr/bin/env python3
"""Replica-fleet orchestrator for dope_hold.py chains.

Launches N independent, fully-instrumented dope_hold workers under one
condition, with distinct recorded RNG seeds, per-replica output files,
memory-gated staggered starts, and a manifest for provenance. The method
rule behind it: low-k / aging claims come from replica ensembles, never
from single chains.

Usage (everything after ``--`` is passed to dope_hold verbatim):

    python scripts/replica_hold.py --name ownpin --replicas 8 \
        --out-dir data/replicas --base-seed 1000 -- \
        --structure c15big --dopant-class Z14 --mu 3.0 --lam 1.0 \
        --host-n6 0 4 --sweeps-total 60000 --chunk 500 --snap-every 10

Per replica k the orchestrator appends:
    --seed <base_seed + k>
    --out-csv <out-dir>/<name>_r<k>.csv
    --save-final <out-dir>/<name>_r<k>_final.mfd
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

DOPE_HOLD = os.path.join(_ROOT, "scripts", "dope_hold.py")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--name", required=True, help="condition name (file prefix)")
    ap.add_argument("--replicas", type=int, default=8)
    ap.add_argument("--out-dir", default="data/replicas")
    ap.add_argument("--base-seed", type=int, default=None,
                    help="replica k gets seed base+k (default: time-derived, "
                         "recorded in the manifest)")
    ap.add_argument("--mem-floor-gb", type=float, default=6.0,
                    help="don't launch the next worker while free memory is "
                         "below this")
    ap.add_argument("--stagger", type=float, default=20.0,
                    help="seconds between worker launches (lets RSS ramp)")
    ap.add_argument("--max-parallel", type=int, default=None)
    ap.add_argument("worker_args", nargs=argparse.REMAINDER,
                    help="args after -- go to dope_hold verbatim")
    args = ap.parse_args()

    wargs = args.worker_args
    if wargs and wargs[0] == "--":
        wargs = wargs[1:]
    if not wargs:
        raise SystemExit("no worker args given (put dope_hold args after --)")
    base_seed = (args.base_seed if args.base_seed is not None
                 else int(time.time()) % 10**6)
    maxp = args.max_parallel or args.replicas
    os.makedirs(args.out_dir, exist_ok=True)

    manifest = dict(
        name=args.name, replicas=args.replicas, base_seed=base_seed,
        worker_args=wargs, git=get_git_info(),
        started=datetime.datetime.now().isoformat(timespec="seconds"),
        workers=[])
    mpath = os.path.join(args.out_dir, f"{args.name}_manifest.json")

    procs = {}
    launched = 0
    while launched < args.replicas or procs:
        # reap finished workers
        for k in list(procs):
            ret = procs[k].poll()
            if ret is not None:
                manifest["workers"][k]["exit"] = ret
                manifest["workers"][k]["finished"] = \
                    datetime.datetime.now().isoformat(timespec="seconds")
                print(f"[replica {k}] exited {ret}", flush=True)
                del procs[k]
        # launch next if allowed
        if launched < args.replicas and len(procs) < maxp:
            free = get_free_memory_gb()
            if free >= args.mem_floor_gb:
                k = launched
                seed = base_seed + k
                prefix = os.path.join(args.out_dir, f"{args.name}_r{k}")
                cmd = [sys.executable, DOPE_HOLD, *wargs,
                       "--seed", str(seed),
                       "--out-csv", f"{prefix}.csv",
                       "--save-final", f"{prefix}_final.mfd"]
                log = open(f"{prefix}.log", "w")
                procs[k] = subprocess.Popen(cmd, stdout=log, stderr=log,
                                            cwd=_ROOT)
                manifest["workers"].append(dict(
                    replica=k, seed=seed, pid=procs[k].pid,
                    csv=f"{prefix}.csv", log=f"{prefix}.log"))
                print(f"[replica {k}] launched pid={procs[k].pid} seed={seed} "
                      f"(free {free:.1f} GB)", flush=True)
                launched += 1
                with open(mpath, "w") as f:
                    json.dump(manifest, f, indent=1)
                time.sleep(args.stagger)
                continue
            else:
                print(f"waiting on memory ({free:.1f} GB free "
                      f"< {args.mem_floor_gb})", flush=True)
        time.sleep(10)

    manifest["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=1)
    bad = [w for w in manifest["workers"] if w.get("exit") != 0]
    print(f"done: {args.replicas - len(bad)}/{args.replicas} workers clean"
          + (f"; FAILED: {[w['replica'] for w in bad]}" if bad else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
