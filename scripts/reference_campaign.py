#!/usr/bin/env python3
"""Reference-state campaign driver: certified, lineage-recorded TCP fleets.

Produces the project's "load-bearing reference states" as replica fleets,
every run a sidecar-tracked dope_hold SEGMENT (see dope_hold.py: the
ensemble is frozen into <base>.run.json and continuations inherit it). A
reference state is defined by a PIPELINE of one or more stages; stage k>0
is a dope_hold ``--from`` continuation of stage k-1's saved final, so the
full melt -> anneal -> hold lineage is explicit rather than inferred.

The whole pipeline runs through dope_hold (no tcp_melt): a "melt" stage is
just a low-lambda hold (weak FK-legality pressure => thermal disorder), an
"anneal" stage a lambda=1 hold at the flat edge target under the protective
tilt. Cocycles initialize on the fresh crystal stage and are auto-carried
across continuations (closed under all Pachner moves), so even the glass
vacuum keeps exact T3 winding classes.

Each replica runs its stages sequentially in one detached shell chain;
replicas are launched memory-gated and staggered. A manifest records the
resolved pipeline, seeds, tiers, pids, and git info.

Usage:
    python scripts/reference_campaign.py --tier small --fleets r_crystal \
        r_flat_vac r_ownpin r_glass --replicas 8 --out-dir data/reference
    python scripts/reference_campaign.py --tier large --fleets r_crystal \
        r_flat_vac r_glass --replicas 4         # after small validates
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
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
from seed_utils import get_free_memory_gb, get_git_info

DOPE = os.path.join(_ROOT, "scripts", "dope_hold.py")
FLAT = "5.1042993"            # arccos(1/3) flat edge-degree target
PROTECT = ["-1", "0", "-4", "0", "-1"]   # FK-stabilizing tilt (stage2b vacuum)

# Structure name per tier. Small = m3 R (V=4293); large = m4 R (V=10176).
STRUCT = {"small": "r", "large": "rbig"}

# A stage is a list of dope_hold flags (WITHOUT --seed/--out-csv/--save-final/
# --from, which the driver fills in). Stage 0 is a fresh --structure start;
# later stages are continuations (the driver inserts --from <prev final> and,
# if any physics flag is present, --override).
COMMON = ["--chunk", "500", "--snap-every", "10", "--flip-log-mb", "32"]


def fleets(struct):
    fresh = ["--structure", struct, "--dopant-class", "Z14"]
    return {
        # perfect R held at its own q-bar: stability baseline + cocycle ref
        "r_crystal": [fresh + ["--mu", "0", "--lam", "1.0", "--cocycle",
                               "--sweeps-total", "20000"] + COMMON],
        # R + dilute Z14 gas driven to flat: the tunable flat vacuum
        "r_flat_vac": [fresh + ["--mu", "1.1", "--lam", "1.0",
                                "--edge-target", FLAT, "--cocycle",
                                "--sweeps-total", "30000"] + COMMON],
        # R dissolving a Z14 gas at its own curvature: own-pin solution
        "r_ownpin": [fresh + ["--mu", "3.0", "--lam", "1.0", "--cocycle",
                              "--sweeps-total", "30000"] + COMMON],
        # quench (lam0.4/4000sw -> pure56~0.2, positionally amorphous) ->
        # anneal at flat under the protective tilt -> hold: the FK glass.
        # A real glass-former recipe: melt past crystalline memory, then
        # recover FK-legality (not position) by annealing.
        "r_glass": [
            fresh + ["--mu", "0", "--lam", "0.4", "--cocycle",
                     "--sweeps-total", "4000"] + COMMON,
            ["--mu", "0", "--lam", "1.0", "--edge-target", FLAT,
             "--tilt"] + PROTECT + ["--sweeps-total", "44000"] + COMMON,
        ],
    }


def build_chain(fleet, stages, k, seed, out_dir, name):
    """Return a bash command string running replica k's stages in sequence."""
    cmds = []
    prev_final = None
    for si, stage in enumerate(stages):
        base = f"{out_dir}/{name}_r{k}" if len(stages) == 1 \
            else f"{out_dir}/{name}_s{si}_r{k}"
        csv = f"{base}.csv"
        final = f"{base}_final.mfd"
        seed_i = seed + 1000 * si          # distinct seed per stage
        args = list(stage)
        if si > 0:
            args = ["--from", prev_final, "--override"] + args
        args += ["--seed", str(seed_i), "--out-csv", csv,
                 "--save-final", final]
        cmds.append("python3 " + DOPE + " " + " ".join(args))
        prev_final = final
    # '&&' so a stage only runs if the previous one succeeded
    return " && ".join(cmds)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tier", choices=list(STRUCT), required=True)
    ap.add_argument("--fleets", nargs="+", required=True)
    ap.add_argument("--replicas", type=int, default=8)
    ap.add_argument("--out-dir", default="data/reference")
    ap.add_argument("--base-seed", type=int, default=6000)
    ap.add_argument("--mem-floor-gb", type=float, default=8.0)
    ap.add_argument("--stagger", type=float, default=25.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    struct = STRUCT[args.tier]
    spec = fleets(struct)
    unknown = [f for f in args.fleets if f not in spec]
    if unknown:
        raise SystemExit(f"unknown fleets {unknown}; have {list(spec)}")
    out_dir = os.path.join(args.out_dir, args.tier)
    os.makedirs(out_dir, exist_ok=True)

    launched = []
    for fi, fleet in enumerate(args.fleets):
        stages = spec[fleet]
        name = f"{fleet}_{args.tier}"
        for k in range(args.replicas):
            seed = args.base_seed + 100 * fi + k
            chain = build_chain(fleet, stages, k, seed, out_dir, name)
            log = f"{out_dir}/{name}_r{k}.log"
            if args.dry_run:
                print(f"# {name}_r{k} (seed {seed}, {len(stages)} stage/s)")
                print(chain + "\n")
                continue
            while get_free_memory_gb() < args.mem_floor_gb:
                print(f"  mem {get_free_memory_gb():.1f}GB < floor; waiting…",
                      flush=True)
                time.sleep(30)
            proc = subprocess.Popen(
                f"cd {_ROOT} && setsid nohup bash -c {json.dumps(chain)} "
                f"> {log} 2>&1 &", shell=True)
            launched.append(dict(fleet=name, replica=k, seed=seed,
                                 stages=len(stages), log=log))
            print(f"launched {name}_r{k} (seed {seed}, {len(stages)} stages)",
                  flush=True)
            time.sleep(args.stagger)

    if args.dry_run:
        return
    manifest = dict(
        campaign="reference_states", tier=args.tier, struct=struct,
        fleets=args.fleets, replicas=args.replicas, flat_target=FLAT,
        protective_tilt=PROTECT, base_seed=args.base_seed,
        started=datetime.datetime.now().isoformat(timespec="seconds"),
        git=get_git_info(), launched=launched)
    mpath = os.path.join(out_dir, f"manifest_{args.tier}.json")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"\n{len(launched)} chains launched; manifest -> {mpath}")


if __name__ == "__main__":
    main()
