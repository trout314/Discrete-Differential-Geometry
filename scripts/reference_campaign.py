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
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import discrete_differential_geometry as ddg
from seed_utils import get_free_memory_gb, get_git_info

DOPE = os.path.join(_ROOT, "scripts", "dope_hold.py")
FLAT = "5.1042993"            # arccos(1/3) flat edge-degree target
COMMON = ["--chunk", "500", "--snap-every", "10", "--flip-log-mb", "32"]

# Crystal families: (small, large) structure names in tcp_melt.CRYSTALS.
# large=None where no big supercell exists yet.
FAMILIES = {
    "r":     ("r", "rbig"),
    "a15":   ("a15", "a15big"),
    "sigma": ("sigma", None),
    "c15":   ("c15", "c15big"),
}

# Per-family flat-vacuum dopant chemical potential (calibrated per crystal;
# the "check the melting/glass range" step sets this). Only meaningful for
# crystals at or below flat that reach flat by adding six-edges (Z14).
FLAT_VAC_MU = {"r": "1.1", "c15": "3.0"}


def census_tilt(struct):
    """Composition-matched FK-stabilizing tilt for `struct`: per-class
    chemical potentials proportional to the crystal's native Z-class census
    (Z12,Z14,Z15,Z16 at n6 indices 0,2,3,4), scaled so the most abundant
    legal class is -2. Returns a list of 5 stringified floats for --tilt."""
    from tcp_melt import CRYSTALS
    from fk_skeleton import edges_from_facets, vertex_class_census
    m = ddg.Manifold.load(CRYSTALS[struct], 3)
    eu, ed, V = edges_from_facets(m.facets())
    fz, _ = vertex_class_census(eu, ed, V)
    frac = {"Z12": fz["Z12"], "Z14": fz["Z14"], "Z15": fz["Z15"],
            "Z16": fz["Z16"]}
    top = max(frac.values())
    idx = {"Z12": 0, "Z14": 2, "Z15": 3, "Z16": 4}
    tilt = [0.0] * 5
    for cls, i in idx.items():
        tilt[i] = -2.0 * frac[cls] / top if top > 0 else 0.0
    return [f"{x:.3f}" for x in tilt]


def fleets(struct, family):
    fresh = ["--structure", struct, "--dopant-class", "Z14"]
    protect = census_tilt(struct)
    fl = {
        # perfect crystal held at its own q-bar: stability + cocycle reference
        "crystal": [fresh + ["--mu", "0", "--lam", "1.0", "--cocycle",
                             "--sweeps-total", "20000"] + COMMON],
        # crystal dissolving a Z14 gas at its own curvature: own-pin solution
        "ownpin": [fresh + ["--mu", "3.0", "--lam", "1.0", "--cocycle",
                            "--sweeps-total", "30000"] + COMMON],
        # LIGHT melt -> anneal at flat under the COMPOSITION-MATCHED tilt ->
        # hold: the FK glass (deep quench kinetically traps; light melt +
        # census-matched tilt recovers pure56 to ~0.9).
        "glass": [
            fresh + ["--mu", "0", "--lam", "0.4", "--cocycle",
                     "--sweeps-total", "600", "--chunk", "100",
                     "--snap-every", "10", "--flip-log-mb", "32"],
            ["--mu", "0", "--lam", "1.0", "--edge-target", FLAT,
             "--tilt"] + protect + ["--sweeps-total", "40600"] + COMMON,
        ],
    }
    # flat vacuum: drive to flat with a dilute Z14 gas (only where calibrated)
    if family in FLAT_VAC_MU:
        fl["flat_vac"] = [fresh + ["--mu", FLAT_VAC_MU[family], "--lam", "1.0",
                                   "--edge-target", FLAT, "--cocycle",
                                   "--sweeps-total", "30000"] + COMMON]
    return fl


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
    ap.add_argument("--struct", choices=list(FAMILIES), required=True,
                    help="crystal family (r/a15/sigma/c15)")
    ap.add_argument("--tier", choices=["small", "large"], required=True)
    ap.add_argument("--fleets", nargs="+", required=True)
    ap.add_argument("--replicas", type=int, default=8)
    ap.add_argument("--out-dir", default="data/reference")
    ap.add_argument("--base-seed", type=int, default=6000)
    ap.add_argument("--mem-floor-gb", type=float, default=8.0)
    ap.add_argument("--stagger", type=float, default=25.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    struct = FAMILIES[args.struct][0 if args.tier == "small" else 1]
    if struct is None:
        raise SystemExit(f"no {args.tier}-tier crystal for family "
                         f"{args.struct} (generate one first)")
    spec = fleets(struct, args.struct)
    unknown = [f for f in args.fleets if f not in spec]
    if unknown:
        raise SystemExit(f"unknown/unavailable fleets {unknown}; "
                         f"have {list(spec)}")
    out_dir = os.path.join(args.out_dir, args.tier)
    os.makedirs(out_dir, exist_ok=True)

    launched = []
    for fi, fleet in enumerate(args.fleets):
        stages = spec[fleet]
        name = f"{args.struct}_{fleet}_{args.tier}"
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
        campaign="reference_states", family=args.struct, tier=args.tier,
        struct=struct, fleets=args.fleets, replicas=args.replicas,
        flat_target=FLAT, census_tilt=census_tilt(struct),
        base_seed=args.base_seed,
        started=datetime.datetime.now().isoformat(timespec="seconds"),
        git=get_git_info(), launched=launched)
    mpath = os.path.join(out_dir, f"manifest_{args.struct}_{args.tier}.json")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"\n{len(launched)} chains launched; manifest -> {mpath}")


if __name__ == "__main__":
    main()
