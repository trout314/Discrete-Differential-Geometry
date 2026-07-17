#!/usr/bin/env python3
"""Quench / relaxation experiments: start from certified equilibrium samples of one
objective and run the sampler under a DIFFERENT (target) objective, recording the
full relaxation trajectory from t=0 (no burn-in skip, no warmup dispersal).

Two-tier recording (see notes/sampler-dynamics-matter-coupling-plan.md):
  Tier 1 (dense, every `--thin` sweeps incl. t=0): observables + degree histograms
          -> chain_NNN.csv (chain,sample,sweeps,num_facets,vdv,edge_deg,hdv,energy,
             accept) + chain_NNN.hist.npz (vhist,ehist). Global/temporal dynamics.
  Tier 2 (sparse, ~`--snapshots` LOG-spaced sweeps): full triangulation .mfd configs
          -> chain_NNN.snap_<sweeps>.mfd. Spatial dynamics (transport, S(k,t)) +
             post-hoc recompute of any observable.
Plus manifest.json: initial + target objectives, IC paths, snapshot schedule, commit.

Each --init file is one independent chain (use the 32 certified _s0NN replicas of a
family as 32 independent initial conditions).

    python scripts/quench.py \
        --init 'seeds/S3_N1e4_1e-1_ED5p1043_2_VDVs_2e-3_s0*.mfd' \
        --n-target 10000 --hinge-target 5.1043 --num-hinges-coef 2 \
        --vdv-coef 80 --hdv-coef 0 --num-facets-coef 0.1 \
        --sweeps 20000 --thin 20 --snapshots 16 \
        --out-dir data/quench/ED5p1043_2_VDVs2e3_to_8e3_N1e4
"""
import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
from seed_utils import load_seed_metadata, get_git_info  # noqa: E402


def snapshot_schedule(total, n_snaps):
    """t=0 plus ~n_snaps log-spaced sweep counts in (0, total]."""
    if total <= 1 or n_snaps <= 1:
        return [0, total]
    pts = np.unique(np.round(np.geomspace(1, total, n_snaps)).astype(int))
    return sorted(set([0] + [int(x) for x in pts] + [total]))


def run_chain(init_path, tgt, total, thin, snaps, out_dir, idx):
    """One relaxation trajectory from init_path under the target objective tgt."""
    from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
    m = Manifold.load(init_path, tgt["dim"])
    params = SamplerParams(
        num_facets_target=tgt["n"], hinge_degree_target=tgt["edge"],
        num_facets_coef=tgt["nfc"], num_hinges_coef=tgt["k"],
        hinge_degree_variance_coef=tgt["hdv"], codim3_degree_variance_coef=tgt["beta"])
    s = ManifoldSampler(m, params)
    v = s.manifold

    rows, vh, eh = [], [], []
    snaps = sorted(set(int(x) for x in snaps))
    si = 0

    def take_snaps(cum):
        nonlocal si
        while si < len(snaps) and cum >= snaps[si]:
            v.save(os.path.join(out_dir, f"chain_{idx:03d}.snap_{cum}.mfd"))
            si += 1

    def record(sample, cum, acc):
        rows.append([idx, sample, cum, v.num_facets, v.degree_variance(0),
                     v.mean_degree(1), v.degree_variance(1),
                     s.current_objective, acc])
        vh.append(np.asarray(v.degree_histogram(0)))
        eh.append(np.asarray(v.degree_histogram(1)))

    # t = 0: the IC evaluated under the TARGET objective (quench starting point)
    record(0, 0, 0.0)
    take_snaps(0)
    cum, sample = 0, 1
    while cum < total:
        step = min(thin, total - cum)
        s.reset_stats(); s.run(sweeps=step); cum += step
        stt = s.get_stats()
        acc = stt.total_accepted / stt.total_tried if stt.total_tried else 0.0
        record(sample, cum, acc)
        take_snaps(cum)
        sample += 1

    base = os.path.join(out_dir, f"chain_{idx:03d}")
    with open(base + ".csv", "w") as f:
        f.write("chain,sample,sweeps,num_facets,vdv,edge_deg,hdv,energy,accept\n")
        for r in rows:
            f.write(",".join(map(repr, r)) + "\n")
    w = max(a.shape[0] for a in vh); e = max(a.shape[0] for a in eh)
    V = np.zeros((len(vh), w), np.int64); E = np.zeros((len(eh), e), np.int64)
    for i, (a, b) in enumerate(zip(vh, eh)):
        V[i, :a.shape[0]] = a; E[i, :b.shape[0]] = b
    np.savez_compressed(base + ".hist.npz", vhist=V, ehist=E)
    return idx, len(rows), si


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--init", required=True,
                   help="Glob of initial-condition .mfd files (one chain each).")
    p.add_argument("--dim", type=int, default=3)
    p.add_argument("--topology", default="S3")
    # target objective to quench TO
    p.add_argument("--n-target", type=int, required=True)
    p.add_argument("--hinge-target", type=float, required=True)
    p.add_argument("--num-facets-coef", type=float, default=0.1)
    p.add_argument("--num-hinges-coef", type=float, required=True)
    p.add_argument("--vdv-coef", type=float, default=0.0, help="codim3 (VDV) coef (raw).")
    p.add_argument("--hdv-coef", type=float, default=0.0, help="hinge (HDV) coef (raw).")
    p.add_argument("--sweeps", type=int, required=True, help="Total relaxation sweeps.")
    p.add_argument("--thin", type=int, default=20, help="Record observables every N sweeps.")
    p.add_argument("--snapshots", type=int, default=16, help="Number of log-spaced full snapshots.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--max-workers", type=int, default=8)
    args = p.parse_args()

    inits = sorted(glob.glob(args.init))
    if not inits:
        sys.exit(f"no init files matched: {args.init}")
    os.makedirs(args.out_dir, exist_ok=True)

    tgt = dict(dim=args.dim, n=args.n_target, edge=args.hinge_target,
               nfc=args.num_facets_coef, k=args.num_hinges_coef,
               hdv=args.hdv_coef, beta=args.vdv_coef)
    snaps = snapshot_schedule(args.sweeps, args.snapshots)

    # manifest / provenance
    ic_md = load_seed_metadata(inits[0])
    commit, dirty = get_git_info()
    manifest = dict(
        topology=args.topology, dim=args.dim, n_chains=len(inits),
        total_sweeps=args.sweeps, thin=args.thin, snapshot_sweeps=snaps,
        target_objective=tgt,
        initial_condition=dict(
            paths=[os.path.basename(x) for x in inits],
            objective_from_metadata={k: ic_md.get(k) for k in (
                "num_facets_target", "num_facets_coef", "hinge_degree_target",
                "num_hinges_coef", "hinge_degree_variance_coef",
                "codim3_degree_variance_coef")}),
        git_commit=commit, git_dirty=dirty)
    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)

    print(f"quench: {len(inits)} chains, target n={args.n_target} edge={args.hinge_target} "
          f"k={args.num_hinges_coef} vdv={args.vdv_coef} hdv={args.hdv_coef}; "
          f"{args.sweeps} sweeps thin {args.thin}; snapshots at {snaps}", flush=True)
    print(f"  IC objective: {manifest['initial_condition']['objective_from_metadata']}", flush=True)

    done = 0
    with ProcessPoolExecutor(max_workers=args.max_workers) as ex:
        futs = {ex.submit(run_chain, ip, tgt, args.sweeps, args.thin, snaps,
                          args.out_dir, i): i for i, ip in enumerate(inits)}
        for fu in as_completed(futs):
            idx, nrec, nsnap = fu.result()
            done += 1
            print(f"  chain {idx:03d} done: {nrec} samples, {nsnap} snapshots "
                  f"({done}/{len(inits)})", flush=True)
    print(f"quench complete -> {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
