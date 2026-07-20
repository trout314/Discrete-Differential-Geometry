#!/usr/bin/env python3
"""Dope-and-hold: time series of dopant population in a TCP crystal.

Holds a crystal (or a continued state) at a fixed physics ENSEMBLE -- volume
pin, edge-degree pin, EDQ coupling lambda, n6 tilt field -- recording the Z
census, dopant count, dopant-cluster statistics, and the disclination-network
census every chunk.

Segments and lineage
--------------------
Every run is a SEGMENT. Its ensemble is resolved by strict precedence

    explicit CLI flag  >  parent sidecar  >  derived from initial state
                                              (fresh starts only)

and frozen into a sidecar ``<base>.run.json`` written at segment start, at
every snapshot, and at the final save (alongside .mfd/.cocycle.npz). The
sidecar records the resolved ensemble + its hash, absolute sweeps done, the
geometry-ledger clock (flip-stream continuity), seed, parent, explicit
overrides, and git rev. Continuations NEVER re-derive physics from the
loaded state -- that rule only applies to fresh starts.

* Fresh start:   --structure <name> (+ --mu/--tilt ...). Unset targets derive
                 from the loaded crystal (own-pin semantics) and are frozen.
* Continuation:  --from <snap.mfd>. Ensemble inherited from the sidecar; any
                 physics flag passed alongside is an ensemble CHANGE and
                 requires --override (the diff is printed and recorded).
                 A sibling .cocycle.npz is picked up automatically.
* Bootstrap:     --from on a pre-sidecar snapshot requires the FULL explicit
                 physics set (segment start parsed from the _snapNNN name).

Sweep numbers are ABSOLUTE (a continuation's CSV starts where the parent
stopped); the CSV carries an ``# ensemble=<hash> segment_start=<n>`` comment
so splice tools can verify segments belong to the same ensemble.

Usage:
    python scripts/dope_hold.py --structure a15big --dopant-class Z16 \
        --mu 2.5 --sweeps-total 12000 --chunk 250 \
        --out-csv data/dope_hold/a15big_z16_mu2.5.csv

    python scripts/dope_hold.py --from data/x_snap30000.mfd --seed 7 \
        --sweeps-total 60000 --out-csv data/x_seg2.csv
"""
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from tcp_melt import CRYSTALS
from dopant_pairs import CLASS_N6, vertex_classes, cluster_census
from fk_skeleton import edges_from_facets, vertex_class_census

COLS = ["sweeps", "acc", "n_dop", "n_companion", "pure56", "mean_edeg",
        "fZ12", "fZ14", "fZ15", "fZ16",
        "n_clusters", "mean_sz", "max_sz", "n_singletons",
        "n_joint", "j_clusters", "j_mean", "j_max", "j_frac_largest",
        # disclination-network census (D-side; see disclination.py)
        "n_six", "net_comps", "net_giant_frac", "cycle_rank",
        "n_seg", "mean_seg", "n_pure_loops", "net_fray",
        "e_dop_dop", "e_dop_host", "e_host_host", "e_imp_any"]
COMPANION = {"Z16": 3, "Z14": 3}   # n6 of the change-species (Z15) for both

#: argparse dests that define the physics ensemble (all default to None so
#: "explicitly given" is distinguishable from "defaulted")
PHYSICS_FLAGS = ["dopant_class", "mu", "lam", "zleg_scale", "cimp_scale",
                 "edge_target", "facet_target", "tilt", "host_n6"]
FRESH_DEFAULTS = dict(lam=1.0, zleg_scale=0.3, cimp_scale=0.3)
NUM_FACETS_COEF = 0.1
NUM_HINGES_COEF = 2.0


def ensemble_hash(ens):
    js = json.dumps(ens, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(js.encode()).hexdigest()[:12]


def git_rev():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=_ROOT, text=True).strip()
    except Exception:
        return None


def build_tilt(mu, dopant_class):
    tilt = [0.0] * 5
    tilt[CLASS_N6[dopant_class]] = -mu
    return tilt


def resolve_ensemble(args, parent, state_qbar, state_facets):
    """Resolve the physics ensemble. Returns (ensemble, overrides_diff)."""
    explicit = {k: getattr(args, k) for k in PHYSICS_FLAGS
                if getattr(args, k) is not None}
    if parent is not None:                        # continuation w/ sidecar
        ens = dict(parent["ensemble"])
        if explicit and not args.override:
            raise SystemExit(
                "ensemble change on --from requires --override; "
                f"explicit physics flags: {sorted(explicit)}")
        diff = {}
        for k, v in explicit.items():
            if ens.get(k) != v:
                diff[k] = [ens.get(k), v]
            ens[k] = v
        if "tilt" not in explicit and ("mu" in explicit
                                       or "dopant_class" in explicit):
            ens["tilt"] = build_tilt(ens["mu"], ens["dopant_class"])
            diff.setdefault("tilt", [parent["ensemble"].get("tilt"),
                                     ens["tilt"]])
    elif args.from_snap:                          # bootstrap (no sidecar)
        need = ["dopant_class", "lam", "zleg_scale", "cimp_scale",
                "edge_target", "facet_target"]
        missing = [k for k in need if k not in explicit]
        if "mu" not in explicit and "tilt" not in explicit:
            missing.append("mu (or tilt)")
        if missing:
            raise SystemExit(
                "BOOTSTRAP continuation (no parent sidecar): the full "
                f"explicit physics set is required; missing: {missing}")
        print("WARNING: bootstrap continuation -- no parent sidecar; "
              "trusting the explicit ensemble on the command line",
              flush=True)
        ens = dict(explicit)
        if "tilt" not in ens:
            ens["tilt"] = build_tilt(ens["mu"], ens["dopant_class"])
        diff = {}
    else:                                         # fresh start
        ens = dict(FRESH_DEFAULTS)
        ens.update(explicit)
        if "dopant_class" not in ens:
            raise SystemExit("--dopant-class is required for a fresh start")
        if "mu" not in ens and "tilt" not in ens:
            raise SystemExit("--mu or --tilt is required for a fresh start")
        # own-pin semantics: unset targets derive from the initial state,
        # then freeze into the sidecar (never re-derived on continuation)
        ens.setdefault("edge_target", state_qbar)
        ens.setdefault("facet_target", state_facets)
        if "tilt" not in ens:
            ens["tilt"] = build_tilt(ens["mu"], ens["dopant_class"])
        diff = {}
    ens.setdefault("mu", None)
    ens.setdefault("host_n6", None)
    ens["num_facets_coef"] = NUM_FACETS_COEF
    ens["num_hinges_coef"] = NUM_HINGES_COEF
    ens["edge_target_coef"] = ens["lam"] * ens["edge_target"] / 6.0
    return ens, diff


def write_sidecar(path, ens, sweeps, start, clock, args, parent_path, diff):
    with open(path, "w") as f:
        json.dump({"ensemble": ens, "ensemble_hash": ensemble_hash(ens),
                   "sweeps": sweeps, "segment_start": start,
                   "ledger_clock": clock, "seed": args.seed,
                   "parent": parent_path, "overrides": diff,
                   "structure": args.structure, "git": git_rev(),
                   "argv": sys.argv[1:]}, f, indent=1)


def census_cols(view, host_n6):
    c = view.disclination_census(host_classes=host_n6)
    return dict(
        n_six=c["n_six_edges"], net_comps=c["n_components"],
        net_giant_frac=round(c["giant_frac"], 4), cycle_rank=c["cycle_rank"],
        n_seg=c["n_segments"], mean_seg=round(c["mean_seg_len"], 3),
        n_pure_loops=c["n_pure_loops"], net_fray=c["n_fray_verts"],
        e_dop_dop=c["e_dop_dop"], e_dop_host=c["e_dop_host"],
        e_host_host=c["e_host_host"], e_imp_any=c["e_imp_any"])


def snapshot_row(view, dop_n6, host_n6=None):
    facets = np.asarray(view.facets())
    n6, imp, adj = vertex_classes(facets)
    eu, edeg, V = edges_from_facets(facets)
    fz, _ = vertex_class_census(eu, edeg, V)
    dop = np.where((n6 == dop_n6) & (imp == 0))[0]
    comp = int(np.sum((n6 == 3) & (imp == 0))) if dop_n6 != 3 else 0
    cc = cluster_census(dop, adj) if len(dop) else {}
    ncl = sum(cc.values())
    sizes = [sz for (sz, ne), c in cc.items() for _ in range(c)]
    row = dict(
        n_dop=len(dop), n_companion=comp, pure56=1.0 - fz["impure"],
        mean_edeg=float(edeg.mean()),
        fZ12=fz["Z12"], fZ14=fz["Z14"], fZ15=fz["Z15"], fZ16=fz["Z16"],
        n_clusters=ncl, mean_sz=(len(dop) / ncl if ncl else 0.0),
        max_sz=(max(sizes) if sizes else 0),
        n_singletons=sum(c for (sz, ne), c in cc.items() if sz == 1),
        n_joint=0, j_clusters=0, j_mean=0.0, j_max=0, j_frac_largest=0.0)
    if host_n6 is not None:
        joint = np.where((imp == 0) & np.isin(n6, [0, 2, 3, 4])
                         & ~np.isin(n6, host_n6))[0]
        jc = cluster_census(joint, adj) if len(joint) else {}
        jn = sum(jc.values())
        jsz = [sz for (sz, ne), c in jc.items() for _ in range(c)]
        row.update(n_joint=len(joint), j_clusters=jn,
                   j_mean=(len(joint) / jn if jn else 0.0),
                   j_max=(max(jsz) if jsz else 0),
                   j_frac_largest=(max(jsz) / max(1, len(joint)) if jsz else 0.0))
    row.update(census_cols(view, host_n6))
    return row


def ledger_clock(s):
    try:
        return int(s.tet_stats()["clock"])
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--structure", choices=list(CRYSTALS),
                     help="fresh start from this perfect crystal")
    src.add_argument("--from", dest="from_snap", metavar="SNAP.mfd",
                     help="continue from a snapshot; ensemble inherited from "
                          "its .run.json sidecar (physics flags then require "
                          "--override)")
    # --- physics flags (None sentinels; see resolve_ensemble) ---
    ap.add_argument("--dopant-class", dest="dopant_class",
                    choices=list(CLASS_N6), default=None)
    ap.add_argument("--mu", type=float, default=None)
    ap.add_argument("--lam", type=float, default=None)
    ap.add_argument("--zleg-scale", dest="zleg_scale", type=float, default=None)
    ap.add_argument("--cimp-scale", dest="cimp_scale", type=float, default=None)
    ap.add_argument("--edge-target", dest="edge_target", type=float,
                    default=None,
                    help="hinge_degree_target (fresh default: the initial "
                         "state's own mean edge degree, frozen at start)")
    ap.add_argument("--facet-target", dest="facet_target", type=int,
                    default=None,
                    help="num_facets_target (fresh default: the initial "
                         "state's facet count, frozen at start)")
    ap.add_argument("--tilt", type=float, nargs=5, default=None,
                    help="ABSOLUTE 5-component tilt override (index=n6); "
                         "replaces the single-class --mu construction")
    ap.add_argument("--host-n6", dest="host_n6", type=int, nargs="+",
                    default=None,
                    help="native n6 classes of the host (C15: 0 4); enables "
                         "joint defect-complex tracking columns")
    ap.add_argument("--override", action="store_true",
                    help="acknowledge an ensemble CHANGE on --from (the diff "
                         "is printed and recorded in the sidecar)")
    # --- bookkeeping flags (freely changeable between segments) ---
    ap.add_argument("--sweeps-total", type=int, default=12000,
                    help="ABSOLUTE end sweep (continuations run "
                         "sweeps-total minus segment start)")
    ap.add_argument("--chunk", type=int, default=250)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--save-final", default=None)
    ap.add_argument("--snap-every", type=int, default=0,
                    help="save a snapshot (.mfd + .cocycle.npz + .run.json) "
                         "every K chunks (0 = off)")
    ap.add_argument("--seed", type=int, default=None,
                    help="D-side RNG seed (required on --from; must differ "
                         "from the parent segment's)")
    ap.add_argument("--cocycle", action="store_true",
                    help="fresh crystal starts only: track T3 winding "
                         "cocycles from the structure's reference coords "
                         "(continuations pick up the sibling .cocycle.npz "
                         "automatically)")
    ap.add_argument("--flip-log-mb", type=float, default=0.0,
                    help="enable the six-edge flip stream with this buffer; "
                         "records append to <out-csv-base>_flips.bin "
                         "(dtype: disclination.SIX_FLIP_DTYPE)")
    args = ap.parse_args()

    # --- locate initial state, parent sidecar, segment start ---
    parent = None
    start = 0
    if args.from_snap:
        init_path = args.from_snap
        base = re.sub(r"\.mfd$", "", init_path)
        sidecar = base + ".run.json"
        if os.path.exists(sidecar):
            with open(sidecar) as f:
                parent = json.load(f)
            start = int(parent["sweeps"])
        else:
            m = re.search(r"_snap(\d+)$", base)
            if not m:
                raise SystemExit("no sidecar and cannot parse segment start "
                                 f"from filename: {init_path}")
            start = int(m.group(1))
        if args.seed is None:
            raise SystemExit("--seed is required on --from")
        if parent is not None and args.seed == parent.get("seed"):
            raise SystemExit(f"--seed {args.seed} equals the parent "
                             "segment's seed; choose a fresh one")
        if args.cocycle:
            raise SystemExit("--cocycle is for fresh crystal starts; "
                             "continuations auto-load the sibling "
                             ".cocycle.npz")
    else:
        init_path = os.path.join(_ROOT, CRYSTALS[args.structure])

    if args.seed is not None:
        ddg.set_random_seed(args.seed)
    if args.sweeps_total <= start:
        raise SystemExit(f"--sweeps-total {args.sweeps_total} <= segment "
                         f"start {start}: nothing to do")

    m = ddg.Manifold.load(init_path, 3)
    eu, edeg, V = edges_from_facets(m.facets())
    ens, diff = resolve_ensemble(args, parent, float(edeg.mean()),
                                 m.num_facets)
    h = ensemble_hash(ens)
    mode = ("continuation" if parent is not None
            else "bootstrap" if args.from_snap else "fresh")
    print(f"[{mode}] ensemble {h}  start={start}  end={args.sweeps_total}\n"
          f"  facet_target={ens['facet_target']} "
          f"edge_target={ens['edge_target']:.7f} lam={ens['lam']} "
          f"tilt={ens['tilt']} zleg={ens['zleg_scale']} "
          f"cimp={ens['cimp_scale']} host_n6={ens['host_n6']}", flush=True)
    for k, (old, new) in diff.items():
        print(f"  OVERRIDE {k}: {old} -> {new}", flush=True)

    params = ddg.SamplerParams(
        num_facets_target=ens["facet_target"],
        num_facets_coef=ens["num_facets_coef"],
        hinge_degree_target=ens["edge_target"],
        num_hinges_coef=ens["num_hinges_coef"],
        hinge_degree_target_coef=ens["edge_target_coef"])
    s = ddg.ManifoldSampler(m, params)
    tilt = list(ens["tilt"])
    if any(tilt) or ens["zleg_scale"] or ens["cimp_scale"]:
        s.set_n6_potential(ens["zleg_scale"] * ens["lam"],
                           ens["cimp_scale"] * ens["lam"], tilt=tilt)
    v = s.manifold
    dop_n6 = CLASS_N6[ens["dopant_class"]]
    host_n6 = ens["host_n6"]

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    snap_base = os.path.splitext(args.out_csv)[0]

    tracking = False
    if args.from_snap:
        cpath = re.sub(r"\.mfd$", "", args.from_snap) + ".cocycle.npz"
        if os.path.exists(cpath):
            e0, w0, _ = coc.load_cocycle(cpath)
            try:
                s.enable_cocycle(e0, w0)
            except RuntimeError:
                # legacy npz saved with raw sampler labels: the .mfd writer
                # renumbers to rank order, so canonicalize and retry
                e0, w0 = coc.canonicalize_labels(e0, w0)
                s.enable_cocycle(e0, w0)
            tracking = True
    elif args.cocycle:
        edges = np.asarray(v.simplices(1))
        from cocycle_check import reference_frac_positions
        mm = re.search(r"T3_([A-Z0-9]+)_m(\d+)_", CRYSTALS[args.structure])
        refname, mcell = mm.group(1).lower(), int(mm.group(2))
        frac = reference_frac_positions(refname, mcell)
        s.enable_cocycle(edges, coc.build_from_positions(edges, frac, mcell))
        tracking = True
    if args.flip_log_mb > 0:
        s.enable_six_flip_log(args.flip_log_mb)
        flip_path = f"{snap_base}_flips.bin"
        open(flip_path, "wb").close()

    def save_state(mfd_path, sweeps):
        v.save(mfd_path)
        stem = re.sub(r"\.mfd$", "", mfd_path)
        if tracking:
            e1, w1 = coc.canonicalize_labels(*s.read_cocycle())
            coc.save_cocycle(stem + ".cocycle.npz", e1, w1, sweeps=sweeps)
        write_sidecar(stem + ".run.json", ens, sweeps, start,
                      ledger_clock(s), args, args.from_snap, diff)

    write_sidecar(f"{snap_base}.run.json", ens, start, start,
                  ledger_clock(s), args, args.from_snap, diff)
    with open(args.out_csv, "w", newline="") as f:
        f.write(f"# ensemble={h} segment_start={start} "
                f"parent={args.from_snap or ''}\n")
        w = csv.writer(f)
        w.writerow(COLS)
        if start == 0:
            r = snapshot_row(v, dop_n6, host_n6)
            w.writerow([0, 0.0] + [r[c] for c in COLS[2:]])
        f.flush()
        done = start
        chunk_i = 0
        while done < args.sweeps_total:
            s.reset_stats()
            s.run(sweeps=min(args.chunk, args.sweeps_total - done))
            done += args.chunk
            chunk_i += 1
            st = s.get_stats()
            acc = st.total_accepted / max(1, st.total_tried)
            r = snapshot_row(v, dop_n6, host_n6)
            w.writerow([done, acc] + [r[c] for c in COLS[2:]])
            f.flush()
            if tracking:
                s.check_cocycle()   # raises on drift -- fail loud, not late
            if args.flip_log_mb > 0:
                ev = s.drain_six_flip_log()
                if s.six_flip_log_overflowed():
                    print(f"WARNING: flip log overflowed in chunk {chunk_i}",
                          flush=True)
                with open(flip_path, "ab") as ff:
                    ev.tofile(ff)
            if args.snap_every and chunk_i % args.snap_every == 0:
                save_state(f"{snap_base}_snap{done}.mfd", done)
            print(f"t={done} n_dop={r['n_dop']} comp={r['n_companion']} "
                  f"clusters={r['n_clusters']} mean_sz={r['mean_sz']:.2f} "
                  f"max={r['max_sz']} jmax={r['j_max']} "
                  f"jfrac={r['j_frac_largest']:.2f} pure56={r['pure56']:.4f}",
                  flush=True)
    if args.save_final:
        save_state(args.save_final, done)


if __name__ == "__main__":
    main()
