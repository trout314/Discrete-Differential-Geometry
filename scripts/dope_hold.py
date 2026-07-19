#!/usr/bin/env python3
"""Dope-and-hold: time series of dopant population in a TCP crystal.

Holds a perfect crystal at fixed EDQ coupling lambda and fixed chemical
potential mu on a dopant class (tcp_melt.py conventions), recording the
Z census, dopant count, and dopant-cluster statistics every chunk. Answers:

  * does n(t) PLATEAU (equilibrium isotherm; the doping measurements are
    thermodynamic) or keep growing (kinetics/nucleation-limited)?
  * birth size vs aggregation: does mean multiplet size grow at fixed n
    (coarsening of bound complexes) or stay at the insertion quantum?

Usage:
    python scripts/dope_hold.py --structure a15big --dopant-class Z16 \
        --mu 2.5 --sweeps-total 12000 --chunk 250 \
        --out-csv data/dope_hold/a15big_z16_mu2.5.csv
"""
import argparse
import csv
import os
import re
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--structure", required=True, choices=list(CRYSTALS))
    ap.add_argument("--dopant-class", required=True, choices=list(CLASS_N6))
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--zleg-scale", type=float, default=0.3)
    ap.add_argument("--cimp-scale", type=float, default=0.3)
    ap.add_argument("--edge-target", type=float, default=None,
                    help="override hinge_degree_target (default: the crystal's "
                         "own mean edge degree). Retargeting to flat 5.1043 "
                         "makes the dopant gas supply the curvature offset — "
                         "the flat-by-doping / lattice-OCP experiment.")
    ap.add_argument("--sweeps-total", type=int, default=12000)
    ap.add_argument("--chunk", type=int, default=250)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--save-final", default=None)
    ap.add_argument("--init", default=None,
                    help="start from this .mfd instead of the perfect crystal "
                         "(extend a previous hold from its saved final state)")
    ap.add_argument("--tilt", type=float, nargs=5, default=None,
                    help="ABSOLUTE 5-component tilt override (index=n6); "
                         "replaces the single-class --mu construction (e.g. "
                         "host-class crystal fields)")
    ap.add_argument("--host-n6", type=int, nargs="+", default=None,
                    help="native n6 classes of the host (C15: 0 4); enables "
                         "joint defect-complex tracking columns")
    ap.add_argument("--snap-every", type=int, default=0,
                    help="save a snapshot .mfd every K chunks (0 = off)")
    ap.add_argument("--seed", type=int, default=None,
                    help="D-side RNG seed (reproducible chain; record it)")
    ap.add_argument("--cocycle", action="store_true",
                    help="track T3 winding cocycles, initialized from the "
                         "structure's reference coordinates (crystal starts "
                         "only); audited every chunk, saved with snapshots")
    ap.add_argument("--cocycle-file", default=None,
                    help="resume cocycle tracking from a .cocycle.npz saved "
                         "by a previous run (for --init starts)")
    ap.add_argument("--flip-log-mb", type=float, default=0.0,
                    help="enable the six-edge flip stream with this buffer; "
                         "records append to <out-csv-base>_flips.bin "
                         "(dtype: disclination.SIX_FLIP_DTYPE)")
    args = ap.parse_args()

    if args.seed is not None:
        ddg.set_random_seed(args.seed)

    init_path = args.init or os.path.join(_ROOT, CRYSTALS[args.structure])
    m = ddg.Manifold.load(init_path, 3)
    eu, edeg, V = edges_from_facets(m.facets())
    qbar = float(edeg.mean())
    target = args.edge_target if args.edge_target is not None else qbar
    params = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=target, num_hinges_coef=2.0,
        hinge_degree_target_coef=args.lam * target / 6.0)
    s = ddg.ManifoldSampler(m, params)
    if args.tilt is not None:
        tilt = list(args.tilt)
    else:
        tilt = [0.0] * 5
        tilt[CLASS_N6[args.dopant_class]] = -args.mu
    if any(tilt) or args.zleg_scale or args.cimp_scale:
        s.set_n6_potential(args.zleg_scale * args.lam,
                           args.cimp_scale * args.lam, tilt=tilt)
    v = s.manifold
    dop_n6 = CLASS_N6[args.dopant_class]
    host_n6 = args.host_n6

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    snap_base = os.path.splitext(args.out_csv)[0]

    tracking = False
    if args.cocycle or args.cocycle_file:
        edges = np.asarray(v.simplices(1))
        if args.cocycle_file:
            e0, w0, _ = coc.load_cocycle(args.cocycle_file)
            s.enable_cocycle(e0, w0)
        else:
            if args.init:
                raise SystemExit("--cocycle needs a crystal start; use "
                                 "--cocycle-file to resume from --init")
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
    with open(args.out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        r = snapshot_row(v, dop_n6, host_n6)
        w.writerow([0, 0.0] + [r[c] for c in COLS[2:]])
        f.flush()
        done = 0
        chunk_i = 0
        while done < args.sweeps_total:
            s.reset_stats()
            s.run(sweeps=args.chunk)
            done += args.chunk
            chunk_i += 1
            st = s.get_stats()
            acc = st.total_accepted / max(1, st.total_tried)
            r = snapshot_row(v, dop_n6, host_n6)
            w.writerow([done, acc] + [r[c] for c in COLS[2:]])
            f.flush()
            if tracking:
                s.check_cocycle()   # raises on drift — fail loud, not late
            if args.flip_log_mb > 0:
                ev = s.drain_six_flip_log()
                if s.six_flip_log_overflowed():
                    print(f"WARNING: flip log overflowed in chunk {chunk_i}",
                          flush=True)
                with open(flip_path, "ab") as ff:
                    ev.tofile(ff)
            if args.snap_every and chunk_i % args.snap_every == 0:
                v.save(f"{snap_base}_snap{done}.mfd")
                if tracking:
                    e1, w1 = s.read_cocycle()
                    coc.save_cocycle(f"{snap_base}_snap{done}.cocycle.npz",
                                     e1, w1, sweeps=done)
            print(f"t={done} n_dop={r['n_dop']} comp={r['n_companion']} "
                  f"clusters={r['n_clusters']} mean_sz={r['mean_sz']:.2f} "
                  f"max={r['max_sz']} jmax={r['j_max']} "
                  f"jfrac={r['j_frac_largest']:.2f} pure56={r['pure56']:.4f}",
                  flush=True)
    if args.save_final:
        v.save(args.save_final)
        if tracking:
            e1, w1 = s.read_cocycle()
            coc.save_cocycle(os.path.splitext(args.save_final)[0]
                             + ".cocycle.npz", e1, w1,
                             sweeps=args.sweeps_total)


if __name__ == "__main__":
    main()
