#!/usr/bin/env python3
"""M1 activation spectroscopy: heating perfect TCP crystals, coupling by coupling.

For each TCP crystal (A15 / C15 / sigma tori from scripts/tcp_reference.py) and
each coupling on a ladder, start a FRESH sampler from the perfect crystal under
a MATCHED objective:

  * num_facets_target = the crystal's own f3 (volume pin exerts no force),
  * hinge_degree_target = the crystal's own mean edge degree (5.1111 A15,
    5.1000 C15/C14, 5.1089 sigma) so the edge pin exerts no composition stress,
  * strictly-local fixed-target penalties (EDQ on edges, VDQ on vertices) at
    per-tet scaled coupling lambda (raw: c_e = lambda*edge/6,
    c_v = lambda/(6/edge - 1), same maps as the seed queue),

run a short window, and record per-move-type acceptance plus the FK census of
the final state. Physics: acceptance ~ exp(-c * dE) per move class, so the
log-slope of acceptance vs coupling measures the ACTIVATION ENERGY of each
elementary excitation of the TCP vacuum; the census distinguishes the
reversible flicker regime (fFK stays ~1: virtual defect pairs) from
proliferation/melting (fFK collapses within the window).

Move-type key: bistellar types indexed by coCenter size: 0: 1->4, 1: 2->3,
2: 3->2, 3: 4->1; hinge = 4-4 hinge moves (dim 3).

Output: printed table + out/tcp_melt_m1.json.
"""
import argparse
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets, vertex_class_census

CRYSTALS = {
    "a15":   "data/tcp_reference/T3_A15_m6_N9936.mfd",
    "c15":   "data/tcp_reference/T3_C15_m4_N8704.mfd",
    "sigma": "data/tcp_reference/T3_SIGMA_m4_N11008.mfd",
    "c36":   "data/tcp_reference/T3_C36_m3_N3672.mfd",
    "r":     "data/tcp_reference/T3_R_m3_N24462.mfd",
    "c14":   "data/tcp_reference/T3_C14_m5_N8500.mfd",
    "z":     "data/tcp_reference/T3_Z_m5_N5000.mfd",
    "mu":    "data/tcp_reference/T3_MU_m3_N5994.mfd",
    "p":     "data/tcp_reference/T3_P_m3_N8640.mfd",
    "delta": "data/tcp_reference/T3_DELTA_m3_N8640.mfd",
    "a15big": "data/tcp_reference/T3_A15_m13_N101062.mfd",
    "c15big": "data/tcp_reference/T3_C15_m9_N99144.mfd",
    "rbig":   "data/tcp_reference/T3_R_m4_N57984.mfd",
    "sigmabig": "data/tcp_reference/T3_SIGMA_m7_N58996.mfd",
}
LADDER = [3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.003]   # per-tet scaled lambda


def census(view):
    eu, edeg, V = edges_from_facets(view.facets())
    fz, _ = vertex_class_census(eu, edeg, V)
    return dict(edv=float(np.var(edeg)),
                p56=float(np.mean((edeg == 5) | (edeg == 6))),
                pure56=1.0 - fz["impure"],
                fFK=fz["Z12"] + fz["Z14"] + fz["Z15"] + fz["Z16"],
                fZ12=fz["Z12"], fZ14=fz["Z14"], fZ15=fz["Z15"],
                fZ16=fz["Z16"])


def run_point(path, lam, sweeps, vdq_scale=1.0, zleg_scale=0.0, cimp_scale=0.0,
              save_path=None, tilt=None):
    m = ddg.Manifold.load(path, 3)
    eu, edeg, V = edges_from_facets(m.facets())
    qbar = float(edeg.mean())
    params = ddg.SamplerParams(
        num_facets_target=m.num_facets,
        num_facets_coef=0.1,
        hinge_degree_target=qbar,
        num_hinges_coef=2.0,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * qbar / 6.0,
        codim3_degree_target_coef=vdq_scale * lam / (6.0 / qbar - 1.0),
        codim3_degree_target=ddg.vertex_degree_target(qbar),
    )
    s = ddg.ManifoldSampler(m, params)
    if zleg_scale or cimp_scale or tilt:
        s.set_n6_potential(zleg_scale * lam, cimp_scale * lam, tilt=tilt)
    s.run(sweeps=sweeps)
    st = s.get_stats()
    v = s.manifold
    if save_path:
        sys.path.insert(0, os.path.join(_ROOT, "tools"))
        from seed_utils import (build_metadata_comments, make_leg, obj_of,
                                 read_history)
        # Provenance: append a `melt` leg to the source crystal's root history so
        # the saved state records exactly which crystal it came from and under
        # what melt objective. obj_of captures the standard couplings; the n6
        # potential (if any) is recorded alongside since it isn't a SamplerParam.
        obj = obj_of(params)
        if zleg_scale or cimp_scale or tilt:
            obj["n6"] = {"zleg_c": zleg_scale * lam, "cimp_c": cimp_scale * lam,
                         "tilt": list(tilt) if tilt else None}
        melt_leg = make_leg("melt", obj, sweeps,
                            tried=st.total_tried, accepted=st.total_accepted)
        prior = read_history(path)
        if prior:                              # normal: source carries a root leg
            legs, root, note = [melt_leg], None, None
        else:                                  # source predates crystal-root stamping
            root = f"crystal:{os.path.splitext(os.path.basename(path))[0]}"
            legs = [make_leg("build", {"src": os.path.basename(path)}, 0,
                             from_=root, tried=0, accepted=0), melt_leg]
            note = ("root synthesized from source path; regenerate the reference "
                    "with tcp_reference.py for the exact Wyckoff root")
        comments = build_metadata_comments(
            topology="T3", dimension=3,
            initial_triangulation=os.path.basename(path),
            num_facets_target=params.num_facets_target,
            num_facets_coef=params.num_facets_coef,
            hinge_degree_target=params.hinge_degree_target,
            num_hinges_coef=params.num_hinges_coef,
            hinge_degree_variance_coef=params.hinge_degree_variance_coef,
            codim3_degree_variance_coef=params.codim3_degree_variance_coef,
            hinge_degree_target_coef=params.hinge_degree_target_coef,
            codim3_degree_target_coef=params.codim3_degree_target_coef,
            codim3_degree_target=params.codim3_degree_target,
            growth_step_size=0, eq_sweeps_per_step=0, equilibration_sweeps=sweeps,
            manifold_view=v, objective=s.current_objective, sampler_stats=st,
            legs=legs, prior_history=(prior or None), root=root,
            history_note=note)
        v.save(save_path, comments=comments)
    row = census(v)
    row.update(
        lam=lam, sweeps=sweeps, qbar=qbar, n=v.num_facets,
        acc_total=st.total_accepted / max(1, st.total_tried),
        tries_14=int(st.bistellar_tries[0]), acc_14=int(st.bistellar_accepts[0]),
        tries_23=int(st.bistellar_tries[1]), acc_23=int(st.bistellar_accepts[1]),
        tries_32=int(st.bistellar_tries[2]), acc_32=int(st.bistellar_accepts[2]),
        tries_41=int(st.bistellar_tries[3]), acc_41=int(st.bistellar_accepts[3]),
        tries_h=int(st.hinge_tries), acc_h=int(st.hinge_accepts),
    )
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sweeps", type=int, default=100)
    ap.add_argument("--vdq-scale", type=float, default=1.0,
                    help="scale factor on the VDQ (vertex) term; 0 = EDQ-only "
                         "(the vertex-uniformity target is frustrated against "
                         "FK order — crystal vertex degrees are 20/24/26/28, "
                         "not the uniform f-vector value)")
    ap.add_argument("--zleg-scale", type=float, default=0.0,
                    help="Z-legality coefficient = this * lambda (vertex "
                         "6-valence potential; frustration-free TCP pressure)")
    ap.add_argument("--cimp-scale", type=float, default=0.0,
                    help="impurity-valence coefficient = this * lambda "
                         "(m^2 anti-clustering term)")
    ap.add_argument("--out", default=os.path.join(_ROOT, "out", "tcp_melt_m1.json"))
    ap.add_argument("--structures", nargs="+", default=["a15", "c15", "sigma"],
                    choices=list(CRYSTALS))
    ap.add_argument("--lams", type=float, nargs="+", default=LADDER,
                    help="per-tet lambda ladder (override for fine scans)")
    ap.add_argument("--save-states", default=None,
                    help="directory to save final .mfd per (structure, lambda) "
                         "for positional analysis with crystal_match.py")
    ap.add_argument("--tilt", type=float, nargs=5, default=None,
                    metavar=("T0", "T1", "T2", "T3", "T4"),
                    help="ABSOLUTE chemical-potential tilts on n6 classes "
                         "(index=n6: [Z12, illegal, Z14, Z15, Z16]; NEGATIVE "
                         "favors; NOT scaled by lambda — doping experiments)")
    args = ap.parse_args()

    if args.save_states:
        os.makedirs(args.save_states, exist_ok=True)
    results = {}
    print(f"{'crystal':>7} {'lam':>7} {'acc_tot':>9} {'a23':>8} {'a32':>8} "
          f"{'a14':>8} {'a44':>8} {'pure56':>7} {'fFK':>6} {'EDV':>6}")
    for name in args.structures:
        path = CRYSTALS[name]
        rows = []
        for lam in args.lams:
            sp = (os.path.join(args.save_states, f"{name}_lam{lam:g}.mfd")
                  if args.save_states else None)
            r = run_point(os.path.join(_ROOT, path), lam, args.sweeps,
                          args.vdq_scale, args.zleg_scale, args.cimp_scale,
                          save_path=sp, tilt=args.tilt)
            rows.append(r)
            def rate(a, t):
                return a / t if t else 0.0
            print(f"{name:>7} {lam:7.3f} {r['acc_total']:9.2e} "
                  f"{rate(r['acc_23'], r['tries_23']):8.2e} "
                  f"{rate(r['acc_32'], r['tries_32']):8.2e} "
                  f"{rate(r['acc_14'], r['tries_14']):8.2e} "
                  f"{rate(r['acc_h'], r['tries_h']):8.2e} "
                  f"{r['pure56']:7.4f} {r['fFK']:6.4f} {r['edv']:6.3f}", flush=True)
        results[name] = rows

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
