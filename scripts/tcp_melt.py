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

CRYSTALS = [
    ("A15",   "data/tcp_reference/T3_A15_m6_N9936.mfd"),
    ("C15",   "data/tcp_reference/T3_C15_m4_N8704.mfd"),
    ("sigma", "data/tcp_reference/T3_SIGMA_m4_N11008.mfd"),
]
LADDER = [3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.003]   # per-tet scaled lambda


def census(view):
    eu, edeg, V = edges_from_facets(view.facets())
    fz, _ = vertex_class_census(eu, edeg, V)
    return dict(edv=float(np.var(edeg)),
                p56=float(np.mean((edeg == 5) | (edeg == 6))),
                pure56=1.0 - fz["impure"],
                fFK=fz["Z12"] + fz["Z14"] + fz["Z15"] + fz["Z16"])


def run_point(path, lam, sweeps, vdq_scale=1.0):
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
    s.run(sweeps=sweeps)
    st = s.get_stats()
    v = s.manifold
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
    ap.add_argument("--out", default=os.path.join(_ROOT, "out", "tcp_melt_m1.json"))
    args = ap.parse_args()

    results = {}
    print(f"{'crystal':>7} {'lam':>7} {'acc_tot':>9} {'a23':>8} {'a32':>8} "
          f"{'a14':>8} {'a44':>8} {'pure56':>7} {'fFK':>6} {'EDV':>6}")
    for name, path in CRYSTALS:
        rows = []
        for lam in LADDER:
            r = run_point(os.path.join(_ROOT, path), lam, args.sweeps, args.vdq_scale)
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
