#!/usr/bin/env python3
"""Lensing strength of a single 2->3 defect, across inequivalent sites.

The defect's boundary distance map is UNIVERSAL (``ball_boundary.py``), so
every 2->3 presents the same optics at dB. All site dependence enters
through d_out -- how the six boundary triangles sit in the surrounding
crystal:

    d(x,y) = min( d_out(x,y),
                  min_{p,q in dB} [ d_out(x,p) + d_B(p,q) + d_out(q,y) ] )

with d_out identical before and after the move. Since min(d_out, d_B) is
just the full-crystal distance, the site-characteristic object is

    A(p,q) = d_R(p,q) - d_D(p,q)   for p, q on dB,

the boundary-to-boundary transit gain AFTER the option of routing around B
has been taken into account. Any distant pair benefits from the defect only
by transiting some (p,q) with A > 0, so A bounds and controls the global
field. A is CANONICAL: its source set is dB itself, so it needs no random
sampling and two translation copies of a site give bitwise identical
answers -- which makes the within-class spread a real zero and any
between-class difference meaningful.

The global field (Delta = d_D - d_R from random sources to all vertices) is
also reported, but it is dominated by source-placement noise: the capture
cone is narrow, so with a few dozen sources the statistic is set by whether
any source happens to fall inside it, not by the site. Read `A`, not `xs`.

Usage:
  python scripts/defect_lens_field.py data/tcp_reference/T3_R_m3_N24462.mfd
  python scripts/defect_lens_field.py <mfd> --order 3 --reps 3 --sources 24
"""
import os, sys, argparse, collections, time
from fractions import Fraction

import numpy as np
from scipy.sparse.csgraph import shortest_path

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_ROOT, "python"))

from discrete_differential_geometry import Manifold
from discrete_differential_geometry.ball_boundary import boundary_nodes
from steiner_geodesic import (steiner_distance, build_steiner_graph,
                              _enumerate, face_interior_bary)
from move_site_census import load, sites, ledger, describe

TOL = 1e-9


def steiner_ids(T, V, n, nodes):
    """Node ids in ``build_steiner_graph(T, V, n)`` for barycentric nodes
    given as {vertex: Fraction} weight dicts. Mirrors that function's
    layout: vertices, then (n-1) points per edge keyed from the global-min
    endpoint, then interior barycentric points per face."""
    E, Fd = _enumerate(np.asarray(T, dtype=np.int64))
    fb = face_interior_bary(n)
    base_edge, base_face = V, V + len(E) * (n - 1)
    out = []
    for w in nodes:
        vs = sorted(w)
        if len(vs) == 1:
            out.append(vs[0])
        elif len(vs) == 2:
            a, b = vs                       # a < b; fraction measured from a
            out.append(base_edge + E[(a, b)] * (n - 1) + int(w[b] * n) - 1)
        else:
            v0, v1, v2 = vs
            ijk = tuple(int(w[v] * n) for v in (v0, v1, v2))
            out.append(base_face + Fd[(v0, v1, v2)] * len(fb) + fb.index(ijk))
    return out


def relabel(order, face, apex):
    """dB nodes of ball_boundary, relabelled onto the host's vertices."""
    lab = {0: face[0], 1: face[1], 2: face[2], 3: apex[0], 4: apex[1]}
    return [{lab[v]: w for v, w in nd.items()} for nd in boundary_nodes(order)]


def transit_gain(F0, FD, V, order, face, apex, GR):
    """A(p,q) = d_R(p,q) - d_D(p,q) on dB, in the FULL crystal both times."""
    nodes = relabel(order, face, apex)
    idR = steiner_ids(F0, V, order, nodes)
    dR = shortest_path(GR, method="D", directed=False, indices=idR)[:, idR]
    GD, _ = build_steiner_graph(FD, V, order)
    idD = steiner_ids(FD, V, order, nodes)
    dD = shortest_path(GD, method="D", directed=False, indices=idD)[:, idD]
    A = dR - dD
    iu = np.triu_indices(len(nodes), 1)
    a = A[iu]
    return dict(A_max=float(a.max()), A_mean=float(a.mean()),
                A_frac=float((a > TOL).mean()), A_min=float(a.min()),
                n_nodes=len(nodes))


def global_field(FD, V, order, src, base):
    dD = steiner_distance(FD, V, order, src)
    ok = np.isfinite(base) & np.isfinite(dD) & (base > 0)
    d = (dD - base)[ok]
    return dict(xs=float((d < -TOL).mean()), min_dd=float(d.min()),
                lens=float(-np.minimum(d, 0).sum() / ok.sum()))


def main(crystal, order, nsrc, reps, seed):
    m, facets, deg, nbr, f2a = load(crystal)
    S = sites(deg, f2a)
    by = collections.defaultdict(list)
    for f, a in S:
        by[ledger(f, a, deg)].append((f, a))
    chosen = []
    for key in sorted(by, key=lambda k: -len(by[k])):
        mem = by[key]
        step = max(1, len(mem) // reps)
        chosen += [(key, f, a) for f, a in mem[::step][:reps]]

    F0 = [tuple(sorted(map(int, f))) for f in m.facets()]
    F0a = np.asarray(F0, np.int64)
    V = int(F0a.max()) + 1
    rng = np.random.default_rng(seed)
    src = sorted(rng.choice(V, min(nsrc, V), replace=False).tolist())

    print(f"{os.path.basename(crystal)}: {V} vertices, {len(F0)} tets")
    print(f"  Steiner order {order}; {len(by)} ledgers x {reps} reps; "
          f"global field from {len(src)} random sources (seed {seed})")
    t0 = time.time()
    GR, _ = build_steiner_graph(F0a, V, order)
    base = steiner_distance(F0a, V, order, src)
    fin = np.isfinite(base) & (base > 0)
    print(f"  reference: mean d = {base[fin].mean():.3f}, "
          f"max = {base[fin].max():.2f}   ({time.time()-t0:.1f}s)\n")

    rows = collections.defaultdict(list)
    for key, face, apex in chosen:
        mm = Manifold(3, [list(t) for t in F0])
        mm.do_bistellar_move(list(face), list(apex))
        FD = np.asarray(mm.facets(), np.int64)
        r = transit_gain(F0a, FD, V, order, face, apex, GR)
        r.update(global_field(FD, V, order, src, base))
        rows[key].append(r)
        d = describe(*key)
        print(f"  new6={d['new_deg6']} new7={d['new_deg7']} {str(key[0]):<10}"
              f"{str(key[1]):<21} face {str(face):<20}"
              f"A_max={r['A_max']:.6f} A_mean={r['A_mean']:.6f} "
              f"A_frac={100*r['A_frac']:5.1f}%   [global xs={100*r['xs']:.2f}%]")

    print(f"\n  CANONICAL boundary-transit gain A (no sampling; identical for "
          f"translation copies)")
    print(f"  {'face edges':<12} {'spokes':<21} {'new6':>4} {'new7':>4} "
          f"{'sites':>7} {'A_max':>18} {'A_mean':>18} {'A_frac %':>16}")
    agg = []
    for key in sorted(rows, key=lambda k: -len(by[k])):
        rs = rows[key]
        d = describe(*key)
        f = lambda nm, s=1: (np.mean([r[nm] for r in rs]) * s,
                             np.std([r[nm] for r in rs]) * s)
        mx, mxd = f("A_max"); mn, mnd = f("A_mean"); fr, frd = f("A_frac", 100)
        agg.append((d["new_deg6"], d["new_deg7"], mx, mn, fr))
        print(f"  {str(key[0]):<12} {str(key[1]):<21} {d['new_deg6']:>4} "
              f"{d['new_deg7']:>4} {len(by[key]):>7} {mx:>11.6f}+/-{mxd:.6f} "
              f"{mn:>11.6f}+/-{mnd:.6f} {fr:>9.2f}+/-{frd:.2f}")

    a = np.array(agg)
    print()
    for j, nm in ((2, "A_max"), (3, "A_mean"), (4, "A_frac")):
        print(f"  corr(new_deg6, {nm}) = {np.corrcoef(a[:,0], a[:,j])[0,1]:+.3f}"
              f"    corr(new_deg7, {nm}) = "
              f"{np.corrcoef(a[:,1], a[:,j])[0,1]:+.3f}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("crystal")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--sources", type=int, default=24)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    p = args.crystal if os.path.isabs(args.crystal) \
        else os.path.join(_ROOT, args.crystal)
    main(p, args.order, args.sources, args.reps, args.seed)
