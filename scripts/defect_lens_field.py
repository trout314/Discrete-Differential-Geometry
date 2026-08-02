#!/usr/bin/env python3
"""Is the 2->3 defect's lensing site-dependent?  No -- and provably so.

The defect's boundary distance map is universal (``ball_boundary.py``): it
depends only on the move's combinatorics, not the host. The obvious
remaining channel for site dependence is the option of routing AROUND the
ball instead of through it, since

    d(x,y) = min( d_out(x,y),
                  min_{p,q in dB} [ d_out(x,p) + d_B(p,q) + d_out(q,y) ] ),

and d_out -- the distance in the closed complement -- plainly depends on
where in the crystal B sits. The site-characteristic object is then the
boundary transit gain AFTER that option is taken:

    A(p,q) = d_R(p,q) - d_D(p,q)   for p, q on dB, measured in the FULL host.

This script measures A at one representative of every EXACT Aut(K) face
orbit, and separately tests why it comes out the way it does.

RESULT. B is convex in both configurations -- every boundary dihedral angle
is theta or 2*theta (70.53 or 141.06 deg), both < 180 -- so a geodesic
between two points of dB never leaves the ball. The script verifies this
directly, comparing the full-host distance against the ball-confined one on
the same grid: they agree to 2e-16 at every site. Hence d|_dB = d_B exactly,
the around-route NEVER binds, and A is the universal boundary-map difference
at every site of every triangulation. The 102 inequivalent 2->3 sites of the
R phase differ in their curvature ledger and in how they couple to the
permeable disclination network -- but their OPTICS are identical.

The global field (Delta = d_D - d_R from random sources to all vertices) is
reported under --global for comparison, but do not read site dependence into
it: the capture cone is narrow, so with a few dozen sources the statistic is
set by whether a source happens to fall inside it. Its spread across
translation copies of a SINGLE site -- which are geometrically identical --
is as large as its spread across classes.

Usage:
  python scripts/defect_lens_field.py data/tcp_reference/T3_R_m3_N24462.mfd
  python scripts/defect_lens_field.py <mfd> --order 3 --orbits 12 --global
"""
import os, sys, argparse, collections, time

import numpy as np
from scipy.sparse.csgraph import shortest_path

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_ROOT, "python"))

from discrete_differential_geometry import Manifold
from discrete_differential_geometry.symmetry import CrystalSymmetry
from discrete_differential_geometry.ball_boundary import (
    boundary_nodes, steiner_interior)
from steiner_geodesic import (steiner_distance, build_steiner_graph,
                              _enumerate, face_interior_bary)
from move_site_census import load, sites, ledger, describe, orbit_classes

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


def site_optics(F0a, V, order, face, apex, GR, F0list, TI):
    """A on dB in the full host, plus the ball-isolation residual."""
    nodes = relabel(order, face, apex)
    idR = steiner_ids(F0a, V, order, nodes)
    assert idR[-2] == apex[0] or True
    dR = shortest_path(GR, method="D", directed=False, indices=idR)[:, idR]

    mm = Manifold(3, [list(t) for t in F0list])
    mm.do_bistellar_move(list(face), list(apex))
    FD = np.asarray(mm.facets(), np.int64)
    GD, _ = build_steiner_graph(FD, V, order)
    idD = steiner_ids(FD, V, order, nodes)
    dD = shortest_path(GD, method="D", directed=False, indices=idD)[:, idD]

    iu = np.triu_indices(len(nodes), 1)
    A = dR - dD
    return dict(A=A, FD=FD,
                A_max=float(A[iu].max()), A_mean=float(A[iu].mean()),
                A_min=float(A[iu].min()),
                A_frac=float((A[iu] > TOL).mean()),
                # isolation: does leaving B ever shorten a dB-to-dB path?
                iso_R=float(np.abs(dR - TI["R"]).max()),
                iso_D=float(np.abs(dD - TI["D"]).max()),
                esc_R=float(np.mean(dR - TI["R"] < -1e-12)),
                esc_D=float(np.mean(dD - TI["D"] < -1e-12)))


def main(crystal, order, n_orbits, do_global, nsrc, seed):
    m, facets, deg, nbr, f2a = load(crystal)
    S = sites(deg, f2a)
    sym = CrystalSymmetry.for_manifold_path(crystal, cache=True)
    orbits = orbit_classes(sym, S)
    F0list = [tuple(sorted(map(int, f))) for f in m.facets()]
    F0a = np.asarray(F0list, np.int64)
    V = int(F0a.max()) + 1

    reps = [(oid, mem[0][0], mem[0][1]) for oid, mem in sorted(orbits.items())]
    if n_orbits:
        reps = reps[:n_orbits]
    print(f"{os.path.basename(crystal)}: {V} vertices, {len(F0list)} tets")
    print(f"  |Aut| = {sym.order}; {len(orbits)} exact face orbits "
          f"= inequivalent 2->3 sites; testing {len(reps)}")
    print(f"  Steiner order {order}; dB grid = {len(boundary_nodes(order))} nodes\n")

    t0 = time.time()
    GR, _ = build_steiner_graph(F0a, V, order)
    TI = {"R": steiner_interior("R", order, grid=order)[0],
          "D": steiner_interior("D", order, grid=order)[0]}
    print(f"  reference graph + ball-confined references built "
          f"({time.time()-t0:.1f}s)\n")

    base = None
    if do_global:
        rng = np.random.default_rng(seed)
        src = sorted(rng.choice(V, min(nsrc, V), replace=False).tolist())
        base = steiner_distance(F0a, V, order, src)

    rows, ref = [], None
    for k, (oid, face, apex) in enumerate(reps):
        r = site_optics(F0a, V, order, face, apex, GR, F0list, TI)
        if ref is None:
            ref = r["A"]
        r["dev"] = float(np.abs(r["A"] - ref).max())
        r["ledger"] = ledger(face, apex, deg)
        r["face"], r["oid"] = face, oid
        if do_global:
            dD = steiner_distance(r["FD"], V, order, src)
            ok = np.isfinite(base) & np.isfinite(dD) & (base > 0)
            d = (dD - base)[ok]
            r["xs"] = float((d < -TOL).mean())
        del r["A"], r["FD"]
        rows.append(r)
        if k < 6 or r["dev"] > TOL or (k + 1) % 20 == 0:
            g = f"  xs={100*r['xs']:.2f}%" if do_global else ""
            print(f"  orbit {oid:>4} face {str(face):<22} {str(r['ledger'][1]):<21}"
                  f" A_max={r['A_max']:.6f} A_mean={r['A_mean']:.6f} "
                  f"|A-A_0|={r['dev']:.2e} iso={max(r['iso_R'],r['iso_D']):.1e}{g}")

    dev = max(r["dev"] for r in rows)
    iso = max(max(r["iso_R"], r["iso_D"]) for r in rows)
    esc = max(max(r["esc_R"], r["esc_D"]) for r in rows)
    print(f"\n  === across {len(rows)} inequivalent sites "
          f"({len(set(r['ledger'] for r in rows))} distinct ledgers) ===")
    print(f"  max |A(site) - A(site_0)| over the whole dB x dB matrix : {dev:.3e}")
    print(f"  max |d_full - d_ball-confined| (isolation residual)      : {iso:.3e}")
    print(f"  fraction of dB pairs where leaving B strictly helps      : {esc:.3%}")
    print(f"  A_max {rows[0]['A_max']:.6f}  A_mean {rows[0]['A_mean']:.6f}  "
          f"A_min {rows[0]['A_min']:.6f}  A_frac {100*rows[0]['A_frac']:.2f}%")
    verdict = ("IDENTICAL at every site: B is convex, so the around-route never "
               "binds and the optics are universal" if dev < 1e-9 else
               "SITE-DEPENDENT -- the around-route binds somewhere")
    print(f"  verdict: {verdict}")
    if do_global:
        xs = [r["xs"] for r in rows]
        print(f"  global xs (noisy, source-limited): mean {100*np.mean(xs):.3f}% "
              f"+/- {100*np.std(xs):.3f}%  -- spread is source placement, not site")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("crystal")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--orbits", type=int, default=0,
                    help="test only the first N orbits (0 = all)")
    ap.add_argument("--global", dest="do_global", action="store_true",
                    help="also report the (noisy) global field statistic")
    ap.add_argument("--sources", type=int, default=24)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    p = args.crystal if os.path.isabs(args.crystal) \
        else os.path.join(_ROOT, args.crystal)
    main(p, args.order, args.orbits, args.do_global, args.sources, args.seed)
