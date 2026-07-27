#!/usr/bin/env python3
"""Ball variance of the FULL vertex curvature field -- the real HU test.

defect_statics.py treats complexes as point charges, which measures the
SOURCE process only. But screening physics says the compensating response
lives in the legal matrix around each defect (elastic distortion = the Z-class
population shifting), and hyperuniformity -- if present -- is a property of
the TOTAL field, sources plus response. This is the Coulomb lesson: ions
cluster, yet the total charge field is hyperuniform because polarization
follows them around.

So here: sigma^2(R) of sum over vertices in a Euclidean ball (lift chart) of
dq(v) = q_R(v) - qbar, on melt snapshots, with two controls:

  pristine   the same measurement on the perfect crystal. Its residual is
             lattice-structure noise (site classes vs ball boundary), which
             is class-I HU by construction (~ R^2 surface scaling). The melt
             curve can only be read against this floor.
  shuffled   same positions, dq values randomly permuted among vertices --
             destroys ALL spatial correlation at fixed charge population:
             the Poisson-like ceiling (~ R^3 volume scaling).

Between floor and ceiling: where the melt's fluctuation field actually sits.

FRAME: registry lift (extrinsic; CONVENTIONS.md sec 6). Ball
boundaries cut the ~25%-anisotropic registry lattice: the
sigma^2(R) scaling CLASS is frame-robust, R-values and
amplitudes are gauge.
"""
import argparse
import glob as globmod
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
import defect_state as ds

CELL = 1.0e6


def minimg(d, box):
    return d - np.round(d / box) * box


def ball_var(pts, w, box, Rgrid, centers, rng, shuffles=0):
    """Variance over ball centers of sum of w inside radius R."""
    cen = rng.uniform(0, box, (centers, 3))
    d = minimg(cen[:, None, :] - pts[None, :, :], box)
    dist = np.sqrt((d ** 2).sum(-1))
    out = {}
    for R in Rgrid:
        mask = dist < R
        out[float(R)] = [float(np.var(mask @ w))]
        for _ in range(shuffles):
            out[float(R)].append(float(np.var(mask @ rng.permutation(w))))
    return out


def vertex_field(mfd):
    st = ds.DefectState(mfd)
    q = st.vertex_charges()
    verts = sorted(st.v2t)
    qa = np.array([q[v] for v in verts])
    return st, verts, qa - qa.mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--glob", required=True, action="append")
    ap.add_argument("--cell", required=True, help="pristine crystal control")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--centers", type=int, default=2000)
    ap.add_argument("--shuffles", type=int, default=6)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    print("[frame] registry lift -- HU class robust, R/amplitudes gauge")

    box = float(args.mcell)
    Rgrid = np.array([0.3, 0.4, 0.5, 0.7, 0.9, 1.1, 1.4, 1.7, 2.0])
    rng = np.random.default_rng(args.seed)

    # pristine control: reference positions ARE the lift
    ref = ddg.Manifold.load(args.cell, 3)
    _, verts0, dq0 = vertex_field(ref)
    # reference_frac_positions returns period-m coordinates, i.e. ALREADY in
    # cell units on the m-cell torus -- do not rescale (doing so aliases the
    # box onto itself and produced a pristine "variance" above the coherent
    # Cauchy-Schwarz bound, which is how the bug was caught).
    frac = reference_frac_positions("r", args.mcell)
    pts0 = np.asarray(frac)[verts0] % box
    prist = ball_var(pts0, dq0, box, Rgrid, args.centers, rng,
                     args.shuffles)
    print(f"pristine {args.cell}: {len(verts0)} vertices, "
          f"var(dq) = {dq0.var():.4f}")

    files = sorted({f for g in args.glob for f in globmod.glob(g)})
    melt = {float(R): [] for R in Rgrid}
    shuf = {float(R): [] for R in Rgrid}
    used = []
    for f in files:
        cocf = f.replace(".mfd", ".cocycle.npz")
        if not os.path.exists(cocf):
            continue
        m = ddg.Manifold.load(f, 3)
        st, verts, dq = vertex_field(m)
        e1, w1, _ = coc.load_cocycle(cocf)
        e1 = np.asarray(e1)
        pos = coc.tree_positions(e1, np.asarray(w1),
                                 int(e1.max()) + 1)[0].astype(float) / CELL
        pts = pos[verts] % box
        bv = ball_var(pts, dq, box, Rgrid, args.centers, rng, args.shuffles)
        for R in Rgrid:
            melt[float(R)].append(bv[float(R)][0])
            shuf[float(R)].extend(bv[float(R)][1:])
        used.append(os.path.basename(f))
        print(f"  {os.path.basename(f)}: {len(verts)} verts, "
              f"var(dq) = {dq.var():.4f}", flush=True)
        del st, m

    print(f"\nsigma^2(R) of total curvature field [rad^2] "
          f"({len(used)} snapshots):")
    print(f"   {'R':>4} {'melt':>10} {'+-':>8} {'shuffled':>10} "
          f"{'pristine':>10} {'melt/shuf':>9}")
    rows = []
    for R in Rgrid:
        mv = np.array(melt[float(R)])
        sv = np.mean(shuf[float(R)])
        pv = prist[float(R)][0]
        rows.append({"R": float(R), "melt": float(mv.mean()),
                     "melt_err": float(mv.std() / np.sqrt(len(mv))),
                     "shuffled": float(sv), "pristine": float(pv)})
        print(f"   {R:4.1f} {mv.mean():10.4f} "
              f"{mv.std()/np.sqrt(len(mv)):8.4f} {sv:10.4f} {pv:10.4f} "
              f"{mv.mean()/sv:9.3f}")

    # scaling exponents by log-log fit over R >= 0.5
    def expo(y):
        m_ = Rgrid >= 0.5
        x = np.log(Rgrid[m_])
        return float(np.polyfit(x, np.log(np.array(y)[m_]), 1)[0])
    e_melt = expo([r["melt"] for r in rows])
    e_shuf = expo([r["shuffled"] for r in rows])
    e_pris = expo([max(r["pristine"], 1e-12) for r in rows])
    print(f"\nscaling exponents (3 = Poisson/volume, 2 = class-I HU/surface):")
    print(f"   melt {e_melt:.2f}   shuffled {e_shuf:.2f}   "
          f"pristine {e_pris:.2f}")

    with open(args.out, "w") as fh:
        json.dump({"files": used, "cell": args.cell, "box": box,
                   "rows": rows, "exponents": {"melt": e_melt,
                                               "shuffled": e_shuf,
                                               "pristine": e_pris}},
                  fh, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
