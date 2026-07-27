#!/usr/bin/env python3
"""INTRINSIC transport constant: Cartan-development MSD of worldlines.

Positions are the exact developed centroids of the knot's site tets
along its ACTUAL worldline (holo_w*.json jump sequences): the
development is the rolling-without-slipping (Cartan) chart, so the
trajectory's flat image is the antidevelopment and its MSD defines the
canonical intrinsic diffusion constant -- unit-edge PL metric, no
embedding gauge, no torus wrapping (the chart is the universal cover).
Exact rational development end to end; floats only for statistics.

Distances in UNIT-EDGE lengths (the seed embedding has edge sqrt(2);
coordinates are rescaled). Approximate cell conversion via the
registry mean edge (0.476 cells/edge) is printed for comparison with
the registry-frame D_slide -- labeled approximate, since "cells" is
itself a registry notion. Time axis: slide hops (jump-chain steps).
Transport channel: slide only, frozen background (v1).
"""
import glob
import json
import os
import sys
from collections import deque

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import development as dev
from discrete_differential_geometry.fpkmc import face_apex_maps

REF = "data/tcp_reference/T3_R_m2_N7248.mfd"
EDGE = np.sqrt(2.0)          # seed-embedding edge length


def main():
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    _, face_of = face_apex_maps(m)
    ctx = dev.TransportContext([tuple(int(x) for x in f) for f in F])

    def site_tet(chord):
        c = tuple(sorted(chord))
        f = face_of[c]
        return frozenset(min(tuple(sorted((c[0],) + tuple(f))),
                             tuple(sorted((c[1],) + tuple(f)))))

    def bfs_path(a, b):
        if a == b:
            return [a]
        prev = {a: None}
        q = deque([a])
        while q:
            t = q.popleft()
            for nb in ctx.dual[t]:
                if nb in prev:
                    continue
                prev[nb] = t
                if nb == b:
                    out = []
                    while nb is not None:
                        out.append(nb)
                        nb = prev[nb]
                    return list(reversed(out))
                q.append(nb)
        raise RuntimeError("no dual path")

    # ---- develop each worldline exactly, collect float centroids ----
    trajs = []
    for fn in sorted(glob.glob("data/fpkmc/holo_w*.json")):
        chords = [tuple(c) for c in json.load(open(fn))["chords"]]
        tets = [site_tet(c) for c in chords]
        cur_tet = tuple(sorted(tets[0]))
        cur_pos = ctx.canon_placement(tets[0])
        pts = [np.array([float(sum(cur_pos[v][k] for v in cur_tet)) / 4
                         for k in range(3)])]
        for k in range(len(tets) - 1):
            path = bfs_path(tets[k], tets[k + 1])
            placements = dev.develop_path(cur_tet, cur_pos,
                                          [tuple(sorted(t))
                                           for t in path])
            cur_pos = placements[-1]
            cur_tet = tuple(sorted(path[-1]))
            pts.append(np.array(
                [float(sum(cur_pos[v][j] for v in cur_tet)) / 4
                 for j in range(3)]))
        trajs.append(np.asarray(pts) / EDGE)     # unit-edge lengths
        print(f"  {os.path.basename(fn)}: {len(pts)} positions, "
              f"|x_end - x_0| = "
              f"{np.linalg.norm(trajs[-1][-1] - trajs[-1][0]):.2f} edges",
              flush=True)

    # ---- MSD over hop lags ----
    Ls = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    msd = {}
    for L in Ls:
        d2 = []
        for x in trajs:
            n = len(x)
            if L >= n:
                continue
            disp = x[L:] - x[:-L]
            d2.extend((disp ** 2).sum(1))
        if len(d2) >= 20:
            msd[L] = (np.mean(d2), np.std(d2) / np.sqrt(len(d2)),
                      len(d2))

    print(f"\n{'l (hops)':>9} {'MSD (edge^2)':>13} {'SE':>8} {'n':>6}")
    for L, (v, se, n) in msd.items():
        print(f"{L:9d} {v:13.3f} {se:8.3f} {n:6d}")
    Lf = np.array([L for L in msd if 4 <= L <= 128], float)
    Mf = np.array([msd[L][0] for L in msd if 4 <= L <= 128])
    alpha, lnA = np.polyfit(np.log(Lf), np.log(Mf), 1)
    Dl = Mf[-1] / (6 * Lf[-1])
    print(f"\nintrinsic MSD exponent alpha = {alpha:.3f} "
          f"(fit l in [4,128])")
    print(f"D_intrinsic = MSD/6l at l={int(Lf[-1])}: {Dl:.4f} "
          f"edge^2/hop")
    print(f"(approx registry conversion: x (0.476 cells/edge)^2 = "
          f"{Dl * 0.476**2:.4f} cells^2/hop -- approximate, 'cells' is "
          f"a registry notion)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5))
    Ls_p = list(msd)
    ax.errorbar(Ls_p, [msd[L][0] for L in Ls_p],
                yerr=[msd[L][1] for L in Ls_p], fmt="o-",
                label=f"intrinsic MSD (alpha = {alpha:.2f})")
    ax.plot(Lf, np.exp(lnA) * Lf ** alpha, "r--", lw=1, label="fit")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("lag (slide hops)")
    ax.set_ylabel("MSD (unit-edge lengths^2)")
    ax.set_title("Intrinsic (Cartan-development) MSD -- lone knot, "
                 "R m2, lam=0.40, slide channel,\nfrozen bg; exact "
                 "development of actual worldlines, universal-cover "
                 "chart", fontsize=9)
    ax.legend()
    fig.tight_layout()
    out = "data/fpkmc/fp_intrinsic_msd.png"
    fig.savefig(out, dpi=140)
    print(f"Saved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
