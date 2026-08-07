#!/usr/bin/env python3
"""Spatial statistics of defect-complex centroids from a --track series.

Sharpens the single-snapshot centroid estimator of defect_density_hu.py by
pooling the per-frame centroid point sets of a crystal_gas_scan --track
worldline record (hundreds of quasi-independent configurations instead of
one or two snapshots). Position-sector only -- number statistics of WHERE
complexes sit:

  S_cen(k)   direct centroid structure factor |sum_i e^{ik r_i}|^2 / N on
             the torus wavevectors k = 2 pi n / L, |n| <= nmax, averaged
             over frames and k-shells; 1 = Poisson, <1 repulsive/HU,
             >1 clustered. Same estimator as defect_density_hu's
             "COMPLEX centroids S(k0)" (its random-site null is 1 + O(1/V)
             here, absorbed into the quoted error).
  g(r)       pair correlation via min-image distances vs the ideal-gas
             shell normalization (valid r < L/2).
  NN         nearest-neighbour distance distribution vs the Poisson
             expectation P(>r) = exp(-rho 4/3 pi r^3).

Errors: moving-block bootstrap over frames (frames are correlated on the
blink timescale; use --every to thin).

Usage:
    python scripts/defect_dynamics/centroid_spatial.py PRE
        [--burn-frames 250] [--every 5] [--nmax 6]
"""
import argparse
import json
import os
import sys

import numpy as np

CELL = 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pre")
    ap.add_argument("--burn-frames", type=int, default=250)
    ap.add_argument("--every", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--nboot", type=int, default=200)
    args = ap.parse_args()

    raw = [json.loads(l) for l in open(args.pre + ".ts.jsonl")]
    hdr = next((r for r in raw if r.get("kind") == "header"), {})
    frames = [r for r in raw if "sweep" in r][args.burn_frames::args.every]
    mcell = int(hdr.get("mcell"))
    L = float(mcell)
    cents = [np.asarray(r["cents"], float) / CELL for r in frames]
    Ns = np.array([len(c) for c in cents])
    rho = Ns.mean() / L ** 3
    print(f"{os.path.basename(args.pre)}: {len(frames)} frames "
          f"(every {args.every}), <N>={Ns.mean():.1f} complexes, "
          f"box {mcell} cells, rho={rho:.2f}/cell^3")

    # ---- S_cen(k) on torus wavevectors, shell-averaged ----
    ns = np.array([(i, j, k) for i in range(-args.nmax, args.nmax + 1)
                   for j in range(-args.nmax, args.nmax + 1)
                   for k in range(-args.nmax, args.nmax + 1)
                   if (i, j, k) != (0, 0, 0)
                   and i * i + j * j + k * k <= args.nmax ** 2])
    kmag = np.linalg.norm(ns, axis=1)          # in units of 2 pi / L
    per_frame = []                             # (frames, nk)
    for c in cents:
        ph = np.exp(2j * np.pi * (c @ ns.T) / L)
        per_frame.append(np.abs(ph.sum(0)) ** 2 / len(c))
    per_frame = np.array(per_frame)

    def boot_mean(x, nb=8):
        blk = len(x) // nb
        rng = np.random.default_rng(0)
        ms = []
        for _ in range(args.nboot):
            idx = np.concatenate([np.arange(s, s + blk) for s in
                                  rng.integers(0, len(x) - blk, nb)])
            ms.append(x[idx].mean(0))
        return x.mean(0), np.std(ms, axis=0)

    Sk, Sk_err = boot_mean(per_frame)
    print("  S_cen by |n| shell (1 = Poisson):")
    for lo, hi in ((0.9, 1.1), (1.3, 1.8), (1.9, 2.5), (2.6, 3.5),
                   (3.6, 5.0), (5.0, 6.1)):
        m = (kmag >= lo) & (kmag < hi)
        if m.sum():
            v = Sk[m].mean()
            e = np.sqrt((Sk_err[m] ** 2).mean() / m.sum())
            print(f"    |n| in [{lo},{hi}): S = {v:.3f} +- {e:.3f} "
                  f"(n_k={m.sum()})")

    # ---- g(r) ----
    edges = np.linspace(0.05, L / 2, 30)
    counts = np.zeros(len(edges) - 1)
    npairs_tot = 0.0
    nn_all = []
    for c in cents:
        d = c[:, None, :] - c[None, :, :]
        d -= L * np.round(d / L)
        dist = np.sqrt((d ** 2).sum(-1))
        iu = np.triu_indices(len(c), 1)
        dv = dist[iu]
        counts += np.histogram(dv, bins=edges)[0]
        npairs_tot += len(dv)
        np.fill_diagonal(dist, np.inf)
        nn_all.append(dist.min(1))
    vol_shell = 4 / 3 * np.pi * np.diff(edges ** 3)
    # ideal-gas expectation per bin, pooled
    g = counts / (npairs_tot * vol_shell / (L ** 3 / 2 * 2))
    g /= (counts.sum() / npairs_tot) / (vol_shell.sum() * 2 / L ** 3 / 2) \
        if False else 1.0
    # normalize simply: expected pairs per bin = npairs_tot * shellvol / Vbox
    g = counts / (npairs_tot * vol_shell / L ** 3)
    mid = 0.5 * (edges[1:] + edges[:-1])
    print("  g(r) (1 = ideal gas):")
    for i in range(0, len(mid), 4):
        print(f"    r={mid[i]:.2f} cells: g={g[i]:.3f}")

    nn = np.concatenate(nn_all)
    print(f"  nearest neighbour: median {np.median(nn):.3f} cells "
          f"(Poisson median {(np.log(2) * 3 / (4 * np.pi * rho)) ** (1/3):.3f}), "
          f"P(NN < 0.25 x median_P) = {np.mean(nn < 0.25 * (np.log(2) * 3 / (4 * np.pi * rho)) ** (1/3)):.3f}")


if __name__ == "__main__":
    main()
