#!/usr/bin/env python3
"""V3b: numerical detailed-balance certification of the HB kernel.

V3's occupancy test is structurally vacuous here: the slide landscape has
flat directions (degenerate knot translates), so the chain glides and
never recurs -- per-microstate occupancy measures translational entropy,
not e^{-S}, at ANY lambda (measured: 682 states/1500 steps at lam=0.4,
757 at lam=1.0). The kernel, however, can be certified DIRECTLY, with no
equilibration at all: for a pair (x, y), y in B(x), every factor of the
detailed-balance identity

    pi(x) q_x(y) alpha(x->y)  =  pi(y) q_y(x) alpha(y->x)

is exactly computable (S from the audited running action, Z from the
scans). Algebraically both sides reduce to min(e^{-d}/Zx, 1/Zy), so the
identity holds IF AND ONLY IF four implementation facts do; those are
what this script measures, pair by pair:

  A1  dS antisymmetry: x found in B(y) at dS = -d (1e-9);
  A2  membership symmetry: y in B(x) => x in B(y), OR the driver's
      safe rejection (alpha=0 when x is not identifiable -- both fluxes
      vanish, DB preserved by rejection; counted, not a failure);
  A3  scan determinism: rescanning x after walk-back reproduces Z(x)
      bitwise (also audits exact state restoration);
  A4  path-validity symmetry: the driver rejects proposals whose scan-
      TREE path has a multi-chord intermediate. The trees at x and y
      differ, so forward-valid/reverse-invalid would be a REAL DB
      violation (nonzero forward flux, zero reverse). Counted exactly.

Where both fluxes are positive, the ln-flux residual is also checked to
machine precision. Pairs are drawn alternately pi-weighted (typical
flux) and uniform over the ball (stresses the dS_max boundary, where
membership asymmetry lives). The wander policy (stay at pi-weighted y)
diversifies x without affecting kernel facts. Running-action audits
against the edge-map oracle as in V3.

Each pair costs 2-3 D scans (~12 MB leaked each, R6); one process does
--pairs pairs (default 60 ~ 2 GB peak) -- run several seeds as separate
processes for more coverage.
"""
import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg  # noqa: F401
from discrete_differential_geometry import fpkmc
import worm_slide as ws
from fpkmc_v3_hb import build


def node_chord(g, j):
    return tuple(sorted((int(g["chord"][j][0]), int(g["chord"][j][1]))))


def tree_path_ok(g, par, t):
    """True iff every INTERMEDIATE node on the scan-tree path root->t is
    single-chord (endpoints are checked by the caller)."""
    j = t
    while j != 0:
        a, _ = par[j]
        j = a
        if j != 0 and int(g["n_chords"][j]) != 1:
            return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--pairs", type=int, default=60)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--seed", type=int, default=201)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    m, s, apx, window = build(args.ref, args.estar, args.lam)
    v = s.manifold
    em0 = ws.edeg_dict(v)
    drv = fpkmc.HBDriver(s, apx, dS_max=args.dS_max, max_depth=args.depth,
                         rng=rng, audit_every=0, oracle=None)

    stats = dict(pairs=0, skipped_multichord=0, anti_max=0.0,
                 not_in_ball=0, path_asym=0, flux_resid_max=0.0,
                 both_zero=0, zx_mismatch=0, stays=0)
    S_rel = 0.0

    for k in range(args.pairs):
        gx = drv._scan()
        wx = np.exp(-np.asarray(gx["dS"]))
        Zx = float(wx.sum())
        n = len(wx)
        if k % 2 == 0:                      # pi-weighted (typical flux)
            j = int(rng.choice(n, p=wx / Zx))
            stay = True
        else:                               # uniform (stresses boundary)
            j = int(rng.integers(1, n)) if n > 1 else 0
            stay = False
        if j == 0:
            continue
        if int(gx["n_chords"][j]) != 1:
            stats["skipped_multichord"] += 1
            continue
        d = float(gx["dS"][j])
        par = drv._parents(gx)
        path = drv._path(gx, par, j)
        fwd_ok = tree_path_ok(gx, par, j)

        committed = drv._replay(gx, path)
        assert abs(committed - d) < 1e-6, "replay dS != node dS"
        old = drv.chord
        drv.chord = node_chord(gx, j)
        gy = drv._scan()
        wy = np.exp(-np.asarray(gy["dS"]))
        Zy = float(wy.sum())

        cand = [jj for jj in range(len(wy))
                if node_chord(gy, jj) == old
                and abs(float(gy["dS"][jj]) + d) < 1e-9]
        if len(cand) > 1:
            raise RuntimeError("ambiguous x in B(y)")
        stats["pairs"] += 1
        if not cand:
            # A2: driver rejects (alpha=0) => both fluxes zero, DB safe.
            stats["not_in_ball"] += 1
            drv._walk_back(gx, path)
            drv.chord = old
            continue
        jj = cand[0]
        d2 = float(gy["dS"][jj])
        stats["anti_max"] = max(stats["anti_max"], abs(d + d2))   # A1
        par_y = drv._parents(gy)
        rev_ok = tree_path_ok(gy, par_y, jj)
        if fwd_ok != rev_ok:                                       # A4
            stats["path_asym"] += 1
            print(f"  PATH ASYMMETRY pair {k}: fwd_ok={fwd_ok} "
                  f"rev_ok={rev_ok} d={d:+.4f}", flush=True)
        afwd = min(1.0, float(np.exp(d)) * Zx / Zy) if fwd_ok else 0.0
        arev = min(1.0, float(np.exp(d2)) * Zy / Zx) if rev_ok else 0.0
        flux_f = (np.exp(-d) / Zx) * afwd            # pi(x)=1 baseline
        flux_r = np.exp(-d) * (np.exp(-d2) / Zy) * arev
        if flux_f > 0 and flux_r > 0:
            resid = abs(np.log(flux_f) - np.log(flux_r))
            stats["flux_resid_max"] = max(stats["flux_resid_max"], resid)
        elif flux_f == 0 and flux_r == 0:
            stats["both_zero"] += 1
        # else: covered by path_asym above (the only one-sided-zero source
        # once x is identified)

        if stay:
            stats["stays"] += 1
            S_rel += d
        else:
            drv._walk_back(gx, path)
            drv.chord = old
            if (k // 2) % 5 == 0:                                  # A3
                # determinism = same NODE SET: equal count, sorted dS
                # equal to 1e-9. NOT bitwise: each node's dS accumulates
                # along its DFS tree path, and the traversal order after
                # restoration differs (D AA iteration), so the same node
                # arrives with ulp-level (~1e-14) path-dependent
                # rounding. Bounds the kernel's DB exactness; reported.
                gx2 = drv._scan()
                d1 = np.sort(np.asarray(gx["dS"]))
                d2 = np.sort(np.asarray(gx2["dS"]))
                if len(d1) != len(d2):
                    stats["zx_mismatch"] += 1
                    print(f"  SCAN NODE-COUNT MISMATCH pair {k}: "
                          f"{len(d1)} vs {len(d2)}", flush=True)
                else:
                    ulp = float(np.max(np.abs(d1 - d2))) if len(d1) else 0.0
                    stats["scan_ulp_drift_max"] = max(
                        stats.get("scan_ulp_drift_max", 0.0), ulp)
                    if ulp > 1e-9:
                        stats["zx_mismatch"] += 1
                        print(f"  SCAN dS MISMATCH pair {k}: {ulp:.3e}",
                              flush=True)
        if (k + 1) % 10 == 0:
            true = args.lam * ws.dS_between(em0, ws.edeg_dict(v),
                                            estar=args.estar)
            assert abs(true - S_rel) < 1e-6, "action audit failed"
            print(f"  pair {k + 1}/{args.pairs}: "
                  f"anti {stats['anti_max']:.2e}, "
                  f"resid {stats['flux_resid_max']:.2e}, "
                  f"not-in-ball {stats['not_in_ball']}, "
                  f"path-asym {stats['path_asym']}", flush=True)

    verdict = ("PASS" if stats["path_asym"] == 0
               and stats["zx_mismatch"] == 0
               and stats["anti_max"] < 1e-9
               and stats["flux_resid_max"] < 1e-12 else "FAIL")
    print(f"\nV3b seed {args.seed}: {verdict}")
    for kk, vv in stats.items():
        print(f"  {kk}: {vv}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"seed": args.seed, "verdict": verdict, **stats}, fh)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
