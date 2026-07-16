#!/usr/bin/env python3
"""Compare two certified seed families at the replica-ensemble level.

Built for the VDQ/EDQ validation: a variance-mode family (e.g. ..._VDVs_2e-3)
and its exact-mapped fixed-target twin (..._VDQ_2e-3) sample measures that are
provably identical conditioned on the f-vector, differing only by a soft
mean-degree pin ~1e-3 of the k=2 edge pin. This script tests that claim
empirically on the stored replicas: equal-time observables (VDV, HDV, edge
degree, f-vector) compared replica-set vs replica-set (Welch z + two-sample
KS), plus pooled vertex/edge degree distributions (total-variation distance
against the resampling null). Works for ANY two families, so it also serves
Stage-1 spot checks and future objective changes.

    python scripts/compare_twins.py \
        S3_N1e3_1e-1_ED5p1043_2_VDVs_2e-3 S3_N1e3_1e-1_ED5p1043_2_VDQ_2e-3
"""

import argparse
import glob
import math
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
from discrete_differential_geometry import Manifold  # noqa: E402
from seed_utils import load_seed_metadata  # noqa: E402


def replica_observables(paths, dim, with_hists):
    """Per-replica equal-time observables (+ pooled degree histograms)."""
    rows = {"vdv": [], "hdv": [], "edge_deg": [], "f0": [], "f1": [], "f3": []}
    vhist_tot, ehist_tot = np.zeros(0), np.zeros(0)
    for p in paths:
        md = load_seed_metadata(p)
        fv = [int(x) for x in md["actual_f_vector"].strip("[]").split(",")]
        rows["f0"].append(fv[0]); rows["f1"].append(fv[1]); rows["f3"].append(fv[3])
        rows["edge_deg"].append(6.0 * fv[3] / fv[1])
        rows["vdv"].append(float(md["actual_vertex_degree_variance"]))
        rows["hdv"].append(float(md["actual_hinge_degree_variance"]))
        if with_hists:
            m = Manifold.load(p, dim)
            vh = np.asarray(m.degree_histogram(0), float)
            eh = np.asarray(m.degree_histogram(1), float)
            if len(vh) > len(vhist_tot):
                vhist_tot = np.pad(vhist_tot, (0, len(vh) - len(vhist_tot)))
            if len(eh) > len(ehist_tot):
                ehist_tot = np.pad(ehist_tot, (0, len(eh) - len(ehist_tot)))
            vhist_tot[:len(vh)] += vh
            ehist_tot[:len(eh)] += eh
            del m
    return {k: np.array(v) for k, v in rows.items()}, vhist_tot, ehist_tot


def ks_two_sample(a, b):
    """Two-sample KS statistic and asymptotic p-value (no scipy)."""
    a, b = np.sort(a), np.sort(b)
    allv = np.concatenate([a, b])
    cdf_a = np.searchsorted(a, allv, side="right") / len(a)
    cdf_b = np.searchsorted(b, allv, side="right") / len(b)
    d = np.abs(cdf_a - cdf_b).max()
    en = math.sqrt(len(a) * len(b) / (len(a) + len(b)))
    t = (en + 0.12 + 0.11 / en) * d
    if t < 0.3:      # series is alternating-degenerate at small t; p -> 1 there
        return d, 1.0
    p = 2 * sum((-1) ** (j - 1) * math.exp(-2 * (j * t) ** 2) for j in range(1, 101))
    return d, min(max(p, 0.0), 1.0)


def tv_distance(h1, h2):
    """Total-variation distance between two count histograms."""
    n = max(len(h1), len(h2))
    p = np.pad(h1, (0, n - len(h1))); q = np.pad(h2, (0, n - len(h2)))
    return 0.5 * np.abs(p / p.sum() - q / q.sum()).sum()


def tv_null(h1, h2, rng, n_boot=200):
    """Null scale for the TV distance: split the POOLED counts into two
    multinomial halves of the observed sizes, n_boot times -> mean +/- sd."""
    n = max(len(h1), len(h2))
    p = np.pad(h1, (0, n - len(h1))); q = np.pad(h2, (0, n - len(h2)))
    pool = (p + q) / (p + q).sum()
    n1, n2 = int(p.sum()), int(q.sum())
    tvs = [tv_distance(rng.multinomial(n1, pool).astype(float),
                       rng.multinomial(n2, pool).astype(float))
           for _ in range(n_boot)]
    return float(np.mean(tvs)), float(np.std(tvs))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stem_a", help="Family stem (filename without _sNNN.mfd)")
    ap.add_argument("stem_b", help="Family stem to compare against")
    ap.add_argument("--seeds-dir", default=os.path.join(_ROOT, "seeds"))
    ap.add_argument("--dim", type=int, default=3)
    ap.add_argument("--no-hists", action="store_true",
                    help="Skip manifold loads / degree-distribution comparison.")
    ap.add_argument("--boot-seed", type=int, default=0)
    args = ap.parse_args()

    fams = {}
    for stem in (args.stem_a, args.stem_b):
        paths = sorted(glob.glob(os.path.join(args.seeds_dir, stem + "_s*.mfd")))
        if not paths:
            sys.exit(f"no replicas match {stem}_s*.mfd in {args.seeds_dir}")
        fams[stem] = paths
    with_hists = not args.no_hists

    obs, hists = {}, {}
    for stem, paths in fams.items():
        print(f"{stem}: {len(paths)} replicas", flush=True)
        obs[stem], vh, eh = replica_observables(paths, args.dim, with_hists)
        hists[stem] = (vh, eh)

    a, b = args.stem_a, args.stem_b
    print(f"\n{'observable':<10} {'A mean+/-sem':>22} {'B mean+/-sem':>22} "
          f"{'z':>7} {'KS D':>7} {'KS p':>7}")
    worst_z, worst_p = 0.0, 1.0
    for key in ("vdv", "hdv", "edge_deg", "f0", "f1", "f3"):
        xa, xb = obs[a][key], obs[b][key]
        ma, mb = xa.mean(), xb.mean()
        sa, sb = xa.std(ddof=1) / math.sqrt(len(xa)), xb.std(ddof=1) / math.sqrt(len(xb))
        denom = math.hypot(sa, sb)
        z = (ma - mb) / denom if denom > 0 else float("nan")
        d, pks = ks_two_sample(xa, xb)
        worst_z = max(worst_z, abs(z)) if not math.isnan(z) else worst_z
        worst_p = min(worst_p, pks)
        print(f"{key:<10} {ma:>13.6g} +/-{sa:<7.2g} {mb:>13.6g} +/-{sb:<7.2g} "
              f"{z:>7.2f} {d:>7.3f} {pks:>7.3f}")

    if with_hists:
        rng = np.random.default_rng(args.boot_seed)
        print()
        for name, idx in (("vertex-degree dist", 0), ("edge-degree dist", 1)):
            tv = tv_distance(hists[a][idx], hists[b][idx])
            mu, sd = tv_null(hists[a][idx], hists[b][idx], rng)
            sig = (tv - mu) / sd if sd > 0 else float("nan")
            print(f"{name}: TV={tv:.5f}  (resampling null {mu:.5f}+/-{sd:.5f}, "
                  f"{sig:+.1f} sigma)")

    print(f"\nsummary: worst |z|={worst_z:.2f}, worst KS p={worst_p:.3f} over "
          f"6 scalar observables ({len(obs[a]['vdv'])}+{len(obs[b]['vdv'])} replicas). "
          "Twins should show |z|<~3 and no tiny KS p beyond look-elsewhere.")


if __name__ == "__main__":
    main()
