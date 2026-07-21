#!/usr/bin/env python3
"""Cross-crystal summary of a reference-state tier: pool each fleet's replica
finals into mean+/-SEM of the observables that tell the story --- FK-legality,
curvature (mean edge degree), disclination-web topology, and S(k)
hyperuniformity (via cocycle harmonic coords where a cocycle exists).

Writes out/reference_summary_<tier>.csv and prints a table.

Usage:
    python scripts/reference_summary.py --tier small
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.structure_factor import structure_factor
from discrete_differential_geometry.vertex_fields import FIELDS
from fk_skeleton import edges_from_facets, vertex_class_census

FAMS = ["a15", "sigma", "r", "c15"]           # ordered above-flat -> below-flat
FLEETS = ["crystal", "flat_vac", "ownpin", "glass"]


def fleet_finals(tier, fam, fleet):
    d = f"data/reference/{tier}"
    # glass is 2-stage: the annealed state is s1
    pats = [f"{d}/{fam}_{fleet}_{tier}_r*_final.mfd",
            f"{d}/{fam}_{fleet}_{tier}_s1_r*_final.mfd"]
    files = []
    for p in pats:
        files += glob.glob(p)
    return sorted(files)


def analyze_final(path):
    v = ddg.Manifold.load(path, 3)
    facets = np.asarray(v.facets())
    eu, edeg, V = edges_from_facets(facets)
    fz, _ = vertex_class_census(eu, edeg, V)
    c = v.disclination_census()
    row = dict(pure56=1.0 - fz["impure"], qbar=float(edeg.mean()),
               n_comp=c["n_components"], giant_frac=c["giant_frac"],
               n_seg=c["n_segments"], mean_seg=c["mean_seg_len"],
               fray=c["n_fray_verts"], n_six=c["n_six_edges"])
    cpath = re.sub(r"\.mfd$", ".cocycle.npz", path)
    if os.path.exists(cpath):
        try:
            edges, omega, _ = coc.load_cocycle(cpath)
            frac, basis = coc.torus_positions(facets, edges, omega)
            kmag, s_obs, s_null = structure_factor(frac, basis, FIELDS["n6"](facets), 6)
            low = kmag <= 2.0 + 1e-9
            row["sk_lowk"] = float(np.mean(s_obs[low] / s_null[low]))
        except Exception as e:
            row["sk_lowk"] = float("nan")
    else:
        row["sk_lowk"] = float("nan")
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tier", default="small")
    args = ap.parse_args()

    cols = ["pure56", "qbar", "sk_lowk", "giant_frac", "n_comp",
            "mean_seg", "n_seg", "fray", "n_six"]
    out_rows = []
    print(f"{'crystal':7} {'state':9} {'nrep':>4} {'pure56':>7} {'qbar':>8} "
          f"{'S(k)lowk':>9} {'giantf':>7} {'ncomp':>6} {'mseg':>6} {'fray':>6}")
    for fam in FAMS:
        for fleet in FLEETS:
            files = fleet_finals(args.tier, fam, fleet)
            if not files:
                continue
            rows = [analyze_final(f) for f in files]
            agg = {}
            for c in cols:
                vals = np.array([r[c] for r in rows], float)
                vals = vals[~np.isnan(vals)]
                agg[c + "_mean"] = float(vals.mean()) if len(vals) else float("nan")
                agg[c + "_sem"] = (float(vals.std(ddof=1) / np.sqrt(len(vals)))
                                   if len(vals) > 1 else 0.0)
            agg.update(crystal=fam, state=fleet, nrep=len(files))
            out_rows.append(agg)
            print(f"{fam:7} {fleet:9} {len(files):>4} "
                  f"{agg['pure56_mean']:7.3f} {agg['qbar_mean']:8.4f} "
                  f"{agg['sk_lowk_mean']:9.4f} {agg['giant_frac_mean']:7.3f} "
                  f"{agg['n_comp_mean']:6.1f} {agg['mean_seg_mean']:6.2f} "
                  f"{agg['fray_mean']:6.1f}")

    os.makedirs("out", exist_ok=True)
    outp = f"out/reference_summary_{args.tier}.csv"
    keys = ["crystal", "state", "nrep"] + [c + s for c in cols
                                           for s in ("_mean", "_sem")]
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nwrote {outp}")


if __name__ == "__main__":
    main()
