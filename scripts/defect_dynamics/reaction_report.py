#!/usr/bin/env python3
"""Analyse the thermal reaction census (reaction_census.py output).

Tests, in order of how much they would change the picture:

 1. IS "REACTION" THE RIGHT LANGUAGE?  Merge is a two-body encounter and
    split a one-body decay, so at steady state
        k_merge = k2 * <n(n-1)/2> / V      k_split = k1 * <n_bound>
    with k1, k2 intensive.  Reduced rates that agree across the density
    ladder (lam=0.40 vs 0.35, differing ~2.5x in complex number) support
    the reaction picture; ones that do not falsify it.

 2. TRANSPORT- OR REACTION-LIMITED?  The slide channel satisfies detailed
    balance, so slide on/off leaves the stationary ensemble alone and
    changes only kinetics.  k2 rising with the slide channel means the
    chemistry is diffusion-limited.

 3. ASSOCIATION CONSTANT K = k_merge/k_split.  The static two-body
    potential is exactly zero beyond contact, so a K that does not match
    the ideal (non-interacting) value is entropic or kinetic in origin.

 4. COMPOUND LIFETIMES.  Survival curves for tracks that absorbed a
    partner vs those that never did -- the quantitative form of the
    run5h observation that fused multimers are long-lived.

Rates carry block-bootstrap errors (CONVENTIONS sec 3): the run is cut
into blocks, a rate computed per block, and the block set resampled.
Burn-in is dropped by --skip-sweeps.
"""
import argparse
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(_HERE, "..", "..")


def block_boot(counts, spans, nboot=4000, rng=None):
    """Rate = sum(counts)/sum(spans) with a block bootstrap over blocks."""
    rng = rng or np.random.default_rng(0)
    counts = np.asarray(counts, float)
    spans = np.asarray(spans, float)
    if len(counts) == 0 or spans.sum() == 0:
        return float("nan"), float("nan")
    est = counts.sum() / spans.sum()
    idx = rng.integers(0, len(counts), size=(nboot, len(counts)))
    b = counts[idx].sum(axis=1) / np.maximum(spans[idx].sum(axis=1), 1e-12)
    return float(est), float(b.std())


def load_run(base):
    """Load one reaction_census run, or None if `base` is some other JSON
    that happens to share the directory."""
    try:
        with open(base + ".json") as fh:
            d = json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(d, dict) or "counters" not in d or "ts" not in d:
        return None
    ev = []
    p = base + ".events.jsonl"
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    try:
                        ev.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass          # truncated final line of a live run
    tr = []
    p = base + ".tracks.jsonl"
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    try:
                        tr.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass
    d["_events"], d["_tracks"] = ev, tr
    d["_name"] = os.path.basename(base)

    # CLEAN WINDOW. DefectState's incremental bookkeeping diverges from a
    # from-scratch rebuild in long runs (first seen 2026-07-27; audits at
    # sw 8k-65k depending on the run). Worldline components are built from
    # that state, so every event after the divergence is suspect. The audit
    # runs every --audit-every chunks, so the last VERIFIED-clean point is
    # one audit interval before the first failure. Truncate there, always:
    # a report must not be able to include post-divergence data by default.
    d["_clean_to"] = d.get("sweeps", 0)
    d["_truncated"] = False
    logp = base + ".log"
    if os.path.exists(logp):
        with open(logp) as fh:
            fails = [int(m.group(1)) for m in
                     re.finditer(r"AUDIT FAILED sw(\d+)", fh.read())]
        if fails:
            d["_clean_to"] = max(0, min(fails) - 1000)
            d["_truncated"] = True
    return d


def analyse(d, skip, nblock):
    """Reduced rates for one run, over the post-burn-in window."""
    cap = d["_clean_to"]
    ts = [r for r in d["ts"] if skip <= r["sw"] <= cap]
    if len(ts) < 2:
        return None
    sw0, sw1 = ts[0]["sw"], ts[-1]["sw"]
    ev = [e for e in d["_events"]
          if skip <= e.get("clock_sw", -1) <= cap]

    # complex-number series, and the pair count that a bimolecular rate
    # should be proportional to
    n = np.array([r["n_components"] for r in ts], float)
    edges = np.linspace(sw0, sw1, nblock + 1)
    mcount, scount, spans, npair, nmean = [], [], [], [], []
    for i in range(nblock):
        lo, hi = edges[i], edges[i + 1]
        sel = [r for r in ts if lo <= r["sw"] < hi]
        if not sel:
            continue
        nn = np.array([r["n_components"] for r in sel], float)
        mcount.append(sum(1 for e in ev if e["k"] == "merge"
                          and lo <= e["clock_sw"] < hi))
        scount.append(sum(1 for e in ev if e["k"] == "split"
                          and lo <= e["clock_sw"] < hi))
        spans.append(hi - lo)
        npair.append(float((nn * (nn - 1) / 2).mean()))
        nmean.append(float(nn.mean()))

    npair = np.array(npair)
    nmean = np.array(nmean)
    km, kme = block_boot(mcount, spans)
    ks, kse = block_boot(scount, spans)
    # reduced (intensive) rates: divide each block's count by that
    # block's mean pair / mean complex count before pooling
    k2, k2e = block_boot(mcount, np.array(spans) * npair)
    k1, k1e = block_boot(scount, np.array(spans) * nmean)
    return {
        "name": d["_name"], "lam": d["lam"], "slide": d["slide_prob"],
        "sweeps": sw1 - sw0, "n_mean": float(n.mean()),
        "n_pair": float(npair.mean()),
        "k_merge": km, "k_merge_e": kme, "k_split": ks, "k_split_e": kse,
        "k2": k2, "k2_e": k2e, "k1": k1, "k1_e": k1e,
        "K": km / ks if ks else float("nan"),
        "n_merge": sum(mcount), "n_split": sum(scount),
        "clean_to": cap, "truncated": d["_truncated"],
        "total_sweeps": d.get("sweeps", 0),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=os.path.join(ROOT, "data", "reactions"))
    ap.add_argument("--skip-sweeps", type=float, default=500.0)
    ap.add_argument("--nblock", type=int, default=20)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    print("[frame] none -- graph identity only (CONVENTIONS sec 6)")

    bases = sorted(p[:-5] for p in glob.glob(os.path.join(args.dir, "*.json")))
    if not bases:
        sys.exit(f"no runs in {args.dir}")
    runs = [r for r in (load_run(b) for b in bases) if r]
    if not runs:
        sys.exit(f'no reaction_census runs in {args.dir}')
    rows = [r for r in (analyse(d, args.skip_sweeps, args.nblock)
                        for d in runs) if r]

    print(f"\n{len(rows)} runs, burn-in {args.skip_sweeps} sweeps dropped")
    tr = [r for r in rows if r["truncated"]]
    if tr:
        print("TRUNCATED at the last verified-clean audit (DefectState\n"
              "incremental divergence; events past it are discarded):")
        for r in tr:
            print(f"    {r['name']:6s} using {r['clean_to']:6d} of "
                  f"{r['total_sweeps']:6d} sweeps "
                  f"({100 * r['clean_to'] / max(r['total_sweeps'], 1):3.0f}%)")
    print()
    print(f"{'run':6s} {'lam':>5s} {'slide':>6s} {'sweeps':>8s} {'<n>':>6s} "
          f"{'merges':>7s} {'splits':>7s} {'k_m/sw':>16s} {'k_s/sw':>16s}")
    for r in rows:
        print(f"{r['name']:6s} {r['lam']:5.2f} {r['slide']:6.2f} "
              f"{r['sweeps']:8.0f} {r['n_mean']:6.2f} "
              f"{r['n_merge']:7d} {r['n_split']:7d} "
              f"{r['k_merge']:8.4f}+-{r['k_merge_e']:<6.4f} "
              f"{r['k_split']:8.4f}+-{r['k_split_e']:<6.4f}")

    print("\nreduced rates -- intensive if the reaction picture holds")
    print(f"{'run':6s} {'lam':>5s} {'slide':>6s} "
          f"{'k2 = k_m/<pairs>':>22s} {'k1 = k_s/<n>':>22s} {'K':>8s}")
    for r in rows:
        print(f"{r['name']:6s} {r['lam']:5.2f} {r['slide']:6.2f} "
              f"{r['k2']:12.6f}+-{r['k2_e']:<8.6f} "
              f"{r['k1']:12.6f}+-{r['k1_e']:<8.6f} {r['K']:8.3f}")

    # --- verdicts -----------------------------------------------------
    def grp(pred):
        return [r for r in rows if pred(r)]
    off = grp(lambda r: r["slide"] == 0)
    on = grp(lambda r: r["slide"] > 0)
    print("\nVERDICTS")
    for lam in sorted({r["lam"] for r in rows}):
        g = [r for r in off if r["lam"] == lam]
        if g:
            k2 = np.mean([r["k2"] for r in g])
            print(f"  lam={lam:.2f} slide off: <n>={np.mean([r['n_mean'] for r in g]):.1f}, "
                  f"k2={k2:.6f}")
    if len(off) >= 2:
        lams = sorted({r["lam"] for r in off})
        if len(lams) >= 2:
            a = np.mean([r["k2"] for r in off if r["lam"] == lams[0]])
            b = np.mean([r["k2"] for r in off if r["lam"] == lams[-1]])
            print(f"  density test: k2 ratio across the ladder = "
                  f"{a / b:.2f} (1.0 = bimolecular picture holds)")
    for lam in sorted({r["lam"] for r in rows}):
        a = [r["k2"] for r in off if r["lam"] == lam]
        b = [r["k2"] for r in on if r["lam"] == lam]
        if a and b:
            print(f"  lam={lam:.2f} transport test: k2(slide on)/k2(off) = "
                  f"{np.mean(b) / np.mean(a):.2f} "
                  f"(>1 = diffusion-limited chemistry)")

    # --- compound lifetimes -------------------------------------------
    print("\ncompound lifetimes (tracks that absorbed >=1 partner vs not)")
    for d in runs:
        tr = [t for t in d["_tracks"]
              if t.get("born_sw", 0) + t.get("life_sw", 0) <= d["_clean_to"]]
        if not tr:
            continue
        comp = [t["life_sw"] for t in tr if t.get("n_merge_in", 0) > 0]
        solo = [t["life_sw"] for t in tr if t.get("n_merge_in", 0) == 0
                and t.get("ended") != "merge"]
        if comp and solo:
            print(f"  {d['_name']:6s}: compounds n={len(comp):5d} "
                  f"median {np.median(comp):8.2f} sw, mean {np.mean(comp):9.2f} | "
                  f"solo n={len(solo):6d} median {np.median(solo):6.2f} sw")

    # --- species table -------------------------------------------------
    print("\nreaction species table (top channels, pooled)")
    shape_of = {}
    for d in runs:
        for t in d["_tracks"]:
            if t.get("born_sw", 0) > d["_clean_to"]:
                continue
            if t.get("shape"):
                shape_of[(d["_name"], t["tid"])] = tuple(t["shape"])
    chan = Counter()
    for d in runs:
        for e in d["_events"]:
            if e["k"] != "merge" or e.get("clock_sw", 0) > d["_clean_to"]:
                continue
            b = tuple(sorted(shape_of.get((d["_name"], f[0]))
                             for f in e["from"]
                             if shape_of.get((d["_name"], f[0]))))
            if b:
                chan[(b, e["size"])] += 1
    for (b, sz), k in chan.most_common(12):
        print(f"  {k:5d}  absorbed {list(b)} -> product size {sz}")

    # --- figure ---------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    ax = axs[0][0]
    for d in runs:
        ts = [r for r in d["ts"]
              if args.skip_sweeps <= r["sw"] <= d["_clean_to"]]
        ax.plot([r["sw"] for r in ts], [r["n_components"] for r in ts],
                lw=0.8, label=d["_name"])
    ax.set_xlabel("sweep"); ax.set_ylabel("n complexes")
    ax.set_title("complex number (stationarity check)")
    ax.legend(fontsize=7)

    ax = axs[0][1]
    x = np.arange(len(rows))
    ax.bar(x - 0.2, [r["k2"] for r in rows], 0.4,
           yerr=[r["k2_e"] for r in rows], capsize=3, label="k2 (merge/pair)")
    ax.bar(x + 0.2, [r["k1"] for r in rows], 0.4,
           yerr=[r["k1_e"] for r in rows], capsize=3, label="k1 (split/complex)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['name']}\nlam{r['lam']:.2f} sl{r['slide']:.2f}"
                        for r in rows], fontsize=7)
    ax.set_ylabel("reduced rate (per sweep)")
    ax.set_title("reduced rates -- flat across a row = intensive")
    ax.legend(fontsize=8)

    ax = axs[1][0]
    for d in runs:
        tr = [t for t in d["_tracks"]
              if t.get("born_sw", 0) + t.get("life_sw", 0) <= d["_clean_to"]]
        life = sorted(t["life_sw"] for t in tr
                      if t.get("n_merge_in", 0) > 0)
        if len(life) > 5:
            ax.step(life, 1 - np.arange(len(life)) / len(life), where="post",
                    lw=1.2, label=f"{d['_name']} compounds")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("lifetime (sweeps)"); ax.set_ylabel("survival")
    ax.set_title("compound survival")
    ax.legend(fontsize=7)

    ax = axs[1][1]
    for d in runs:
        h = d.get("life_hist") or {}
        if h:
            k = np.array(sorted(int(x) for x in h))
            vv = np.array([h[str(x)] for x in k], float)
            ax.step(np.maximum(k, 0.5), vv / vv.sum(), where="post", lw=1.0,
                    label=d["_name"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("lifetime (sweeps)"); ax.set_ylabel("fraction of tracks")
    ax.set_title("all-track lifetime distribution (flicker-dominated)")
    ax.legend(fontsize=7)

    fig.suptitle(
        "Thermal reaction census -- R m4 N57984, e*=5.105025, full Pachner "
        "dynamics; merge/split of defect complexes at accepted-move "
        "resolution", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = args.out or os.path.join(args.dir, "reaction_report.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
