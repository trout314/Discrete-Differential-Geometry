#!/usr/bin/env python3
"""Pool several contact_census.py chain outputs and report by rung.

The rung Q is conserved for a whole flight episode (no Q-refresh), and the
per-rung free network is NOT uniform -- some rungs carry a crystal-spanning
free web, others fragment or do not span at all. Pooling rungs averages a
mobile sector together with trapped ones and hides exactly the sector
dependence that decides whether a contact handoff has anything to hand to.

Usage:
    python scripts/defect_dynamics/contact_census_merge.py \
        data/ecmc_ab/contact_census_*.json --out data/ecmc_ab/contact_census.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.normpath(os.path.join(              # archived:
    os.path.dirname(os.path.abspath(__file__)), "..", "..",   # siblings
    "scripts", "defect_dynamics")))                            # stayed in dd/
from contact_census import report


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--min-n", type=int, default=100,
                    help="minimum rejected contacts to report a rung")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows, episodes = [], []
    for f in args.jsons:
        d = json.load(open(f))
        rows += d["rows"]
        episodes += d.get("episodes", [])
    print(f"pooled {len(args.jsons)} chains: {len(rows):,} rejected contacts, "
          f"{len(episodes):,} episodes")

    summ = {"all": report(rows, "ALL rejected contacts (pooled)")}

    qs = collections.Counter(r["Q"] for r in rows)
    print("\n--- rung occupancy (share of rejected contacts): "
          + ", ".join(f"Q={q}: {100*n/len(rows):.1f}%"
                      for q, n in sorted(qs.items())))
    summ["by_Q"] = {}
    for q, n in sorted(qs.items()):
        if n >= args.min_n:
            summ["by_Q"][str(q)] = report([r for r in rows if r["Q"] == q],
                                          f"RUNG Q={q}")

    # Is a trapped sector INTRINSIC (the rung's free web does not span) or an
    # artefact of blob CROWDING (every site's support brushes a neighbouring
    # flier while the blob is still packed)? Cross rung against time: crowding
    # relaxes as the blob disperses, web topology does not.
    half = max(r["sw"] for r in rows) // 2
    print(f"\n--- free-continuation rate by rung x time  (split at sw={half})")
    print(f"    {'Q':>4} {'n_early':>8} {'resume':>8}   {'n_late':>8} {'resume':>8}"
          f"   {'verdict':>12}")
    xtab = {}
    for q in sorted(qs):
        E = [r for r in rows if r["Q"] == q and r["sw"] <= half]
        L = [r for r in rows if r["Q"] == q and r["sw"] > half]
        if len(E) < 30 or len(L) < 30:
            continue
        re_ = np.mean([r["resume_k"] is not None for r in E])
        rl = np.mean([r["resume_k"] is not None for r in L])
        v = ("crowding" if rl - re_ > 0.15 else
             "intrinsic" if max(re_, rl) < 0.25 else "spans")
        xtab[str(q)] = {"early": re_, "late": rl, "verdict": v}
        print(f"    {q:>4} {len(E):>8,} {re_:>7.1%}   {len(L):>8,} {rl:>7.1%}"
              f"   {v:>12}")
    summ["rung_x_time"] = xtab

    per_q = collections.defaultdict(collections.Counter)
    nep = collections.Counter()
    for e in episodes:
        per_q[e["Q"]].update(e["events"])
        nep[e["Q"]] += 1
    print("\n--- per-rung kernel-step budget (share of that rung's steps)")
    print(f"    {'Q':>4} {'eps':>6} {'steps':>9} {'free':>8} {'rejctc':>8} "
          f"{'accctc':>8} {'accdeep':>8} {'noslot':>8}")
    tab = {}
    for q in sorted(per_q):
        c = per_q[q]
        st = sum(v for k, v in c.items()
                 if k.startswith(("accept_", "reject_"))
                 or k in ("noslot_flip", "stuck_flip", "illegal_flip"))
        if not st:
            continue
        tab[str(q)] = {"episodes": nep[q], "steps": st,
                       "free": c["accept_free"] / st,
                       "reject_contact": c["reject_contact_flip"] / st,
                       "accept_contact": c["accept_contact"] / st,
                       "accept_deeper": c["accept_deeper"] / st,
                       "noslot": c["noslot_flip"] / st}
        t = tab[str(q)]
        print(f"    {q:>4} {nep[q]:>6} {st:>9,} {t['free']:>7.1%} "
              f"{t['reject_contact']:>7.1%} {t['accept_contact']:>7.1%} "
              f"{t['accept_deeper']:>7.1%} {t['noslot']:>7.1%}")
    summ["per_rung_steps"] = tab

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"summary": summ, "n_chains": len(args.jsons)}, f,
                      indent=1, default=float)
        print(f"\nwrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
