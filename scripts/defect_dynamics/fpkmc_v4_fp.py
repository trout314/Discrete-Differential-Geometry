#!/usr/bin/env python3
"""V4: FP frozen-driver kinetics vs brute-force slide dynamics (M3).

Setup: two (3,4,4) knots on pristine R m2 -- knot A (mobile) and knot B
(FROZEN target) `--sep-sites` chain sites apart on the same BC orbit.
B's five window vertices are the blocked set; the FP domain is A's scan
ball with dock absorption one tet-layer outside B (plus the dS/depth
frontier and multichord absorbers).

FP side (one scan): exact splitting probabilities + exact mean
absorption time (fpkmc.FPFlight dense solves) and `--fp-flights` exact
jump-chain samples (exit histogram + FPT sample).

Brute side (`--brute-runs` independent runs, exact same law): the
sampler's pure slide channel with B frozen, in attempted-move units.
Per attempt, the probability of proposing A's current chord (with some
slot) is 3/(6*N3); attempts between such hits are geometric, at a hit
the slot is uniform over 12 and Metropolis runs on the trial dS
(slide_at2). This is an exact marginalization of the uniform
(facet, pair, slot) proposal -- proposals of non-deg-3 pairs reject, and
proposals of B's chord are frozen out of the model (v1 FROZEN
hierarchy, design section 5). A run absorbs when its state
(chord, round(S_rel, 6)) leaves the scan's interior set; the landing
state is matched to an absorbing node (interior completeness guarantees
one exists -- unmatched landings are counted loudly as bugs). Runs are
walked back to the start state through exact inverse slides (arrival-
chord matching, as in the HB driver) and audited.

PASS: chi^2 p-value of brute exit counts vs FP exact splitting > 0.01,
and two-sample KS p-value of FP vs brute FPT samples > 0.01.
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import fpkmc
from discrete_differential_geometry._sampler import SLIDE_SLOTS
import worm_slide as ws
from worm_helix import bc_orbit


def make_knot(m, orb, w):
    face = sorted(int(x) for x in orb[w + 1:w + 4])
    apx = tuple(sorted((int(orb[w]), int(orb[w + 4]))))
    m.do_bistellar_move(face, list(apx))
    return apx, set(face) | set(apx)


def build(args):
    m = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(m.facets())
    orb = bc_orbit(m, [int(x) for x in F[0]])
    wA = args.window
    wB = wA + 4 * args.sep_sites
    apxB, vertsB = make_knot(m, orb, wB)     # frozen target first
    apxA, vertsA = make_knot(m, orb, wA)
    if vertsA & vertsB:
        raise SystemExit("knot windows overlap; increase --sep-sites")
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=args.estar, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * args.estar / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
    return m, s, apxA, apxB, sorted(vertsB)


def walk_back(s, recs):
    """Invert committed slides in reverse, identifying each exact inverse
    by arrival chord + negated dS before committing (no impostors)."""
    for src_ch, arr_ch, dS in reversed(recs):
        want = -dS
        done = False
        for slot in range(SLIDE_SLOTS):
            try:
                d2, arr2 = s.slide_at2(arr_ch[0], arr_ch[1], slot,
                                       commit=False)
            except RuntimeError:
                continue
            if d2 is None or arr2 is None:
                continue
            if tuple(sorted(arr2)) != src_ch or abs(d2 - want) > 1e-9:
                continue
            s.slide_at(arr_ch[0], arr_ch[1], slot, commit=True)
            done = True
            break
        if not done:
            raise RuntimeError("V4 walk-back: inverse not found")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--window", type=int, default=40)
    ap.add_argument("--sep-sites", type=int, default=5)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--fp-flights", type=int, default=20000)
    ap.add_argument("--brute-runs", type=int, default=150)
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    m, s, apxA, apxB, blocked = build(args)
    v = s.manifold
    em0 = ws.edeg_dict(v)
    n3 = v.num_facets
    nu = fpkmc.nu_per_attempt(1.0, n3)
    print(f"A chord {apxA}, B chord {apxB} (frozen), sep "
          f"{args.sep_sites} sites; N3={n3}, nu={nu:.3e}", flush=True)

    # ---------------- FP side: one scan ----------------
    g = s.slide_graph_scan(apxA, dS_max=args.dS_max, max_depth=args.depth,
                           blocked_verts=blocked)
    fl = fpkmc.FPFlight(g, args.dS_max, args.depth, nu)
    n_int = int(fl.interior.sum())
    n_dock = int((fl.reason == fpkmc.ABSORB_DOCK).sum())
    print(f"scan: {len(g['dS'])} nodes ({n_int} interior, {n_dock} dock, "
          f"{len(g['edge_dS'])} edges)", flush=True)
    if not fl.interior[0]:
        raise SystemExit(f"start state is absorbing ({fl.reason[0]}); "
                         f"increase --sep-sites")
    split = fl.splitting_exact(0)
    t_mean = fl.mean_time_exact(0)
    rng_fp = np.random.default_rng(args.seed)
    fp_exits = Counter()
    fp_times = np.zeros(args.fp_flights, dtype=np.int64)
    for i in range(args.fp_flights):
        j, t, _ = fl.sample(rng_fp)
        fp_exits[j] += 1
        fp_times[i] = t
    print(f"FP: mean time exact {t_mean:.3e} attempts, sampled "
          f"{np.mean(fp_times):.3e}; P(dock) exact "
          f"{sum(p for a, p in split.items() if fl.reason[a] == 'dock'):.4f}",
          flush=True)
    # internal consistency: sampled exit law vs exact splitting
    for a, p in sorted(split.items(), key=lambda kv: -kv[1])[:5]:
        print(f"   exit {a:4d} [{fl.reason[a]:14s}] exact {p:.4f} "
              f"sampled {fp_exits[a] / args.fp_flights:.4f}", flush=True)

    # node lookup tables for the brute chain
    def ckey(j):
        return (tuple(sorted((int(g["chord"][j][0]), int(g["chord"][j][1])))),
                round(float(g["dS"][j]), 6))
    interior_key = {}
    absorb_key = {}
    for j in range(len(g["dS"])):
        k = ckey(j)
        tbl = interior_key if fl.interior[j] else absorb_key
        if k in tbl:
            raise SystemExit(f"state key collision at node {j} vs {tbl[k]}")
        tbl[k] = j
    absorb_by_S = {}
    for j in np.nonzero(~fl.interior)[0]:
        absorb_by_S.setdefault(round(float(g["dS"][int(j)]), 6),
                               []).append(int(j))

    # ---------------- brute side ----------------
    rng = np.random.default_rng(args.seed + 1)
    p_hit = 3.0 / (6.0 * n3)
    br_exits = Counter()
    br_times = []
    unmatched = 0
    for run in range(args.brute_runs):
        chord = apxA
        S_rel = 0.0
        t = 0
        recs = []
        while True:
            t += int(rng.geometric(p_hit))
            slot = int(rng.integers(SLIDE_SLOTS))
            try:
                dS, arr = s.slide_at2(chord[0], chord[1], slot,
                                      commit=False)
            except RuntimeError:
                dS, arr = None, None
            if dS is None:
                continue
            if not (dS <= 0 or rng.random() < np.exp(-dS)):
                continue
            s.slide_at(chord[0], chord[1], slot, commit=True)
            new_chord = tuple(sorted(int(x) for x in arr))
            recs.append((chord, new_chord, float(dS)))
            chord = new_chord
            S_rel += float(dS)
            k = (chord, round(S_rel, 6))
            if k in interior_key:
                continue
            if k in absorb_key:
                br_exits[absorb_key[k]] += 1
            else:
                # multichord absorbers: rep chord != arrival chord;
                # match by action level if unique
                cands = absorb_by_S.get(round(S_rel, 6), [])
                if len(cands) == 1:
                    br_exits[cands[0]] += 1
                else:
                    unmatched += 1
                    print(f"  UNMATCHED landing run {run}: {k} "
                          f"({len(cands)} S-candidates)", flush=True)
            br_times.append(t)
            break
        walk_back(s, recs)
        if run % 10 == 9:
            true = args.lam * ws.dS_between(em0, ws.edeg_dict(v),
                                            estar=args.estar)
            assert abs(true) < 1e-6, "restore audit failed"
            print(f"  brute {run + 1}/{args.brute_runs}: "
                  f"median t {int(np.median(br_times))}", flush=True)

    # ---------------- compare ----------------
    print(f"\nbrute: {args.brute_runs} runs, {unmatched} unmatched")
    nb = sum(br_exits.values())
    # chi^2 vs exact splitting, pooling exits with expected < 5
    exp, obs, pool_e, pool_o = [], [], 0.0, 0
    for a, p in sorted(split.items(), key=lambda kv: -kv[1]):
        e = p * nb
        o = br_exits.get(a, 0)
        if e < 5:
            pool_e += e
            pool_o += o
        else:
            exp.append(e)
            obs.append(o)
    if pool_e > 0:
        exp.append(pool_e)
        obs.append(pool_o)
    exp, obs = np.array(exp), np.array(obs, dtype=float)
    chi2 = float(((obs - exp) ** 2 / exp).sum())
    dof = max(1, len(exp) - 1)
    try:
        from scipy import stats as sps
        p_chi = float(sps.chi2.sf(chi2, dof))
        ks = sps.ks_2samp(fp_times, np.array(br_times))
        p_ks, ks_stat = float(ks.pvalue), float(ks.statistic)
    except ImportError:
        p_chi = p_ks = float("nan")
        ks_stat = float("nan")
        print("(scipy unavailable: reporting statistics without p-values)")
    print(f"exit law: chi2 {chi2:.2f} / dof {dof}  p {p_chi:.3f}")
    print(f"FPT: KS {ks_stat:.4f}  p {p_ks:.3f}  "
          f"(FP mean {np.mean(fp_times):.3e}, brute mean "
          f"{np.mean(br_times):.3e}, exact {t_mean:.3e})")
    verdict = "PASS" if (p_chi > 0.01 and p_ks > 0.01 and unmatched == 0) \
        else "FAIL"
    print(f"\nV4 seed {args.seed}: {verdict}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"verdict": verdict, "chi2": chi2, "dof": dof,
                       "p_chi": p_chi, "ks": ks_stat, "p_ks": p_ks,
                       "t_mean_exact": t_mean,
                       "fp_mean": float(np.mean(fp_times)),
                       "brute_mean": float(np.mean(br_times)),
                       "unmatched": unmatched,
                       "split": {str(a): p for a, p in split.items()},
                       "brute_exits": {str(a): c
                                       for a, c in br_exits.items()},
                       "brute_times": br_times}, fh)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
