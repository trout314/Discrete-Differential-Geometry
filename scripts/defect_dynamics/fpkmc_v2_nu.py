#!/usr/bin/env python3
"""V2: measure the slide-proposal rate constant against the nu formula.

The FPKMC time bookkeeping (and nothing else -- M0 result E) depends on
nu: the per-attempt probability of proposing one specific (chord, slot),

    nu = slide_prob * (3/N3) * (1/6) * (1/12)

read off sampler.d's slide branch. A wrong constant silently rescales all
kinetics (design R2), so it is MEASURED here, not trusted:

  tries rate    slideCfg.tries increments exactly when a proposal forms a
                legal slide, so with one knot present (n_legal = 12
                measured at every site) the prediction is
                    tries/attempt = slide_prob * (3/N3) * (1/6)
                -- independent of Metropolis, dS, or lambda.
  accept ratio  accepts/tries ~ <min(1, e^-dS)> over the menu of the
                sites the knot actually visits; the survey gives the
                starting site's value as a sanity band.

Protocol: create the knot, run in short chunks, count (attempts, tries,
accepts) only over chunks where the knot's chord is still degree 3
(thermal churn eventually morphs or kills it; those chunks are dropped).
Flicker complexes carry their own degree-3 chords and inflate the tries
rate slightly; the contamination is estimated and reported rather than
ignored (flicker lives ~2 sweeps and holds ~1 deg-3 chord when present).
"""
import argparse
import sys
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry._dlang import _lib
import ctypes
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit
from cocycle_check import reference_frac_positions


def chord_degree(s, a, b):
    arr = (ctypes.c_int * 2)(a, b)
    try:
        return int(_lib.ddg_sampler_degree(s._handle, arr, 2))
    except RuntimeError:
        return 0            # edge no longer exists: the knot died


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--slide-prob", type=float, default=0.5)
    ap.add_argument("--window", type=int, default=1624)
    ap.add_argument("--lives", type=int, default=120)
    ap.add_argument("--chunk-moves", type=int, default=4000)
    ap.add_argument("--freeze-radius", type=float, default=1.2,
                    help="cells: freeze every vertex farther than this from "
                         "the knot chord. Kills flicker everywhere (no "
                         "contaminating deg-3 chords) and pins the tries "
                         "prediction exactly; frozen-rejected slides do NOT "
                         "count as tries (anyFrozen precedes valid=true).")
    ap.add_argument("--seed", type=int, default=7)
    add_chain_args(ap, default=None)
    args = ap.parse_args()

    ddg.set_random_seed(args.seed)
    m0 = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(m0.facets())
    _cc, _kcls, _seq, chain_prov = chain_for_run(
        args.ref, F, args.chain_class, seed_tet=0)
    orb = _seq
    k = args.window
    face = sorted(int(x) for x in orb[k + 1:k + 4])
    apx = sorted((int(orb[k]), int(orb[(k + 4) % len(orb)])))
    # survey the start site on a PRISTINE sampler (the survey creates the
    # knot itself; running it on an already-knotted state is invalid)
    params0 = ddg.SamplerParams(
        num_facets_target=m0.num_facets, num_facets_coef=0.1,
        hinge_degree_target=args.estar, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * args.estar / 6.0)
    s0 = ddg.ManifoldSampler(m0, params0)
    s0.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
    sv = s0.site_survey([apx[0]] + face + [apx[1]])
    dS0 = sv["slot_dS"][0]
    legal0 = ~np.isnan(dS0)
    acc0 = float(np.minimum(1.0, np.exp(-dS0[legal0])).mean())
    print(f"start site: n_legal = {int(legal0.sum())}, "
          f"<min(1,e^-dS)> = {acc0:.4f}")
    del s0, m0

    n3 = None
    att_tot = tries_tot = acc_tot = 0
    lives = 0
    while lives < args.lives:
        m = ddg.Manifold.load(args.ref, 3)
        m.do_bistellar_move(face, apx)
        params = ddg.SamplerParams(
            num_facets_target=m.num_facets, num_facets_coef=0.1,
            hinge_degree_target=args.estar, num_hinges_coef=0.0,
            hinge_degree_variance_coef=0.0,
            codim3_degree_variance_coef=0.0,
            hinge_degree_target_coef=args.lam * args.estar / 6.0)
        s = ddg.ManifoldSampler(m, params)
        s.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
        s.set_slide_prob(args.slide_prob)
        if n3 is None:
            n3 = s.manifold.num_facets
        # freeze the complement of the knot patch ON THE SAMPLER'S copy
        if args.freeze_radius > 0:
            rp = np.asarray(reference_frac_positions("r", 4))
            mid = (rp[apx[0]] + rp[apx[1]]) / 2.0
            d = rp - mid
            d -= np.round(d / 4.0) * 4.0
            far = np.where((d ** 2).sum(1) > args.freeze_radius ** 2)[0]
            far = [int(v) for v in far]
            mh = _lib.ddg_sampler_get_manifold(s._handle)
            arr = (ctypes.c_int * len(far))(*far)
            rc = _lib.ddg_manifold_freeze_vertices(mh, arr, len(far), 1)
            assert rc == 0
        alive = True
        while alive:
            t0, a0 = s.slide_stats()
            s.run(moves=args.chunk_moves)
            t1, a1 = s.slide_stats()
            if chord_degree(s, apx[0], apx[1]) == 3:
                att_tot += args.chunk_moves
                tries_tot += t1 - t0
                acc_tot += a1 - a0
            else:
                alive = False
        lives += 1
        if lives % 20 == 0:
            r = tries_tot / max(att_tot, 1)
            print(f"  {lives} lives: {att_tot} attempts, {tries_tot} tries "
                  f"({r:.3e}/attempt), {acc_tot} accepts", flush=True)
        del s, m
        _lib.ddg_gc_collect(); _lib.ddg_gc_minimize()

    pred = args.slide_prob * (3.0 / n3) / 6.0
    r = tries_tot / att_tot
    se = np.sqrt(tries_tot) / att_tot
    print(f"\nN3 = {n3}")
    print(f"tries/attempt: measured {r:.4e} +- {se:.1e}")
    print(f"               predicted {pred:.4e}  "
          f"(slide_prob*(3/N3)*(1/6), n_legal = 12/12)")
    print(f"ratio measured/predicted = {r / pred:.3f} "
          f"+- {se / pred:.3f}")
    if acc_tot:
        print(f"accepts/tries: {acc_tot}/{tries_tot} = "
              f"{acc_tot / max(tries_tot, 1):.3f} "
              f"(start-site prediction {acc0:.3f}; drifts as the knot "
              f"visits other sites)")
    if args.freeze_radius > 0:
        print("\n(frozen-complement mode: no flicker, prediction exact -- "
              "the ratio should be 1.000 within errors)")
    else:
        print("\nNOTE: flicker chords inflate the tries rate; a ratio "
              "modestly above 1 is expected.")


if __name__ == "__main__":
    main()
