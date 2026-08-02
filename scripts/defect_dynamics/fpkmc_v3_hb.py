#!/usr/bin/env python3
"""V3: HB equilibrium vs brute-force slide-channel equilibrium vs Boltzmann.

Three ways to the single-defect equilibrium on R m2, which must agree:

  HB       the Metropolized-ball driver (fpkmc.HBDriver), audited;
  brute    an exact Python emulation of the sampler's slide channel:
           uniform facet, uniform vertex pair, uniform slot, Metropolis
           accept via slide_at trial + force. (The p=1.0 trick in-process:
           no thermal channel exists here at all, so the defect is
           immortal and the chain is exactly the slide-channel dynamics.)
  exact    per-state Boltzmann weights e^{-S} with S from the running
           committed dS (both chains carry an audited running action).

Comparison: per-state occupancy (dwell-weighted for brute, visit-weighted
for HB -- both chains are reversible w.r.t. pi, so both histograms
estimate pi) against e^{-S}, as ln(occupancy ratio) vs -(S_i - S_j) for
all well-visited state pairs. Slope 1, intercept 0 within errors = PASS.

States are keyed by (chord, round(S_rel, 6)): the chord alone can recur
with different decorations; adding the running action disambiguates
(collision would need two states with equal chord AND equal action --
and any keying error would surface in the audits, which recompute the
true action from the edge map).
"""
import argparse
import json
import os
import sys
from collections import Counter

# The D-side graph scan leaks ~12 MB per call under the conservative GC
# (never reclaimed; GC.collect measured ineffective -- same class as the
# Manifold-load leak, see FPKMC_DESIGN R-notes). Until the scan gets
# allocation-free scratch buffers, the HB phase EXEC-RESTARTS itself every
# --hb-batch steps: exact state (manifold file + chain bookkeeping + RNG
# state) is saved and reloaded, memory resets to zero, statistics are
# unaffected.

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import fpkmc
import worm_slide as ws
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit


def build(ref, estar, lam, window=None):
    m = ddg.Manifold.load(ref, 3)
    F = np.asarray(m.facets())
    _cc, _kcls, _seq, chain_prov = chain_for_run(
        ref, F, None, seed_tet=0)
    orb = _seq
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=estar, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * estar / 6.0)
    if window is None:
        s0 = ddg.ManifoldSampler(m, p)
        s0.set_n6_potential(0.6 * lam, 1.0 * lam, tilt=[0.0] * 5)
        sv = s0.site_survey(orb[:200])
        window = int(np.nanargmin(sv["dS_create"]))
        del s0
    face = sorted(int(x) for x in orb[window + 1:window + 4])
    apx = tuple(sorted((int(orb[window]), int(orb[window + 4]))))
    m.do_bistellar_move(face, list(apx))
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.6 * lam, 1.0 * lam, tilt=[0.0] * 5)
    return m, s, apx, window


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--hb-steps", type=int, default=1500)
    ap.add_argument("--hb-batch", type=int, default=250,
                    help="HB steps per process (exec-restart to shed the "
                         "D scan leak)")
    ap.add_argument("--state", default=None,
                    help="internal: resume-state prefix")
    ap.add_argument("--brute-props", type=int, default=4_000_000,
                    help="brute slide-channel proposals (only ~1/(2 N3) hit "
                         "the chord, so this is ~270 tries)")
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    state_prefix = args.state or (args.out or "v3_hb") + ".state"
    st_json = state_prefix + ".json"
    st_mfd = state_prefix + ".mfd"

    # ---------------- HB chain (exec-restart batches) ----------------
    if os.path.exists(st_json):
        with open(st_json) as fh:
            stj = json.load(fh)
        m = ddg.Manifold.load(st_mfd, 3)
        pp = ddg.SamplerParams(
            num_facets_target=stj["n3_target"], num_facets_coef=0.1,
            hinge_degree_target=args.estar, num_hinges_coef=0.0,
            hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
            hinge_degree_target_coef=args.lam * args.estar / 6.0)
        s = ddg.ManifoldSampler(m, pp)
        s.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
        apx = tuple(stj["chord"])
        window = stj["window"]
        hb = Counter({(tuple(json.loads(k)[0]), json.loads(k)[1]): v
                      for k, v in stj["hb"].items()})
        done = stj["done"]
        rng.bit_generator.state = stj["rng"]
        S0 = stj["S_rel"]
        prev_acc = stj["naccept"]
    else:
        m, s, apx, window = build(args.ref, args.estar, args.lam)
        hb = Counter()
        done = 0
        S0 = 0.0
        prev_acc = 0
    v = s.manifold
    em0 = ws.edeg_dict(v)

    def oracle():
        return S0 + args.lam * ws.dS_between(em0, ws.edeg_dict(v),
                                             estar=args.estar) - S0

    drv = fpkmc.HBDriver(s, apx, dS_max=args.dS_max, max_depth=args.depth,
                         rng=rng, audit_every=100, oracle=oracle)
    drv.S_rel = 0.0            # audits are per-batch; totals carry S0
    batch = min(args.hb_batch, args.hb_steps - done)
    for step in range(batch):
        drv.step()
        hb[(drv.chord, round(S0 + drv.S_rel, 6))] += 1
        if (done + step + 1) % 200 == 0:
            print(f"  HB {done + step + 1}/{args.hb_steps}: "
                  f"{len(hb)} states, acc {prev_acc + drv.naccept}",
                  flush=True)
    done += batch
    if done < args.hb_steps:
        # save exact state and exec-restart to shed the leak
        v.save(st_mfd)
        with open(st_json, "w") as fh:
            json.dump({"chord": list(drv.chord),
                       "S_rel": S0 + drv.S_rel,
                       "hb": {json.dumps([list(k[0]), k[1]]): c
                              for k, c in hb.items()},
                       "done": done, "window": window,
                       "n3_target": v.num_facets,
                       "naccept": prev_acc + drv.naccept,
                       "rng": rng.bit_generator.state}, fh)
        print(f"  [batch done: {done}/{args.hb_steps}; exec-restart]",
              flush=True)
        argv_clean, skip = [], False
        for a in sys.argv[1:]:
            if skip:
                skip = False
                continue
            if a == "--state":
                skip = True
                continue
            argv_clean.append(a)
        os.execv(sys.executable,
                 [sys.executable, os.path.abspath(__file__)]
                 + argv_clean + ["--state", state_prefix])
    print(f"HB done: {len(hb)} states, "
          f"{prev_acc + drv.naccept}/{args.hb_steps} accepted, "
          f"audits passed", flush=True)
    # persist the HB histogram NOW -- the brute phase must not be able to
    # take these 1500 steps down with it. State files stay on disk until
    # the final output is written.
    if args.out:
        with open(args.out + ".hb.json", "w") as fh:
            json.dump({"hb": {json.dumps([list(k[0]), k[1]]): c
                              for k, c in hb.items()},
                       "window": window,
                       "naccept": prev_acc + drv.naccept}, fh)
    del drv, s, m

    # ---------------- brute chain ----------------
    m2, s2, apx2, _ = build(args.ref, args.estar, args.lam, window=window)
    v2 = s2.manifold
    em02 = ws.edeg_dict(v2)
    S_rel = 0.0
    facets = np.asarray(v2.facets())
    n3 = len(facets)
    pairs6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    brute = Counter()
    cur = tuple(apx2)
    naccept = 0
    em = ws.edeg_dict(v2)
    for it in range(args.brute_props):
        f = facets[rng.integers(n3)]
        i, j = pairs6[rng.integers(6)]
        a, b = int(f[i]), int(f[j])
        slot = int(rng.integers(12))
        brute[(cur, round(S_rel, 6))] += 1        # dwell time in attempts
        # Python-side prefilter: only deg-3 edges can slide, so anything
        # else is a rejection of the emulated channel -- decided HERE, not
        # by letting slide_at2 throw. Each D throw allocates (exception +
        # trace) under the never-shrinking conservative GC; 4M of them was
        # the third jetsam kill. Same proposal distribution, same rng
        # stream (the old exception path also skipped the accept draw).
        if em.get((a, b) if a < b else (b, a)) != 3:
            continue
        try:
            dS, arr = s2.slide_at2(a, b, slot, commit=False)
        except RuntimeError:
            dS, arr = None, None                  # e.g. illegal slot
        if dS is None:
            continue
        if dS <= 0 or rng.random() < np.exp(-dS):
            s2.slide_at(a, b, slot, commit=True)
            S_rel += dS
            cur = tuple(sorted(arr))
            naccept += 1
            # accepts are rare (~1e-4 of proposals): refresh everything
            facets = np.asarray(v2.facets())
            n3 = len(facets)
            em = ws.edeg_dict(v2)
            if naccept % 100 == 0:
                true = args.lam * ws.dS_between(em02, em, estar=args.estar)
                assert abs(true - S_rel) < 1e-6, "brute audit failed"
        if (it + 1) % 1_000_000 == 0:
            print(f"  brute {it + 1}/{args.brute_props}: "
                  f"{naccept} accepts, {len(brute)} states", flush=True)
    print(f"brute done: {naccept} accepted slides, {len(brute)} states, "
          f"audits passed")

    # ---------------- compare ----------------
    def top_states(c, n=12):
        return [k for k, _ in c.most_common(n)]

    common = [k for k in top_states(hb, 20) if k in brute
              and brute[k] > 50 and hb[k] > 10]
    print(f"\n{len(common)} well-visited common states")
    print(f"{'state (chord, S)':>34} {'HB occ':>8} {'brute occ':>10}")
    rows = []
    for k in common:
        rows.append((k, hb[k] / sum(hb.values()),
                     brute[k] / sum(brute.values())))
        print(f"{str(k):>34} {rows[-1][1]:8.4f} {rows[-1][2]:10.4f}")

    # pairwise Boltzmann test: ln(occ_i/occ_j) vs -(S_i - S_j)
    print("\npairwise ln-ratio vs -dS (slope 1 = Boltzmann):")
    xs, ys_hb, ys_br = [], [], []
    for ii in range(len(common)):
        for jj in range(ii + 1, len(common)):
            Si, Sj = common[ii][1], common[jj][1]
            if abs(Si - Sj) < 1e-9:
                continue
            xs.append(-(Si - Sj))
            ys_hb.append(np.log(hb[common[ii]] / hb[common[jj]]))
            ys_br.append(np.log(brute[common[ii]] / brute[common[jj]]))
    xs = np.array(xs)
    if len(xs) >= 3:
        for lab, ys in (("HB", np.array(ys_hb)), ("brute", np.array(ys_br))):
            A = np.vstack([xs, np.ones_like(xs)]).T
            (a, b0), *_ = np.linalg.lstsq(A, ys, rcond=None)
            res = ys - (a * xs + b0)
            print(f"   {lab:6s} slope {a:+.3f}  intercept {b0:+.3f}  "
                  f"rms residual {res.std():.3f}  (n={len(xs)} pairs)")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"hb": {str(k): v for k, v in hb.items()},
                       "brute": {str(k): v for k, v in brute.items()},
                       "window": window}, fh)
        print(f"wrote {args.out}")
    for f in (st_json, st_mfd):
        if os.path.exists(f):
            os.remove(f)


if __name__ == "__main__":
    main()
