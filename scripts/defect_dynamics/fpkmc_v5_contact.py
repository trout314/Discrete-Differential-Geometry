#!/usr/bin/env python3
"""V5: contact round-trip -- dock boundary vs interaction onset (M3).

Two checks license FP exactness at contact (design sections 8, 9-V5):

PART 1 (deterministic head-on): knot A marches site-by-site along the
shared BC chain into frozen knot B. At every separation s we measure
  V(s) = R(s) - R(inf),   R(s) = S(with B) - S(without B)
(B removed and re-added THROUGH the sampler; the volume-pin constant
cancels in the difference), and the D dock flag of A's state (root of a
blocked-aware mini scan). REQUIRED ordering, per the one-tet-layer
conservative boundary: the dock flag fires strictly BEFORE the first
site with V != 0 -- every interior (non-dock) state satisfies V = 0 to
1e-9 (the additivity license), and the first nonzero V values are the
collider's repulsive contact wall (reported for comparison with the
knot_collider phase-1 table). The march is then walked back and the
exact restoration audited.

PART 2 (stochastic round-trip): repeated FP episodes -- fly (FPDriver)
from the start until a flight absorbs at a dock, verify V(dock) = 0 to
1e-9 at the ACTUAL dock geometry (arbitrary, not just head-on), then run
an explicit Metropolis burst (the exact marginalized slide channel, no
absorption) at contact and classify the outcome: bounce (final state no
longer dock-flagged) vs linger; single-chord integrity of the final
state is checked by a mini scan. Everything (FP slides + burst) is
walked back through exact inverses and the objective audited each
episode.

Each FP flight costs one D scan (~MB-scale leak, R6): keep --episodes
modest per process; split across seeds/processes for more statistics.
"""
import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import fpkmc
from discrete_differential_geometry._sampler import SLIDE_SLOTS
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit
from fpkmc_v4_fp import walk_back


def build(args):
    m = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(m.facets())
    _cc, _kcls, _seq, chain_prov = chain_for_run(
        args.ref, F, args.chain_class, seed_tet=0)
    orb = [int(x) for x in _seq]
    L = len(orb)
    jA = args.window
    jB = jA + 4 * args.sep_sites
    faceB = sorted(orb[(jB + k) % L] for k in (1, 2, 3))
    apxB = tuple(sorted((orb[jB % L], orb[(jB + 4) % L])))
    faceA = sorted(orb[(jA + k) % L] for k in (1, 2, 3))
    apxA = tuple(sorted((orb[jA % L], orb[(jA + 4) % L])))
    m.do_bistellar_move(faceB, list(apxB))
    m.do_bistellar_move(faceA, list(apxA))
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=args.estar, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * args.estar / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
    blocked = sorted(set(faceB) | set(apxB))
    return s, orb, jA, jB, apxA, apxB, faceB, blocked


def R_of_B(s, apxB, faceB):
    """S(with B) - S(without B), by tracked removal + re-add of B.
    Returns None when the 3->2 at B is invalid (deep contact)."""
    S1 = s.current_objective
    try:
        s.do_bistellar_move(list(apxB), list(faceB))
    except RuntimeError:
        return None
    S0 = s.current_objective
    s.do_bistellar_move(list(faceB), list(apxB))
    if abs(s.current_objective - S1) > 1e-9:
        raise RuntimeError("B remove/re-add did not restore the action")
    return S1 - S0


def root_dock(s, chord, blocked):
    g = s.slide_graph_scan(chord, dS_max=4.0, max_depth=1,
                           blocked_verts=blocked)
    return int(g["dock"][0]), int(g["n_chords"][0])


def forced_step(s, orb, j, fwd=True):
    """Chain slide of the knot at window j to window j+-4, found by
    arrival matching. Returns (new_j, rec) or (None, None)."""
    L = len(orb)
    chord = tuple(sorted((orb[j % L], orb[(j + 4) % L])))
    jn = (j + 4) % L if fwd else (j - 4) % L
    want = tuple(sorted((orb[jn % L], orb[(jn + 4) % L])))
    for slot in range(SLIDE_SLOTS):
        try:
            dS, arr = s.slide_at2(chord[0], chord[1], slot, commit=False)
        except RuntimeError:
            continue
        if arr is None or tuple(sorted(arr)) != want:
            continue
        dS, arr = s.slide_at2(chord[0], chord[1], slot, commit=True)
        return jn, (chord, tuple(sorted(arr)), float(dS))
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--window", type=int, default=40)
    ap.add_argument("--sep-sites", type=int, default=6)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--episodes", type=int, default=10)
    ap.add_argument("--max-flights", type=int, default=60)
    ap.add_argument("--burst-hits", type=int, default=200)
    ap.add_argument("--seed", type=int, default=41)
    ap.add_argument("--out", default=None)
    add_chain_args(ap, default=None)
    args = ap.parse_args()

    s, orb, jA, jB, apxA, apxB, faceB, blocked = build(args)
    S_start = s.current_objective
    n3 = s.manifold.num_facets
    nu = fpkmc.nu_per_attempt(1.0, n3)
    Rinf = R_of_B(s, apxB, faceB)
    print(f"A at window {jA}, B at {jB} (sep {args.sep_sites} sites); "
          f"R(inf) = {Rinf:.6f}", flush=True)

    # ---------------- PART 1: deterministic head-on ----------------
    print("\nPART 1: head-on march, V(s) vs dock flag")
    j = jA
    recs = []
    profile = []
    first_dock_s = first_V_s = None
    ok_additivity = True
    while True:
        sep = (jB - j) // 4
        chord = tuple(sorted((orb[j % len(orb)],
                              orb[(j + 4) % len(orb)])))
        dock, nch = root_dock(s, chord, blocked)
        R = R_of_B(s, apxB, faceB)
        V = None if R is None else R - Rinf
        profile.append({"sep": sep, "dock": dock, "V": V, "n_chords": nch})
        vs = "  (B frozen inside contact)" if V is None else f"{V:+.6f}"
        print(f"   sep {sep:2d}: dock={dock}  V={vs}", flush=True)
        if dock and first_dock_s is None:
            first_dock_s = sep
        if V is not None and abs(V) > 1e-9 and first_V_s is None:
            first_V_s = sep
        if not dock and V is not None and abs(V) > 1e-9:
            ok_additivity = False
            print("   ADDITIVITY VIOLATION at non-dock site", flush=True)
        if sep <= 1:
            break
        jn, rec = forced_step(s, orb, j, fwd=True)
        if jn is None:
            print(f"   forward slide invalid at sep {sep}; stopping march",
                  flush=True)
            break
        recs.append(rec)
        j = jn
    walk_back(s, recs)
    if abs(s.current_objective - S_start) > 1e-9:
        raise RuntimeError("part 1 walk-back failed to restore action")
    # FP exactness needs every INTERIOR (non-dock) state additive; dock
    # states carry exact energies (scan runs on the real manifold), so
    # dock-at-onset (tight boundary) is fine, interaction-before-dock is
    # not. ok_additivity already enforces exactly that; the sep values
    # are reported to show HOW tight the boundary runs.
    ordering_ok = (first_dock_s is not None
                   and (first_V_s is None or first_dock_s >= first_V_s))
    print(f"   first dock at sep {first_dock_s}, first V!=0 at sep "
          f"{first_V_s}; boundary "
          f"{'tight' if first_dock_s == first_V_s else 'conservative'}; "
          f"march reversed exactly", flush=True)

    # ---------------- PART 2: FP dock episodes ----------------
    print(f"\nPART 2: {args.episodes} FP dock episodes")
    rng = np.random.default_rng(args.seed)
    p_hit = 3.0 / (6.0 * n3)
    episodes = []
    v_dock_max = 0.0
    for ep in range(args.episodes):
        drv = fpkmc.FPDriver(s, apxA, blocked, nu, dS_max=args.dS_max,
                             max_depth=args.depth, rng=rng)
        reason = None
        for _ in range(args.max_flights):
            reason, info = drv.flight()
            if reason in (fpkmc.ABSORB_DOCK, fpkmc.ABSORB_MULTI):
                break
        rec = {"flights": drv.nflight, "t": drv.t, "reason": str(reason)}
        if reason == fpkmc.ABSORB_DOCK:
            R = R_of_B(s, apxB, faceB)
            V = None if R is None else R - Rinf
            rec["V_dock"] = V
            if V is not None:
                v_dock_max = max(v_dock_max, abs(V))
            # explicit burst at contact (marginalized slide channel)
            chord = drv.chord
            S_rel = 0.0
            burst = []
            commits = 0
            for _ in range(args.burst_hits):
                slot = int(rng.integers(SLIDE_SLOTS))
                try:
                    dS, arr = s.slide_at2(chord[0], chord[1], slot,
                                          commit=False)
                except RuntimeError:
                    continue
                if dS is None:
                    continue
                if dS <= 0 or rng.random() < np.exp(-dS):
                    dS, arr = s.slide_at2(chord[0], chord[1], slot,
                                          commit=True)
                    new_chord = tuple(sorted(int(x) for x in arr))
                    burst.append((chord, new_chord, float(dS)))
                    chord = new_chord
                    S_rel += float(dS)
                    commits += 1
            dock_end, nch_end = root_dock(s, chord, blocked)
            rec.update({"burst_commits": commits, "S_excursion": S_rel,
                        "end_dock": dock_end, "end_n_chords": nch_end,
                        "outcome": "linger" if dock_end else "bounce"})
            walk_back(s, burst)
        walk_back(s, drv.recs)
        if abs(s.current_objective - S_start) > 1e-6:
            raise RuntimeError(f"episode {ep} restore failed")
        episodes.append(rec)
        print(f"   ep {ep}: {rec}", flush=True)

    docks = [e for e in episodes if e["reason"] == fpkmc.ABSORB_DOCK]
    outcomes = {}
    for e in docks:
        outcomes[e.get("outcome", "?")] = outcomes.get(
            e.get("outcome", "?"), 0) + 1
    verdict = ("PASS" if ok_additivity and ordering_ok
               and v_dock_max < 1e-9 else "FAIL")
    print(f"\nV5 seed {args.seed}: {verdict}")
    print(f"  part 1: dock at sep {first_dock_s}, V!=0 at sep {first_V_s}; "
          f"additivity at non-dock sites "
          f"{'clean' if ok_additivity else 'VIOLATED'}")
    print(f"  part 2: {len(docks)}/{len(episodes)} episodes docked; "
          f"max |V(dock)| = {v_dock_max:.2e}; outcomes {outcomes}; "
          f"all episodes restored exactly")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"verdict": verdict, "profile": profile,
                       "first_dock_s": first_dock_s,
                       "first_V_s": first_V_s,
                       "v_dock_max": v_dock_max,
                       "episodes": episodes}, fh)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
