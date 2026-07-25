#!/usr/bin/env python3
"""The KNOT-SLIDE move: translate a (3,4,4) knot 4 positions along its local
BC (Boerdijk-Coxeter) chain in 4 bistellar moves -- the elementary transport
operator, packaged as a proposable sampler move.

Frame (from worm_helix motif, verified move-by-move 2026-07-25): the knot's
deg-3 edge is the CHORD (v_b, v_b+4); its link is {v_b+1, v_b+2, v_b+3}. The
slide template, with c_i = chain vertex v_b+i derived locally:

    M1  3-2 on chord (c0, c4), link {l1,l2,l3}
    M2  2-3 on face (c3, c4, c5), apexes (c2, c6)
    M3  2-3 on face (c6, c5, c7), apexes (c4, c8)
    M4  3-2 on chord (c2, c6)

after which the knot chord is (c4, c8): the knot moved +4 chain steps. The
chain continuation derives LOCALLY by the sliding-window rule (face_apexes
with exclusion) -- no global chain needed:

    c5 = apex of (c2,c3,c4) != c0,  c6 = apex of (c3,c4,c5) != c2,
    c7 = apex of (c4,c5,c6) != c3,  c8 = apex of (c5,c6,c7) != c4.

A direction candidate = choice of chord orientation (2) x ordered pair
(c2,c3) from the link (6) = up to 12 per knot; each validates (or not) by
existence + apex checks + has_bistellar_move. Invalid in a thermal
environment => proposal rejected (correct MH behavior). The inverse slide is
the mirrored template (alphabet inverse-closed).

Modes:
  --crystal-test         create a knot on the reference crystal, slide +N,
                         slide -N, close; assert EXACT facet restoration.
  --thermal SNAP.mfd     count, for every deg-3 edge of the state, the valid
                         slide directions + exact dS -- the mobilization
                         measurement.
"""
import argparse
import os
import sys
from collections import defaultdict
from itertools import permutations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import worm_moves as wm
from crystal_grains import REF_GLOB, best_refs

ESTAR = 5.105025


def edge_link_verts(m, a, b):
    try:
        lk = m.edge_link(a, b).tolist()
    except RuntimeError:
        return []
    return sorted({x for pr in lk for x in pr})


def apex_excluding(m, face, excl):
    """The apex of `face` that is not `excl`; None if face missing or excl
    is not one of its apexes (frame assumption broken)."""
    try:
        a1, a2 = m.face_apexes(*face)
    except RuntimeError:
        return None
    if a1 == excl:
        return a2
    if a2 == excl:
        return a1
    return None


def derive_frame(m, c0, c4, c2, c3):
    """Given chord (c0,c4) and ordered link pick (c2,c3), derive c5..c8 by
    the sliding-window rule. Returns dict or None."""
    c5 = apex_excluding(m, (c2, c3, c4), c0)
    if c5 is None:
        return None
    c6 = apex_excluding(m, (c3, c4, c5), c2)
    if c6 is None:
        return None
    c7 = apex_excluding(m, (c4, c5, c6), c3)
    if c7 is None:
        return None
    c8 = apex_excluding(m, (c5, c6, c7), c4)
    if c8 is None:
        return None
    cs = dict(c0=c0, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8)
    if len(set(cs.values())) != len(cs):
        return None                          # degenerate frame
    return cs


def slide_moves(m, cs):
    """The 4-move template for frame cs, with full validity checking.
    Returns list of ('23'|'32', center, cocenter) or None."""
    mv = []
    l1 = edge_link_verts(m, cs["c0"], cs["c4"])
    if len(l1) != 3 or not m.has_bistellar_move([cs["c0"], cs["c4"]], l1):
        return None
    mv.append(("32", [cs["c0"], cs["c4"]], l1))
    f2 = (cs["c3"], cs["c4"], cs["c5"])
    ap2 = (cs["c2"], cs["c6"])
    f3 = (cs["c6"], cs["c5"], cs["c7"])
    ap3 = (cs["c4"], cs["c8"])
    mv.append(("23", sorted(f2), list(ap2)))
    mv.append(("23", sorted(f3), list(ap3)))
    mv.append(("32", [cs["c2"], cs["c6"]], None))   # link resolved at apply
    return mv


def apply_slide(m, mv, dv=None):
    """Apply the 4 moves, checking validity stepwise; returns applied records
    for undo, or None (with automatic rollback) if any step invalid."""
    recs = []

    def undo_all():
        for kind, cen, coc in reversed(recs):
            m.do_bistellar_move(coc, cen)   # inverse: center/cocenter swap
        return None

    for i, (kind, cen, coc) in enumerate(mv):
        if kind == "32":
            if coc is None:
                coc = edge_link_verts(m, cen[0], cen[1])
                if len(coc) != 3:
                    return undo_all()
            if not m.has_bistellar_move(cen, coc):
                return undo_all()
            m.do_bistellar_move(cen, coc)
            recs.append(("32", cen, coc))
        else:
            try:
                a1, a2 = m.face_apexes(*cen)
            except RuntimeError:
                return undo_all()
            if {a1, a2} != set(coc):
                return undo_all()
            if not m.has_bistellar_move(cen, coc):
                return undo_all()
            m.do_bistellar_move(cen, coc)
            recs.append(("23", cen, coc))
    return recs


def undo_slide(m, recs):
    for kind, cen, coc in reversed(recs):
        m.do_bistellar_move(coc, cen)   # inverse: center/cocenter swap


def knot_chords(m):
    """All deg-3 edges with their complex signature sizes."""
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    e3 = [e for e, d in ill.items() if d == 3]
    return e3, ill


def candidate_slides(m, chord):
    """All valid slide frames from this chord (both orientations x 6 link
    orderings), deduplicated by resulting move list."""
    out = []
    link = edge_link_verts(m, *chord)
    if len(link) != 3:
        return out
    for c0, c4 in (chord, chord[::-1]):
        for c2, c3 in permutations(link, 2):
            cs = derive_frame(m, c0, c4, c2, c3)
            if cs is None:
                continue
            mv = slide_moves(m, cs)
            if mv is not None:
                out.append((cs, mv))
    return out


def facset(m):
    A = np.sort(np.asarray(m.facets()), axis=1)
    return frozenset(map(tuple, A.tolist()))


def dS_between(em_before, em_after, zleg=0.6, cimp=1.0, estar=ESTAR):
    """FULL shape-action difference (lam=1 units): edge-target term PLUS the
    n6 potential (zleg*U(n6) + cimp*m^2) — the terms that feel the halo.
    `estar` is both the edge-degree target and the per-edge coupling scale
    (coef estar/6); pass the run's actual target when matching a sampler."""
    from collections import defaultdict
    x = estar - int(estar)
    dS = 0.0
    changed = set()
    for k in set(em_before) | set(em_after):
        d0, d1 = em_before.get(k), em_after.get(k)
        if d0 == d1:
            continue
        changed |= set(k)
        if d0 is not None:
            dS -= (estar / 6.0) * ((d0 - estar) ** 2 - x * (1 - x))
        if d1 is not None:
            dS += (estar / 6.0) * ((d1 - estar) ** 2 - x * (1 - x))

    def counters(em):
        n6 = defaultdict(int)
        mm = defaultdict(int)
        for (a, b), d in em.items():
            if a not in changed and b not in changed:
                continue
            if d >= 6:
                n6[a] += 1
                n6[b] += 1
            if d not in (5, 6):
                mm[a] += 1
                mm[b] += 1
        return n6, mm

    n60, mm0 = counters(em_before)
    n61, mm1 = counters(em_after)
    for v in changed:
        dS += zleg * (wm.U_zleg(n61.get(v, 0)) - wm.U_zleg(n60.get(v, 0)))
        dS += cimp * (mm1.get(v, 0) ** 2 - mm0.get(v, 0) ** 2)
    return dS


def edeg_dict(m):
    eu, cnt = wm.edge_degrees_np(np.asarray(m.facets()))
    return {(int(a), int(b)): int(c) for (a, b), c in zip(eu, cnt)}


# ---------------------------------------------------------------------------
# MH proposal integration
# ---------------------------------------------------------------------------

def slide_net_diff(recs):
    """Net facet multiset change of a move-record list, as a canonical
    hashable key. From a fixed start state, equal keys <=> equal end state,
    so this identifies proposal multiplicities without facet-set hashing."""
    from collections import Counter
    net = Counter()
    for kind, cen, coc in recs:
        if kind == "32":
            a, b = cen
            x, y, z = coc
            for p, q in ((x, y), (y, z), (x, z)):
                net[tuple(sorted((a, b, p, q)))] -= 1
            net[tuple(sorted((a, x, y, z)))] += 1
            net[tuple(sorted((b, x, y, z)))] += 1
        else:
            f1, f2, f3 = cen
            d, e = coc
            net[tuple(sorted((d, f1, f2, f3)))] -= 1
            net[tuple(sorted((e, f1, f2, f3)))] -= 1
            for p, q in ((f1, f2), (f2, f3), (f1, f3)):
                net[tuple(sorted((d, e, p, q)))] += 1
    return frozenset((t, c) for t, c in net.items() if c != 0)


def enumerate_slides(m):
    """All slide candidates of the state, over every deg-3 chord:
    [(chord, cs, mv), ...]. Template-valid only; apply-time validity of the
    later moves is checked by apply_slide."""
    e3, _ = knot_chords(m)
    out = []
    for ch in e3:
        out += [(ch, cs, mv) for cs, mv in candidate_slides(m, ch)]
    return out


def mh_slide_step(sampler, lam, rng, estar=ESTAR, zleg=0.6, cimp=1.0):
    """One Metropolis-Hastings knot-slide proposal through a sampler.

    REFERENCE IMPLEMENTATION / ORACLE. This is the readable, checkable
    version of the slide as a sampler move; it is ~4s and ~1.8GB of D-side
    garbage per proposal (144 apply/undo ctypes round trips plus a manifold
    copy), i.e. ~6 orders of magnitude slower per move than the D-side
    thermal channel, so it is NOT a production sampling path. Production
    slides belong in sampler.d as a genuine move type; this stays as the
    crossval oracle for that implementation.

    Proposal: uniform over all template-valid slide candidates of the
    current state (2 orientations x 6 link orderings per deg-3 chord); a
    candidate that fails apply-time validation is a rejected proposal.
    Hastings correction: q(x'|x) = k_f/n_f and q(x|x') = k_r/n_r, where n =
    candidate count and k = multiplicity of the transition among candidates
    (by exact end-state key), so

        alpha = min(1, exp(-lam*dS) * (k_r * n_f) / (k_f * n_r)).

    dS is the FULL shape-action difference in lam=1 units (dS_between:
    edge-target with `estar` as target and coupling scale, zleg*U(n6),
    cimp*m^2). The composite kernel (thermal sweeps + slide proposals)
    preserves the production ensemble exactly when estar/zleg/cimp match the
    sampler couplings, the variance (VDV) terms are off, and the volume /
    global terms cancel (a slide preserves N3 and E exactly).

    Enumeration and trial applications run on a scratch copy; only an
    accepted slide replays through sampler.do_bistellar_move, which keeps
    cocycle tracking and the tracked objective coherent. Frozen-vertex
    constraints are NOT checked here.

    Returns a diagnostics dict: status (accepted / rejected / invalid /
    no-candidates), dS, alpha, n_f, n_r, k_f, k_r, chord, chord_after.
    """
    W = sampler.manifold.dup()
    fwd = enumerate_slides(W)
    n_f = len(fwd)
    if n_f == 0:
        return {"status": "no-candidates", "n_f": 0}
    keys = []
    for _, _, mv in fwd:
        recs = apply_slide(W, mv)
        if recs is None:
            keys.append(None)
            continue
        keys.append(slide_net_diff(recs))
        undo_slide(W, recs)
    i = int(rng.integers(n_f))
    ch, cs, mv = fwd[i]
    if keys[i] is None:
        return {"status": "invalid", "n_f": n_f}
    k_f = sum(1 for k in keys if k == keys[i])
    em0 = edeg_dict(W)
    recs = apply_slide(W, mv)
    em1 = edeg_dict(W)
    dS = dS_between(em0, em1, zleg=zleg, cimp=cimp, estar=estar)
    inv_key = frozenset((t, -c) for t, c in keys[i])
    rev = enumerate_slides(W)
    n_r = len(rev)
    k_r = 0
    for _, _, mv_r in rev:
        rr = apply_slide(W, mv_r)
        if rr is None:
            continue
        if slide_net_diff(rr) == inv_key:
            k_r += 1
        undo_slide(W, rr)
    out = {"status": "rejected", "dS": dS, "n_f": n_f, "n_r": n_r,
           "k_f": k_f, "k_r": k_r,
           "chord": (int(cs["c0"]), int(cs["c4"])),
           "chord_after": (int(cs["c4"]), int(cs["c8"]))}
    if k_r == 0:
        # the mirrored template should always be a reverse candidate; a hit
        # here means the frame derivation lost it -- reject (q(x|x')=0) and
        # flag for investigation.
        out["alpha"] = 0.0
        out["warn"] = "no reverse candidate"
        return out
    alpha = float(np.exp(-lam * dS)) * (k_r * n_f) / (k_f * n_r)
    out["alpha"] = min(1.0, alpha)
    if rng.random() < alpha:
        for _, cen, coc in recs:
            sampler.do_bistellar_move(cen, coc)
        out["status"] = "accepted"
    return out


def crystal_test(args):
    path = args.ref or best_refs(REF_GLOB)["r"]
    m = ddg.Manifold.load(path, 3)
    fs_crystal = facset(m)
    F = np.asarray(m.facets())
    faces0, edeg0, vedges0 = wm.build_tables(F)
    # create a knot on a (5,5,6) face (clean (3,4,4) creation)
    site = None
    for face, d, e, valid in wm.two_three_sites(F, faces0, edeg0):
        if not valid:
            continue
        degs = sorted(edeg0[frozenset(p)]
                      for p in [(sorted(face)[0], sorted(face)[1]),
                                (sorted(face)[1], sorted(face)[2]),
                                (sorted(face)[2], sorted(face)[0])])
        if degs == [5, 5, 6]:
            deltas = wm.delta_two_three(face, d, e, edeg0)
            if sorted(nw for _, nw in deltas.values()
                      if nw is not None and nw not in (5, 6)) == [3, 4, 4]:
                site = (face, d, e)
                break
    assert site is not None
    face, d, e = site
    m.do_bistellar_move(sorted(face), [d, e])
    e3, ill = knot_chords(m)
    print(f"created knot: chord {e3[0]}, illegal {sorted(ill.values())}")

    chord = e3[0]
    trail = []
    clean_frac = []
    for step in range(args.steps):
        cands = candidate_slides(m, chord)
        done = None
        nclean = 0
        for cs, mv in cands:
            recs = apply_slide(m, mv)
            if recs is None:
                continue
            want = tuple(sorted((cs["c4"], cs["c8"])))
            e3n, illn = knot_chords(m)
            clean = (sorted(illn.values()) == [3, 4, 4]
                     and tuple(sorted(e3n[0])) == want)
            if clean:
                done = (recs, want)
                nclean += 1
                break                        # keep it applied; stop scanning
            undo_slide(m, recs)
        if done is None:
            print(f"stuck at step {step} ({len(cands)} template-valid, "
                  f"0 clean)")
            return
        # re-apply the chosen one if we undid it while counting
        clean_frac.append((nclean, len(cands)))
        recs, chord = done
        trail.append(recs)
    cf = np.array([c for c, _ in clean_frac])
    tv = np.array([t for _, t in clean_frac])
    print(f"slid +{args.steps} steps (chord now {chord}); species (3,4,4) "
          f"preserved every step")
    print(f"  per step: template-valid {tv.mean():.1f}, "
          f"clean {cf.mean():.1f} directions")
    print("sliding back...")
    for recs in reversed(trail):
        undo_slide(m, recs)
    e3, ill = knot_chords(m)
    print(f"returned: chord {e3[0]}")
    # close and verify exact restoration
    link = edge_link_verts(m, *e3[0])
    m.do_bistellar_move(list(e3[0]), link)
    pairs, _ = m.illegal_edges()
    assert len(pairs) == 0
    assert facset(m) == fs_crystal, "crystal NOT exactly restored"
    print("closed: crystal restored EXACTLY (facet-set equality). "
          f"Slide is invertible transport over {args.steps} steps.")


def complex_census(m):
    """Complexes of illegal vertices with their illegal-edge signatures."""
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    illv = sorted({v for e in ill for v in e})
    adj = {v: set() for v in illv}
    for (a, b) in ill:
        adj[a].add(b)
        adj[b].add(a)
    seen = set()
    comps = []
    for v0 in illv:
        if v0 in seen:
            continue
        st, comp = [v0], set()
        seen.add(v0)
        while st:
            u = st.pop()
            comp.add(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
        sig = tuple(sorted(d for e, d in ill.items() if set(e) & comp))
        comps.append((sorted(comp), sig))
    return comps, ill


def thermal_test(args):
    m = ddg.Manifold.load(args.thermal, 3)
    em = edeg_dict(m)
    comps, ill = complex_census(m)
    ill_multiset = sorted(ill.values())
    e3 = [e for e, d in ill.items() if d == 3]
    from collections import Counter
    print(f"state: {os.path.basename(args.thermal)}  "
          f"illegal edges {len(ill)}  complexes {len(comps)}")
    print("  complex signatures: "
          + ", ".join(f"{sig}x{n}" for sig, n in
                      Counter(s for _, s in comps).most_common()))
    print(f"  deg-3 chords (slide handles): {len(e3)}")
    n_mobile = 0
    clean_dS = []
    per_clean = []
    per_valid = []
    for chord in e3:
        cands = candidate_slides(m, chord)
        nclean = nvalid = 0
        for cs, mv in cands:
            recs = apply_slide(m, mv)
            if recs is None:
                continue
            nvalid += 1
            em2 = edeg_dict(m)
            want = tuple(sorted((cs["c4"], cs["c8"])))
            clean = (sorted(d for d in em2.values() if d not in (5, 6))
                     == ill_multiset and em2.get(want) == 3)
            if clean:
                nclean += 1
                clean_dS.append(dS_between(em, em2))
            undo_slide(m, recs)
        per_valid.append(nvalid)
        per_clean.append(nclean)
        if nclean:
            n_mobile += 1
    print(f"chords with >=1 CLEAN slide (species-preserving): "
          f"{n_mobile}/{len(e3)}")
    print(f"  per chord: template-valid {per_valid}, clean {per_clean}")
    if clean_dS:
        d = np.array(clean_dS)
        print(f"CLEAN-slide FULL dS_shape (lam=1 units): mean {d.mean():+.3f} "
              f"med {np.median(d):+.3f}  min {d.min():+.3f}  max {d.max():+.3f}"
              f"  (n={len(d)})")
        acc = np.exp(-0.4 * d)
        print(f"  acceptance exp(-0.4 dS): med {np.median(acc):.3f}  "
              f"mean {acc.mean():.3f}")


def mh_test(args):
    """Integration test of mh_slide_step against a live sampler: cocycle
    stays attached across accepted slides, the tracked objective moves by
    exactly lam*dS, and the Hastings counts are sane."""
    from discrete_differential_geometry import cocycle as coc
    snap = args.mh_test
    lam = args.lam
    et = args.etarget
    ddg.set_random_seed(args.seed)
    m = ddg.Manifold.load(snap, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, hinge_degree_target=et,
        num_facets_coef=0.1, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * et / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(args.zleg * lam, args.cimp * lam, tilt=[0.0] * 5)
    cocpath = snap.replace(".mfd", ".cocycle.npz")
    if os.path.exists(cocpath):
        e0, w0, _ = coc.load_cocycle(cocpath)
        s.enable_cocycle(np.asarray(e0), np.asarray(w0))
        s.check_cocycle()
        print(f"state: {os.path.basename(snap)}  N3={m.num_facets}  "
              f"cocycle attached")
    else:
        cocpath = None
        print(f"state: {os.path.basename(snap)}  N3={m.num_facets}  "
              f"(no cocycle file; cocycle checks skipped)")
    v = s.manifold
    rng = np.random.default_rng(args.seed)
    stats = defaultdict(int)
    nobj = 0
    for k in range(args.props):
        obj0 = s.current_objective
        r = mh_slide_step(s, lam, rng, estar=et, zleg=args.zleg,
                          cimp=args.cimp)
        stats[r["status"]] += 1
        if r.get("warn"):
            stats["warn"] += 1
            print(f"  prop {k:2d}: WARN {r['warn']} chord {r.get('chord')}")
        if r["status"] in ("accepted", "rejected"):
            assert r["k_f"] >= 1, f"k_f=0 on a valid proposal: {r}"
        if r["status"] != "accepted":
            continue
        err = abs((s.current_objective - obj0) - lam * r["dS"])
        assert err < 1e-9, f"objective moved by != lam*dS (err {err})"
        nobj += 1
        if cocpath:
            s.check_cocycle()
            e1, _ = s.read_cocycle()
            eset = {tuple(sorted(e))
                    for e in np.asarray(v.simplices(1)).tolist()}
            cset = {tuple(sorted(e)) for e in np.asarray(e1).tolist()}
            assert eset == cset, f"cocycle DETACHED ({len(eset ^ cset)} edges)"
        print(f"  prop {k:2d}: ACCEPT chord {r['chord']} -> "
              f"{r['chord_after']}  dS={r['dS']:+.3f}  alpha={r['alpha']:.3f}"
              f"  (kf={r['k_f']} kr={r['k_r']} nf={r['n_f']} nr={r['n_r']})"
              f"  obj-check {err:.1e}"
              + ("  cocycle OK" if cocpath else ""))
        ddg.gc_collect()
    print(f"\n{args.props} proposals: {dict(stats)}")
    print(f"objective-vs-lam*dS verified on {nobj} accepts")
    # composite kernel: thermal sweeps after slides must leave everything
    # consistent (this is the production interleaving)
    s.run(sweeps=2)
    if cocpath:
        s.check_cocycle()
    m2 = v.dup()
    s2 = ddg.ManifoldSampler(m2, p)
    s2.set_n6_potential(args.zleg * lam, args.cimp * lam, tilt=[0.0] * 5)
    drift = abs(s.current_objective - s2.current_objective)
    assert drift < 1e-6, f"objective drift after composite kernel: {drift}"
    print(f"composite kernel (slides + 2 sweeps): objective drift {drift:.2e}"
          + (", cocycle OK" if cocpath else ""))
    print("PASS")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--crystal-test", action="store_true")
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--thermal", default=None)
    ap.add_argument("--mh-test", default=None,
                    help="snapshot .mfd: run MH slide proposals against a "
                         "live sampler and verify bookkeeping")
    ap.add_argument("--props", type=int, default=8,
                    help="--mh-test: number of proposals (each is ~4s / "
                         "~1.8GB of D-side garbage in this reference "
                         "implementation; keep it small)")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--etarget", type=float, default=ESTAR)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=4242)
    args = ap.parse_args()
    if args.crystal_test:
        crystal_test(args)
    if args.thermal:
        thermal_test(args)
    if args.mh_test:
        mh_test(args)


if __name__ == "__main__":
    main()
