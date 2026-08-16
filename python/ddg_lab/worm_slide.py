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

import discrete_differential_geometry as ddg
from . import worm_moves as wm
from discrete_differential_geometry.crystal_grains import REF_GLOB, best_refs

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


N_SLOTS = 12   # 2 chord orientations x 6 ordered (c2,c3) picks from the link


def slide_slots(m, chord):
    """The FIXED slot family at a chord: 2 orientations x 6 ordered link
    picks = N_SLOTS entries, each (cs, mv) or None where the frame or the
    move template fails to validate.

    The slot count is deliberately constant (not a count of valid slides):
    the MH proposal draws a slot uniformly from N_SLOTS and treats an
    invalid slot as a rejected proposal, so the denominator is the same in
    both states and cancels from the Hastings ratio."""
    link = edge_link_verts(m, *chord)
    if len(link) != 3:
        return [None] * N_SLOTS
    out = []
    for c0, c4 in (tuple(chord), tuple(chord)[::-1]):
        for c2, c3 in permutations(link, 2):
            cs = derive_frame(m, c0, c4, c2, c3)
            mv = slide_moves(m, cs) if cs is not None else None
            out.append((cs, mv) if mv is not None else None)
    return out


def candidate_slides(m, chord):
    """The valid slots at this chord, as [(cs, mv), ...]."""
    return [s for s in slide_slots(m, chord) if s is not None]


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


def clean_verdict(em_before, em_after, arrival_chord):
    """Is this slide species-preserving ("clean")?

    Returns (global_clean, local_clean). The GLOBAL test compares the whole
    illegal-degree multiset; the LOCAL test compares it only over edges
    whose degree changed. The two are equivalent -- unchanged edges
    contribute identically to both multisets -- and the local form is what
    the D implementation uses, since a slide's changed edges are exactly the
    template's touched edges (O(1) off the incremental degree map, no global
    scan). closure_test asserts the equivalence on real states.

    Both also require the arrival chord to be a degree-3 edge, i.e. the knot
    actually landed with a usable handle."""
    def ill(d):
        return sorted(x for x in d if x not in (5, 6))
    landed = em_after.get(arrival_chord) == 3
    gclean = ill(em_before.values()) == ill(em_after.values()) and landed
    changed = [k for k in set(em_before) | set(em_after)
               if em_before.get(k) != em_after.get(k)]
    lb = ill([em_before[k] for k in changed if k in em_before])
    la = ill([em_after[k] for k in changed if k in em_after])
    return gclean, (lb == la and landed)


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


def clean_slots_at(m, chord, em0=None):
    """Every CLEAN slide at `chord`, as [(cs, mv, key, em_after), ...]:
    frame, move template, end-state key (net facet diff), post-state edge
    degrees. Also returns the count of slots where the local and global
    cleanliness verdicts disagreed (must be 0)."""
    if em0 is None:
        em0 = edeg_dict(m)
    out = []
    ndisagree = 0
    for slot in slide_slots(m, chord):
        if slot is None:
            continue
        cs, mv = slot
        recs = apply_slide(m, mv)
        if recs is None:
            continue
        em1 = edeg_dict(m)
        arrival = tuple(sorted((cs["c4"], cs["c8"])))
        gclean, lclean = clean_verdict(em0, em1, arrival)
        if gclean != lclean:
            ndisagree += 1
        key = slide_net_diff(recs)
        undo_slide(m, recs)
        if gclean:
            out.append((cs, mv, key, em1))
    return out, ndisagree


def closure_test(args):
    """EXHAUSTIVE inverse-closure check for the clean-slide move class.

    For every clean slide of the state: apply it, then enumerate the slot
    family at the ARRIVAL chord and count the clean slots there whose net
    facet diff is the exact inverse. That count is k_r; the number of clean
    slots at the departure chord mapping to the same end state is k_f.

    Option B (clean-only move class) requires k_r >= 1 for every clean
    slide -- otherwise the class is not inverse-closed and the Hastings
    ratio has a zero denominator branch that is a BUG rather than a
    physical rejection. Also verifies local == global cleanliness.
    """
    if args.crystal_test or args.ref:
        path = args.ref or best_refs(REF_GLOB)["r"]
        m = ddg.Manifold.load(path, 3)
        F = np.asarray(m.facets())
        faces0, edeg0, vedges0 = wm.build_tables(F)
        site = None
        for face, d, e, valid in wm.two_three_sites(F, faces0, edeg0):
            if not valid:
                continue
            sf = sorted(face)
            degs = sorted(edeg0[frozenset(p)] for p in
                          [(sf[0], sf[1]), (sf[1], sf[2]), (sf[2], sf[0])])
            if degs == [5, 5, 6]:
                deltas = wm.delta_two_three(face, d, e, edeg0)
                if sorted(nw for _, nw in deltas.values()
                          if nw is not None and nw not in (5, 6)) == [3, 4, 4]:
                    site = (face, d, e)
                    break
        assert site is not None
        m.do_bistellar_move(sorted(site[0]), [site[1], site[2]])
        label = f"crystal {os.path.basename(path)} + one knot"
    else:
        m = ddg.Manifold.load(args.closure_test, 3)
        label = os.path.basename(args.closure_test)
    e3, ill = knot_chords(m)
    print(f"state: {label}  N3={m.num_facets}  "
          f"illegal {len(ill)}  deg-3 chords {len(e3)}")

    em0 = edeg_dict(m)
    n_clean = n_closed = 0
    kf_hist = defaultdict(int)
    kr_hist = defaultdict(int)
    disagree = 0
    failures = []
    for chord in e3:
        fwd, nd = clean_slots_at(m, chord, em0)
        disagree += nd
        bykey = defaultdict(list)
        for cs, mv, key, em1 in fwd:
            bykey[key].append((cs, mv, em1))
        for key, group in bykey.items():
            k_f = len(group)
            cs, mv, em1 = group[0]
            n_clean += 1
            recs = apply_slide(m, mv)
            assert recs is not None
            arrival = tuple(sorted((cs["c4"], cs["c8"])))
            inv = frozenset((t, -c) for t, c in key)
            rev, nd2 = clean_slots_at(m, arrival, em1)
            disagree += nd2
            k_r = sum(1 for _, _, rk, _ in rev if rk == inv)
            del em1
            undo_slide(m, recs)
            kf_hist[k_f] += 1
            kr_hist[k_r] += 1
            if k_r >= 1:
                n_closed += 1
            else:
                failures.append((chord, arrival, k_f))
                print(f"  NOT CLOSED: {chord} -> {arrival} (k_f={k_f}, "
                      f"{len(rev)} clean slots at arrival)")
        ddg.gc_collect()
    print(f"\nclean slides (distinct transitions): {n_clean}")
    print(f"  inverse-closed (k_r >= 1): {n_closed}/{n_clean}")
    print(f"  k_f histogram: {dict(sorted(kf_hist.items()))}")
    print(f"  k_r histogram: {dict(sorted(kr_hist.items()))}")
    print(f"  local-vs-global cleanliness disagreements: {disagree}")
    ok = (n_closed == n_clean and disagree == 0 and n_clean > 0)
    if ok and set(kf_hist) == {1} and set(kr_hist) == {1}:
        print("PASS: clean-slide class is inverse-closed, and k_f = k_r = 1 "
              "on every transition (Hastings ratio reduces to exp(-lam dS))")
    elif ok:
        print("PASS: clean-slide class is inverse-closed and the local "
              "cleanliness test is exact; multiplicities are NOT all 1, so "
              "the k_r/k_f factor must be kept")
    else:
        print("FAIL")
        sys.exit(1)


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


def dcross_test(args):
    """CROSSVAL: the D-side slide move against this oracle, slot for slot.

    For every degree-3 chord of the state and every one of the SLIDE_SLOTS
    slots, compare:

      * the CLEAN / not-clean verdict -- D's O(1) local changed-edge test
        against the oracle's global illegal-multiset comparison;
      * the action change -- D's speculative-delta accumulation through the
        real tracked objective against the oracle's independent recount (in
        lam = 1 units, so D's value must equal lam * the oracle's);
      * the resulting state, facet set for facet set, on a committed slide;
      * that the tracked objective moves by exactly that dS.

    This is the regression pinning the D move type to the reference
    implementation the inverse-closure proof was established on. The sampler
    state is restored after each committed slide by replaying the composite's
    inverse through the bookkeeping-safe targeted-move API.
    """
    snap = args.dcross
    lam, et = args.lam, args.etarget
    ddg.set_random_seed(args.seed)
    mo = ddg.Manifold.load(snap, 3)            # oracle playground
    m = ddg.Manifold.load(snap, 3)             # sampler manifold
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, hinge_degree_target=et,
        num_facets_coef=0.1, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * et / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(args.zleg * lam, args.cimp * lam, tilt=[0.0] * 5)
    v = s.manifold
    e3, _ = knot_chords(mo)
    print(f"state: {os.path.basename(snap)}  N3={mo.num_facets}  "
          f"deg-3 chords {len(e3)}  lam={lam}  e*={et:.6f}  "
          f"(zleg={args.zleg} cimp={args.cimp})")

    em0 = edeg_dict(mo)
    obj_base = s.current_objective
    nslot = nclean = nclean_seen = 0
    worst_dS = worst_obj = 0.0
    verdict_bad, state_bad, restore_bad = [], [], []

    for chord in e3:
        for j, slot in enumerate(slide_slots(mo, chord)):
            nslot += 1

            # --- oracle verdict for this slot
            py_ds, py_fac, py_recs = None, None, None
            if slot is not None:
                cs, mv = slot
                recs = apply_slide(mo, mv)
                if recs is not None:
                    em1 = edeg_dict(mo)
                    arrival = tuple(sorted((cs["c4"], cs["c8"])))
                    gclean, lclean = clean_verdict(em0, em1, arrival)
                    if gclean != lclean:
                        verdict_bad.append((chord, j, "oracle local != global"))
                    py_ds = dS_between(em0, em1, zleg=args.zleg,
                                       cimp=args.cimp, estar=et)
                    py_fac = facset(mo)
                    py_recs = recs
                    nclean_seen += 1 if gclean else 0
                    undo_slide(mo, recs)

            # --- D verdict for the same slot (trial: no state change)
            d_ds = s.slide_at(chord[0], chord[1], j, commit=False)
            if (py_ds is None) != (d_ds is None):
                verdict_bad.append((chord, j, f"valid: oracle="
                                    f"{py_ds is not None} D={d_ds is not None}"))
                continue
            if py_ds is None:
                continue
            nclean += 1
            err = abs(d_ds - lam * py_ds)
            worst_dS = max(worst_dS, err)
            assert err < 1e-9, (f"dS mismatch at chord {chord} slot {j}: "
                                f"D={d_ds!r} vs lam*oracle={lam * py_ds!r}")

            # --- commit through D and compare the end state
            obj0 = s.current_objective
            d_commit = s.slide_at(chord[0], chord[1], j, commit=True)
            assert d_commit is not None and abs(d_commit - d_ds) < 1e-12
            worst_obj = max(worst_obj,
                            abs((s.current_objective - obj0) - d_ds))
            if facset(v) != py_fac:
                state_bad.append((chord, j))

            # --- restore: replay the composite's inverse through the
            # bookkeeping-safe targeted-move API (reverse order, swap
            # center/coCenter).
            for _kind, cen, coc in reversed(py_recs):
                s.do_bistellar_move(coc, cen)
            if abs(s.current_objective - obj_base) > 1e-9:
                restore_bad.append((chord, j, s.current_objective - obj_base))
            ddg.gc_collect()

    print(f"\nslots examined: {nslot}    valid (verdicts agree): {nclean}"
          f"    of which clean: {nclean_seen}")
    print(f"  clean/not-clean verdict mismatches: {len(verdict_bad)}")
    for x in verdict_bad[:10]:
        print(f"    {x}")
    print(f"  end-state (facet set) mismatches:   {len(state_bad)}")
    for x in state_bad[:10]:
        print(f"    {x}")
    print(f"  restore-objective failures:         {len(restore_bad)}")
    for x in restore_bad[:10]:
        print(f"    {x}")
    print(f"  worst |dS_D - lam*dS_oracle|:       {worst_dS:.3e}")
    print(f"  worst |tracked objective step - dS|:{worst_obj:.3e}")
    ok = not (verdict_bad or state_bad or restore_bad) and nclean > 0
    print("PASS: the D slide move reproduces the oracle exactly "
          "(verdict, dS and end state)" if ok else "FAIL")
    if not ok:
        sys.exit(1)


def dsampler_test(args):
    """Integration test of the D-side slide MOVE TYPE inside a live sampler.

    Runs the ordinary sampler with slides enabled and verifies that after
    thousands of interleaved thermal moves and slides:

      * slides actually fire (tries/accepts > 0) at the expected rate;
      * the incremental degree/ridge maps still validate;
      * N_3 is unchanged by the slide channel (slides are f-vector neutral);
      * the cocycle is still attached and closed;
      * the tracked objective agrees with a from-scratch recompute.

    Also times the same run with slides off, to price the channel."""
    from discrete_differential_geometry import cocycle as coc
    import time
    snap = args.dsampler
    lam, et = args.lam, args.etarget

    def build(seed):
        ddg.set_random_seed(seed)
        m = ddg.Manifold.load(snap, 3)
        p = ddg.SamplerParams(
            num_facets_target=m.num_facets, hinge_degree_target=et,
            num_facets_coef=0.1, num_hinges_coef=0.0,
            hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
            hinge_degree_target_coef=lam * et / 6.0)
        s = ddg.ManifoldSampler(m, p)
        s.set_n6_potential(args.zleg * lam, args.cimp * lam, tilt=[0.0] * 5)
        return m, p, s

    m0, _, s0 = build(args.seed)
    n3_0 = m0.num_facets
    e3_0, _ = knot_chords(m0)
    print(f"state: {os.path.basename(snap)}  N3={n3_0}  "
          f"deg-3 chords {len(e3_0)}  lam={lam}  e*={et:.6f}")
    t0 = time.time()
    s0.run(sweeps=args.sweeps)
    t_off = time.time() - t0
    del m0, s0

    m, p, s = build(args.seed)
    cocpath = snap.replace(".mfd", ".cocycle.npz")
    have_coc = os.path.exists(cocpath)
    if have_coc:
        e0, w0, _ = coc.load_cocycle(cocpath)
        s.enable_cocycle(np.asarray(e0), np.asarray(w0))
        s.check_cocycle()
    v = s.manifold
    s.set_slide_prob(args.slide_prob)
    t0 = time.time()
    s.run(sweeps=args.sweeps)
    t_on = time.time() - t0
    tries, accepts = s.slide_stats()

    print(f"\n{args.sweeps} sweeps, slide_prob={args.slide_prob}")
    print(f"  slides: {tries} tries, {accepts} accepts "
          f"({100.0 * accepts / max(tries, 1):.1f}% acceptance)")
    print(f"  wall clock: {t_on:.1f}s with slides vs {t_off:.1f}s without "
          f"({100.0 * (t_on - t_off) / t_off:+.1f}%)")
    m.validate_maps()
    print(f"  incremental maps: OK")
    assert m.num_facets == n3_0, f"N3 changed: {n3_0} -> {m.num_facets}"
    print(f"  N3 preserved: {m.num_facets}")
    if have_coc:
        s.check_cocycle()
        e1, _ = s.read_cocycle()
        eset = {tuple(sorted(e)) for e in np.asarray(v.simplices(1)).tolist()}
        cset = {tuple(sorted(e)) for e in np.asarray(e1).tolist()}
        assert eset == cset, f"cocycle DETACHED ({len(eset ^ cset)} edges)"
        print(f"  cocycle: closed and attached")
    m2 = v.dup()
    s2 = ddg.ManifoldSampler(m2, p)
    s2.set_n6_potential(args.zleg * lam, args.cimp * lam, tilt=[0.0] * 5)
    drift = abs(s.current_objective - s2.current_objective)
    print(f"  tracked objective vs from-scratch: {drift:.3e}")
    assert drift < 1e-6, f"objective drift {drift}"
    assert tries > 0, "no slides fired -- the move type is not wired in"
    print("PASS: the D slide move type runs clean inside the sampler")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", nargs="?", default=None)
    ap.add_argument("--crystal-test", action="store_true")
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--thermal", default=None)
    ap.add_argument("--mh-test", default=None,
                    help="snapshot .mfd: run MH slide proposals against a "
                         "live sampler and verify bookkeeping")
    ap.add_argument("--closure-test", nargs="?", const="", default=None,
                    help="exhaustive inverse-closure + k_f/k_r check for the "
                         "clean-slide class; give a snapshot .mfd, or pass "
                         "--crystal-test to run on the reference crystal")
    ap.add_argument("--dcross", default=None,
                    help="snapshot .mfd: crossval the D-side slide move type "
                         "against this oracle, slot for slot (verdict, dS, "
                         "end state)")
    ap.add_argument("--dsampler", default=None,
                    help="snapshot .mfd: run the live sampler with the D-side "
                         "slide move type enabled and verify bookkeeping")
    ap.add_argument("--sweeps", type=int, default=200,
                    help="--dsampler: sweeps per run")
    ap.add_argument("--slide-prob", type=float, default=1.0,
                    help="--dsampler: probability of proposing a slide once "
                         "the unified proposal lands on a degree-3 edge")
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
    if args.dsampler:
        dsampler_test(args)
        return
    if args.dcross:
        dcross_test(args)
        return
    if args.closure_test is not None:
        closure_test(args)
        return
    if args.crystal_test:
        crystal_test(args)
    if args.thermal:
        thermal_test(args)
    if args.mh_test:
        mh_test(args)


if __name__ == "__main__":
    main()
