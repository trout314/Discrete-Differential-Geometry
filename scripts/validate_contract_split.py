#!/usr/bin/env python3
"""Detailed-balance validation of the D-side contract/split channel.

THE CRITERION is the labeled pair test (--pair-test): from a fixed state x
(labels compacted to 0..n-1) with ONLY the channel active, the empirical
single-transition ratio P(x->y)/P(y->x) must equal exp(S(x)-S(y)). The
construction makes transitions exactly identifiable: a split creates
fresh = n deterministically and the reverse contraction of (w, n) keeps
w < n, restoring x bit-for-bit, so sorted-facet-bytes keys count single
labeled transitions. Forensics recomputes every Hastings ingredient and
prints predicted absolute rates for the correct formula and candidate
mis-formulas, plus the labeled path multiplicity (a 1->4-type split is
realized by 4 (w, gamma, side) paths, each in bijection with a reverse
contract-edge path; per-path acceptance pairs exactly its partner, so each
sub-kernel satisfies detailed balance individually -- single-path rate
predictions are underestimates by the multiplicity, ratios are exact).

Historical traps this design evolved through (kept as warnings):
  * The A/B equilibrium comparison (default mode, retained as a
    DIAGNOSTIC ONLY) is NOT a correctness criterion: the baseline
    sampler's own 1<->4 bistellar pair carries an uncorrected O(1)
    proposal asymmetry (hence run_exact / importance_weight), so the
    baseline does not sample exp(-S), and adding a CORRECT f0 channel
    shifts the mixture equilibrium anyway.
  * A pair test keyed on a weak isomorphism invariant aggregates classes
    of distinct targets and fakes a detailed-balance violation of size
    ~ class multiplicity (measured forward rates exceeding the
    single-transition proposal bound are the tell).

THE SENSITIVE DIAGNOSTIC is --circulation: with BOTH f0-changing families
live (the channel and the bistellar 1<->4 pair, the latter independently
Hastings-corrected), stationarity forces
    net channel f0 flux  +  net bistellar f0 flux  =  df0/dt = 0,
and net flux through EITHER family alone must vanish -- a family is closed
under reversal, so a biased-but-stationary chain still shows zero net flux
through it. A nonzero circulation is therefore two families DISAGREEING
about the same distribution, which is exactly what a wrong Hastings factor
in one of them produces. It is far more sensitive than any equilibrium
comparison (an O(1) rate imbalance moves <f0> by well under one sigma but
shows up as many sigma of circulation), and --max-ring scans it over the
cycle-length cap, which BISECTS the suspect machinery: maxRing = 3 exercises
only triangle cycles (a facial one is literally the 1->4 move), while 4..8
pull in the bounded-DFS / FK-catalog cycle counting that supplies q_rev.

Usage:
    caffeinate -i python scripts/validate_contract_split.py \
        --pair-test [--trials 150000] [--pairs 4]     # the criterion
    caffeinate -i python scripts/validate_contract_split.py \
        --circulation [--max-ring 3,4,5,6] [--sweeps 20000]
    caffeinate -i python scripts/validate_contract_split.py \
        [--moves 4000000]                             # A/B diagnostic
"""
import argparse
import math
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))

import discrete_differential_geometry as ddg


def block_err(x, nblock=20):
    x = np.asarray(x, float)
    k = len(x) // nblock
    if k == 0:
        return float("nan")
    means = x[: k * nblock].reshape(nblock, k).mean(axis=1)
    return float(means.std(ddof=1) / math.sqrt(nblock))


def make_sphere():
    m = ddg.Manifold.standard_sphere(3)
    return m


def run_chain(cs_prob, cimp, moves, sample_every, burn_frac=0.3, seed=0):
    m = make_sphere()
    params = ddg.SamplerParams(
        num_facets_target=120, num_facets_coef=0.05,
        hinge_degree_target=4.8, num_hinges_coef=0.02,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, params)
    if cimp > 0:
        s.set_n6_potential(0.0, cimp, tilt=[0.0] * 5)
    if cs_prob > 0:
        s.set_contract_split(cs_prob, 6)

    f0s, f3s = [], []
    n_chunks = moves // sample_every
    for i in range(n_chunks):
        s.run(moves=sample_every)
        v = s.manifold
        f0s.append(int(v.f_vector[0]))
        f3s.append(int(v.num_facets))
    burn = int(len(f0s) * burn_frac)
    stats = s.contract_split_stats() if cs_prob > 0 else (0, 0, 0, 0, 0)
    return (np.array(f0s[burn:]), np.array(f3s[burn:]), stats, s)


def compare(tag, moves, sample_every, cimp):
    print(f"\n=== {tag} (cimp={cimp}) ===")
    f0A, f3A, _, sA = run_chain(0.0, cimp, moves, sample_every)
    f0B, f3B, st, sB = run_chain(0.3, cimp, moves, sample_every)
    ct, ca, sp, sa, nv = st
    print(f"  channel: contract {ca}/{ct}, split {sa}/{sp}, noValid {nv}")
    ok = True
    for name, a, b in (("f0", f0A, f0B), ("f3", f3A, f3B)):
        ea, eb = block_err(a), block_err(b)
        diff = a.mean() - b.mean()
        err = math.hypot(ea, eb)
        z = diff / err if err > 0 else float("inf")
        line = (f"  <{name}>: A {a.mean():.3f}+-{ea:.3f}  "
                f"B {b.mean():.3f}+-{eb:.3f}  z={z:+.2f}")
        print(line)
        if abs(z) > 4.0:
            ok = False
    # histogram overlap for f0
    lo = min(f0A.min(), f0B.min())
    hi = max(f0A.max(), f0B.max())
    bins = np.arange(lo, hi + 2)
    hA, _ = np.histogram(f0A, bins=bins, density=True)
    hB, _ = np.histogram(f0B, bins=bins, density=True)
    overlap = float(np.minimum(hA, hB).sum())
    print(f"  f0 histogram overlap: {overlap:.3f}"
          f"  (A range {f0A.min()}-{f0A.max()}, B {f0B.min()}-{f0B.max()})")
    if overlap < 0.9:
        ok = False
    if ca + sa == 0:
        print("  WARNING: channel never accepted -- test has no power")
        ok = False
    print(f"  {'PASS' if ok else 'FAIL'}")
    return ok


def _base_action(F, params_tuple):
    """Mirror of the sampler's base objective for a facet array."""
    import itertools
    target, fc, hdt, hc = params_tuple
    F = np.asarray(F)
    f3 = len(F)
    E = np.vstack([F[:, [i, j]]
                   for i, j in itertools.combinations(range(4), 2)])
    E.sort(axis=1)
    f1 = len(np.unique(E, axis=0))
    return fc * (f3 - target) ** 2 + hc * (f1 - 6.0 * f3 / hdt) ** 2


def _state_key(v):
    """Label-invariant state fingerprint: sorted per-vertex profiles of
    (tet degree, sorted incident edge degrees). Needed because the reverse
    split draws a different fresh label than the contract retired."""
    import hashlib
    import itertools
    from collections import Counter

    F = np.asarray(v.facets())
    E = np.vstack([F[:, [i, j]]
                   for i, j in itertools.combinations(range(4), 2)])
    E.sort(axis=1)
    eu, edeg = np.unique(E, axis=0, return_counts=True)
    prof = {}
    for (a, b), d in zip(eu, edeg):
        prof.setdefault(int(a), []).append(int(d))
        prof.setdefault(int(b), []).append(int(d))
    vdeg = Counter(F.flatten().tolist())
    sig = sorted((vdeg[x], tuple(sorted(p))) for x, p in prof.items())
    return hashlib.blake2b(repr((len(F), sig)).encode(),
                           digest_size=12).hexdigest()


def pair_test(trials, seed_moves=150):
    """Direct detailed-balance check of the channel ALONE (cs.prob = 1).

    From a fixed labeled state x, estimate P(x -> y) over single-move
    trials for the first accepted contract's target y; then estimate
    P(y -> x). Detailed balance demands
        P(x->y) / P(y->x) = exp(S(x) - S(y)).
    This isolates the channel from the bistellar baseline's own
    (uncorrected, O(1)) 1<->4 proposal asymmetry, which invalidates naive
    A/B equilibrium comparisons.
    """
    pt = (120, 0.05, 4.8, 0.02)

    def make_params():
        return ddg.SamplerParams(
            num_facets_target=pt[0], num_facets_coef=pt[1],
            hinge_degree_target=pt[2], num_hinges_coef=pt[3],
            hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
            hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)

    # churned start; compact labels so the fresh-vertex draw (fVector[0])
    # is deterministic and collision-free
    m = make_sphere()
    s = ddg.ManifoldSampler(m, make_params())
    s.run(moves=seed_moves)
    FX0 = np.asarray(s.manifold.dup().facets())
    lab0, inv0 = np.unique(FX0, return_inverse=True)
    X = inv0.reshape(FX0.shape).tolist()

    # build each start manifold ONCE; ManifoldSampler copies it internally,
    # so per-trial cost is the sampler copy + one move, not a rebuild
    _mcache = {}

    def one_move_from(F):
        k = id(F)
        if k not in _mcache:
            _mcache[k] = ddg.Manifold(3, F)
        sm = ddg.ManifoldSampler(_mcache[k], make_params())
        sm.set_contract_split(1.0, 6)
        sm.run(moves=1)
        return sm.manifold

    # probe phase: find the most-frequented SPLIT neighbour (f0-increasing).
    # Labeled-exact keys demand split-forward pairs: from x (labels 0..n-1)
    # a split creates fresh = n deterministically, and the reverse
    # contraction of (w, n) keeps w < n, restoring x BIT-FOR-BIT -- so
    # sorted-facet-bytes keys identify single labeled transitions with no
    # isomorphism-class aggregation (a weak invariant key aggregates all
    # gamma's with the same degree profile and fakes a detailed-balance
    # violation of size ~ class multiplicity).
    from collections import Counter

    def labeled_key(v):
        F = np.asarray(v.facets())
        F = np.sort(F, axis=1)
        F = F[np.lexsort(F.T[::-1])]
        return F.tobytes()

    mX = ddg.Manifold(3, X)
    kx = labeled_key(mX)
    f0x = int(mX.f_vector[0])
    seen = Counter()
    example = {}
    for _ in range(6000):
        v = one_move_from(X)
        k = labeled_key(v)
        if k != kx and int(v.f_vector[0]) == f0x + 1:
            seen[k] += 1
            if k not in example:
                example[k] = np.asarray(v.dup().facets()).tolist()
    assert seen, "channel never accepted a split from x"
    # among the well-trafficked targets, pick the smallest |dS| (both
    # directions then have measurable rates -- an argmax-count pick can
    # land on a steeply uphill pair with no forward statistics)
    cands = [k for k, c in seen.most_common(10) if c >= 2] or \
        [seen.most_common(1)[0][0]]
    SX = _base_action(X, pt)
    ky = min(cands, key=lambda k: abs(_base_action(example[k], pt) - SX))
    Y = example[ky]
    dS = _base_action(Y, pt) - _base_action(X, pt)
    print(f"  pair: f0 {len(np.unique(np.asarray(X)))} -> "
          f"{len(np.unique(np.asarray(Y)))}, dS = {dS:+.4f}")

    def count_hits(F, target_key, n):
        hits = 0
        for _ in range(n):
            if labeled_key(one_move_from(F)) == target_key:
                hits += 1
        return hits

    n_xy = count_hits(X, ky, trials)
    n_yx = count_hits(Y, kx, trials)
    print(f"  P(x->y) = {n_xy}/{trials}, P(y->x) = {n_yx}/{trials}")
    assert n_xy > 10 and n_yx > 10, "too few transitions -- raise --trials"
    ratio = (n_xy / trials) / (n_yx / trials)
    expect = math.exp(-dS)
    # binomial errors on the log-ratio
    err = math.sqrt(1.0 / n_xy + 1.0 / n_yx)
    z = (math.log(ratio) - math.log(expect)) / err
    print(f"  ratio {ratio:.4f} vs exp(-dS) {expect:.4f}  "
          f"(log-z = {z:+.2f})")
    ok = abs(z) < 4.0
    print(f"  {'PASS' if ok else 'FAIL'}")
    forensics(X, Y, n_xy, n_yx, trials, pt)
    return ok


def _count_cycles_le(adj_sets, L=6):
    """Simple-cycle count (length <= L) in an adjacency-set graph."""
    verts = sorted(adj_sets)
    count = 0

    def dfs(start, path, seen):
        nonlocal count
        u = path[-1]
        for w2 in adj_sets[u]:
            if w2 == start and len(path) >= 3:
                if path[1] < path[-1]:
                    count += 1
            elif w2 > start and w2 not in seen and len(path) < L:
                seen.add(w2)
                path.append(w2)
                dfs(start, path, seen)
                path.pop()
                seen.remove(w2)

    for s in verts:
        dfs(s, [s], {s})
    return count


def forensics(X, Y, n_xy, n_yx, trials, pt):
    """Recompute every Hastings ingredient for the pair and print predicted
    absolute transition probabilities for the correct formula and for
    candidate mis-formulas -- the best match to the measured rates
    fingerprints the bug."""
    import itertools
    sys.path.insert(0, os.path.join(_ROOT, "scripts"))
    from graft_signature import CrystalContext
    from contract_relax import contraction_delta, apply_contraction

    FX, FY = np.asarray(X), np.asarray(Y)
    f0X, f0Y = len(np.unique(FX)), len(np.unique(FY))
    # orient: L = lower-f0 state, H = higher; the split direction is L->H
    if f0X < f0Y:
        FL, FH, measured_LH, measured_HL = FX, FY, n_xy, n_yx
    else:
        FL, FH, measured_LH, measured_HL = FY, FX, n_yx, n_xy

    from validate_contract_split import _state_key  # self-import ok
    kL = _state_key(ddg.Manifold(3, FL.tolist()))

    # find the contractible edge of H whose contraction is L
    ctxH = CrystalContext(FH)
    found = None
    for (a, b) in ctxH.edge_deg:
        d = contraction_delta(ctxH, int(a), int(b))
        if d is None:
            continue
        newF = apply_contraction(FH, d, int(FH.max()) + 1)
        if _state_key(ddg.Manifold(3, newF.tolist())) == kL:
            found = (int(a), int(b), d, newF)
            break
    assert found is not None, "could not identify the pair's move"
    a, b, d, FLp = found
    rl = ctxH.edge_deg[(a, b)]
    f3H, f3L = len(FH), len(FLp)
    w = min(a, b)   # surviving label convention

    # merged-link quantities in the contracted state
    ctxL = CrystalContext(FLp)
    starW = ctxL.star_of_vertex(w)
    deg3 = len(starW)
    adj = {}
    for t in FLp[starW]:
        rest = [int(x) for x in t if x != w]
        for p, q in itertools.combinations(rest, 2):
            adj.setdefault(p, set()).add(q)
            adj.setdefault(q, set()).add(p)
    N = _count_cycles_le(adj, 6)

    SL = _base_action(FLp, pt)
    SH = _base_action(FH, pt)
    dS_LH = SH - SL

    q_split = (deg3 / (4.0 * f3L)) * (1.0 / N) * 0.5
    q_con = rl / (6.0 * f3H)

    def pred(qs, qc):
        aS = min(1.0, math.exp(-dS_LH) * qc / qs)
        aC = min(1.0, math.exp(+dS_LH) * qs / qc)
        return qs * aS, qc * aC

    # PROPOSAL MULTIPLICITY: a labeled transition can be realized by several
    # (w, gamma, side) paths -- e.g. a 1->4-type split by 4 (one per vertex
    # of the subdivided tet), each in bijection with a reverse contract-edge
    # path (contract (a, fresh) keeps a and restores L bit-for-bit). The
    # per-path acceptance pairs exactly the corresponding reverse path, so
    # each sub-kernel -- and hence the total -- satisfies detailed balance
    # EXACTLY; but the single-path theory rows below underpredict the
    # measured TOTAL rates by roughly this multiplicity:
    def sorted_bytes(F):
        F = np.sort(np.asarray(F), axis=1)
        return F[np.lexsort(F.T[::-1])].tobytes()

    labs_fresh = sorted(set(int(x) for x in np.unique(FH))
                        - set(int(x) for x in np.unique(FL)))
    mult = 0
    if len(labs_fresh) == 1:
        nf = labs_fresh[0]
        target = sorted_bytes(FL)
        for (p, q) in list(ctxH.edge_deg):
            if nf not in (p, q):
                continue
            keep = p if q == nf else q
            dd = contraction_delta(ctxH, int(p), int(q))
            if dd is None:
                continue
            if sorted_bytes(apply_contraction(FH, dd, keep)) == target:
                mult += 1
    print(f"  forensics: rl={rl} deg3(w)={deg3} N={N} "
          f"f3L={f3L} f3H={f3H} dS(L->H)={dS_LH:+.4f} "
          f"labeled path multiplicity={mult}")
    print(f"  measured:  P(L->H)={measured_LH/trials:.3e}  "
          f"P(H->L)={measured_HL/trials:.3e}")
    variants = {
        "correct":        (q_split, q_con),
        "no side coin":   (q_split * 2, q_con),
        "no N":           (deg3 / (4.0 * f3L) * 0.5, q_con),
        "N=len<=8":       None,   # filled below
        "swap 4/6":       ((deg3 / (6.0 * f3L)) * (1.0 / N) * 0.5,
                           rl / (4.0 * f3H)),
        "no deg3":        ((1.0 / (4.0 * f3L)) * (1.0 / N) * 0.5, q_con),
        "no rl":          (q_split, 1.0 / (6.0 * f3H)),
    }
    N8 = _count_cycles_le(adj, 8)
    variants["N=len<=8"] = ((deg3 / (4.0 * f3L)) * (1.0 / N8) * 0.5, q_con)
    for name, qq in variants.items():
        pLH, pHL = pred(*qq)
        print(f"    {name:14s} P(L->H)={pLH:.3e}  P(H->L)={pHL:.3e}")


def circulation(max_ring, cs_prob, sweeps, chunks, target, cimp,
                burn_frac=0.25, nh=0.1, seed_moves=20000):
    """Net f0 flux carried by each family in a stationary mixed chain.

    Returns a dict of per-sweep rates with block errors. The exact identity
    d f0 = (splits - contracts) + (n_{1->4} - n_{4->1}) is asserted per
    chunk -- those are the only moves that change f0, so any mismatch is a
    counter bug, not statistics.
    """
    params = ddg.SamplerParams(
        num_facets_target=target, num_facets_coef=0.05,
        hinge_degree_target=4.8, num_hinges_coef=nh,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    m = make_sphere()
    s = ddg.ManifoldSampler(m, params)
    if cimp > 0:
        s.set_n6_potential(0.0, cimp, tilt=[0.0] * 5)
    s.set_contract_split(cs_prob, max_ring)

    # grow + equilibrate (channel live throughout, so f0 starts stationary)
    s.run(moves=seed_moves)
    s.run(sweeps=sweeps * burn_frac)

    def snap():
        ct, ca, sp, sa, nv = s.contract_split_stats()
        st = s.get_stats()
        v = s.manifold
        return dict(ca=ca, sa=sa, ct=ct, sp=sp, nv=nv,
                    # index = len(coCenter) - 1: 1->4 puts a single new
                    # vertex in coCenter (index 0), 4->1 a whole facet (3)
                    b14=int(st.bistellar_accepts[0]),
                    b41=int(st.bistellar_accepts[3]),
                    f0=int(v.f_vector[0]), f3=int(v.num_facets),
                    tried=int(st.total_tried))

    prev = snap()
    chan, bist, f0s, f3s, msweeps = [], [], [], [], []
    per = sweeps * (1.0 - burn_frac) / chunks
    for _ in range(chunks):
        s.run(sweeps=per)
        cur = snap()
        dchan = (cur["sa"] - prev["sa"]) - (cur["ca"] - prev["ca"])
        dbist = (cur["b14"] - prev["b14"]) - (cur["b41"] - prev["b41"])
        df0 = cur["f0"] - prev["f0"]
        assert df0 == dchan + dbist, (
            f"f0 bookkeeping broken: df0={df0} chan={dchan} bist={dbist}")
        nsw = (cur["tried"] - prev["tried"]) / max(1.0, 0.5 *
                                                   (cur["f3"] + prev["f3"]))
        chan.append(dchan / nsw)
        bist.append(dbist / nsw)
        f0s.append(cur["f0"])
        f3s.append(cur["f3"])
        msweeps.append(nsw)
        prev = cur

    chan = np.array(chan, float)
    bist = np.array(bist, float)
    nb = min(20, len(chan))
    ec = block_err(chan, nb)
    eb = block_err(bist, nb)
    tot = prev["ca"] + prev["sa"]
    return dict(max_ring=max_ring, circ=chan.mean(), circ_err=ec,
                bist=bist.mean(), bist_err=eb,
                z=chan.mean() / ec if ec > 0 else float("nan"),
                f0=float(np.mean(f0s)), f0_sd=float(np.std(f0s)),
                f3=float(np.mean(f3s)), accepts=tot,
                acc_rate=tot / max(1, prev["ct"] + prev["sp"]),
                noValid=prev["nv"], sweeps=float(np.sum(msweeps)))


def fugacity(cs_prob, max_ring, sweeps, chunks, target, cimp,
             burn_frac=0.25, nh=0.1, hdt=4.8, seed_moves=20000,
             bist_hastings=True):
    """Spurious vertex fugacity of the channel, read off the equilibrium.

    The base action depends on the f-vector only through (f1, f3), and a
    3-manifold has f1 = f0 + f3, so
        mu := dS/df0 |_f3 = 2 nh (f1 - 6 f3 / hdt)
    is the chemical potential the chain equilibrates against. Two kernels
    that both target exp(-S) must settle at the same <mu>; a kernel carrying
    a spurious fugacity z^f0 settles where <mu> is displaced by exactly
    ln z. So ln z = <mu>(cs_prob) - <mu>(0) reads the missing Hastings
    factor directly, in nats -- this is the same meter that identified the
    bistellar 1->4 clip as ln 2.
    """
    params = ddg.SamplerParams(
        num_facets_target=target, num_facets_coef=0.05,
        hinge_degree_target=hdt, num_hinges_coef=nh,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(make_sphere(), params)
    s.set_bistellar_hastings(bist_hastings)
    if cimp > 0:
        s.set_n6_potential(0.0, cimp, tilt=[0.0] * 5)
    if cs_prob > 0:
        s.set_contract_split(cs_prob, max_ring)
    s.run(moves=seed_moves)
    s.run(sweeps=sweeps * burn_frac)

    mus, f0s, f3s = [], [], []
    per = sweeps * (1.0 - burn_frac) / chunks
    for _ in range(chunks):
        s.run(sweeps=per)
        v = s.manifold
        f0 = int(v.f_vector[0])
        f3 = int(v.num_facets)
        f1 = f0 + f3                      # Euler, closed 3-manifold
        mus.append(2.0 * nh * (f1 - 6.0 * f3 / hdt))
        f0s.append(f0)
        f3s.append(f3)
    mus = np.array(mus, float)
    return dict(cs=cs_prob, mu=mus.mean(), mu_err=block_err(mus, 20),
                f0=float(np.mean(f0s)), f3=float(np.mean(f3s)),
                bh=bist_hastings)


def run_fugacity(probs, max_ring, sweeps, chunks, target, cimp):
    print(f"\n=== fugacity scan (ring {max_ring}, target f3 {target}, "
          f"{sweeps} sweeps/setting) ===")
    print("  ln z = <mu> - <mu>(cs=0); a correct channel gives ln z = 0.")
    print(f"  {'cs_prob':>8} {'<f0>':>8} {'<f3>':>7} {'<mu> = dS/df0':>22} "
          f"{'ln z':>18} {'z':>8}")
    base = None
    rows = []
    for p in probs:
        d = fugacity(p, max_ring, sweeps, chunks, target, cimp)
        rows.append(d)
        if base is None:
            base = d
            lnz, lnz_err = 0.0, 0.0
        else:
            lnz = d["mu"] - base["mu"]
            lnz_err = math.hypot(d["mu_err"], base["mu_err"])
        print(f"  {d['cs']:>8.3f} {d['f0']:>8.2f} {d['f3']:>7.1f} "
              f"{d['mu']:>+14.4f} +- {d['mu_err']:.4f} "
              f"{lnz:>+11.4f} +- {lnz_err:.4f} {math.exp(lnz):>8.4f}")
    print("  reference constants: ln 2 = 0.6931, ln 3 = 1.0986, "
          "ln(3/2) = 0.4055, ln 4 = 1.3863")
    return rows


def run_circulation(rings, cs_prob, sweeps, chunks, target, cimp):
    print(f"\n=== circulation scan (target f3 {target}, cs_prob {cs_prob}, "
          f"cimp {cimp}, {sweeps} sweeps/setting) ===")
    print("  net f0 flux per sweep through each family; both must be 0 in a\n"
          "  correct stationary chain, and they must cancel.")
    print(f"  {'ring':>4} {'<f0>':>8} {'<f3>':>7} {'acc':>6} "
          f"{'channel circ/sweep':>24} {'bistellar/sweep':>20} {'z':>7}")
    rows = []
    ok = True
    for r in rings:
        d = circulation(r, cs_prob, sweeps, chunks, target, cimp)
        rows.append(d)
        print(f"  {d['max_ring']:>4} {d['f0']:>8.2f} {d['f3']:>7.1f} "
              f"{d['acc_rate']:>6.3f} "
              f"{d['circ']:>+14.4f} +- {d['circ_err']:.4f} "
              f"{d['bist']:>+12.4f} +- {d['bist_err']:.4f} "
              f"{d['z']:>+7.2f}")
        if abs(d["z"]) > 4.0:
            ok = False
    print(f"  {'PASS' if ok else 'FAIL'} "
          f"(|z| < 4 on every ring cap)")
    return ok, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--moves", type=int, default=4_000_000)
    ap.add_argument("--sample-every", type=int, default=400)
    ap.add_argument("--pair-test", action="store_true")
    ap.add_argument("--trials", type=int, default=60000)
    ap.add_argument("--pairs", type=int, default=3)
    ap.add_argument("--fugacity", action="store_true",
                    help="read the channel's spurious vertex fugacity in nats")
    ap.add_argument("--cs-scan", default="0.0,0.3,0.9,0.99",
                    help="channel probabilities for --fugacity")
    ap.add_argument("--circulation", action="store_true",
                    help="net f0 flux per family (the sensitive diagnostic)")
    ap.add_argument("--max-ring", default="3,4,5,6",
                    help="comma-separated ring caps to scan")
    ap.add_argument("--sweeps", type=int, default=20000)
    ap.add_argument("--chunks", type=int, default=200)
    ap.add_argument("--target", type=int, default=240)
    ap.add_argument("--cs-prob", type=float, default=0.3)
    ap.add_argument("--cimp", type=float, default=0.0)
    args = ap.parse_args()

    if args.fugacity:
        probs = [float(x) for x in args.cs_scan.split(",") if x.strip()]
        rings = [int(x) for x in args.max_ring.split(",") if x.strip()]
        for r in rings:
            run_fugacity(probs, r, args.sweeps, args.chunks,
                         args.target, args.cimp)
        sys.exit(0)

    if args.circulation:
        rings = [int(x) for x in args.max_ring.split(",") if x.strip()]
        ok, _ = run_circulation(rings, args.cs_prob, args.sweeps,
                                args.chunks, args.target, args.cimp)
        sys.exit(0 if ok else 1)

    if args.pair_test:
        ok = True
        for i in range(args.pairs):
            print(f"\n=== pair test {i} ===")
            ok &= pair_test(args.trials, seed_moves=150 + 37 * i)
        print(f"\nOVERALL: {'PASS' if ok else 'FAIL'}")
        sys.exit(0 if ok else 1)

    ok = compare("base action", args.moves, args.sample_every, cimp=0.0)
    ok &= compare("with n6 potential", args.moves, args.sample_every,
                  cimp=0.3)
    print(f"\nOVERALL: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
