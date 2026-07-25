#!/usr/bin/env python3
"""Migrating (3,4,4)+halo species in a quiet EDQ-only R-crystal host.

Action is the volume pin plus the EDQ (flat-edge-degree) term ONLY -- no n6
potential. Under that action a CLEAN knot slide has dS identically 0:

  * a slide is f-vector neutral, so E and N3 are preserved;
  * cleanliness preserves the whole illegal-degree multiset, hence n_ill and
    the illegal degree sum;
  * legal degrees are only 5 and 6, so n5 + n6 and 5*n5 + 6*n6 are both
    preserved => n5 and n6 individually are => the FULL degree multiset is.

The EDQ term is a function of that multiset alone, and the volume pin of N3
alone, so both cancel exactly. Every clean slide is therefore accepted with
probability 1: the knot free-diffuses on the constraint surface and the ONLY
control on its speed is `--slide-prob`. Nothing energetic binds a knot to its
halo here, so how far a knot outruns its halo is purely kinetic -- the point
this script measures.

Starting from a pristine crystal at stiff lam_EDQ nothing nucleates (n_ill
stays 0), so knots are IMPLANTED: a 2-3 move on a (5,5,6) face whose delta
signature is exactly [3,4,4]. Implanted knots are excitations above the
crystal ground state, so at stiff lam_EDQ a thermal move can annihilate one;
measuring that lifetime against the slide rate is the other half of the job.

Per TS-sweep block, records to <out>.ts.jsonl:
  n_illegal, complex sizes/members, per-complex tree-lift centroids, the
  knot census (how many complexes still carry a [3,4,4] signature), the
  cumulative slide tries/accepts, and <edeg>.
"""
import argparse
import json
import os
import sys
import time
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts", "../../tools"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes
import worm_moves as wm

ESTAR_R = 5.104225352112676     # native mean edge degree of the R crystal


def implant_knots(m, n_want, min_sep, rng, verbose=True):
    """Implant up to `n_want` (3,4,4) knots by 2-3 moves on (5,5,6) faces
    whose delta signature is exactly [3,4,4]. Sites are kept `min_sep` apart
    in vertex-label-free terms (no shared vertices within a graph ball), so
    the knots start as isolated monomers rather than a fused multimer."""
    placed = []
    used = set()
    if n_want <= 0:
        return placed          # naturally-nucleated defects only
    F = np.asarray(m.facets())
    faces0, edeg0, vedges0 = wm.build_tables(F)
    sites = list(wm.two_three_sites(F, faces0, edeg0))
    rng.shuffle(sites)
    for face, d, e, valid in sites:
        if len(placed) >= n_want:
            break
        if not valid:
            continue
        sf = sorted(face)
        degs = sorted(edeg0[frozenset(p)] for p in
                      [(sf[0], sf[1]), (sf[1], sf[2]), (sf[2], sf[0])])
        if degs != [5, 5, 6]:
            continue
        deltas = wm.delta_two_three(face, d, e, edeg0)
        sig = sorted(nw for _, nw in deltas.values()
                     if nw is not None and nw not in (5, 6))
        if sig != [3, 4, 4]:
            continue
        support = set(sf) | {d, e}
        if support & used:
            continue
        m.do_bistellar_move(sf, sorted((d, e)))
        placed.append(tuple(sf))
        # crude separation: forbid reuse of any vertex within min_sep hops
        ball = set(support)
        for _ in range(min_sep):
            ball = set(int(v) for f in np.asarray(m.facets())
                       if ball & set(int(x) for x in f) for v in f)
        used |= ball
        # tables are stale after the move; rebuild for the next site
        F = np.asarray(m.facets())
        faces0, edeg0, vedges0 = wm.build_tables(F)
    if verbose:
        print(f"implanted {len(placed)}/{n_want} knots", flush=True)
    return placed


def complexes(fac):
    """Connected components of illegal-degree vertices, with their edge-degree
    signature -- the object we call a knot/species."""
    n6, imp, adj = vertex_classes(fac)
    lab = np.unique(fac)
    V = len(lab)
    illv = [i for i in range(V) if imp[i] > 0]
    seen, comps = set(), []
    for s0 in illv:
        if s0 in seen:
            continue
        st, comp = [s0], []
        seen.add(s0)
        while st:
            u = st.pop()
            comp.append(u)
            for w in adj[u]:
                if imp[w] > 0 and w not in seen:
                    seen.add(w)
                    st.append(w)
        comps.append(sorted(int(lab[x]) for x in comp))
    return comps, imp


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", default=os.path.join(
        _HERE, "..", "..", "data", "tcp_reference", "T3_R_m3_N24462.mfd"))
    ap.add_argument("--mcell", type=int, default=3)
    ap.add_argument("--lam", type=float, required=True,
                    help="lam_EDQ (the ONLY shape coupling; n6 is off)")
    ap.add_argument("--etarget", type=float, default=None,
                    help="edge-degree target (default: the cell's native mean)")
    ap.add_argument("--slide-prob", type=float, default=0.0)
    ap.add_argument("--zleg", type=float, default=0.0,
                    help="n6 leg-potential coupling (0 = EDQ-only action)")
    ap.add_argument("--cimp", type=float, default=0.0,
                    help="n6 impurity coupling (0 = EDQ-only action)")
    ap.add_argument("--start", default=None,
                    help="resume from this .mfd instead of the pristine cell")
    ap.add_argument("--startcoc", default=None,
                    help="cocycle .npz matching --start (else rebuilt from "
                         "reference positions, which needs the crystal labels)")
    ap.add_argument("--knots", type=int, default=8)
    ap.add_argument("--min-sep", type=int, default=2)
    ap.add_argument("--burn", type=int, default=0)
    ap.add_argument("--span", type=int, default=6000)
    ap.add_argument("--ts", type=int, default=150)
    ap.add_argument("--snap", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    assert args.snap % args.ts == 0, "SNAP must be a multiple of TS"

    ddg.set_random_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    ref = ddg.Manifold.load(args.cell, 3)
    native = float(edges_from_facets(ref.facets())[1].mean())
    et = args.etarget if args.etarget is not None else native

    m = ddg.Manifold.load(args.start or args.cell, 3)
    placed = implant_knots(m, args.knots, args.min_sep, rng)
    n3_0 = m.num_facets

    params = ddg.SamplerParams(
        num_facets_target=ref.num_facets, num_facets_coef=0.1,
        hinge_degree_target=et, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * et / 6.0)
    s = ddg.ManifoldSampler(m, params)
    if args.zleg or args.cimp:
        s.set_n6_potential(args.zleg * args.lam, args.cimp * args.lam,
                           tilt=[0.0] * 5)
    if args.slide_prob:
        s.set_slide_prob(args.slide_prob)
    v = s.manifold
    if args.startcoc:
        e0, w0, _ = coc.load_cocycle(args.startcoc)
        s.enable_cocycle(np.asarray(e0), np.asarray(w0))
    else:
        edges = np.asarray(v.simplices(1))
        s.enable_cocycle(edges, coc.build_from_positions(
            edges, reference_frac_positions("r", args.mcell), args.mcell))
    s.check_cocycle()

    log = open(args.out + ".ts.jsonl", "w")
    t0 = time.time()
    if args.burn:
        s.run(sweeps=args.burn)
    action = ("EDQ-only" if not (args.zleg or args.cimp)
              else f"EDQ+n6(zleg={args.zleg},cimp={args.cimp})")
    print(f"[{os.path.basename(args.out)}] start lam={args.lam} [{action}] "
          f"e*={et:.9f} slide_prob={args.slide_prob} N3={n3_0} "
          f"implanted={len(placed)}", flush=True)

    done = args.burn
    nsnap = 0
    prev_tries = prev_acc = 0
    while done - args.burn < args.span:
        s.run(sweeps=args.ts)
        done += args.ts
        fac = np.asarray(v.facets())
        eu, edeg, V = edges_from_facets(fac)
        comps, imp = complexes(fac)
        sizes = sorted((len(c) for c in comps), reverse=True)
        # per-complex tree-lift centroids, min-imaged about each complex's
        # first member so a complex straddling a lift seam cannot split
        cents = []
        if comps:
            e1, w1 = s.read_cocycle()
            e1 = np.asarray(e1)
            pos = coc.tree_positions(e1, np.asarray(w1), int(e1.max()) + 1)[0]
            pbox = 1.0e6 * args.mcell
            for c in comps:
                p0 = pos[c[0]].astype(float)
                rel = (pos[c].astype(float) - p0 + pbox / 2) % pbox - pbox / 2
                cents.append([round(float(x), 1)
                              for x in (p0 + rel.mean(0)) % pbox])
        # per-complex illegal edge-degree signature: a live (3,4,4) knot
        edeg_of = {frozenset(int(x) for x in e): int(d)
                   for e, d in zip(eu, edeg) if d not in (5, 6)}
        sigs = []
        for c in comps:
            cs = set(c)
            sig = sorted(d for e, d in edeg_of.items() if e <= cs)
            sigs.append(sig)
        n_knot = sum(1 for sg in sigs if sg == [3, 4, 4])
        tries, acc = s.slide_stats() if args.slide_prob else (0, 0)
        rec = {"sweep": done, "t": round(time.time() - t0, 1),
               "n_illegal": int(sum(sizes)), "sizes": sizes,
               "members": comps, "cents": cents, "sigs": sigs,
               "n_knot344": n_knot,
               "slide_tries": tries, "slide_accepts": acc,
               "slide_tries_blk": tries - prev_tries,
               "slide_accepts_blk": acc - prev_acc,
               "mean_edeg": float(edeg.mean()),
               "legalvert": float(np.mean(imp == 0))}
        prev_tries, prev_acc = tries, acc
        log.write(json.dumps(rec) + "\n")
        log.flush()
        if (done - args.burn) % args.snap == 0:
            nsnap += 1
            stem = f"{args.out}_snap{done}"
            v.save(stem + ".mfd")
            e1, w1 = s.read_cocycle()
            e1s, w1s = coc.canonicalize_labels(np.asarray(e1), np.asarray(w1))
            coc.save_cocycle(stem + ".cocycle.npz", e1s, w1s, sweeps=done)
            drift = ""
            try:
                s.check_cocycle()
            except Exception as e:
                drift = f" COCYCLE-DRIFT {e}"
            eset = {tuple(sorted(e))
                    for e in np.asarray(v.simplices(1)).tolist()}
            cset = {tuple(sorted(e)) for e in np.asarray(e1).tolist()}
            if eset != cset:
                drift += f" COCYCLE-DETACHED({len(eset ^ cset)} edges)"
            print(f"[{os.path.basename(args.out)}] sw{done} "
                  f"({time.time()-t0:.0f}s): nill={sum(sizes)} "
                  f"knots344={n_knot} ncomp={len(sizes)} "
                  f"slides={acc}/{tries} <edeg>={edeg.mean():.9f}"
                  f" snap{nsnap}{drift}", flush=True)
    log.close()
    # The volume pin is soft (a quadratic penalty, not a constraint), and
    # implantation itself puts N3 one tet above target per knot, so N3 is
    # expected to relax back toward the target rather than sit still.
    print(f"[{os.path.basename(args.out)}] N3 {n3_0} -> {v.num_facets} "
          f"(target {ref.num_facets})")
    print(f"[{os.path.basename(args.out)}] DONE {done} sweeps, "
          f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
