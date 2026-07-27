#!/usr/bin/env python3
"""Phase-1 knot collider: two knots on one BC helix, exact V(s), contact
outcome. Deterministic -- no sampling anywhere.

Preparation exploits the no-halo theorem: a knot leaves EXACTLY pristine
crystal behind it, so creating knot A on a chain window, sliding it +k
slides down the helix, and creating knot B at the original window yields a
controlled two-knot state on one helix at separation 4k chain steps.

The BC chain is mapped once with worm_helix.bc_orbit (sliding-window walk;
cyclic vertex list v[0..L-1]); a knot on the chain has chord
(v_j, v_j+4) and the face between windows k, k+1 is (v_k+1, v_k+2, v_k+3)
with apexes (v_k, v_k+4), so creation and steering are exact chain
arithmetic, not heuristics.

Measured while walking A back toward B, separation s = chain steps between
chords:

  V(s)      full shape action (worm_slide.dS_between vs the crystal, lam=1
            units) minus two isolated-knot creation costs. The no-halo
            theorem predicts EXACTLY 0 until skins touch; any nonzero value
            at range would falsify it dynamically.
  mobility  number of template-valid and clean slide directions for A --
            where the slide channel jams as the knots approach.
  census    components and signatures at each s (when do "two knots" stop
            being two (3,4,4)s?).

Endgame at contact: classify the outcome (merge -> species f/sig/sumZ;
jam -> closest approach), then test reversibility (slide A back out one
step -- does it un-fuse into two clean knots?).

FRAME NOTE: registry coordinates are used for TOPOLOGICAL
WINDING only (integer-valued, frame-robust; CONVENTIONS.md
sec 6).
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
import worm_slide as ws
import worm_moves as wm
from worm_helix import bc_orbit, orbit_winding
from cocycle_check import reference_frac_positions
import defect_state as ds

ESTAR = 5.105025


def chain_face(chain, j):
    """Face between windows j, j+1 (indices mod L)."""
    L = len(chain)
    return [chain[(j + 1) % L], chain[(j + 2) % L], chain[(j + 3) % L]]


def make_knot(m, chain, j):
    """2->3 on the window-j face; returns the chord (v_j, v_j+4) sorted."""
    L = len(chain)
    face = chain_face(chain, j)
    a = m.face_apexes(*face)
    want = {chain[j % L], chain[(j + 4) % L]}
    assert set(a) == want, f"chain face apexes {a} != {want}"
    m.do_bistellar_move(sorted(face), sorted(want))
    return tuple(sorted(want))


def slide_along(m, chain, j, fwd, recs_out=None):
    """One slide of the knot whose chord is (chain[j], chain[j+4]), in the
    chain direction `fwd`. Deterministic: the only accepted candidate is the
    one whose new chord matches the chain arithmetic. Returns
    (new_j, n_template_valid, n_clean_chain) or (None, nv, nc).

    recs_out: optional list; the applied slide's undo records are appended,
    so a caller can unwind with ws.undo_slide (memory-flat sweeps restore
    state by inverse moves instead of reloading manifolds)."""
    L = len(chain)
    chord = tuple(sorted((chain[j % L], chain[(j + 4) % L])))
    jn = (j + 4) % L if fwd else (j - 4) % L
    want = tuple(sorted((chain[jn % L], chain[(jn + 4) % L])))
    cands = ws.candidate_slides(m, chord)
    nclean = 0
    for cs, mv in cands:
        newc = tuple(sorted((cs["c4"], cs["c8"])))
        if newc != want:
            continue
        recs = ws.apply_slide(m, mv)
        if recs is None:
            continue
        nclean += 1
        if recs_out is not None:
            recs_out.append(recs)
        return jn, len(cands), nclean
    return None, len(cands), nclean


def separation(L, jA, jB):
    d = (jA - jB) % L
    return min(d, L - d)


def census_str(m):
    st = ds.DefectState(m)
    out = []
    for cx in st.components(broad=True):
        f = st.induced_shape(cx.verts)["f"]
        out.append((tuple(f), tuple(cx.sig),
                    st.total_coordination(cx.verts)))
    del st
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref",
                    default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--kmax", type=int, default=10,
                    help="initial separation in slides (4 chain steps each)")
    ap.add_argument("--sites", type=int, default=3,
                    help="number of distinct chain windows to repeat at")
    ap.add_argument("--seed-tet", type=int, default=0)
    ap.add_argument("--lam", type=float, default=0.40,
                    help="reporting only: V is stored in lam=1 units")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ref = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(ref.facets())
    chain = bc_orbit(ref, [int(x) for x in F[args.seed_tet]])
    L = len(chain)
    rp = np.asarray(reference_frac_positions("r", args.mcell))
    wind = orbit_winding(chain, rp, args.mcell)
    print(f"orbit: L = {L} chain steps, winding {tuple(wind)} "
          f"(box periods); {L // 4} slide positions")
    nuniq = len(set(chain))
    print(f"  ({nuniq} distinct vertices; chain "
          f"{'revisits vertices' if nuniq < L else 'is simple'})")

    em0 = ws.edeg_dict(ref)
    fs0 = ws.facset(ref)
    del ref

    results = []
    step = max(1, (L // 4) // max(args.sites, 1))
    for si in range(args.sites):
        b0 = (si * step * 4) % L
        m = ddg.Manifold.load(args.ref, 3)
        # -- knot A + its isolated creation cost
        chordA = make_knot(m, chain, b0)
        S1 = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
        cen = census_str(m)
        rungA = cen[0][2]
        print(f"\n=== site {si}: window {b0}, knot A chord {chordA}, "
              f"sig {cen[0][1]}, sumZ {rungA}, S1 = {S1:.4f} (lam=1)")
        # -- walk A out kmax slides, recording the single-knot energy
        # landscape S_single(j): the WASHBOARD. The knot's decorated species
        # (rung) varies with chain site, so its isolated energy does too;
        # the interaction V(s) is defined against the measured single-knot
        # energy AT THE SAME SITE, not against the creation cost.
        jA = b0
        Ssingle = {b0: S1}
        stuck = False
        for k in range(args.kmax):
            jA2, nv, nc = slide_along(m, chain, jA, fwd=True)
            if jA2 is None:
                print(f"  A stuck after {k} outbound slides "
                      f"({nv} template-valid)")
                stuck = True
                break
            jA = jA2
            Ssingle[jA] = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
        if stuck:
            continue
        wb = [Ssingle[b0 + 4 * k] for k in range(args.kmax + 1)
              if (b0 + 4 * k) in Ssingle]
        print(f"  washboard S_single along walk: "
              f"{', '.join(f'{x:.2f}' for x in wb)}")
        # -- knot B at the original window (no-halo: site is pristine again)
        chordB = make_knot(m, chain, b0)
        jB = b0
        S2c = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
        print(f"  created B at window {b0}; S(A,B) = {S2c:.4f}, "
              f"2*S1 = {2 * S1:.4f}, V = {S2c - 2 * S1:+.6f}")
        # -- walk A back toward B, recording V(s)
        rows = []
        outcome = None
        while True:
            s = separation(L, jA, jB)
            S = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
            cen = census_str(m)
            ncomp = len(cen)
            # mobility of A toward B before the next step
            base = Ssingle.get(jA)
            V = S - base - S1 if base is not None else None
            rows.append({"sep_chain_steps": s, "V": V,
                         "S_total": S, "S_single_A": base,
                         "n_components": ncomp,
                         "census": [list(map(list, c[:2])) + [c[2]]
                                    for c in cen]})
            print(f"  s={s:3d} chain steps: V = "
                  f"{(f'{V:+.6f}' if V is not None else '   n/a  ')}  "
                  f"components = {ncomp} "
                  f"{[c[1] for c in cen] if ncomp <= 3 else ''}")
            if ncomp < 2:
                outcome = {"type": "merged", "at_sep": s,
                           "species_f": list(cen[0][0]),
                           "sig": list(cen[0][1]), "sumZ": cen[0][2],
                           "V": V, "S_total": S}
                break
            jA2, nv, nc = slide_along(m, chain, jA, fwd=False)
            if jA2 is None:
                outcome = {"type": "jammed", "at_sep": s,
                           "template_valid": nv, "clean": nc,
                           "V": S - 2 * S1}
                print(f"  JAMMED at s={s}: {nv} template-valid, "
                      f"{nc} clean toward B")
                break
            jA = jA2
        # -- reversibility: does any single slide un-fuse the product
        # into two clean knots?
        rev = None
        if outcome and outcome["type"] == "merged":
            e3, _ = ws.knot_chords(m)
            rev = {"chords_after_merge": [list(c) for c in e3],
                   "unfuses": False}
            for c in e3:
                for cs, mv in ws.candidate_slides(m, tuple(c)):
                    recs = ws.apply_slide(m, mv)
                    if recs is None:
                        continue
                    cen2 = census_str(m)
                    # two components, each a decorated knot (exactly one
                    # degree-3 chord) -- rung decorations vary by site, so
                    # demanding bare (3,4,4) would wrongly reject
                    if (len(cen2) == 2
                            and all(list(x[1]).count(3) == 1
                                    for x in cen2)):
                        rev["unfuses"] = True
                        rev["via_chord"] = list(c)
                        break
                    ws.undo_slide(m, recs)
                if rev["unfuses"]:
                    break
        results.append({"site": si, "window": b0, "rungA": rungA,
                        "S1": S1, "rows": rows, "outcome": outcome,
                        "reversal": rev})
        del m

    print(f"\n{'='*66}\nSUMMARY over {len(results)} sites "
          f"(V in lam=1 shape units; x{args.lam} for the run action):")
    for r in results:
        far = [x for x in r["rows"] if x["sep_chain_steps"]
               and x["sep_chain_steps"] > 8 and x["V"] is not None]
        vmax = max((abs(x["V"]) for x in far), default=0.0)
        o = r["outcome"]
        print(f"  site {r['site']} (rung {r['rungA']}): "
              f"max|V| at s>8: {vmax:.2e}  outcome: {o['type']} at s="
              f"{o['at_sep']}"
              + (f" -> f={tuple(o['species_f'])} sig={tuple(o['sig'])} "
                 f"sumZ={o['sumZ']}"
                 + (f" V={o['V']:+.3f}" if o.get("V") is not None else "")
                 if o["type"] == "merged" else "")
              + (f"  unfuses={r['reversal']['unfuses']}"
                 if r.get("reversal") else ""))

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"ref": args.ref, "orbit_len": L,
                       "winding": [int(x) for x in wind],
                       "kmax": args.kmax, "results": results}, fh, indent=1)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
