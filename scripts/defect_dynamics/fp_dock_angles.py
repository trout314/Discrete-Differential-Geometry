#!/usr/bin/env python3
"""Segregate the FP production docks by BC-chain crossing angle.

A knot's chord spans one helix period, so the minimal-image displacement
between its two vertex labels in reference crystal coordinates is the
local BC-chain axis at that site. For every dock in the production data
(encounter + recombination episodes both record dock_chord), compute the
undirected angle in [0, 90] deg between A's chain axis at dock and B's
chain axis. Expect quantized angle classes (discrete crystallographic
chain directions); on-chain docks must land at 0 deg (consistency).

Outputs: angle histogram, dock census by angle class, and recombination
outcome (recombine vs escape) by the angle of the ORIGINAL dock.
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from cocycle_check import reference_frac_positions
from worm_helix import bc_orbit

ROOT = os.path.join(_HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")
REF = os.path.join(ROOT, "data", "tcp_reference", "T3_R_m2_N7248.mfd")
MCELL = 2


def main():
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    orb = [int(x) for x in bc_orbit(m, [int(x) for x in F[0]])]
    L = len(orb)
    frac = np.asarray(reference_frac_positions("r", MCELL)) * MCELL
    box = float(MCELL)

    # chord -> creation face (pristine crystal registry)
    from discrete_differential_geometry.fpkmc import face_apex_maps
    _, face_of = face_apex_maps(m)

    def unwrapped(vs):
        p = frac[vs[0]].copy()
        out = [p.copy()]
        for k in range(1, len(vs)):
            d = frac[vs[k]] - frac[vs[k - 1]]
            d -= box * np.round(d / box)
            p = p + d
            out.append(p.copy())
        return np.asarray(out)

    def axis(chord, rod_steps=8):
        """Local rod axis at the chord: walk the BC chain rod_steps
        each way from the chord's window and take the end-to-end
        displacement. rod_steps=8 spans two full periods of the BC
        walk's direction precession (period 8 windows, measured), so
        the winding cancels exactly. rod_steps=0
        falls back to the bare chord vector."""
        if rod_steps == 0:
            d = frac[chord[1]] - frac[chord[0]]
            d -= box * np.round(d / box)
            return d / np.linalg.norm(d)
        face = face_of[tuple(sorted(chord))]
        fwd = m.chain_walk([chord[0]] + list(face), rod_steps)
        bwd = m.chain_walk([chord[1]] + list(face), rod_steps)
        pf = unwrapped([int(x) for x in fwd])[-1]
        pb = unwrapped([int(x) for x in bwd])[-1]
        d0 = frac[int(bwd[0])] - frac[int(fwd[0])]
        d0 -= box * np.round(d0 / box)
        d = (pf - unwrapped([int(x) for x in fwd])[0]) \
            - (pb - unwrapped([int(x) for x in bwd])[0]) - d0
        return d / np.linalg.norm(d)

    def angle_deg(c1, c2):
        c = abs(float(np.dot(axis(c1), axis(c2))))
        return float(np.degrees(np.arccos(min(1.0, c))))

    def bchord(jB):
        return tuple(sorted((orb[jB % L], orb[(jB + 4) % L])))

    # ---------------- collect docks ----------------
    enc = []           # (angle, sep, on_chain, V_dock)
    rec = []           # (angle, outcome)
    for f in sorted(glob.glob(os.path.join(DATA, "prodA_s*_p*.json"))):
        d = json.load(open(f))
        for e in d["episodes"]:
            if e.get("dock_chord"):
                a = angle_deg(tuple(e["dock_chord"]), bchord(e["jB"]))
                enc.append((a, d["sep_sites"],
                            e.get("on_chain_site") is not None,
                            e.get("V_dock")))
    for f in sorted(glob.glob(os.path.join(DATA, "prodB_p*.json"))) + \
            sorted(glob.glob(os.path.join(DATA, "prodB2_p*.json"))):
        d = json.load(open(f))
        for e in d["episodes"]:
            if e.get("dock_chord") and e.get("freed") and "outcome" in e:
                a = angle_deg(tuple(e["dock_chord"]), bchord(e["jB"]))
                rec.append((a, e["outcome"]))

    angles = np.array([x[0] for x in enc])
    print(f"{len(enc)} encounter docks, {len(rec)} recombination docks "
          f"with outcomes")

    # BC-walk direction persistence along the orbit: u_k . u_{k+lag},
    # u_k the unit step between window heads k, k+4. Decides whether
    # chains are straight rods (plateau) or curve (decay).
    P = unwrapped(orb + orb[:5])
    U = P[4:] - P[:-4]
    U /= np.linalg.norm(U, axis=1)[:, None]
    print("\nBC-walk direction autocorrelation <u_k . u_{k+lag}>:")
    for lag in (1, 2, 4, 8, 12, 16, 24, 32, 48):
        c = float(np.mean(np.sum(U[:-lag] * U[lag:], axis=1)))
        print(f"  lag {lag:3d} windows: {c:+.3f}")

    # continuous spectrum -> coarse 15-deg bins
    classes = [7, 22, 37, 52, 67, 82]        # bin centers, 15-deg wide
    def cls_of(a):
        return classes[min(len(classes) - 1, int(a // 15))]
    print("\nangle distribution (15-deg bins, encounter docks):")
    for c in classes:
        n = sum(1 for a in angles if cls_of(a) == c)
        print(f"  {c - 7:2d}-{c + 8:2d} deg: {n}")

    print("\nrecombination by dock angle:")
    stats = {}
    for a, o in rec:
        c = cls_of(a)
        stats.setdefault(c, {}).setdefault(o, 0)
        stats[c][o] += 1
    for c in sorted(stats):
        s = stats[c]
        nr, ne = s.get("recombine", 0), s.get("escape", 0)
        tot = nr + ne
        print(f"  {c:3d} deg: {tot:3d} freed docks, P(rec) = "
              f"{nr}/{tot} = {nr / max(1, tot):.2f}")

    # V_dock by class (contact docks live where?)
    print("\ncontact docks (V > 1e-6) by angle:")
    for a, s0, onc, v in enc:
        if v is not None and abs(v) > 1e-6:
            print(f"  angle {a:.1f} deg, s0={s0}, on_chain={onc}, "
                  f"V={v:+.3f}")

    # ---------------- figure ----------------
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    bins = np.arange(-2.5, 95, 5)
    ax.hist(angles, bins=bins, color="tab:blue", alpha=0.8,
            label="docks")
    th = np.radians(np.clip(bins, 0, 90))
    iso = len(angles) * (np.cos(th[:-1]) - np.cos(th[1:]))
    ax.plot(0.5 * (bins[:-1] + bins[1:]), iso, "k--", lw=1.5,
            label="isotropic (sin)")
    ax.set_xlabel("rod crossing angle (deg)")
    ax.set_ylabel("dock count")
    ax.set_title(f"dock angles ({len(enc)} encounter docks)")
    ax.legend(fontsize=8)

    ax = axs[1]
    seps = sorted({x[1] for x in enc})
    bot = np.zeros(len(classes))
    for s0 in seps:
        cnt = [sum(1 for x in enc if x[1] == s0 and cls_of(x[0]) == c)
               for c in classes]
        ax.bar([str(c) for c in classes], cnt, bottom=bot,
               label=f"s0={s0}")
        bot += np.array(cnt, float)
    ax.set_xlabel("angle class (deg)")
    ax.set_ylabel("dock count")
    ax.set_title("dock census by angle class (stacked by s0)")
    ax.legend(fontsize=7)

    ax = axs[2]
    cs = sorted(stats)
    pr = [stats[c].get("recombine", 0) /
          max(1, sum(stats[c].values())) for c in cs]
    n = [sum(stats[c].values()) for c in cs]
    err = [np.sqrt(p * (1 - p) / max(1, k)) for p, k in zip(pr, n)]
    ax.bar([str(c) for c in cs], pr, yerr=err, capsize=4,
           color="tab:orange", alpha=0.8)
    for i, k in enumerate(n):
        ax.text(i, 0.03, f"n={k}", ha="center", fontsize=8)
    ax.axhline(61 / 87, color="k", ls="--", lw=1, label="pooled 0.70")
    ax.set_ylim(0, 1)
    ax.set_xlabel("dock angle class (deg)")
    ax.set_ylabel("P(recombine | freed)")
    ax.set_title("recombination by dock angle")
    ax.legend(fontsize=8)

    fig.suptitle(
        "FP production docks by BC-chain crossing angle -- R m2 N7248, "
        "lam=0.40, slide channel, frozen background; chord vector = local "
        "chain axis (reference registry)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(DATA, "fp_dock_angles.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
