#!/usr/bin/env python3
"""Fit the small-k exponent alpha of the curvature-charge structure factor
S(k) ~ k^alpha, with jackknife-over-chains error bars.

S_obs is mean-subtracted (isolates the CHARGE arrangement, not the point
pattern). Sub-Bragg window k in [KMIN, KCUT] (charge Bragg peaks start ~k=5).
Pool modes over all snapshots in a group; bin in log-k; the binned ensemble
power <|A(k)|^2> is averaged then logged; weighted LS gives alpha; jackknife
over independent chains gives the error bar.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "python")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from discrete_differential_geometry.structure_factor import structure_factor

SP = sys.argv[1]
NMAX = 8            # k<5 fit window needs |n|<~5; nmax=8 covers it + Bragg onset
KMIN, KCUT = 1.0, 5.0
NSHELL = 9


def sk_of(p):
    fac = np.asarray(Manifold.load(p, 3).facets())
    edges, omega, _ = coc.load_cocycle(os.path.splitext(p)[0] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    q = FIELDS["curvature_charge"](fac)
    kmag, s_obs, s_null = structure_factor(frac, basis, q, NMAX)
    return kmag, s_obs


def bin_curve(kmag, s_obs):
    """Sub-Bragg binned S(k): geometric-mean k, mean power, SEM-of-power per shell."""
    m = (kmag >= KMIN) & (kmag <= KCUT)
    k, s = kmag[m], s_obs[m]
    edg = np.geomspace(KMIN, KCUT, NSHELL + 1)
    kc, sm, se, nn = [], [], [], []
    for lo, hi in zip(edg[:-1], edg[1:]):
        sel = (k >= lo) & (k < hi)
        if sel.sum() >= 3:
            kc.append(np.exp(np.log(k[sel]).mean()))
            sm.append(s[sel].mean())
            se.append(s[sel].std(ddof=1) / np.sqrt(sel.sum()))
            nn.append(sel.sum())
    return map(np.array, (kc, sm, se, nn))


def fit_alpha(kc, sm, se):
    """Weighted LS slope of log S vs log k; returns (alpha, intercept)."""
    x, y = np.log(kc), np.log(sm)
    w = (sm / np.maximum(se, 1e-30)) ** 2               # weight = 1/var(log S)
    a, b = np.polyfit(x, y, 1, w=np.sqrt(w))
    return a, b


def alpha_of_snaps(snaps):
    kmag = np.concatenate([S[0] for S in snaps])
    s_obs = np.concatenate([S[1] for S in snaps])
    kc, sm, se, nn = bin_curve(kmag, s_obs)
    a, b = fit_alpha(kc, sm, se)
    return a, (kc, sm, se, nn)


def jackknife_alpha(chains):
    """chains: list of lists-of-snapshot-Sk. Jackknife alpha over chains."""
    n = len(chains)
    full = [s for ch in chains for s in ch]
    a_full, curve = alpha_of_snaps(full)
    a_jk = []
    for i in range(n):
        loo = [s for j, ch in enumerate(chains) if j != i for s in ch]
        a_jk.append(alpha_of_snaps(loo)[0])
    a_jk = np.array(a_jk)
    se = np.sqrt((n - 1) / n * np.sum((a_jk - a_jk.mean()) ** 2))
    return a_full, se, curve


def chains_for(prefix, ids):
    out = []
    for i in ids:
        snaps = [sk_of(p) for p in sorted(glob.glob(f"{SP}/run5h/{prefix}{i}_snap*.mfd"))]
        if snaps:
            out.append(snaps)
    return out


print("computing S(k) for all snapshots (nmax=%d) ..." % NMAX, flush=True)
above = chains_for("above1", [0, 1, 2, 3])            # above10..above13
below = chains_for("below", [0, 1, 2, 3])
constrained = above + below
crystal = sk_of(f"{SP}/r_m4.mfd")

# diagnostic: print the pooled binned curve for crystal vs constrained
for tag, snaps in [("perfect crystal", [crystal]),
                   ("all constrained", [s for ch in constrained for s in ch])]:
    km = np.concatenate([S[0] for S in snaps]); so = np.concatenate([S[1] for S in snaps])
    kc, sm, se, nn = bin_curve(km, so)
    print(f"\n-- {tag}: sub-Bragg binned S(k) --")
    for a_, b_, c_ in zip(kc, sm, nn):
        print(f"   k={a_:6.3f}  S={b_:10.4e}  n={int(c_)}")
    print(f"   naive slope log-log: {np.polyfit(np.log(kc), np.log(sm), 1)[0]:+.3f}")
print()

fig, ax = plt.subplots(figsize=(7.2, 5.6))
results = {}
for name, chains, col in [("above-native (4 chains)", above, "C2"),
                          ("below-native (4 chains)", below, "C3"),
                          ("all constrained (8 chains)", constrained, "C0")]:
    a, se, (kc, sm, sem, nn) = jackknife_alpha(chains)
    results[name] = (a, se)
    print(f"{name:28s}: alpha = {a:.3f} +/- {se:.3f}   "
          f"(shells n_modes {nn.min()}-{nn.max()})")
    if name.startswith("all"):
        ax.errorbar(kc, sm, yerr=sem, fmt="o", color=col, ms=6, capsize=3,
                    label=f"S(k), all constrained")
        kk = np.array([kc.min(), kc.max()])
        ax.plot(kk, np.exp(np.polyfit(np.log(kc), np.log(sm), 1, w=sm / sem)[1])
                * kk ** a, "-", color=col, lw=1.8,
                label=fr"fit $\alpha={a:.2f}\pm{se:.2f}$")
    else:
        ax.plot(kc, sm, "o--", color=col, ms=3.5, lw=0.9, alpha=0.75,
                label=f"{name.split(' (')[0]}: "+fr"$\alpha={a:.2f}\pm{se:.2f}$")

# reference slopes
kk = np.geomspace(KMIN, KCUT, 2)
anchor = results["all constrained (8 chains)"]
for al, ls, lab in [(2.0, "--", r"$k^{2}$ (Gauss / Hamiltonian constraint)"),
                    (1.0, "-.", r"$k^{1}$")]:
    ax.plot(kk, 0.006 * (kk / KMIN) ** al, ls, color="k", lw=1.1, label=lab)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$|k|$ (units of smallest reciprocal length)")
ax.set_ylabel(r"charge structure factor  $S(k)$")
ax.set_title("Small-k exponent of the curvature-charge structure factor\n"
             "(mean-subtracted; sub-Bragg window k<5; jackknife-over-chains errors)")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.25, which="both")
fig.tight_layout()
out = f"{SP}/sk_exponent.png"
fig.savefig(out, dpi=130)
print(f"\nSaved to: {out}")
