#!/usr/bin/env python3
"""Is the DEFECT DENSITY hyperuniform? -- and why the curvature-charge S(k)
cannot answer that question.

Motivation (2026-08-05). [[statics-hu-verdict]] reported the melt's total
curvature field as "crystal-grade class-I hyperuniform" (sigma^2_Q ~ R^1.98 vs
pristine 2.01; sk_torus low-k ratio ~2e-3). That is a true measurement of the
wrong thing. Two independent confounds make it blind to the defect
arrangement, and this script measures both instead of assuming them:

 1. NORMALIZATION. S(k) is normalised by S2 = sum_v dq_v^2, and in a dilute
    gas the periodic CRYSTAL carries essentially all of it (measured: 99.5%
    at lam=0.40) while contributing nothing at sub-Bragg k. So the low-k
    ratio is bounded by roughly the defect share of S2 no matter how the
    defects are arranged. Reported here as `s2_frac_defect`.
 2. THE NULL. sk_torus normalises against the full charge-permutation null,
    which scrambles the crystal itself -- so the comparison is
    melt-vs-amorphous, and any crystal-like state passes trivially. The null
    that answers the question RELOCATES the defect charges and leaves the
    crystal alone (`reloc_null_power` below -- and see its docstring for the
    biased first attempt that the --control run caught).

The fields, and what each is good for:

  defect charge     the anomalous dq on impurity vertices, zero elsewhere,
                    vs the matched relocation null: does the actual
                    arrangement suppress low-k curvature power relative to
                    the same charges placed at random? (Source charge only --
                    the legal-vertex screening halo is not separable from the
                    crystal without a background model.)
  defect_indicator  1 on impurity vertices. Its analytic permutation null
                    coincides with "same number of defect vertices placed at
                    random on the same skeleton" to O(n_d/V), so
                    S_obs/S_null IS the defect point-process structure
                    factor. But it is form-factor-loaded: defects come in
                    complexes of mean size <m>, which inflates S(k->0) by
                    <m^2>/<m> even for perfectly Poisson complex centres.
  centroids         one point per connected defect complex, against the
                    random-crystal-site null. Form-factor free -- this is the
                    inter-complex arrangement, the physical quantity.

Consistency check the script reports: vertex_ratio ~= centroid_ratio *
form_factor. If that fails, one of the three estimators is wrong.

VALIDATION: `--control` re-runs the whole pipeline on a synthetic state whose
defect vertices are relocated uniformly at random over the skeleton. Every
ratio must come back 1.0 within error; that is the end-to-end estimator test.

FRAME: harmonic (periodic Tutte) cocycle embedding, as sk_torus
(CONVENTIONS.md sec 6). Small-k exponents and RATIOS are
frame-robust; k-values and bare amplitudes are gauge.

Usage:
    python scripts/defect_dynamics/defect_density_hu.py \
        data/mgas/lam40_snap*.mfd data/mgas/lam40x_snap*.mfd \
        --group '_snap\\d+$' --out data/defect_hu/lam40.json \
        --plot data/figs/defect_density_hu_lam40.png
"""
import argparse
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS, edges_and_degrees

# The k -> 0 window, in units of the smallest reciprocal length. Since
# torus_positions whitens (2026-08-13) this is a PHYSICAL |k| cut, not the
# |n| <= 2 mode count it used to coincide with: on a non-cubic host the two
# differ, and modes along the long cell axis now correctly land at smaller k.
LOWK = 2.0 + 1e-9


# ---------------------------------------------------------------------------
# k-space primitives
# ---------------------------------------------------------------------------

def kvectors(nmax):
    """Half-space of nonzero commensurate reciprocal indices (real field)."""
    r = np.arange(-nmax, nmax + 1)
    n = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    n = n[np.any(n != 0, axis=1)]
    keep = ((n[:, 0] > 0)
            | ((n[:, 0] == 0) & (n[:, 1] > 0))
            | ((n[:, 0] == 0) & (n[:, 1] == 0) & (n[:, 2] > 0)))
    return n[keep]


def phase_matrix(frac, nvec):
    """exp(2 pi i frac . n) as an explicit (V, K) matrix.

    Built ONCE per snapshot and reused for the observed power, the skeleton
    |F|^2, and every relocation-null draw -- the transcendentals dominate, so
    recomputing them per shuffle made the null ~30x more expensive than the
    measurement it normalises."""
    return np.exp(2j * np.pi * frac @ nvec.T)


def power(PH, w):
    """|sum_v w_v exp(i k . x_v)|^2 at each k (w need not be centred)."""
    return np.abs(w @ PH) ** 2


def skeleton_F2(PH):
    """|F(k)|^2 of the bare vertex set -- the skeleton's own structure, which
    every random-placement null must carry."""
    return np.abs(PH.sum(axis=0)) ** 2


def point_process_sk(PH, sel, V, F2):
    """Structure factor of the point subset `sel` of the skeleton `frac_all`,
    against the EXACT random-placement null.

    Choosing n of V sites uniformly without replacement gives
        E|A(k)|^2 = n + n(n-1)(|F|^2 - V) / (V(V-1)),
    so S_null(k) = E|A|^2 / n is analytic -- no shuffling. At sub-Bragg k,
    |F|^2 -> 0 and S_null -> 1 (the Poisson reference).
    """
    n = int(sel.sum())
    if n < 2:
        return None
    s_obs = power(PH, sel.astype(float)) / n
    s_null = 1.0 + (n - 1) * (F2 - V) / (V * (V - 1))
    return s_obs, s_null


def points_sk(frac_pts, nvec, V, F2):
    """Same, for points NOT drawn from the skeleton (complex centroids): the
    observed power uses the true centroid positions, the null still places the
    same number of points on random crystal sites (centroids can only sit
    near lattice sites, so that -- not a continuum Poisson -- is the honest
    reference)."""
    n = len(frac_pts)
    if n < 2:
        return None
    s_obs = power(phase_matrix(frac_pts, nvec), np.ones(n)) / n
    s_null = 1.0 + (n - 1) * (F2 - V) / (V * (V - 1))
    return s_obs, s_null


# ---------------------------------------------------------------------------
# geometry / combinatorics
# ---------------------------------------------------------------------------

def complexes(eu, defect, V, ill_mask=None):
    """Connected components -> label array over defect vertices (-1
    elsewhere) and component sizes.

    Default: the historical defect SUBGRAPH (all edges among defect
    vertices). With ill_mask (boolean over eu rows): components of the
    ILLEGAL-EDGE GRAPH -- required at high defect-vertex density, where the
    subgraph definition percolates (a strain gas has ~44% defect vertices
    and the subgraph fuses everything into one complex, killing the
    centroid estimator)."""
    if ill_mask is not None:
        e = eu[ill_mask]
    else:
        idx = np.nonzero(defect)[0]
        sub = np.isin(eu[:, 0], idx) & np.isin(eu[:, 1], idx)
        e = eu[sub]
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(V, V))
    _, lab = connected_components(A + A.T, directed=False)
    lab = np.where(defect, lab, -1)
    keep = np.unique(lab[lab >= 0])
    remap = -np.ones(lab.max() + 2, np.int64)
    remap[keep] = np.arange(len(keep))
    lab = np.where(lab >= 0, remap[lab], -1)
    sizes = np.bincount(lab[lab >= 0], minlength=len(keep))
    return lab, sizes


def centroids(frac, lab, ncomp):
    """Circular (torus) mean position of each complex -- a plain mean would be
    wrong for a complex straddling the periodic boundary."""
    out = np.empty((ncomp, 3))
    z = np.exp(2j * np.pi * frac)
    for c in range(ncomp):
        m = lab == c
        out[c] = (np.angle(z[m].mean(axis=0)) / (2 * np.pi)) % 1.0
    return out


def defect_charge_field(dq, defect):
    """The defect-ATTRIBUTABLE charge: the anomalous dq on impurity vertices,
    zero on the crystal. Excludes the legal-vertex screening halo, which cannot
    be separated from the crystal without a background model -- so this is the
    SOURCE charge, not the total field."""
    w = np.zeros_like(dq)
    w[defect] = dq[defect]
    return w


def rigid_complex_sk(PH, w, lab, ncomp, rng=None, ncheck=0):
    """S(k) of the defect population under RIGID-COMPLEX relocation -- each
    complex moved as a whole to a random place, internal structure intact.

    This is the null that answers 'are the complexes ARRANGED with
    correlation?'. Relocating individual vertex charges (reloc_null_power)
    instead destroys the intra-complex +/- cancellation, so it measures how
    neutral a complex is, not where the complexes are.

    With F_i(k) = sum_{v in i} w_v exp(i k . x_v) the complex's own amplitude
    and independent uniform torus translations t_i, every cross term carries
    E[exp(i k . (t_i - t_j))] = 0 at nonzero commensurate k, so the null is
    ANALYTIC:
        obs(k)  = |sum_i F_i(k)|^2
        null(k) = sum_i |F_i(k)|^2
    Each complex's form factor AND net charge divide out; only the relative
    PHASES of the complexes -- i.e. their arrangement -- survive. As k -> 0,
    F_i -> Q_i and the ratio becomes the charge-weighted centroid S(k).

    Complexes are translated, NOT reoriented: the crystal locks defect
    orientation (the P2 null is +0.118, not 0), so this asks about positional
    arrangement alone.

    ncheck > 0 runs the Monte-Carlo self-test -- apply actual random phase
    shifts and confirm the mean recovers the analytic null.
    """
    Fi = np.empty((ncomp, PH.shape[1]), complex)
    for i in range(ncomp):
        m = lab == i
        Fi[i] = w[m] @ PH[m]
    obs = np.abs(Fi.sum(axis=0)) ** 2
    null = (np.abs(Fi) ** 2).sum(axis=0)
    mc = None
    if ncheck:
        acc = np.zeros(PH.shape[1])
        for _ in range(ncheck):
            ph = np.exp(2j * np.pi * rng.random((ncomp, 1)))   # per-complex phase
            acc += np.abs((Fi * ph).sum(axis=0)) ** 2
        mc = acc / ncheck
    return obs, null, mc, Fi


def complex_charges(dq, lab, ncomp, defect):
    """Net charge Q_i = sum_{v in i} dq_v per complex, and the neutrality
    fraction |Q_i| / sum_{v in i} |dq_v| (0 = perfectly neutral multipole,
    1 = all vertices the same sign)."""
    Q = np.array([dq[lab == i].sum() for i in range(ncomp)])
    A = np.array([np.abs(dq[lab == i]).sum() for i in range(ncomp)])
    return Q, np.where(A > 0, np.abs(Q) / np.maximum(A, 1e-300), 0.0)


def reloc_null_power(PH, vals, rng, nshuf):
    """Relocation null for the defect charge field: the SAME multiset of
    anomalous values placed on uniformly random sites, zero elsewhere.

    NOTE (the control caught this): the obvious null -- transpose the defect
    charges with those of random crystal sites, keeping the crystal pattern --
    is BIASED LOW. A transposition perturbs both ends, so the null state
    carries strictly more perturbation than the observed one, and the ratio
    comes back 0.61 +/- 0.06 on randomised defects where it must be 1.0. That
    would have read as a 40% low-k suppression that is pure construction.
    Relocating a zero-background field is exactly matched: observed and null
    differ only in the POSITIONS of an identical value multiset."""
    V, n = PH.shape[0], len(vals)
    acc = np.zeros(PH.shape[1])
    for _ in range(nshuf):
        w = np.zeros(V)
        w[rng.choice(V, size=n, replace=False)] = vals
        acc += power(PH, w - w.mean())
    return acc / nshuf


# ---------------------------------------------------------------------------
# per-snapshot measurement
# ---------------------------------------------------------------------------

def measure(mfd, nmax, rng, nshuf, control=False, whiten=True,
            complex_def="illegal"):
    facets = np.asarray(ddg.Manifold.load(mfd, 3).facets())
    edges, omega, _ = coc.load_cocycle(os.path.splitext(mfd)[0] + ".cocycle.npz")
    frac, basis = coc.torus_positions(facets, edges, omega, whiten=whiten)
    eu, ecnt, deg, V = edges_and_degrees(facets)

    defect = FIELDS["defect_indicator"](facets) > 0
    qR = FIELDS["curvature_charge"](facets)
    dq = qR - qR.mean()

    if control:
        # End-to-end estimator test: block-swap the whole defect sector (both
        # the indicator and its anomalous charges) onto random crystal sites.
        # The relocated set is scattered, so it also has no complexes -- the
        # form factor must collapse to 1 along with every ratio.
        src = np.nonzero(defect)[0]
        tgt = rng.choice(np.nonzero(~defect)[0], size=len(src), replace=False)
        dq2, newd = dq.copy(), defect.copy()
        dq2[tgt], dq2[src] = dq[src], dq[tgt]
        newd[tgt], newd[src] = True, False
        dq, defect = dq2, newd

    nvec = kvectors(nmax)
    kmag = np.linalg.norm(nvec @ np.linalg.inv(basis).T, axis=1) \
        / np.abs(np.linalg.inv(basis)).max()
    low = kmag <= LOWK
    PH = phase_matrix(frac, nvec)
    F2 = skeleton_F2(PH)

    # Which graph defines a "complex" -- this is NOT a detail, it decides
    # what the arrangement estimators measure (see --complex-def).  The
    # --control state's relocated scatter has no illegal edges, so it always
    # uses the subgraph; at control densities the two coincide anyway.
    ill_mask = None if (control or complex_def == "subgraph") \
        else (ecnt < 5) | (ecnt > 6)
    lab, sizes = complexes(eu, defect, V, ill_mask)
    ncomp = len(sizes)
    cen = centroids(frac, lab, ncomp) if ncomp >= 2 else None

    out = {"path": os.path.basename(mfd), "V": V, "n_defect": int(defect.sum()),
           "n_complex": ncomp,
           "mean_size": float(sizes.mean()) if ncomp else 0.0,
           "form_factor": float((sizes ** 2).mean() / sizes.mean()) if ncomp else 0.0,
           "s2_frac_defect": float(dq[defect] @ dq[defect] / (dq @ dq)),
           # principal cell lengths (singular values) -- the diagonal is not
           # the period once the basis is whitened / Euclidean-reduced
           "M": np.linalg.svd(basis, compute_uv=False).astype(float).tolist()}

    # --- 1. charge field, sk_torus's own (permutation) null, for reference
    S2 = float(dq @ dq)
    p_obs = power(PH, dq)
    s_null_perm = 1.0 - (F2 - V) / (V * (V - 1))
    out["charge_permnull_lowk"] = float(np.mean(p_obs[low] / S2 / s_null_perm[low]))

    # --- 2. defect SOURCE charge vs the matched RELOCATION null: does the
    #        actual arrangement suppress low-k curvature power relative to the
    #        same charges placed at random?
    wdef = defect_charge_field(dq, defect)
    p_def = power(PH, wdef - wdef.mean())
    p_rel = reloc_null_power(PH, dq[defect], rng, nshuf)
    out["charge_relocnull_lowk"] = float(np.mean(p_def[low] / p_rel[low]))

    # --- 3. defect vertex point process
    r = point_process_sk(PH, defect, V, F2)
    out["vertex_lowk"] = float(np.mean(r[0][low] / r[1][low])) if r else None

    # --- 4b. RIGID-COMPLEX relocation: the arrangement question, charge- and
    #         size-weighted. Internal structure of each complex is preserved.
    if ncomp >= 2:
        wdq = defect_charge_field(dq, defect)
        obs, nul, mc, _ = rigid_complex_sk(PH, wdq, lab, ncomp, rng, nshuf)
        out["rigid_charge_lowk"] = float(np.mean(obs[low] / nul[low]))
        out["rigid_charge_shells"] = shells(kmag, obs / nul, nmax)
        # self-test of the analytic null against actual random translations
        out["rigid_mc_check"] = float(np.mean(mc[low] / nul[low]))

        ind = defect.astype(float)
        obs, nul, _, _ = rigid_complex_sk(PH, ind, lab, ncomp)
        out["rigid_count_lowk"] = float(np.mean(obs[low] / nul[low]))

        Q, neut = complex_charges(dq, lab, ncomp, defect)
        out["Q_mean"] = float(Q.mean())
        out["Q_std"] = float(Q.std())
        out["Q_frac_negative"] = float((Q < 0).mean())
        out["neutrality_frac"] = float(neut.mean())

        # --- 5. real-space charge-charge correlator of complex centroids:
        #        C_QQ(r) = <dQ_i dQ_j>(r) / <dQ^2> (defect_statics'
        #        estimator, on the dq-based complex charges) plus the
        #        running screened fraction
        #        B(r) = -(2 / (N <dQ^2>)) sum_{pairs < r} dQ_i dQ_j,
        #        which climbs toward 1 under perfect screening
        #        (Stillinger-Lovett) and toward 1 - S_Q(0) generally.
        dQ = Q - Q.mean()
        d = cen[:, None, :] - cen[None, :, :]
        d -= np.round(d)                       # min-image, box units
        dist = np.sqrt((d ** 2).sum(-1))
        iu = np.triu_indices(ncomp, 1)
        rr, qq = dist[iu], (dQ[:, None] * dQ[None, :])[iu]
        edges_r = np.linspace(0.0, 0.5, 21)
        cqq, nb = [], []
        for a, b in zip(edges_r[:-1], edges_r[1:]):
            m = (rr >= a) & (rr < b)
            nb.append(int(m.sum()))
            cqq.append(float(qq[m].mean() / (dQ @ dQ / ncomp))
                       if m.sum() else None)
        order = np.argsort(rr)
        brun = -2.0 * np.cumsum(qq[order]) / (ncomp * (dQ @ dQ / ncomp))
        out["cqq_r_edges"] = edges_r.tolist()
        out["cqq_r"] = cqq
        out["cqq_npairs"] = nb
        out["B_of_r"] = [[float(rr[order][i]), float(brun[i])]
                         for i in range(0, len(rr), max(1, len(rr) // 40))]
    else:
        for k in ("rigid_charge_lowk", "rigid_count_lowk", "rigid_mc_check",
                  "Q_mean", "Q_std", "Q_frac_negative", "neutrality_frac"):
            out[k] = None
        out["rigid_charge_shells"] = None

    # --- 4. complex-centroid point process (form-factor free)
    if cen is not None:
        r = points_sk(cen, nvec, V, F2)
        out["centroid_lowk"] = float(np.mean(r[0][low] / r[1][low]))
        out["centroid_shells"] = shells(kmag, r[0] / r[1], nmax)
    else:
        out["centroid_lowk"] = None
        out["centroid_shells"] = None
    return out


def shells(kmag, ratio, nmax):
    edges = np.arange(0.5, nmax + 0.5)
    return [{"k": float((lo + hi) / 2), "n": int(m.sum()),
             "ratio": float(ratio[m].mean())}
            for lo, hi in zip(edges[:-1], edges[1:])
            for m in [(kmag > lo) & (kmag <= hi)] if m.any()]


# ---------------------------------------------------------------------------

def pooled(rows, key):
    v = np.array([r[key] for r in rows if r.get(key) is not None], float)
    if not len(v):
        return None, None
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+", help=".mfd files (need sibling .cocycle.npz)")
    ap.add_argument("--group", default=r"_snap\d+$", help="regex stripped to pool")
    ap.add_argument("--label", default=None,
                    help="pool EVERY input under this one label (overrides --group)")
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--nshuf", type=int, default=64, help="relocation-null draws")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--complex-def", choices=("illegal", "subgraph"),
                    default="illegal",
                    help="what counts as one defect complex. 'subgraph' = all "
                         "edges among defect vertices (the definition the "
                         "2026-08-05 dilute-gas numbers were measured with; "
                         "correct when defects are dilute, but PERCOLATES at "
                         "high defect-vertex density). 'illegal' = components "
                         "of the illegal-edge graph (added 2026-08-06 for "
                         "strain gases; it FRAGMENTS a dilute complex into "
                         "~2.6 pieces, and the pieces of one physical defect "
                         "then read as clustering -- the form-factor trap one "
                         "level up). Pick by density, and say which you used.")
    ap.add_argument("--no-whiten", dest="whiten", action="store_false",
                    help="reproduce the pre-2026-08-13 unfixed GL(3) frame "
                         "(audit only: it mislabels |k| on a non-cubic host, "
                         "so the LOWK window selects the wrong modes)")
    ap.add_argument("--control", action="store_true",
                    help="relocate defects at random first; every ratio must come back 1")
    ap.add_argument("--out", default=None)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    print("[frame] harmonic (periodic Tutte) -- ratios robust, k/amplitudes gauge")
    if args.control:
        print("[CONTROL] defects relocated at random: all ratios must be 1.0\n")

    groups = defaultdict(list)
    for p in args.snapshots:
        if not os.path.exists(os.path.splitext(p)[0] + ".cocycle.npz"):
            print(f"  SKIP {os.path.basename(p)}: no cocycle")
            continue
        r = measure(p, args.nmax, rng, args.nshuf, args.control,
                    whiten=args.whiten, complex_def=args.complex_def)
        label = args.label or re.sub(args.group, "",
                                     os.path.splitext(os.path.basename(p))[0])
        groups[label].append(r)
        print(f"  {r['path']:<26} n_def={r['n_defect']:>4} n_cplx={r['n_complex']:>3} "
              f"<m>={r['mean_size']:>5.2f}  vertex={r['vertex_lowk']:>7.2f} "
              f"centroid={r['centroid_lowk'] if r['centroid_lowk'] is None else round(r['centroid_lowk'],2)}")

    summary = {}
    for label, rows in sorted(groups.items()):
        s = {"n_snapshots": len(rows)}
        for k in ("n_defect", "n_complex", "mean_size", "form_factor",
                  "s2_frac_defect", "charge_permnull_lowk",
                  "charge_relocnull_lowk", "vertex_lowk", "centroid_lowk",
                  "rigid_charge_lowk", "rigid_count_lowk", "rigid_mc_check",
                  "Q_mean", "Q_std", "Q_frac_negative", "neutrality_frac"):
            s[k], s[k + "_sem"] = pooled(rows, k)
        summary[label] = s
        print(f"\n=== {label}  ({len(rows)} snapshots) ===")
        print(f"  defect vertices        {s['n_defect']:.1f}  "
              f"({100*s['n_defect']/rows[0]['V']:.2f}% of V={rows[0]['V']})")
        print(f"  complexes              {s['n_complex']:.1f}   mean size {s['mean_size']:.2f}")
        print(f"  share of sum dq^2      {100*s['s2_frac_defect']:.2f}%  "
              f"<- the crystal carries the rest, and only at Bragg k")
        print(f"  CHARGE  vs perm null   {s['charge_permnull_lowk']:.4f}"
              f" +/- {s['charge_permnull_lowk_sem']:.4f}   (the old 'crystal-grade HU')")
        print(f"  DEFECT charge vs reloc  {s['charge_relocnull_lowk']:.3f}"
              f" +/- {s['charge_relocnull_lowk_sem']:.3f}   (1 = arrangement does nothing)")
        print(f"  DEFECT vertices S(k0)  {s['vertex_lowk']:.2f}"
              f" +/- {s['vertex_lowk_sem']:.2f}   (1 = Poisson)")
        if s["rigid_charge_lowk"] is not None:
            print(f"  --- rigid-complex relocation (whole defects shuffled) ---")
            print(f"  net charge Q_i         {s['Q_mean']:+.3f} +/- {s['Q_std']:.3f}"
                  f"   ({100*s['Q_frac_negative']:.0f}% negative)")
            print(f"  neutrality |Q|/sum|dq| {s['neutrality_frac']:.3f}"
                  f"   (0 = perfect multipole, 1 = one-signed)")
            print(f"  RIGID charge-weighted  {s['rigid_charge_lowk']:.3f}"
                  f" +/- {s['rigid_charge_lowk_sem']:.3f}   (1 = Poisson arrangement)")
            print(f"  RIGID count-weighted   {s['rigid_count_lowk']:.3f}"
                  f" +/- {s['rigid_count_lowk_sem']:.3f}")
            print(f"  [self-test] MC vs analytic null = {s['rigid_mc_check']:.3f}"
                  f" (must be 1.000)")
        if s["centroid_lowk"] is not None:
            print(f"  COMPLEX centroids S(k0){s['centroid_lowk']:.2f}"
                  f" +/- {s['centroid_lowk_sem']:.2f}   (form-factor free)")
            print(f"  consistency: centroid x form_factor = "
                  f"{s['centroid_lowk']*s['form_factor']:.2f}  vs vertex "
                  f"{s['vertex_lowk']:.2f}")

        # pooled C_QQ(r) + screened fraction across the group's snapshots
        cq_rows = [r for r in rows if r.get("cqq_r")]
        if cq_rows:
            edges_r = np.array(cq_rows[0]["cqq_r_edges"])
            mids = 0.5 * (edges_r[1:] + edges_r[:-1])
            print("  C_QQ(r) (box units; <0 = opposite-charge neighbours):")
            for i in range(len(mids)):
                vals = [r["cqq_r"][i] for r in cq_rows
                        if r["cqq_r"][i] is not None]
                npr = sum(r["cqq_npairs"][i] for r in cq_rows)
                if vals and npr > 20 and (i < 6 or i % 3 == 0):
                    print(f"    r={mids[i]:.3f}: C_QQ = {np.mean(vals):+.3f}"
                          + (f" +- {np.std(vals):.3f}" if len(vals) > 1
                             else "") + f"  (pairs {npr})")
            print("  screened fraction B(r) (-> 1 = perfect screening):")
            for rq in (0.1, 0.2, 0.3, 0.45):
                vals = []
                for r in cq_rows:
                    br = np.array(r["B_of_r"])
                    j = np.searchsorted(br[:, 0], rq)
                    if 0 < j <= len(br):
                        vals.append(br[min(j, len(br) - 1), 1])
                if vals:
                    print(f"    B({rq:.2f}) = {np.mean(vals):+.3f}"
                          + (f" +- {np.std(vals):.3f}" if len(vals) > 1
                             else ""))

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"summary": summary,
                       "per_snapshot": {k: v for k, v in groups.items()},
                       "args": vars(args)}, f, indent=1, default=float)
        print(f"\nwrote {args.out}")

    if args.plot:
        make_plot(groups, summary, args.plot, args.control)


def make_plot(groups, summary, path, control):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))

    for label, rows in sorted(groups.items()):
        sh = [r["centroid_shells"] for r in rows if r["centroid_shells"]]
        if not sh:
            continue
        k = [d["k"] for d in sh[0]]
        m = np.array([[d["ratio"] for d in s] for s in sh])
        a1.errorbar(k, m.mean(0), yerr=m.std(0, ddof=1) / np.sqrt(len(m)),
                    marker="o", ms=4, lw=1, capsize=2, label=label)
    a1.axhline(1.0, color="k", ls="--", lw=0.8)
    a1.set_xlabel("|k|  (units of 2π/L)")
    a1.set_ylabel("S(k) / S_random   (complex centroids)")
    a1.set_title("Defect complex arrangement", fontsize=10)
    a1.legend(fontsize=8); a1.grid(alpha=0.3)

    labels = sorted(summary)
    keys = [("charge_permnull_lowk", "curvature charge\nvs permutation null"),
            ("charge_relocnull_lowk", "curvature charge\nvs relocation null"),
            ("vertex_lowk", "defect vertices\n(form-factor loaded)"),
            ("centroid_lowk", "complex centroids\n(form-factor free)")]
    x = np.arange(len(keys))
    w = 0.8 / len(labels)
    for i, lab in enumerate(labels):
        vals = [summary[lab][k] or np.nan for k, _ in keys]
        errs = [summary[lab][k + "_sem"] or 0 for k, _ in keys]
        a2.bar(x + i * w, vals, w, yerr=errs, capsize=3, label=lab)
    a2.axhline(1.0, color="k", ls="--", lw=0.8)
    a2.set_yscale("log")
    a2.set_xticks(x + 0.4 - w / 2)
    a2.set_xticklabels([n for _, n in keys], fontsize=7)
    a2.set_ylabel("low-k ratio  (|n| ≤ 2)")
    a2.set_title("Same snapshots, four estimators", fontsize=10)
    a2.legend(fontsize=8); a2.grid(alpha=0.3, axis="y")

    fig.suptitle(("CONTROL (defects randomised): every bar must be 1  |  "
                  if control else "") +
                 "Defect-density vs curvature-charge hyperuniformity  "
                 "[harmonic frame; PROVISIONAL, uncertified]", fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
