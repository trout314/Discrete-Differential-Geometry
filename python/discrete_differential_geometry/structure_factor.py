"""Real-position structure factor S(k) of a per-vertex field on the metric torus.

Given fractional torus coordinates (from :func:`cocycle.torus_positions`) and a
per-vertex charge (from :mod:`vertex_fields`), evaluate the charge-weighted S(k)
at torus-commensurate wavevectors and normalize by the EXACT charge-permutation
null (charges shuffled over the fixed positions, so the point pattern itself is
factored out). Ratio S_obs/S_null < 1 at small k => hyperuniform, ~1 => Poisson,
> 1 => clustered. Deterministic (analytic null, no shuffling).
"""
import numpy as np


def structure_factor(frac, basis, field, nmax):
    """Structure factor of ``field`` at the positions ``frac``.

    Parameters
    ----------
    frac : (V, 3) fractional torus coordinates in [0, 1)^3.
    basis : (3, 3) winding-lattice basis (rows = lattice vectors).
    field : (V,) per-vertex charge.
    nmax : commensurate reciprocal indices n span [-nmax, nmax]^3 (half-space,
        since the field is real).

    Returns
    -------
    (kmag (M,), s_obs (M,), s_null (M,)) -- physical |k| (units of the smallest
    reciprocal length), observed normalized S(k), and the exact permutation null.
    """
    dq = field - field.mean()
    S2 = float(dq @ dq)
    N = len(field)
    Binv = np.linalg.inv(basis)

    rng = np.arange(-nmax, nmax + 1)
    nvec = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"),
                    axis=-1).reshape(-1, 3)
    nvec = nvec[np.any(nvec != 0, axis=1)]
    keep = ((nvec[:, 0] > 0)
            | ((nvec[:, 0] == 0) & (nvec[:, 1] > 0))
            | ((nvec[:, 0] == 0) & (nvec[:, 1] == 0) & (nvec[:, 2] > 0)))
    nvec = nvec[keep]

    s_obs = np.empty(len(nvec))
    s_null = np.empty(len(nvec))
    for i0 in range(0, len(nvec), 256):
        blk = nvec[i0:i0 + 256]
        phase = np.exp(2j * np.pi * frac @ blk.T)          # (V, blk)
        A = dq @ phase
        F = phase.sum(axis=0)
        s_obs[i0:i0 + len(blk)] = np.abs(A) ** 2 / S2
        s_null[i0:i0 + len(blk)] = 1.0 - (np.abs(F) ** 2 - N) / (N * (N - 1))
    kphys = nvec @ Binv.T
    kmag = np.linalg.norm(kphys, axis=1) / np.abs(Binv).max()
    return kmag, s_obs, s_null
