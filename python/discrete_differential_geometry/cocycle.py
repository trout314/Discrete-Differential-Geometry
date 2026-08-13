"""Integer 1-cocycles on T^3: construction, harmonic gauge, winding tools.

Companion to the D-side tracking (``sampler.CocycleState``): the sampler
maintains an EXACT integer representative of each of the three generators of
H^1(T^3; Z) across Pachner moves; this module builds the initial assignment
from reference coordinates, computes the HARMONIC representative (the gauge
in which edges acquire direction vectors — the discrete analog of dx, dy, dz
on the flat torus, equivalently the periodic Tutte embedding), and provides
the gauge-invariant integer readouts (loop windings, spanning, class rank).

Conventions: an edge (u, v) with u < v stores omega(u->v) in Z^3; a loop's
winding vector is the signed sum of omega along it, and depends only on the
loop's homology class. Windings are in scaled units: a loop winding once
around direction i pairs to M_i = scale * box_i.
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Construction from reference coordinates
# ---------------------------------------------------------------------------

def build_from_positions(edges: np.ndarray, frac: np.ndarray, box,
                         scale: int = 10**6) -> np.ndarray:
    """Initial integer cocycle from torus coordinates.

    edges: (n, 2) vertex pairs; frac: (n_verts, 3) coordinates with period
    ``box`` (scalar or length-3) per direction. Returns omega (n, 3) int32
    with omega[i] = winding values of edges[i][0] -> edges[i][1]:

        omega(u->v) = p(v) - p(u) - M * k(u, v),

    p = round(scale * frac) (an exact coboundary — closed regardless of
    rounding) and k(u, v) = round((frac(v) - frac(u)) / box) the wrap
    indicator (closed because a triangle's boundary is null-homotopic and
    every edge is short). Exactly closed by construction; edge displacements
    must satisfy |wrapped displacement| < box/3 (asserted).
    """
    frac = np.asarray(frac, dtype=float)
    box = np.broadcast_to(np.asarray(box, dtype=float), (3,))
    edges = np.asarray(edges)
    M = np.round(scale * box).astype(np.int64)

    d = frac[edges[:, 1]] - frac[edges[:, 0]]
    k = np.round(d / box).astype(np.int64)
    wrapped = d - k * box
    if not np.all(np.abs(wrapped) < box / 3):
        raise ValueError("edge displacement >= box/3: closedness of the wrap "
                         "cochain is not guaranteed (increase supercell m)")
    p = np.round(scale * frac).astype(np.int64)
    omega = p[edges[:, 1]] - p[edges[:, 0]] - M * k
    out = omega.astype(np.int32)
    if not np.array_equal(out.astype(np.int64), omega):
        raise OverflowError("cocycle values exceed int32; reduce scale")
    return out


# ---------------------------------------------------------------------------
# Persistence (checkpoint/resume companion to .mfd snapshots)
# ---------------------------------------------------------------------------

def canonicalize_labels(edges: np.ndarray, omega: np.ndarray):
    """Remap vertex labels to their rank in sorted order (0..V-1).

    The .mfd writer renumbers vertices exactly this way, so a cocycle saved
    beside a snapshot must be canonicalized to match what Manifold.load will
    return. The map is monotone, so edge orientation (u < v) and omega signs
    are unchanged. Live sampler labels can have holes (1->4/4->1 churn);
    always canonicalize before persisting.
    """
    edges = np.asarray(edges)
    labels = np.unique(edges)
    lookup = np.full(int(labels.max()) + 1, -1, dtype=np.int64)
    lookup[labels] = np.arange(len(labels))
    return lookup[edges], np.asarray(omega)


def save_cocycle(path: str, edges: np.ndarray, omega: np.ndarray,
                 **meta) -> None:
    """Save a cocycle next to its .mfd snapshot (compressed npz). Extra
    keyword scalars/strings are stored as metadata (e.g. scale=, box=,
    sweeps=). The pair (snapshot.mfd, snapshot.cocycle.npz) is the complete
    resumable state of a tracked chain."""
    np.savez_compressed(path, edges=np.asarray(edges, dtype=np.int32),
                        omega=np.asarray(omega, dtype=np.int32),
                        **{k: np.asarray(v) for k, v in meta.items()})


def load_cocycle(path: str):
    """Load (edges, omega, meta_dict) saved by save_cocycle."""
    z = np.load(path, allow_pickle=False)
    meta = {k: z[k] for k in z.files if k not in ("edges", "omega")}
    return z["edges"], z["omega"], meta


# ---------------------------------------------------------------------------
# Harmonic gauge (graph-Laplacian projection)
# ---------------------------------------------------------------------------

def harmonic_gauge(edges: np.ndarray, omega: np.ndarray, n_verts: int,
                   tol: float = 1e-10):
    """Project to the harmonic representative: omega_h = omega - delta(phi),
    phi = argmin sum_e ||omega - delta(phi)||^2 (one CG solve of the graph
    Laplacian per direction).

    Returns (omega_h (n, 3) float, phi (n_verts, 3) float, gram (3, 3)):
    omega_h are the edge DIRECTION VECTORS (scaled units); gram is the period
    Gram matrix G_ij = <omega_h^i, omega_h^j> — the emergent metric on H^1
    (cell-shape modulus: for a cubic box, G = M^2 * (f1/3-ish) * I up to
    normalization; track its anisotropy, not its trace).

    GAUGE WARNING. The Laplacian is solved per COLUMN, so this projection is
    GL(3)-equivariant on the Z^3 value space: omega -> omega A gives
    omega_h -> omega_h A for any invertible A. The Euclidean structure on that
    space is therefore a free modulus and no length or angle read off omega_h
    means anything until it is fixed -- see :func:`whitening_transform`.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import cg

    edges = np.asarray(edges)
    omega = np.asarray(omega, dtype=float)
    n = len(edges)
    u, v = edges[:, 0], edges[:, 1]

    rows = np.concatenate([u, v, u, v])
    cols = np.concatenate([u, v, v, u])
    vals = np.concatenate([np.ones(n), np.ones(n), -np.ones(n), -np.ones(n)])
    L = coo_matrix((vals, (rows, cols)), shape=(n_verts, n_verts)).tocsr()

    phi = np.zeros((n_verts, 3))
    for i in range(3):
        div = np.zeros(n_verts)
        np.add.at(div, u, omega[:, i])          # div(x) = sum_w omega(x->w)
        np.add.at(div, v, -omega[:, i])
        # stationarity of sum_e (omega - delta phi)^2 gives L phi = -div
        try:
            x, info = cg(L, -div, rtol=tol, atol=0, maxiter=10000)
        except TypeError:   # scipy < 1.12 spells it tol=
            x, info = cg(L, -div, tol=tol, atol=0, maxiter=10000)
        if info != 0:
            raise RuntimeError(f"CG failed to converge (info={info})")
        phi[:, i] = x - x.mean()
    omega_h = omega - (phi[v] - phi[u])
    gram = omega_h.T @ omega_h
    return omega_h, phi, gram


def whitening_transform(gram: np.ndarray, unimodular: bool = True):
    """Fix the Euclidean structure the harmonic gauge leaves free.

    ``harmonic_gauge`` determines the embedding only up to GL(3) acting on the
    value space (see its gauge warning), so the canonical choice is the one
    that makes the edge second moment isotropic: W = G^{-1/2}, after which
    <(omega_h W)^T (omega_h W)> = I. Returns W (3, 3), symmetric positive
    definite; apply on the RIGHT (``omega_h @ W``, ``basis @ W``).

    ``unimodular`` (default) rescales to det(W) = 1, so the transform corrects
    the cell SHAPE without touching the fundamental-domain volume -- the scale
    is meaningful (it is ``scale * box``), the anisotropy is the gauge. With
    that normalization W is invariant under an overall rescaling of gram, so
    it does not matter whether gram was divided by the edge count.

    Validated against known lattices: on a pristine non-cubic crystal, whitened
    harmonic coordinates reproduce the true Cartesian embedding (frac @ L) to
    ~1% in per-tet edge ratios WITHOUT being given the lattice matrix L; on a
    cubic host W is the identity to a few parts in 1e3.
    """
    G = np.asarray(gram, dtype=float)
    G = 0.5 * (G + G.T)
    w, V = np.linalg.eigh(G)
    if w.min() <= 0:
        raise RuntimeError(f"Gram matrix is not positive definite: {w}")
    W = (V * w ** -0.5) @ V.T
    if unimodular:
        W = W * float(np.prod(w)) ** (1.0 / 6.0)
    return W


def min_image(dfrac, basis, shells: int = 1):
    """Shortest Cartesian displacement(s) for FRACTIONAL differences on the
    torus with lattice ``basis`` (rows = lattice vectors).

    Two traps this exists to close, both of which bit consumers of
    ``torus_positions``:

    1. ``lattice_basis`` returns a Euclidean-reduced basis, which is NOT
       diagonal in general (R m4 gives [[4,4,0],[0,4,0],[0,0,4]]). So
       ``frac * np.diag(basis)`` is not the Cartesian position and
       ``d -= round(d / diag) * diag`` is not the min-image rule -- both shear
       every distance. Go through the full matrix: ``dfrac @ basis``.
    2. On a non-orthogonal cell, wrapping to the primary cell
       (``dfrac -= round(dfrac)``) does not always land on the NEAREST image;
       a further lattice translation can be closer. This minimizes over the
       (2*shells+1)^3 neighbouring images, which is exact whenever the true
       separation is below the cell's inradius (ample at shells=1 for the
       hosts in use).

    dfrac: (..., 3) fractional differences. Returns (vec (..., 3) Cartesian,
    dist (...)). Note a per-axis scaling cancels in a centroid computed as
    ``frac[0] + mean(wrap(frac - frac[0]))``, so code shaped that way was
    already right; it is DISTANCES that need this.
    """
    import itertools
    d = np.asarray(dfrac, dtype=float)
    d = d - np.round(d)
    offs = np.array(list(itertools.product(range(-shells, shells + 1),
                                           repeat=3)), dtype=float)
    cand = (d[..., None, :] + offs) @ np.asarray(basis, dtype=float)
    nrm = np.linalg.norm(cand, axis=-1)
    best = np.argmin(nrm, axis=-1)
    vec = np.take_along_axis(cand, best[..., None, None], axis=-2)[..., 0, :]
    return vec, np.take_along_axis(nrm, best[..., None], axis=-1)[..., 0]


# ---------------------------------------------------------------------------
# Gauge-invariant integer readouts
# ---------------------------------------------------------------------------

def tree_positions(edges: np.ndarray, omega: np.ndarray, n_verts: int):
    """Integrate omega over a BFS spanning tree from vertex 0.

    Returns (pos (n_verts, 3) int64 — a lift of each vertex to the cover,
    cycle_windings (m, 3) int64 — winding vectors of the fundamental cycles
    of the non-tree edges, tree_mask (n,) bool). The winding lattice of the
    fundamental cycles generates the full image of the pairing H_1 -> Z^3.
    """
    import collections

    edges = np.asarray(edges)
    omega = np.asarray(omega, dtype=np.int64)
    adj = [[] for _ in range(n_verts)]
    for i, (a, b) in enumerate(edges):
        adj[a].append((b, i, 1))
        adj[b].append((a, i, -1))
    pos = np.zeros((n_verts, 3), dtype=np.int64)
    seen = np.zeros(n_verts, bool)
    tree = np.zeros(len(edges), bool)
    for root in range(n_verts):
        if seen[root]:
            continue
        seen[root] = True
        q = collections.deque([root])
        while q:
            x = q.popleft()
            for (w, ei, sgn) in adj[x]:
                if not seen[w]:
                    seen[w] = True
                    tree[ei] = True
                    pos[w] = pos[x] + sgn * omega[ei]
                    q.append(w)
    nt = ~tree
    cyc = pos[edges[nt, 0]] + omega[nt] - pos[edges[nt, 1]]
    return pos, cyc, tree


def lattice_basis(W: np.ndarray) -> np.ndarray:
    """Triangular integer basis (rows) of the lattice generated by the rows
    of W (m, 3), via Euclidean row reduction. Empty rows dropped; the lattice
    equals scale*box*Z^3 iff the basis is diag(M1, M2, M3) up to signs.
    """
    W = [list(map(int, w)) for w in np.asarray(W, dtype=np.int64) if np.any(w)]
    basis: list[list[int]] = []
    for w in W:
        w = list(w)
        for col in range(3):
            if w[col] == 0:
                continue
            piv = next((b for b in basis
                        if b[col] != 0 and all(x == 0 for x in b[:col])), None)
            if piv is None:
                basis.append([0] * col + w[col:])
                w = [0, 0, 0]
                break
            while w[col] != 0:                  # Euclidean step
                qt = w[col] // piv[col]
                for j in range(3):
                    w[j] -= qt * piv[j]
                if w[col] != 0:
                    piv[:], w[:] = list(w), list(piv)
        basis.sort(key=lambda b: next((i for i, x in enumerate(b) if x), 3))
    return np.array(basis, dtype=np.int64) if basis else np.empty((0, 3), np.int64)


def loop_winding(loop_vertices, edges: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """Winding vector of a closed vertex loop [v0, v1, ..., v0] — the pairing
    with the three cohomology classes; depends only on the loop's homology
    class (in scaled units M per full wrap)."""
    idx = {(min(a, b), max(a, b)): i for i, (a, b) in enumerate(np.asarray(edges))}
    w = np.zeros(3, dtype=np.int64)
    for a, b in zip(loop_vertices[:-1], loop_vertices[1:]):
        i = idx[(min(a, b), max(a, b))]
        w += omega[i] if a < b else -omega[i]
    return w


def torus_positions(facets, edges, omega, whiten: bool = True):
    """Harmonic-gauge (periodic-Tutte) torus coordinates for a snapshot + its
    tracked cocycle -- the one entry point turning (state, cocycle) into a point
    set on the metric torus for structure-factor analysis.

    Composes: canonicalize labels -> integrate the cocycle over a spanning tree
    (lift to the cover) -> read the winding lattice off the fundamental cycles ->
    harmonic gauge X_v = lift_v - phi_v -> fractional coordinates s = X B^{-1}
    -> metric gauge fix B -> B W (``whiten``, default on).

    The whitening is what makes ``basis`` an actual METRIC rather than an
    arbitrary GL(3) frame (:func:`whitening_transform`). It is unimodular and
    acts only on the value space, so ``frac`` is UNCHANGED by it -- consumers
    that use fractional coordinates alone are unaffected, and consumers that
    build lengths or wavevectors out of ``basis`` (``frac @ basis``, k = n
    B^{-T}) get the corrected cell shape. On a cubic host the correction is
    the identity to ~1e-3; on the hexagonal / tetragonal hosts (R, sigma, mu,
    Z, C14, C36, P, delta) it is a real anisotropy of up to ~8x in G's
    eigenvalues, previously carried into every |k| bin. Pass ``whiten=False``
    to recover the pre-2026-08-13 unfixed frame.

    facets: (F, 4) of the .mfd; edges, omega: as from load_cocycle (raw or
    canonical -- canonicalized here). Returns (frac (V, 3) in [0, 1)^3, basis
    (3, 3) winding-lattice, rows = lattice vectors). Raises on a vertex-count
    mismatch or a degenerate winding lattice."""
    facets = np.asarray(facets)
    edges, omega = canonicalize_labels(edges, omega)
    n_verts = int(facets.max()) + 1
    if edges.max() + 1 != n_verts:
        raise RuntimeError(f"vertex count mismatch: facets {n_verts} vs "
                           f"cocycle {edges.max() + 1}")
    pos, cyc, _ = tree_positions(edges, omega, n_verts)
    basis = lattice_basis(cyc).astype(float)
    if basis.shape != (3, 3) or abs(np.linalg.det(basis)) < 1:
        raise RuntimeError(f"winding lattice degenerate: {basis}")
    _, phi, gram = harmonic_gauge(edges, omega, n_verts)
    frac = ((pos - phi) @ np.linalg.inv(basis)) % 1.0
    if whiten:
        basis = basis @ whitening_transform(gram)
    return frac, basis


# ---------------------------------------------------------------------------
# Consumers
# ---------------------------------------------------------------------------

def nematic_q(directions: np.ndarray, weights=None):
    """Nematic order tensor of a set of (unoriented) direction vectors:
    Q = <d^ d^T> - I/3 over unit-normalized rows. Returns (Q, eigenvalues
    ascending, S) with S = 1.5 * lambda_max the uniaxial order parameter
    (S = 0 isotropic, S = 1 perfectly aligned)."""
    d = np.asarray(directions, dtype=float)
    norms = np.linalg.norm(d, axis=1)
    d = d[norms > 0] / norms[norms > 0, None]
    w = np.ones(len(d)) if weights is None else np.asarray(weights)[norms > 0]
    Q = (w[:, None, None] * d[:, :, None] * d[:, None, :]).sum(0) / w.sum()
    Q -= np.eye(3) / 3
    ev = np.linalg.eigvalsh(Q)
    return Q, ev, 1.5 * ev[-1]
