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
