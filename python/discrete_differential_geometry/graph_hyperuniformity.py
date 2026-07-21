"""Graph-proxy hyperuniformity estimators for a per-vertex field -- no
coordinates, topology only.

Two estimators, each observed / charge-permutation-shuffle (the same null idea as
:mod:`structure_factor`, here evaluated by Monte-Carlo shuffles rather than
analytically):

* window variance -- Var over random centers of the charge summed in BFS balls
  grown to exactly m vertices, vs m. Hyperuniform <=> ratio decays with m;
  ~1 = Poisson; > 1 = clustered.
* spectral       -- graph-Laplacian mode power |phi_n . dq|^2 (lambda ~ k^2).
  Hyperuniform <=> ratio -> 0 as lambda -> 0. :func:`spectral_ratio` uses an
  eigensolver; :func:`lowpass_ratio` is a matrix-free surrogate for highly
  symmetric states (perfect crystals), whose degenerate Laplacian spectra make
  ARPACK shift-invert thrash.

These are STOCHASTIC (shuffle nulls / shuffled BFS shells): pass an explicit
``rng`` for reproducibility; results match the pre-refactor estimators up to
Monte-Carlo noise, not bit-for-bit.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def adjacency(eu, V):
    """Neighbour-list adjacency (list of int arrays) from unique edges eu (E,2)."""
    adj = [[] for _ in range(V)]
    for a, b in eu:
        adj[a].append(b)
        adj[b].append(a)
    return [np.array(x) for x in adj]


def bfs_order(adj, V, src, mmax, rng):
    """Vertices in BFS order from src, each shell shuffled (rng), length mmax."""
    seen = np.zeros(V, bool)
    seen[src] = True
    order = [src]
    frontier = [src]
    while len(order) < mmax and frontier:
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if not seen[w]:
                    seen[w] = True
                    nxt.append(w)
        rng.shuffle(nxt)
        order.extend(nxt)
        frontier = nxt
    return np.array(order[:mmax])


def window_variances(q, orders, mgrid):
    """Var over centers of the charge summed in each BFS window of size m in
    mgrid (orders: list of BFS-orders, one per center)."""
    sums = np.array([np.cumsum(q[o])[mgrid - 1] for o in orders])
    return sums.var(axis=0)


def window_ratio(q, orders, mgrid, nshuf, rng):
    """Window-variance ratio Var(q) / mean_shuffle Var(perm q), one value per m."""
    vr = window_variances(q, orders, mgrid)
    vs = np.mean([window_variances(rng.permutation(q), orders, mgrid)
                  for _ in range(nshuf)], axis=0)
    return vr / vs


def _laplacian(eu, V):
    ii = np.r_[eu[:, 0], eu[:, 1]]
    jj = np.r_[eu[:, 1], eu[:, 0]]
    A = coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(V, V)).tocsr()
    return A, A.sum(1).A1


def spectral_ratio(eu, V, dq, k=41, nshuf=32, rng=None):
    """(lambda, S_obs/S_shuffle) over the k smallest nonzero Laplacian modes.
    dq must be mean-subtracted. Uses ARPACK shift-invert; see lowpass_ratio for
    crystals."""
    rng = np.random.default_rng(0) if rng is None else rng
    A, d = _laplacian(eu, V)
    L = coo_matrix((d, (range(V), range(V))), shape=(V, V)).tocsr() - A
    lam, phi = eigsh(L, k=min(k, V - 2), sigma=0, which="LM")
    o = np.argsort(lam)
    lam, phi = lam[o][1:], phi[:, o][:, 1:]            # drop the constant mode
    Sr = (phi.T @ dq) ** 2
    Ss = np.mean([(phi.T @ rng.permutation(dq)) ** 2 for _ in range(nshuf)], axis=0)
    return lam, Sr / Ss


def lowpass_ratio(eu, V, dq, n_iter=3000, n_shuf=16, rng=None):
    """Matrix-free low-k spectral power ratio: power of dq surviving n_iter steps
    of the low-pass filter y -> (I - L/lmax) y (effective cutoff lambda ~
    lmax/n_iter), observed / shuffled. Use INSTEAD of spectral_ratio on highly
    symmetric states (perfect crystals) whose degenerate spectra thrash ARPACK;
    exact crystal answer is 0 (a unit-cell-periodic field has no sub-BZ power)."""
    rng = np.random.default_rng(0) if rng is None else rng
    A, d = _laplacian(eu, V)
    lmax = 2 * d.max()

    def low_power(x):
        y = x - x.mean()
        for _ in range(n_iter):
            y = y - (d * y - A @ y) / lmax
            y = y - y.mean()
        return float(y @ y)

    obs = low_power(dq)
    shuf = np.mean([low_power(rng.permutation(dq)) for _ in range(n_shuf)])
    return obs / shuf
