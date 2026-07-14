#!/usr/bin/env python3
"""Heat-method geodesic distance on an equilateral 3-manifold triangulation
(Crane-Weischedel-Wardetzky 2013), prototype.

Every tet is a congruent regular simplex, so the P1 FEM operators are UNIVERSAL:
one reference regular tet gives the local stiffness K_loc (4x4) and the 4 constant
basis gradients used for all tets. Everything is intrinsic -- each tet is worked in
its own local frame; the divergence sums frame-independent scalars (grad_i . X_tet).

Steps (source s):
  1. heat:   (M + t K) u = 1_s          t = m*h^2, h=edge len=1, K=stiffness(=-cotanLap)
  2. field:  X_tet = -grad(u)/|grad(u)|  per tet
  3. div:    b_i = sum_tets V (grad_i . X_tet)
  4. Poisson:(K + eps I) phi = b ; phi -= phi[s]   -> geodesic distance from s
Compare to edge-graph BFS distance (a PL upper bound) on the same source/targets.

NOTE: the heat method carries NO error guarantee. Validated essentially exact on the
uniform 600-cell (see sixhundred_cell.py) but it OVERSHOOTS and violates the crude edge
upper bound on disordered / large-diameter seeds. Use steiner_geodesic.py when a
rigorous bound is required.
"""
import os, sys, glob
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, identity
from scipy.sparse.linalg import factorized
from scipy.sparse.csgraph import shortest_path

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold

# ---- reference regular tetrahedron (edge length 1) ----
P = np.array([[0,0,0],[1,0,0],[0.5,np.sqrt(3)/2,0],
              [0.5,np.sqrt(3)/6,np.sqrt(6)/3]])
Em = (P[1:]-P[0]).T                      # 3x3, columns p1-p0,p2-p0,p3-p0
VOL = abs(np.linalg.det(Em))/6.0         # = sqrt(2)/12
Einv = np.linalg.inv(Em)
GRAD = np.zeros((4,3))
GRAD[1:] = Einv                          # grad(lambda_1,2,3) = rows of Em^{-1}
GRAD[0] = -GRAD[1:].sum(0)               # grad(lambda_0) = -sum
KLOC = VOL * (GRAD @ GRAD.T)             # 4x4 universal local stiffness
# sanity: |grad|^2=3/2, off=-1/2 ; KLOC diag=VOL*1.5, off=VOL*-0.5
assert np.allclose(np.diag(GRAD@GRAD.T), 1.5) and np.allclose((GRAD@GRAD.T)[0,1], -0.5)

def load_mesh(path):
    m=Manifold.load(path,3)
    F=np.asarray(m.facets(),np.int64)
    labels,inv=np.unique(F,return_inverse=True)
    return inv.reshape(F.shape), len(labels)   # tets (nt,4) in 0..V-1, V

def build_operators(T, V):
    nt=len(T)
    # stiffness K (V x V): add KLOC over each tet's 4 global vertices
    rows=np.repeat(T,4,axis=1).reshape(nt,4,4)          # a-index
    cols=np.tile(T,(1,4)).reshape(nt,4,4)               # b-index
    data=np.broadcast_to(KLOC,(nt,4,4))
    K=coo_matrix((data.ravel(),(rows.ravel(),cols.ravel())),shape=(V,V)).tocsr()
    # lumped mass (diagonal): each tet gives VOL/4 to each of its verts
    Mdiag=np.bincount(T.ravel(),weights=np.full(T.size,VOL/4.0),minlength=V)
    return K, Mdiag

def edge_csr(T, V):
    pr=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    e=np.unique(np.sort(np.vstack([T[:,[i,j]] for i,j in pr]),axis=1),axis=0)
    d=np.ones(len(e))
    return coo_matrix((np.r_[d,d],(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),
                      shape=(V,V)).tocsr()

def heat_distance(T, K, Mdiag, V, sources, m=1.0, eps=1e-8):
    """Return (nsrc, V) geodesic distances via the heat method."""
    t=m*1.0**2                                   # h=1 (equilateral)
    A=(K + coo_matrix((Mdiag/t,(range(V),range(V))),shape=(V,V)).tocsr())  # (K + M/t)
    solveA=factorized(A.tocsc())
    solveK=factorized((K+eps*identity(V,format='csr')).tocsc())
    Tg=T                                          # (nt,4)
    out=np.zeros((len(sources),V))
    for si,s in enumerate(sources):
        rhs=np.zeros(V); rhs[s]=1.0/t*Mdiag[s]    # (M/t) delta_s
        u=solveA(rhs)
        # per-tet gradient of u: sum_a u[g_a] GRAD[a]
        gu=np.einsum('ta,ak->tk', u[Tg], GRAD)    # (nt,3)
        nrm=np.linalg.norm(gu,axis=1,keepdims=True); nrm[nrm<1e-14]=1
        X=-gu/nrm                                  # (nt,3) unit, toward source
        # divergence at vertices: b_i = sum_tets VOL*(GRAD[a].X_tet)
        contrib=VOL*np.einsum('ak,tk->ta', GRAD, X)   # (nt,4) per local vertex
        b=np.bincount(Tg.ravel(), weights=contrib.ravel(), minlength=V)
        phi=solveK(b)
        phi-=phi[s]
        # heat method yields distance up to sign; enforce phi(source)=0 & phi>=0-ish
        if np.median(phi)<0: phi=-phi; phi-=phi[s]
        out[si]=phi
    return out

if __name__=="__main__":
    seed=sys.argv[1] if len(sys.argv)>1 else \
        "seeds/S3_N1e4_1e-1_ED5p0043_1_VDVs_2e-3_HDVs_0p032_s000.mfd"
    nsrc=int(sys.argv[2]) if len(sys.argv)>2 else 20
    T,V=load_mesh(os.path.join(_ROOT,seed))
    print(f"seed {os.path.basename(seed)}: V={V} vertices, {len(T)} tets")
    K,Md=build_operators(T,V)
    # operator sanity
    rowsum=np.abs(np.asarray(K.sum(1)).ravel()).max()
    print(f"  K row-sum max |.|={rowsum:.2e} (should be ~0: constant in nullspace)")
    Ecsr=edge_csr(T,V)
    rng=np.random.default_rng(0); src=rng.choice(V,min(nsrc,V),replace=False)
    ed=shortest_path(Ecsr,method="D",unweighted=True,directed=False,indices=src)
    mask=(ed>0)&np.isfinite(ed)
    print(f"\n  edge dist: mean {ed[mask].mean():.2f}  max {ed[mask].max():.0f}")
    print(f"  sweep diffusion time t = m*h^2  (h=1):")
    print(f"  {'m':>5} {'heat_mean':>9} {'heat_max':>8} {'ratio_med':>9} {'ratio5-95':>16} {'%heat<=edge':>11}")
    best=None
    for m in [1,4,16,64,256]:
        hd=heat_distance(T,K,Md,V,src,m=m)
        h=hd[mask]; e=ed[mask]; ratio=h/e
        ub=(h<=e+1e-6).mean()*100
        print(f"  {m:>5} {h.mean():>9.2f} {h.max():>8.2f} {np.median(ratio):>9.3f} "
              f"[{np.percentile(ratio,5):>5.3f},{np.percentile(ratio,95):>5.3f}] {ub:>11.1f}")
        if best is None or ub>best[1]: best=(m,ub,ratio)
    print(f"\n  best upper-bound respect at m={best[0]} ({best[1]:.0f}%); "
          f"median ratio {np.median(best[2]):.3f}")
    h=heat_distance(T,K,Md,V,src,m=best[0])[mask]; e=ed[mask]; ratio=h/e
    # scatter
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        SP=os.path.dirname(os.path.abspath(__file__))
        fig,ax=plt.subplots(1,2,figsize=(12,5))
        ax[0].hexbin(e,h,gridsize=40,mincnt=1,cmap="viridis")
        mx=e.max(); ax[0].plot([0,mx],[0,mx],"r--",label="heat=edge")
        ax[0].set_xlabel("edge-graph distance"); ax[0].set_ylabel("heat geodesic distance")
        ax[0].set_title("heat vs edge distance"); ax[0].legend()
        ax[1].hist(ratio,bins=60,color="#2471a3"); ax[1].axvline(np.median(ratio),color="r")
        ax[1].set_xlabel("heat / edge"); ax[1].set_title(f"ratio (median {np.median(ratio):.3f})")
        fig.tight_layout(); fig.savefig(os.path.join(SP,"heat_vs_edge.png"),dpi=130)
        print("  wrote heat_vs_edge.png")
    except Exception as ex:
        print("  (plot skipped:",ex,")")
