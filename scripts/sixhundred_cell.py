#!/usr/bin/env python3
"""Construct the 600-cell as an equilateral triangulation of S^3 and test the
heat-method geodesic distance on it (the most uniform / round best case:
vertex-deg 20, edge-deg 5, VDV=HDV=0). Vertices = 120 unit icosians; cells =
nearest-neighbor 4-cliques. Ground truth = round-S^3 angular distance arccos<qi,qj>.

Importable: the constructors vertices_600cell()/build_cells()/round_ground_truth()
are used by steiner_geodesic.py for its rigorous-bound validation; the self-test
(heat method vs analytic ground truth) runs only under __main__."""
import os, sys, itertools
import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # sibling: heat_geodesic

phi=(1+5**0.5)/2

def vertices_600cell():
    V=set()
    # 8: (+-1,0,0,0) and permutations
    for i in range(4):
        for s in (1.0,-1.0):
            v=[0.0]*4; v[i]=s; V.add(tuple(v))
    # 16: (+-1/2,+-1/2,+-1/2,+-1/2)
    for s in itertools.product((0.5,-0.5),repeat=4):
        V.add(tuple(s))
    # 96: even perms of (phi,1,1/phi,0)/2 with independent signs
    base=[phi/2, 0.5, 1.0/(2*phi), 0.0]
    def even(p):
        inv=sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        return inv%2==0
    for p in itertools.permutations(range(4)):
        if not even(p): continue
        w=[base[p[k]] for k in range(4)]
        for sg in itertools.product((1,-1),repeat=4):
            V.add(tuple(round(sg[k]*w[k],12)+0.0 for k in range(4)))
    return np.array(sorted(V))

def build_cells(Q):
    n=len(Q)
    D=np.sqrt(((Q[:,None,:]-Q[None,:,:])**2).sum(-1))
    dmin=D[D>1e-6].min()
    A=(np.abs(D-dmin)<1e-6)                       # nearest-neighbor adjacency
    adj=[set(np.where(A[i])[0].tolist()) for i in range(n)]
    cells=set()
    for a in range(n):
        na=sorted(adj[a])
        for b in na:
            if b<=a: continue
            for c in adj[a]&adj[b]:
                if c<=b: continue
                for d in adj[a]&adj[b]&adj[c]:
                    if d<=c: continue
                    cells.add((a,b,c,d))
    return sorted(cells), dmin, A

def round_ground_truth(Q):
    """Round-S^3 angular distance rescaled so nearest neighbor -> 1."""
    G=np.clip(Q@Q.T,-1,1); ang=np.arccos(G)
    theta_min=ang[ang>1e-6].min()
    return ang/theta_min

if __name__=="__main__":
  Q=vertices_600cell()
  cells,dmin,A=build_cells(Q)
  nv=len(Q); ne=int(A.sum()//2); nc=len(cells)
  # f-vector checks
  from collections import Counter
  vdeg=Counter(); edeg=Counter()
  faces=set()
  for (a,b,c,d) in cells:
      for t in itertools.combinations((a,b,c,d),3): faces.add(tuple(sorted(t)))
      for v in (a,b,c,d): vdeg[v]+=1
      for e in itertools.combinations((a,b,c,d),2): edeg[tuple(sorted(e))]+=1
  nf=len(faces)
  print(f"600-cell: V={nv} E={ne} F={nf} cells={nc}  (want 120,720,1200,600)")
  print(f"  Euler V-E+F-C = {nv-ne+nf-nc} (want 0)")
  print(f"  vertex-degree: {set(vdeg.values())} (want {{20}}); edge-degree: {set(edeg.values())} (want {{5}})")
  print(f"  min chord (edge) = {dmin:.6f}, nearest-neighbor angle = {np.degrees(np.arccos(1-dmin**2/2)):.2f} deg")

  # --- run heat method directly on T=cells, V=120 ---
  from heat_geodesic import build_operators, heat_distance
  from scipy.sparse import coo_matrix
  from scipy.sparse.csgraph import shortest_path
  T=np.array(cells); K,Md=build_operators(T,nv)
  print(f"  K row-sum max|.| = {np.abs(np.asarray(K.sum(1)).ravel()).max():.1e}")

  # edge-graph distance
  pr=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
  e=np.unique(np.sort(np.vstack([T[:,[i,j]] for i,j in pr]),axis=1),axis=0)
  Ecsr=coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(nv,nv)).tocsr()

  round_d=round_ground_truth(Q)                     # rescaled: nearest neighbor = 1

  src=[0, 30, 77]
  for m in [0.5,1,2,4]:
      hd=heat_distance(T,K,Md,nv,src,m=m)
      ed=shortest_path(Ecsr,method="D",unweighted=True,directed=False,indices=src)
      rd=round_d[src]
      mask=ed>0
      h=hd[mask]; ee=ed[mask]; rr=rd[mask]
      # neighbor calibration
      nb=[hd[si,np.where(A[s])[0]].mean() for si,s in enumerate(src)]
      # correlation of heat with the round ground truth + slope (heat = slope*round)
      slope=np.polyfit(rr,h,1)[0]
      resid=np.std(h-np.polyval(np.polyfit(rr,h,1),rr))
      print(f"  m={m}: nbr-dist={np.mean(nb):.3f}(want~1)  heat/round slope={slope:.3f} "
            f"scatter={resid:.3f}  heat/edge_med={np.median(h/ee):.3f}  %heat<=edge={100*(h<=ee+1e-6).mean():.0f}")
