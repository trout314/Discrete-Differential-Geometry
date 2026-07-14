#!/usr/bin/env python3
"""Steiner-point Dijkstra geodesic distance on an equilateral 3-cone-manifold.

Every tet is a congruent regular Euclidean simplex (edge length 1); the metric is
flat except for cone singularities concentrated ON THE EDGES. A geodesic is a
straight chord inside each tet, gluing across shared faces/edges. So if we sprinkle
Steiner points on the tet boundaries (edges + face interiors), connect every pair of
boundary points OF THE SAME TET by their straight-line distance in that tet's regular
embedding, and run Dijkstra, any shortest graph path is a REAL piecewise-straight path
in the manifold ==> the graph distance is a rigorous UPPER BOUND on the true PL
geodesic. Refining the Steiner grid (raising subdivision order n) makes the bound
converge DOWN to the true geodesic. (Contrast the heat method in heat_geodesic.py:
no guarantee, and it silently violates even the crude edge bound on disordered seeds.)

Subdivision order n (single knob):
  n=1  -> no Steiner points; each tet's 4 vertices are mutually distance 1.
          The graph is exactly the unit-weight 1-skeleton  ==>  recovers EDGE distance.
  n>=2 -> (n-1) interior points per edge (fractions k/n) SHARED by every tet on the edge,
          plus (n-1)(n-2)/2 interior barycentric points per face SHARED by the 2 tets
          on the face. Sharing is what lets a path cross tet boundaries.

Gluing correctness: an edge point is keyed by its fraction from the edge's global-min
endpoint; a face point by barycentric weights over the face's sorted global vertices --
so the same physical point gets one global id regardless of which tet places it, and
its two intra-tet embeddings coincide (faces/edges are congruent equilateral pieces).

Usage:
  python scripts/steiner_geodesic.py 600 --orders 1,2,3,4,6,8      # 600-cell validation
  python scripts/steiner_geodesic.py seeds/<seed>.mfd --orders 1,2,3,4 --sources 12
"""
import os, sys, glob, itertools, argparse
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import shortest_path

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)                          # siblings: heat_geodesic, sixhundred_cell
sys.path.insert(0, os.path.join(_ROOT, "python"))
from heat_geodesic import P, load_mesh             # P = reference regular tet (edge len 1), 4x3

EDGE_PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
FACE_TRIS =[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]

def face_interior_bary(n):
    """Interior barycentric weights (i,j,k), i+j+k=n, all >=1. Count (n-1)(n-2)/2."""
    return [(i,j,n-i-j) for i in range(1,n) for j in range(1,n-i) if n-i-j>=1]

def _enumerate(T):
    """Assign a global id-block to every unique edge and face of the tet mesh."""
    E={}; Fd={}
    for tet in T:
        for i,j in EDGE_PAIRS:
            e=(tet[i],tet[j]) if tet[i]<tet[j] else (tet[j],tet[i])
            if e not in E: E[e]=len(E)
        for i,j,k in FACE_TRIS:
            f=tuple(sorted((tet[i],tet[j],tet[k])))
            if f not in Fd: Fd[f]=len(Fd)
    return E,Fd

def build_steiner_graph(T, V, n):
    """Return (csr graph over Ntot nodes, Ntot). Nodes 0..V-1 are the mesh vertices."""
    T=np.asarray(T,dtype=np.int64)
    E,Fd=_enumerate(T)
    nE,nF=len(E),len(Fd)
    fbary=face_interior_bary(n); nfi=len(fbary)
    base_edge=V
    base_face=V + nE*(n-1)
    Ntot=base_face + nF*nfi

    II=[]; JJ=[]; WW=[]
    tri=np.triu_indices(4 + 6*(n-1) + 4*nfi, k=1)   # upper-tri template (fixed per-tet size)
    for tet in T:
        ids=[]; pos=[]
        # 4 vertices
        for l in range(4):
            ids.append(int(tet[l])); pos.append(P[l])
        # 6 edges: (n-1) interior points each, fraction from global-min endpoint
        for li,lj in EDGE_PAIRS:
            a,b=tet[li],tet[lj]
            e=(a,b) if a<b else (b,a)
            blk=base_edge + E[e]*(n-1)
            pa,pb=(P[li],P[lj]) if a<b else (P[lj],P[li])   # pa at the global-min end
            for kk in range(1,n):
                f=kk/n
                ids.append(blk+(kk-1)); pos.append((1-f)*pa+f*pb)
        # 4 faces: interior barycentric points over the face's sorted global vertices
        for li,lj,lk in FACE_TRIS:
            loc=[li,lj,lk]; gl=[tet[li],tet[lj],tet[lk]]
            order=np.argsort(gl)                    # local indices sorted by global id
            p0,p1,p2=P[loc[order[0]]],P[loc[order[1]]],P[loc[order[2]]]
            fkey=tuple(int(gl[o]) for o in order)
            blk=base_face + Fd[fkey]*nfi
            for idx,(i,j,k) in enumerate(fbary):
                ids.append(blk+idx); pos.append((i*p0+j*p1+k*p2)/n)
        ids=np.asarray(ids); pos=np.asarray(pos)
        d=np.linalg.norm(pos[tri[0]]-pos[tri[1]],axis=1)
        II.append(ids[tri[0]]); JJ.append(ids[tri[1]]); WW.append(d)

    ii=np.concatenate(II); jj=np.concatenate(JJ); ww=np.concatenate(WW)
    # order endpoints so i<j, then keep the MIN weight per undirected pair
    # (coo_matrix SUMS duplicates -- shared boundary points recur across tets -- so we
    #  must reduce by min first or the glued weights get inflated).
    lo=np.minimum(ii,jj); hi=np.maximum(ii,jj)
    key=lo.astype(np.int64)*Ntot + hi
    o=np.argsort(key); key=key[o]; ww=ww[o]; lo=lo[o]; hi=hi[o]
    first=np.ones(len(key),bool); first[1:]=key[1:]!=key[:-1]
    starts=np.flatnonzero(first)
    wmin=np.minimum.reduceat(ww,starts)
    lo=lo[starts]; hi=hi[starts]
    # symmetric csr
    r=np.concatenate([lo,hi]); c=np.concatenate([hi,lo]); w=np.concatenate([wmin,wmin])
    G=coo_matrix((w,(r,c)),shape=(Ntot,Ntot)).tocsr()
    return G, Ntot

def steiner_distance(T, V, n, sources):
    """(nsrc, V) upper-bound geodesic distances from `sources` to all mesh vertices."""
    G,Ntot=build_steiner_graph(T,V,n)
    d=shortest_path(G,method="D",directed=False,indices=list(sources))
    return d[:, :V]

# --------------------------------------------------------------------------------------
def _run_600cell(orders):
    from sixhundred_cell import vertices_600cell, build_cells, round_ground_truth
    Q=vertices_600cell(); cells,dmin,A=build_cells(Q)
    V=len(Q); T=np.array(cells)
    round_d=round_ground_truth(Q)
    src=[0,30,77]
    ed=steiner_distance(T,V,1,src)                  # n=1 == edge distance (baseline)
    mask=ed>0
    rd=np.vstack([round_d[s] for s in src])
    print(f"600-cell: V={V} tets={len(T)}   (Steiner upper bound; n=1 == edge distance)")
    print(f"  ground truth = round-S^3 angular distance (nearest neighbor rescaled to 1)")
    print(f"  {'n':>3} {'edgepts':>7} {'facepts':>7} {'mean_d':>8} {'/edge':>7} "
          f"{'/round':>8} {'slope_v_round':>13} {'scatter':>8} {'%<=edge':>8} {'%<=prev':>8}")
    prev=None; ed_v=ed[mask]
    for n in orders:
        sd=steiner_distance(T,V,n,src)
        s=sd[mask]; rr=rd[mask]
        slope=np.polyfit(rr,s,1)[0]
        scatter=np.std(s-np.polyval(np.polyfit(rr,s,1),rr))
        le=100*(s<=ed_v+1e-9).mean()
        lp=100*(s<=prev+1e-9).mean() if prev is not None else 100.0
        ne=(n-1)                                    # interior pts per edge
        nf=(n-1)*(n-2)//2                           # interior pts per face
        print(f"  {n:>3} {ne:>7} {nf:>7} {s.mean():>8.3f} {np.median(s/ed_v):>7.3f} "
              f"{np.median(s/rr):>8.3f} {slope:>13.3f} {scatter:>8.4f} {le:>8.1f} {lp:>8.1f}")
        prev=s
    print("\n  Expect: monotone decrease in the mean (converging toward the heat-method /"
          " round value, heat/edge_med ~0.80); %<=edge=100 at every order (rigorous bound).")

def _run_seed(seed, orders, nsrc):
    T,V=load_mesh(seed)
    rng=np.random.default_rng(0); src=rng.choice(V,min(nsrc,V),replace=False).tolist()
    ed=steiner_distance(T,V,1,src); mask=(ed>0)&np.isfinite(ed); ed_v=ed[mask]
    print(f"seed {os.path.basename(seed)}: V={V} tets={len(T)}  sources={len(src)}")
    print(f"  edge dist (n=1): mean {ed_v.mean():.2f} max {ed_v.max():.0f}")
    print(f"  {'n':>3} {'mean_d':>8} {'/edge_med':>10} {'%<=edge':>8} {'%<=prev':>8}")
    prev=None
    for n in orders:
        sd=steiner_distance(T,V,n,src); s=sd[mask]
        le=100*(s<=ed_v+1e-9).mean()
        lp=100*(s<=prev+1e-9).mean() if prev is not None else 100.0
        print(f"  {n:>3} {s.mean():>8.3f} {np.median(s/ed_v):>10.3f} {le:>8.1f} {lp:>8.1f}")
        prev=s
    print("  (Steiner is a rigorous upper bound: %<=edge must be 100 at every order.)")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("target",nargs="?",default="600",
                    help="'600' for the 600-cell validation, or a seed .mfd path")
    ap.add_argument("--orders",default="1,2,3,4,6",help="comma list of subdivision orders n")
    ap.add_argument("--sources",type=int,default=12,help="random sources (seed mode)")
    args=ap.parse_args()
    orders=[int(x) for x in args.orders.split(",")]
    if args.target=="600":
        _run_600cell(orders)
    else:
        seed=args.target
        if not os.path.isabs(seed):
            seed=os.path.join(_ROOT,seed)
        _run_seed(seed,orders,args.sources)
