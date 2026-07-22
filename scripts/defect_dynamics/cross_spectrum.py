#!/usr/bin/env python3
"""Screening-mechanism decomposition: who cancels the carrier monopole power?

Partition each snapshot's mean-subtracted charge field dq by vertex class:
  carrier = illegal vertices (complex members)
  halo    = legal but registry-defect (defect_cores minus carriers)
  bulk    = everything else
so dq = dq_c + dq_h + dq_b exactly. Sub-Bragg:
  S_tot = S_cc + S_hh + S_bb + 2(S_ch + S_cb + S_hb)
Negative cross terms = the medium's screening cloud, measured directly.
S_cc vs the point-monopole floor tests the form-factor alternative.

Usage: cross_spectrum.py PREFIX...   (each PREFIX has _snap*.mfd + .cocycle.npz)
"""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from dopant_pairs import vertex_classes
from defect_census import defect_cores
from tcp_melt import CRYSTALS

NMAX, KCUT = 8, 5.0
REF_FACETS = np.asarray(ddg.Manifold.load(CRYSTALS["r"], 3).facets())

PARTS = ["cc", "hh", "bb", "ch", "cb", "hb"]


def kgrid():
    r = np.arange(-NMAX, NMAX + 1)
    nv = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    nv = nv[np.any(nv != 0, 1)]
    keep = ((nv[:, 0] > 0) | ((nv[:, 0] == 0) & (nv[:, 1] > 0))
            | ((nv[:, 0] == 0) & (nv[:, 1] == 0) & (nv[:, 2] > 0)))
    nv = nv[keep]
    km = np.linalg.norm(nv, axis=1)
    m = km <= KCUT
    return nv[m], km[m]


NV, KM = kgrid()
SHELLS = [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0)]

acc = defaultdict(list)          # (part, shell) -> values;  also totals
floors = []
for pre in sys.argv[1:]:
    for snap in sorted(glob.glob(pre + "_snap*.mfd")):
        fac = np.asarray(ddg.Manifold.load(snap, 3).facets())
        qR = FIELDS["curvature_charge"](fac)
        dq = qR - qR.mean()
        S2 = float(dq @ dq)
        n6, imp, adj = vertex_classes(fac)
        V = len(n6)
        dv, _, _, _, _ = defect_cores(fac, REF_FACETS, "r", 30)
        carrier = imp > 0
        halo = np.zeros(V, bool)
        halo[list(dv)] = True
        halo &= ~carrier
        bulk = ~carrier & ~halo
        edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
        frac, basis = coc.torus_positions(fac, edges, omega)
        ph = np.exp(2j * np.pi * (frac @ NV.T))
        A = {"c": (dq * carrier) @ ph, "h": (dq * halo) @ ph,
             "b": (dq * bulk) @ ph}
        # point-monopole floor from carrier complexes (for S_cc comparison)
        illv = np.nonzero(carrier)[0]
        seen, Q2 = set(), 0.0
        for s0 in illv:
            if s0 in seen:
                continue
            st, comp = [int(s0)], []
            seen.add(int(s0))
            while st:
                u = st.pop(); comp.append(u)
                for w in adj[u]:
                    if imp[w] > 0 and w not in seen:
                        seen.add(int(w)); st.append(int(w))
            Q2 += float(qR[comp].sum() - len(comp) * qR.mean()) ** 2
        floors.append(Q2 / S2)
        for lo, hi in SHELLS + [(1.0, KCUT)]:
            m = (KM >= lo) & (KM < hi)
            tot = 0.0
            for part in PARTS:
                a, b = part
                v = float(np.mean((A[a][m] * np.conj(A[b][m])).real)) / S2
                w = v if a == b else 2 * v
                acc[(part, (lo, hi))].append(w)
                tot += w
            acc[("tot", (lo, hi))].append(tot)

print(f"snapshots: {len(floors)}   point-monopole floor <sum Q^2/S2> = "
      f"{np.mean(floors):.4f}")
print(f"\n{'shell':>9s} " + " ".join(f"{p:>8s}" for p in PARTS)
      + f" {'total':>8s}")
for sh in SHELLS + [(1.0, KCUT)]:
    lab = f"{sh[0]:.0f}-{sh[1]:.0f}"
    print(f"{lab:>9s} " + " ".join(
        f"{np.mean(acc[(p, sh)]):8.4f}" for p in PARTS)
        + f" {np.mean(acc[('tot', sh)]):8.4f}")
print("\n(cc=carrier auto, hh=halo, bb=bulk; ch/cb/hb = 2x cross terms; "
      "negative cross = screening cloud)")
