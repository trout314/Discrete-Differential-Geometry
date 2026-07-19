"""Disclination-network censuses and the six-edge flip stream (dim=3).

Python mirror of the disclination-network observables in ``source/sampler.d``:
the network is the graph of degree>=6 edges ("six-edges"); legal Z14/Z15/Z16
vertices carry exactly 2/3/4 six-edge ends, so this graph is the Frank-Kasper
skeleton and its components/segments/loops are the physical disclination
lines. See ``DisclinationCensus`` in the D source for field semantics.
"""
from __future__ import annotations

import ctypes

import numpy as np

from . import _dlang

_lib = _dlang._lib

#: One record per edge crossing the degree 5<->6 threshold in an accepted
#: move; dir = +1 (six-edge born) / -1 (six-edge died). Cumulative sum of dir
#: tracks E6 exactly; (u, v) is the sorted edge label; clock is the shared
#: geometry-ledger clock (attempted moves), so records correlate with the
#: move event log.
SIX_FLIP_DTYPE = np.dtype([("clock", "<u8"), ("u", "<i4"), ("v", "<i4"),
                           ("dir", "<i4")], align=False)
assert SIX_FLIP_DTYPE.itemsize == 20

#: Scalar slots 0..20 of the flattened census, in D declaration order.
CENSUS_SCALARS = [
    "n_net_verts", "n_six_edges", "n_components", "giant_size", "second_size",
    "giant_diameter", "cycle_rank", "n_segments", "sum_seg_len",
    "n_pure_loops", "sum_loop_len", "n_endpoints", "n_fray_verts",
    "n_imp_end_edges", "n_z14", "n_z15", "n_z16",
    "e_dop_dop", "e_dop_host", "e_host_host", "e_imp_any",
]
CENSUS_SLOTS = 24 + 8 + 64 + 32


def host_mask(host_classes) -> int:
    """Bitmask from an iterable of native host n6 classes (C15: [0, 4])."""
    mask = 0
    for k in host_classes or ():
        if not 0 <= int(k) <= 4:
            raise ValueError(f"host class must be an n6 value in 0..4, got {k}")
        mask |= 1 << int(k)
    return mask


def valence_census_from_handle(handle, n6_cap: int, m_cap: int) -> np.ndarray:
    """Joint (n6, m) vertex census as an (n6_cap+1, m_cap+1) int64 array."""
    n = (n6_cap + 1) * (m_cap + 1)
    buf = (ctypes.c_long * n)()
    _lib.ddg_manifold_valence_census(handle, buf, n6_cap, m_cap)
    return np.array(buf[:n], dtype=np.int64).reshape(n6_cap + 1, m_cap + 1)


def disclination_census_from_handle(handle, host_classes=None) -> dict:
    """Disclination-network census as a dict (see CENSUS_SCALARS plus
    ``net_deg_census`` [degree 1..7, 8+], ``seg_len_hist`` [len clamped at 63],
    ``comp_size_hist`` [log2 bins], and the derived ``mean_seg_len``,
    ``giant_frac``)."""
    buf = (ctypes.c_long * CENSUS_SLOTS)()
    _lib.ddg_manifold_disclination_census(
        handle, host_mask(host_classes), buf, CENSUS_SLOTS)
    flat = np.array(buf[:CENSUS_SLOTS], dtype=np.int64)
    out = {name: int(flat[i]) for i, name in enumerate(CENSUS_SCALARS)}
    out["net_deg_census"] = flat[24:32].copy()
    out["seg_len_hist"] = flat[32:96].copy()
    out["comp_size_hist"] = flat[96:128].copy()
    out["mean_seg_len"] = (out["sum_seg_len"] / out["n_segments"]
                           if out["n_segments"] else 0.0)
    out["giant_frac"] = (out["giant_size"] / out["n_net_verts"]
                         if out["n_net_verts"] else 0.0)
    return out
