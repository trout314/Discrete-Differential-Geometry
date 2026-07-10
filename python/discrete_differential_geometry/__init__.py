"""Discrete Differential Geometry — Python bindings."""

from ._simplicial_complex import SimplicialComplex, join
from ._manifold import Manifold
from ._manifold_view import ManifoldView
from ._sampler import ManifoldSampler, SamplerParams, SamplerStats
from .convergence import (
    split_rhat,
    rank_normalized_rhat,
    quantized_split_rhat,
    effective_sample_size,
    integrated_autocorrelation_time,
    weighted_ess,
)
from ._dlang import _lib as _lib


def gc_collect():
    """Trigger a D garbage collection cycle to reclaim temporary allocations."""
    _lib.ddg_gc_collect()


def gc_minimize():
    """Minimize the D GC heap, returning free pages to the OS."""
    _lib.ddg_gc_minimize()


def gc_stats():
    """Return (used_bytes, free_bytes) of the D GC heap.

    `used` is live (reachable) memory; `free` is reclaimable pool space that
    gc_minimize() can return to the OS. Live memory that climbs while the
    manifold's f-vector stays flat indicates a leak (see project history).
    """
    import ctypes
    used, free = ctypes.c_long(0), ctypes.c_long(0)
    _lib.ddg_gc_stats(ctypes.byref(used), ctypes.byref(free))
    return used.value, free.value


__all__ = [
    "SimplicialComplex",
    "Manifold",
    "ManifoldView",
    "ManifoldSampler",
    "SamplerParams",
    "SamplerStats",
    "join",
    "gc_collect",
    "gc_minimize",
    "gc_stats",
    "split_rhat",
    "rank_normalized_rhat",
    "quantized_split_rhat",
    "effective_sample_size",
    "integrated_autocorrelation_time",
    "weighted_ess",
]
