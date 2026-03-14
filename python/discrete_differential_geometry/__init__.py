"""Discrete Differential Geometry — Python bindings."""

from ._simplicial_complex import SimplicialComplex, join
from ._manifold import Manifold
from ._sampler import ManifoldSampler, SamplerParams
from ._dlang import _lib as _lib


def gc_collect():
    """Trigger a D garbage collection cycle to reclaim temporary allocations."""
    _lib.ddg_gc_collect()


def gc_minimize():
    """Minimize the D GC heap, returning free pages to the OS."""
    _lib.ddg_gc_minimize()


__all__ = [
    "SimplicialComplex",
    "Manifold",
    "ManifoldSampler",
    "SamplerParams",
    "join",
    "gc_collect",
    "gc_minimize",
]
