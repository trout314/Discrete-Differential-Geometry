"""Diagnose D GC allocation sources."""

import ctypes
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams, gc_collect
from discrete_differential_geometry._dlang import _lib

_lib.ddg_gc_stats.argtypes = [ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
_lib.ddg_gc_stats.restype = None


def gc_used_mb():
    used = ctypes.c_long()
    _lib.ddg_gc_stats(ctypes.byref(used), None)
    return used.value / 1024 / 1024


def main():
    params = SamplerParams(
        num_facets_target=30, hinge_degree_target=4.5,
        num_facets_coef=0.5, num_hinges_coef=0.1,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    )

    seed = Manifold.standard_sphere(3)
    sampler = ManifoldSampler(seed, params)
    sampler.run(moves=2000, exact=False)

    gc_collect()
    print(f"Baseline GC used: {gc_used_mb():.2f} MB\n")

    # Test A: One big run of 5000 moves
    gc_collect()
    before = gc_used_mb()
    sampler.run(moves=5000, exact=False)
    gc_collect()
    after = gc_used_mb()
    print(f"1x run(5000 moves):  GC delta = {after - before:+.2f} MB  ({(after-before)/5000*1000:.1f} KB/move)")

    # Test B: 100x run of 50 moves each
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        sampler.run(moves=50, exact=False)
    gc_collect()
    after = gc_used_mb()
    print(f"100x run(50 moves):  GC delta = {after - before:+.2f} MB  ({(after-before)/5000*1000:.1f} KB/move)")

    # Test C: 5000x run of 1 move each
    gc_collect()
    before = gc_used_mb()
    for _ in range(5000):
        sampler.run(moves=1, exact=False)
    gc_collect()
    after = gc_used_mb()
    print(f"5000x run(1 move):   GC delta = {after - before:+.2f} MB  ({(after-before)/5000*1000:.1f} KB/move)")


if __name__ == "__main__":
    main()
