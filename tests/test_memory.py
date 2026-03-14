"""Diagnose D GC allocation sources."""

import ctypes
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams, gc_collect, gc_minimize
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
    base = gc_used_mb()
    print(f"Baseline GC used: {base:.2f} MB")

    # Test 1: Just sampler.run()
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        sampler.run(moves=50, exact=False)
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x run(moves=50):    GC delta = {after - before:+.2f} MB")

    # Test 2: Just simplices queries
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        verts = sampler.simplices(0)
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x simplices(0):     GC delta = {after - before:+.2f} MB")

    # Test 3: Just degree queries
    gc_collect()
    before = gc_used_mb()
    verts = sampler.simplices(0)
    for _ in range(100):
        for v in verts:
            sampler.degree(v)
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x degree loops:     GC delta = {after - before:+.2f} MB")

    # Test 4: Just f_vector
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        fv = sampler.f_vector
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x f_vector:         GC delta = {after - before:+.2f} MB")

    # Test 5: Just importance_weight
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        w = sampler.importance_weight()
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x importance_wt:    GC delta = {after - before:+.2f} MB")

    # Test 6: Combined (like the real test)
    gc_collect()
    before = gc_used_mb()
    for _ in range(100):
        sampler.run(moves=50, exact=False)
        fv = sampler.f_vector
        verts = sampler.simplices(0)
        for v in verts:
            sampler.degree(v)
        w = sampler.importance_weight()
    gc_collect()
    after = gc_used_mb()
    print(f"After 100x combined:         GC delta = {after - before:+.2f} MB")


if __name__ == "__main__":
    main()
