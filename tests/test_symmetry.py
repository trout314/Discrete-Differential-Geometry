"""Ground-truth tests for the exact combinatorial symmetry group (symmetry.py).

Crystals are rebuilt in-memory from Wyckoff positions (tcp_reference), so these
tests need no data files -- the same convention as test_crystal_grains.py.

Three kinds of assertion:

  * STRUCTURAL, which must hold for any triangulation: generators really are
    automorphisms, orbit-stabilizer, chain lengths partition the frames, and
    Aut acts FREELY on frames (the fact that makes the whole computation cheap
    -- an automorphism fixing one ordered tet is forced to be the identity).
  * SUPERCELL-INDEPENDENCE: orbit counts are properties of the crystal, not of
    the m x m x m box, so a15 at m=2 and m=3 must agree. Chain LENGTHS scale
    with the box; the number of chain classes does not.
  * KNOWN VALUES: |Aut| factorises as (m^3 translations) x (centering) x
    (point group), with no accidental automorphisms, for every crystal tested.
    See the module docstring of symmetry.py for why that equality is an
    empirical result rather than something guaranteed -- Aut(K) is a property
    of the abstract complex and is only guaranteed to CONTAIN the crystal's
    space group.
"""
import itertools
import os
import sys
from functools import lru_cache

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
for _p in ("scripts", "tools"):
    sys.path.insert(0, os.path.join(_ROOT, _p))

import tcp_reference as tr
from discrete_differential_geometry.symmetry import (
    CrystalSymmetry, TriView, develop_total)

# (phase, supercell, |Aut|, point-group x centering) -- smallest quotients that
# validate, kept tiny for speed.
CRYSTALS = [
    ("a15", 2, 384, 48),        # Pm-3n,  |P| = 48
    ("a15", 3, 1296, 48),
    ("c15", 2, 1536, 192),      # Fd-3m,  |P| = 48, fcc centering x4
    ("r", 2, 144, 18),          # R-3,    |P| =  6, rhombohedral centering x3
]


@lru_cache(maxsize=None)
def _sym(name, m):
    facets = np.asarray(tr.build_t3_triangulation(name, m)[0])
    return facets, CrystalSymmetry.compute(facets)


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_order_factorises_as_space_group(name, m, order, point):
    _, sym = _sym(name, m)
    assert sym.order == order
    assert sym.order == m ** 3 * point, (
        "|Aut| is not (supercell translations) x (point group x centering); "
        "either the crystal has accidental automorphisms or this is a bug")


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_generators_are_automorphisms(name, m, order, point):
    """Each generator must carry the facet SET onto itself."""
    facets, sym = _sym(name, m)
    view = sym.view
    ref = {tuple(sorted(int(x) for x in tv)) for tv in view.tets}
    assert len(ref) == view.nT
    for g in sym.generators:
        moved = {tuple(sorted(int(g[int(x)]) for x in tv)) for tv in view.tets}
        assert moved == ref


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_free_action_on_frames(name, m, order, point):
    """An automorphism fixing an ordered tet is the identity, so every frame
    orbit has size exactly |Aut| and 24*nT is divisible by it."""
    _, sym = _sym(name, m)
    assert (24 * sym.view.nT) % sym.order == 0
    base = tuple(int(x) for x in sym.view.tets[0])
    assert len({tuple(int(g[v]) for v in base) for g in sym.elements}) == sym.order
    ident = develop_total(sym.view, base, base)
    assert ident is not None and np.array_equal(ident, np.arange(sym.view.V))


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
@pytest.mark.parametrize("kind", ["vertex", "edge", "face", "tet"])
def test_orbit_stabilizer(name, m, order, point, kind):
    _, sym = _sym(name, m)
    totals = {"vertex": sym.view.V, "edge": len(sym.view.edges),
              "face": len(sym.view.faces), "tet": sym.view.nT}
    sizes = sym.orbit_sizes(kind)
    assert sum(sizes) == totals[kind]
    for i, rep in enumerate(sym.orbit_representatives(kind)):
        assert sym.orbit_id(kind, rep) == i
        assert sizes[i] * sym.stabilizer_order(kind, rep) == sym.order
        assert len(sym.stabilizer(kind, rep)) == sym.stabilizer_order(kind, rep)


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_orbits_are_closed_under_the_group(name, m, order, point):
    """Every element must map each orbit representative back into its orbit."""
    _, sym = _sym(name, m)
    for kind in ("vertex", "edge", "face", "tet"):
        for i, rep in enumerate(sym.orbit_representatives(kind)):
            members = set(sym.orbit_members(kind, i))
            for g in sym.elements[:16]:
                img = sym.act(g, rep)
                assert (img if kind == "vertex" else tuple(sorted(img))) in members


def test_orbit_counts_are_supercell_independent():
    """Orbit counts describe the crystal, not the box, so m=2 and m=3 agree."""
    _, s2 = _sym("a15", 2)
    _, s3 = _sym("a15", 3)
    for kind in ("vertex", "edge", "face", "tet"):
        assert s2.n_orbits(kind) == s3.n_orbits(kind)
    assert len(s2.chain_orbits()[1]) == len(s3.chain_orbits()[1])
    assert 24 * s2.view.nT // s2.order == 24 * s3.view.nT // s3.order


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_chains_partition_the_frames(name, m, order, point):
    _, sym = _sym(name, m)
    chains = sym.chains
    assert sum(len(c) for c in chains) == 24 * sym.view.nT
    assert len({int(f) for c in chains for f in c}) == 24 * sym.view.nT


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_chain_walk_matches_independent_implementation(name, m, order, point):
    """The chain through a frame must agree with a direct sliding-window walk
    (the `worm_helix.bc_orbit` recurrence, rebuilt here from face lookups)."""
    _, sym = _sym(name, m)
    view = sym.view
    for fid in (0, 7, 24 * view.nT // 2 + 3):
        w = list(view.frame_window(fid))
        walk = []
        while True:
            if tuple(w) in walk:
                assert walk.index(tuple(w)) == 0
                break
            walk.append(tuple(w))
            face = tuple(sorted(w[1:]))
            a, b = (x[1] for x in view.face2[face])
            w = w[1:] + [b if a == w[0] else a]
        cid = int(sym._chain_tables()[1][fid])
        chain = [int(f) for f in sym.chains[cid]]
        assert len(walk) == len(chain)
        k = chain.index(fid)                        # cycles are stored from
        chain = chain[k:] + chain[:k]               # wherever enumeration
        assert [view.frame_window(f) for f in chain] == \
            [tuple(x) for x in walk]                # first entered them


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_chain_orbits_are_consistent(name, m, order, point):
    """Chains in one orbit must have equal length, and orbit sizes must divide
    |Aut| (the group acts on the chain set)."""
    _, sym = _sym(name, m)
    _, members, _ = sym.chain_orbits()
    assert sum(len(mm) for mm in members) == len(sym.chains)
    for mm in members:
        assert len({len(sym.chains[i]) for i in mm}) == 1
        assert sym.order % len(mm) == 0


def test_melt_has_trivial_symmetry():
    """A crystal with one Pachner move applied has no automorphisms left, and
    the computation must return the trivial group rather than fail."""
    facets = np.asarray(tr.build_t3_triangulation("a15", 2)[0]).tolist()
    view = TriView(np.asarray(facets))
    face = sorted(view.face2)[0]
    (t1, a1), (t2, a2) = view.face2[face]
    lab = view.labels
    keep = [f for f in facets
            if tuple(sorted(f)) not in
            (tuple(sorted(int(lab[x]) for x in view.tets[t1])),
             tuple(sorted(int(lab[x]) for x in view.tets[t2])))]
    p, q = int(lab[a1]), int(lab[a2])
    for e in itertools.combinations((int(lab[x]) for x in face), 2):
        keep.append(sorted([e[0], e[1], p, q]))
    sym = CrystalSymmetry.compute(np.asarray(keep))
    assert sym.order == 1
    assert sym.generators == []
    assert sym.n_orbits("vertex") == sym.view.V


def test_disconnected_input_is_refused():
    """Development from one frame cannot see automorphisms that permute
    components, so it would report a subgroup (typically trivial). That is a
    silently wrong answer, so it must raise instead."""
    a = np.asarray(tr.build_t3_triangulation("a15", 2)[0])
    two = np.vstack([a, a + (int(a.max()) + 1)])       # two disjoint copies
    with pytest.raises(ValueError, match="disconnected"):
        CrystalSymmetry.compute(two)


def test_cache_is_verified_not_trusted(tmp_path):
    """A sidecar whose digest matches but whose generators are wrong must be
    rejected, not loaded."""
    from discrete_differential_geometry import Manifold
    facets = np.asarray(tr.build_t3_triangulation("a15", 2)[0])
    path = str(tmp_path / "a15_m2.mfd")
    Manifold(3, facets.tolist()).save(path)
    good = CrystalSymmetry.for_manifold_path(path)
    side = path + ".sym.npz"
    digest = good.view.digest()

    # a permutation that is NOT an automorphism, under the right digest
    bogus = np.arange(good.view.V, dtype=np.int32)
    bogus[0], bogus[1] = bogus[1], bogus[0]
    np.savez_compressed(side, digest=digest, order=2,
                        gens=np.stack([bogus]))
    with pytest.warns(UserWarning, match="not an automorphism"):
        assert CrystalSymmetry.for_manifold_path(path).order == good.order

    # a real generator, but an inflated order claim
    np.savez_compressed(side, digest=digest, order=999999,
                        gens=np.stack(good.generators))
    with pytest.warns(UserWarning, match="close to"):
        assert CrystalSymmetry.for_manifold_path(path).order == good.order

    # not a permutation at all
    np.savez_compressed(side, digest=digest, order=1,
                        gens=np.zeros((1, good.view.V), np.int32))
    with pytest.warns(UserWarning, match="non-permutation"):
        assert CrystalSymmetry.for_manifold_path(path).order == good.order


def test_cache_round_trip(tmp_path):
    from discrete_differential_geometry import Manifold
    facets = np.asarray(tr.build_t3_triangulation("a15", 2)[0])
    path = str(tmp_path / "a15_m2.mfd")
    Manifold(3, facets.tolist()).save(path)

    first = CrystalSymmetry.for_manifold_path(path)
    assert os.path.exists(path + ".sym.npz")
    second = CrystalSymmetry.for_manifold_path(path)
    assert second.order == first.order
    assert second.n_orbits("face") == first.n_orbits("face")

    # a stale sidecar (wrong content hash) must be ignored, not trusted
    np.savez_compressed(path + ".sym.npz", digest="deadbeef", order=999,
                        gens=np.zeros((0, first.view.V), np.int32))
    third = CrystalSymmetry.for_manifold_path(path)
    assert third.order == first.order
