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

from discrete_differential_geometry import tcp_reference as tr
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


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_orientation_sign_is_a_homomorphism(name, m, order, point):
    """The sign must be well defined (same at every tet) and multiplicative."""
    _, sym = _sym(name, m)
    view = sym.view
    for g in sym.elements[:8]:
        seen = set()
        for t in range(min(view.nT, 40)):
            o = view.oriented_order(t)
            img = tuple(int(g[x]) for x in o)
            o2 = view.oriented_order(view.tetid[tuple(sorted(img))])
            from discrete_differential_geometry.symmetry import _perm_sign
            seen.add(_perm_sign(tuple(o2.index(x) for x in img)))
        assert len(seen) == 1, "orientation sign differs between tets"
    from discrete_differential_geometry.symmetry import _compose
    for g in sym.elements[:6]:
        for h in sym.elements[:6]:
            assert sym.orientation_sign(_compose(g, h)) == \
                sym.orientation_sign(g) * sym.orientation_sign(h)


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_orientation_preserving_subgroup(name, m, order, point):
    """Aut+ must be an honest index-1-or-2 subgroup of orientation-preserving
    elements, usable as a group in its own right."""
    _, sym = _sym(name, m)
    plus = sym.orientation_preserving
    assert plus.order in (sym.order, sym.order // 2)
    assert plus.order == (sym.order // 2 if sym.has_orientation_reversing
                          else sym.order)
    assert all(sym.orientation_sign(g) > 0 for g in plus.elements)
    full = {g.tobytes() for g in sym.elements}
    assert all(g.tobytes() in full for g in plus.elements)
    # every inherited method works against the subgroup
    for kind in ("vertex", "edge", "face", "tet"):
        assert sum(plus.orbit_sizes(kind)) == sum(sym.orbit_sizes(kind))
        assert plus.n_orbits(kind) >= sym.n_orbits(kind)
        rep = plus.orbit_representatives(kind)[0]
        assert len(plus.stabilizer(kind, rep)) == \
            plus.stabilizer_order(kind, rep)


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_chirality_matches_orbit_splitting(name, m, order, point):
    """is_chiral must agree with the stabilizer definition: a class is achiral
    exactly when an orientation-reversing element fixes it."""
    _, sym = _sym(name, m)
    if not sym.has_orientation_reversing:
        assert sym.is_chiral("tet", sym.orbit_representatives("tet")[0]) is None
        return
    for kind in ("vertex", "edge", "face", "tet"):
        for rep in sym.orbit_representatives(kind)[:6]:
            by_stab = all(sym.orientation_sign(g) > 0
                          for g in sym.stabilizer(kind, rep))
            assert sym.is_chiral(kind, rep) == by_stab
    # orbit counts: each chiral orbit splits in two, achiral ones do not
    plus = sym.orientation_preserving
    for kind in ("vertex", "edge", "face", "tet"):
        n_chiral = sum(1 for r in sym.orbit_representatives(kind)
                       if sym.is_chiral(kind, r))
        assert plus.n_orbits(kind) == sym.n_orbits(kind) + n_chiral


def test_bc_chains_are_always_chiral():
    """A tetrahelix is intrinsically handed, so no orientation-reversing
    automorphism can ever fix a BC chain -- every chain class splits."""
    for name, m in (("a15", 2), ("c15", 2), ("r", 2)):
        _, sym = _sym(name, m)
        assert sym.has_orientation_reversing
        n = len(sym.chain_orbits()[1])
        assert all(sym.is_chiral("chain", i) for i in sym.chain_orbits()[2])
        assert len(sym.orientation_preserving.chain_orbits()[1]) == 2 * n


def test_sohncke_group_has_no_mirror():
    """The delta phase is P2_1 2_1 2_1, a chiral (Sohncke) space group: the
    combinatorics alone must show zero orientation-reversing automorphisms,
    and then chirality questions have no answer."""
    facets = np.asarray(tr.build_t3_triangulation("delta", 2)[0])
    sym = CrystalSymmetry.compute(facets)
    assert not sym.has_orientation_reversing
    assert sym.orientation_preserving is sym
    assert sym.is_chiral("chain", 0) is None
    assert sym.is_chiral("tet", sym.orbit_representatives("tet")[0]) is None


def test_chain_tables_shared_between_group_and_subgroup():
    """Chains belong to the triangulation, so Aut and Aut+ must not each pay
    for the 24*nT enumeration."""
    _, sym = _sym("a15", 2)
    assert sym.chains is sym.orientation_preserving.chains


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


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_config_key_matches_canonical_frame_key(name, m, order, point):
    """The O(1)-per-query table form must induce the SAME partition as the
    min-over-group form -- it is an optimisation, not a different notion."""
    _, sym = _sym(name, m)
    nF = 24 * sym.view.nT
    frames = list(range(0, nF, max(1, nF // 60)))[:60]
    a, b = {}, {}
    for f in frames:
        a.setdefault(sym.config_key([f]), []).append(f)
        b.setdefault(sym.canonical_frame_key([f]), []).append(f)
    assert sorted(map(sorted, a.values())) == sorted(map(sorted, b.values()))
    pairs = [(frames[i], frames[(i + 3) % len(frames)])
             for i in range(len(frames))]
    a, b = {}, {}
    for p in pairs:
        a.setdefault(sym.config_key(list(p)), []).append(p)
        b.setdefault(sym.canonical_frame_key(list(p)), []).append(p)
    assert sorted(map(sorted, a.values())) == sorted(map(sorted, b.values()))


@pytest.mark.parametrize("name,m,order,point", CRYSTALS)
def test_config_key_is_an_equivalence_certificate(name, m, order, point):
    """Equal keys must mean an automorphism really carries one to the other,
    and unequal keys that no automorphism does."""
    _, sym = _sym(name, m)
    base = 0
    imgs = [sym.view.frame_id(tuple(int(g[x]) for x in sym.view.frame_window(base)))
            for g in sym.elements[:12]]
    for f in imgs:                       # in the orbit by construction
        assert sym.config_key([f]) == sym.config_key([base])
        assert sym.frames_equivalent([f], [base])
    # a frame in a DIFFERENT orbit must differ
    orb, _, _ = sym.frame_canonicalizer()
    other = int(np.nonzero(orb != orb[base])[0][0])
    assert sym.config_key([other]) != sym.config_key([base])
    assert not sym.frames_equivalent([other], [base])


def test_anchored_key_has_no_residual_freedom():
    """With the anchor pinned there is no symmetry left -- Aut acts freely on
    frames -- so two configurations sharing an anchor are equivalent only if
    their other frames are literally equal. This is why pinning knot A leaves
    every distinct B geometry distinct."""
    _, sym = _sym("a15", 2)
    nF = 24 * sym.view.nT
    a = 0
    seen = {}
    for f in range(0, nF, max(1, nF // 40)):
        seen.setdefault(sym.config_key([a, f]), []).append(f)
    assert all(len(v) == 1 for v in seen.values())


# ---------------------------------------------------------------------------
# the 2->3 ball's own symmetry (ball_boundary), distinct from a host's Aut
# ---------------------------------------------------------------------------


def test_ball_config_symmetry_preserves_the_configuration():
    """S3 on the face x S2 on the apexes must fix TETS_R, TETS_D and BOUNDARY
    setwise -- otherwise it is not a symmetry of the ball at all."""
    from discrete_differential_geometry import ball_boundary as bb
    assert len(bb.CONFIG_SYMMETRY) == 12
    for name in ("TETS_R", "TETS_D", "BOUNDARY"):
        ref = {tuple(sorted(t)) for t in getattr(bb, name)}
        for g in bb.CONFIG_SYMMETRY:
            assert {tuple(sorted(g[v] for v in t))
                    for t in getattr(bb, name)} == ref, name


@pytest.mark.parametrize("config", ["R", "D"])
def test_ball_boundary_map_is_equivariant(config):
    """d(gx, gy) == d(x, y) EXACTLY. The ball is the whole world here, so this
    is an exact symmetry of the object, and it is independent of the other
    checks: a wrong strip placement, chord test or node index would break it."""
    from discrete_differential_geometry import ball_boundary as bb
    m = bb.BallBoundaryMap(config, 3)
    M, missing = m.matrix()
    assert not missing
    for g in bb.CONFIG_SYMMETRY:
        perm = m.node_permutation(g)
        assert np.array_equal(M, M[np.ix_(perm, perm)])


def test_ball_node_permutation_is_a_permutation():
    from discrete_differential_geometry import ball_boundary as bb
    m = bb.BallBoundaryMap("D", 3)
    for g in bb.CONFIG_SYMMETRY:
        p = m.node_permutation(g)
        assert sorted(p) == list(range(len(m.nodes)))
