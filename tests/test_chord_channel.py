"""Unit-testable invariants of the strict chord channel (ddg_lab.f0_worm).

The campaign validators (twosided_chord.py, agg_knobs_test.py) certify the
channel statistically on equilibrated production states; these tests pin the
exact invariants they rely on, on a synthesized state:

  * a chord episode conserves f0 EXACTLY (the channel's defining property --
    two-sided tests are only valid because both starts share an f0 sector);
  * the manifold stays a closed 3-manifold (Euler 0) through episodes;
  * fixed seed => bit-identical episode stream;
  * the aggregation knobs at (region_max=0, agg_beta=0.0) reproduce the
    knobless channel bit-for-bit (the agg_knobs_test regression).
"""
import numpy as np
import pytest

import discrete_differential_geometry as ddg
from discrete_differential_geometry.tcp_reference import build_t3_triangulation
from ddg_lab import f0_worm as F

NPREP = 6          # forced 2->3 flips so the channel has chords to work with
NEP = 150


def _prepped_sampler(seed):
    """Fresh c15 m=2 state with NPREP forced flips + the strict chord channel
    armed exactly as in twosided_chord.py."""
    fac, _ = build_t3_triangulation("c15", 2)
    m = ddg.Manifold(3, np.asarray(fac).tolist())
    F.F3T = m.num_facets                    # pin volume at the actual state
    ddg.set_random_seed(seed)
    s, L = F.fresh(m)
    done = 0
    tets = [t for ts in L.v2t.values() for t in ts]
    for tt in tets:
        if done >= NPREP:
            break
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(tt) if j != i))
            fs = set(face)
            sh = [q for q in L.v2t[face[0]] if fs <= set(q)]
            if len(sh) != 2:
                continue
            ap = tuple(sorted({x for q in sh for x in q} - fs))
            if len(ap) != 2 or ap in L.edeg:
                continue
            try:
                L.do(face, ap)
            except Exception:
                continue
            done += 1
            break
    assert done == NPREP
    s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, F.Z0, lmax=F.LMAX,
                  zeta=1.0, aof=0.5, ph=0.5, pg=0.3, bcf=1e-4,
                  bc4=1.0 - 0.3 - 0.5 - 1e-4, maxstep=800,
                  ucap_hi=50.0, ucap_lo=-50.0, mu=F.MU)
    s.set_worm_pair(zeta2=float("nan"), bcp=0.05, chain_k=20)
    return s


def _run(seed, nep=NEP, agg=None):
    s = _prepped_sampler(seed)
    if agg is not None:
        s.set_worm_chord_agg(region_max=agg[0], agg_beta=agg[1])
    f0_start = int(s.manifold.f_vector[0])
    stream = []
    for _ in range(nep):
        r = s.worm_chord_strict_episode()
        if r["changed"]:
            stream.append(tuple(r["df"]))
    fv = [int(x) for x in s.manifold.f_vector]
    return f0_start, fv, stream


def test_chord_episodes_conserve_f0_and_manifoldness():
    f0_start, fv, stream = _run(11)
    assert len(stream) > 0, "no episode committed; test has no power"
    assert fv[0] == f0_start                       # f0 conserved EXACTLY
    assert fv[0] - fv[1] + fv[2] - fv[3] == 0      # still a closed 3-manifold
    assert all(df[0] == 0 for df in stream)        # per-episode df0 == 0


def test_chord_episode_stream_is_deterministic():
    a = _run(23)
    b = _run(23)
    assert a == b


def test_agg_knobs_zero_reproduce_knobless_channel():
    base = _run(37)
    knobs = _run(37, agg=(0, 0.0))
    assert base == knobs
