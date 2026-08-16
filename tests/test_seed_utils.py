"""Tests for seed metadata / filename / provenance helpers (ddg.seed_utils)."""
import numpy as np
import pytest

from discrete_differential_geometry import seed_utils as su


def test_encode_float():
    assert su.encode_float(0) == "0"
    assert su.encode_float(1000) == "1e3"
    assert su.encode_float(0.1) == "1e-1"
    assert su.encode_float(0.002) == "2e-3"
    assert su.encode_float(2000) == "2e3"
    assert su.encode_float(5.1) == "5p1"
    assert su.encode_float(7) == "7"


@pytest.mark.parametrize("x", [0, 1000, 0.1, 0.002, 5.1, 7, 5.1043, 0.032])
def test_encode_decode_roundtrip(x):
    assert su.decode_float(su.encode_float(x)) == pytest.approx(x)


def test_round_sig():
    assert su.round_sig(0.1) == 0.1
    assert su.round_sig(0.101266) == 0.101
    assert su.round_sig(0) == 0.0


class _Params:
    num_facets_target = 10000
    num_facets_coef = 0.1
    hinge_degree_target = 5.1043
    num_hinges_coef = 2.0
    hinge_degree_variance_coef = 0.0
    codim3_degree_variance_coef = 20.0     # beta/N = 2e-3
    codim3_degree_target_coef = 0.0
    hinge_degree_target_coef = 0.0


def test_build_seed_filename():
    name = su.build_seed_filename("S3", _Params(), seed_index=7)
    assert name == "S3_N1e4_1e-1_ED5p1043_2_VDVs_2e-3_s007.mfd"


def test_header_roundtrip(tmp_path):
    p = tmp_path / "x.mfd"
    p.write_text("# a = 1\n# b = two\n0 1 2 3\n")
    su.set_header_field(str(p), "b", "three")
    su.set_header_field(str(p), "c", "4")
    md = su.load_seed_metadata(str(p))
    assert md == {"a": "1", "b": "three", "c": "4"}
    assert p.read_text().endswith("0 1 2 3\n")   # body untouched


def test_history_legs():
    assert su.is_root_leg_from("sphere")
    assert su.is_root_leg_from("crystal:a15@wyckoff")
    assert not su.is_root_leg_from("prev")
    leg = su.make_leg("build", {"struct": "a15"}, sweeps=5, from_="sphere",
                      commit="deadbeefcafe", dirty=False, tried=10, accepted=4)
    assert leg["commit"] == "deadbee" and not leg["dirty"]
    fields = su.history_fields([leg])
    assert fields[0] == "history_root = sphere"
    assert "history_total_sweeps = 5" in fields
    assert fields[-1].startswith("history = [")


def test_build_metadata_requires_lineage():
    common = dict(
        topology="S3", dimension=3, initial_triangulation="sphere",
        num_facets_target=100, num_facets_coef=0.1, hinge_degree_target=5.1,
        num_hinges_coef=2.0, hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.0, growth_step_size=0,
        eq_sweeps_per_step=0, equilibration_sweeps=0,
        manifold_view=None, objective=0.0)
    with pytest.raises(ValueError, match="legs"):
        su.build_metadata_comments(**common)
    leg = su.make_leg("eq", {}, sweeps=1, from_="sphere",
                      commit="x" * 8, dirty=False)
    with pytest.raises(ValueError, match="exactly one"):
        su.build_metadata_comments(**common, legs=[leg])          # neither
    with pytest.raises(ValueError, match="exactly one"):
        su.build_metadata_comments(**common, legs=[leg],
                                   prior_history=[], root="sphere")  # both
