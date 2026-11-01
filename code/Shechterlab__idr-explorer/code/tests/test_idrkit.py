import math

import numpy as np
import pytest

import idrkit as K


def test_charge_and_his_scales():
    assert K.charge_vector("DEKRHA").tolist() == [-1, -1, 1, 1, 0, 0]
    assert K.charge_vector("H", his=K.HIS["emboss"]).tolist() == [0.5]
    assert math.isclose(K.ncpr("KKDA"), 0.25)
    assert math.isclose(K.fcr("KKDA"), 0.75)


def test_scd_reference_values():
    # (K8E8)3 is -4.290 in the literature (arXiv:2201.11834); alternation is far less negative, a polyelectrolyte positive
    assert round(K.scd("KKKKKKKKEEEEEEEE" * 3), 3) == -4.290
    assert K.scd("KE" * 24) > -1.0
    assert K.scd("K" * 20) > 0
    assert K.scd("AAAA") == 0.0


def test_scd_matches_the_double_sum():
    rng = np.random.default_rng(0)
    for _ in range(20):
        seq = "".join(rng.choice(list("KRDEAGSH"), size=int(rng.integers(2, 80))))
        q = K.charge_vector(seq)
        brute = sum(q[i] * q[j] * math.sqrt(j - i) for i in range(len(q)) for j in range(i + 1, len(q))) / len(q)
        assert math.isclose(K.scd(seq), brute, abs_tol=1e-9)


def test_acidic_stretches_as_glu_tests_them():
    seq = "MKTAYIAKQRQ" + "GSDEEDEEEDSEEG" + "KKLLVV"
    (s, e), = K.acidic_stretches(seq)
    assert seq[s:e] == "DEEDEEEDSEE"
    assert K.acidic_stretches("AAAA" + "EEKEEKEEKEEK" + "AAAA") == []
    assert len(K.acidic_stretches("A" * 10 + "EEEEEDDE" + "GG" + "EEEDDEEE" + "A" * 10)) == 1


def test_runs_rg_and_motifs():
    assert K.residue_runs("AEEEEAEEE", "E", 4) == [(1, 5)]
    assert K.acidic_runs("AADEEDEEAEEE") == [(2, 8)]
    assert K.rg_clusters("AARGGRGGSRGGYAA") == [(2, 11)]
    assert K.rg_motif_count("RGRGR") == 4
    assert K.e_tract_count("EEDEE") == 6
    assert K.ck2_motif_count("SAAEAATDDE") == 2 and K.ck2_motif_sites("SAAEAATDDE") == [0, 6]
    assert K.longest_acidic_run("ADEEDAEE") == 4


def test_basic_patches():
    assert K.basic_patches("A" * 5 + "KRKAAAAAAA" + "A" * 20) == [(0, 15)]
    assert K.basic_patches("KDKDKDKDKD") == []


def test_scd_z_separates_blocks_from_mixing():
    blocky = "A" * 10 + "KKKKKKKKEEEEEEEE" * 3 + "A" * 10
    mixed = "A" * 10 + "KE" * 24 + "A" * 10
    assert K.scd_z(blocky)[1] < -1.5
    assert K.scd_z(mixed)[1] > 1
    assert K.scd_z("AAAA") == (0.0, 0.0)


def test_idr_intervals():
    assert K.merge_intervals([(5, 9), (1, 4), (8, 20)]) == [[1, 4], [5, 20]]
    assert K.idr_intervals([(1, 39), (30, 60), (70, 80), (90, 140)], 120) == [(1, 60), (90, 120)]


def _localcider_delta(seq):
    # localCIDER backend/sequence.py deltaForm() and delta(), as a plain loop
    cp = np.array([1 if c in "KR" else -1 if c in "DE" else 0 for c in seq])
    n = len(seq)
    fp, fn = (cp > 0).sum() / n, (cp < 0).sum() / n
    sigma = (fp - fn) ** 2 / (fp + fn) if fp + fn else 0

    def form(g):
        nb, ans = n - g + 1, 0.0
        for i in range(nb):
            b = cp[i:i + g]
            bp, bn = (b > 0).sum(), (b < 0).sum()
            ans += (sigma - (((bp - bn) / g) ** 2 / ((bp + bn) / g) if bp + bn else 0)) ** 2 / nb
        return ans
    return (form(5) + form(6)) / 2


def test_delta_matches_localcider():
    rng = np.random.default_rng(1)
    for _ in range(50):
        seq = "".join(rng.choice(list("KRDEAGSHP"), size=int(rng.integers(10, 200))))
        assert math.isclose(K.delta(seq), _localcider_delta(seq), abs_tol=1e-12)


def test_delta_z_sees_local_blocks():
    assert K.delta_z("KKKKKEEEEE" * 6)[1] > 2
    assert K.delta_z("KEKEKEKEKE" * 6)[1] < -1


def test_modified_charge_and_shift():
    seq = "KKKKKSSSSSEEEEE"
    q = K.modified_charge(seq, {5: "phospho", 0: "acetyl"})
    assert q[5] == -2 and q[0] == 0 and q[10] == -1
    sh = K.charge_shift(seq, {5: "phospho", 6: "phospho"})
    assert sh["ncpr0"] == 0 and abs(sh["ncpr1"] - (-4 / 15)) < 1e-12
    assert sh["fcr1"] > sh["fcr0"]
    # phosphorylating the serines next to the acidic block lengthens it: charge blocks become more segregated
    assert K.charge_shift(seq, {i: "phospho" for i in range(5, 10)})["delta1"] > sh["delta0"]
    # methylation keeps the charge; a modification the residue cannot carry is refused
    assert (K.modified_charge(seq, {0: "methyl"}) == K.charge_vector(seq)).all()
    with pytest.raises(ValueError):
        K.modified_charge(seq, {0: "phospho"})


def test_charge_type_in_words():
    from idrkit import charge_type, fcr, ncpr
    for seq, want in (("EEEEEDDDDDAAAAAAAAAA", "acidic"), ("KKKKKRRRRRAAAAAAAAAA", "basic"), ("EKEKEKEKAAAAAAAAAAAA", "mixed"),
                      ("EKAAAAAAAAAAAAAAAAAA", "few"), ("DDDDAAAAAAAAAAAAAAAA", "acidic")):
        assert charge_type(fcr(seq), ncpr(seq)) == want, seq
