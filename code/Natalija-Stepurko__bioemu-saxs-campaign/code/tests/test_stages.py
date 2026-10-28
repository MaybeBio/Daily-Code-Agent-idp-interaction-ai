"""Small units of the coordrg and envelopes stages that need no data."""
import numpy as np
import pandas as pd

from bsc import coordrg, envelopes


def test_coordrg_selection_is_seeded_and_covers_the_examples():
    n = 120
    rng = np.random.default_rng(0)
    cls = np.array(["folded"] * 70 + ["partly disordered"] * 30 + ["disordered"] * 20)
    d = pd.DataFrame({"model": "bioemu", "clean": True, "label": [f"S{i:03d}" for i in range(n)],
                      "disorder_class": cls, "chi2_raw": rng.lognormal(1, 0.5, n),
                      "length": rng.integers(50, 400, n), "rg_ratio": rng.lognormal(0, 0.1, n)})
    d.loc[5, "clean"] = False
    sel = coordrg.selection(d)
    assert sel == coordrg.selection(d)                      # deterministic
    assert sum(v == "disordered" for v in sel.values()) == 20
    assert sum(v == "folded sample" for v in sel.values()) == coordrg.N_FOLDED_SAMPLE
    assert "S005" not in sel                                 # flagged entries are not analysed
    members = coordrg.archive_members(sorted(sel), ["S001"])
    assert len(members) == 2 * len(sel) + 1 and members[-1].endswith("alphafold/S001.pdb")


def test_denss_logs_are_parsed(tmp_path):
    log = tmp_path / "x_0.log"
    log.write_text("... - Final Chi2: 1.529e+00\n... - Final Rg: 21.353\n"
                   "... - Final Support Volume: 55507.452\n")
    rec = envelopes.parse_final_log(log)
    assert rec == {"final_chi2": 1.529, "final_rg": 21.353, "support_volume": 55507.452}
    avg = tmp_path / "x_final.log"
    avg.write_text("Mean of correlation scores: 0.858\nNumber of aligned maps accepted: 9\n"
                   "Resolution = 25.3 +- 3.7 A\n")
    rec = envelopes.parse_average_log(avg)
    assert rec == {"resolution_A": 25.3, "resolution_sd_A": 3.7, "mean_correlation": 0.858,
                   "maps_accepted": 9}


def test_operating_point_is_the_first_point_reaching_the_target():
    assert envelopes.operating_point([5.0, 3.0, 1.9, 1.2], raw_chi2=6.0) == 2
    assert envelopes.operating_point([1.5, 1.2], raw_chi2=1.8) == -1          # the raw ensemble already fits
    assert envelopes.operating_point([5.0, 4.0, 3.0], raw_chi2=6.0) is None   # never reached


def test_density_levels_are_the_particle_and_the_protein_volume():
    assert envelopes.density_volumes(100.0, 120.0) == [150.0, 100.0]       # particle at least 1.5x
    assert envelopes.density_volumes(100.0, 400.0) == [400.0, 100.0]


def test_examples_are_nearest_the_class_median_raw_chi2():
    from bsc import analyse
    rows = []
    for cls in ("folded", "partly disordered", "disordered"):
        for i, (chi2, length) in enumerate([(2.0, 120), (2.0, 90), (8.0, 100), (1.0, 80), (4.0, 110)]):
            rows.append({"model": "bioemu", "clean": True, "label": f"{cls[:3]}{i}", "disorder_class": cls,
                         "chi2_raw": chi2, "length": length})
    ex = analyse.choose_examples(pd.DataFrame(rows))
    # median 2.0 is shared by two entries; the shorter chain is taken
    assert ex == {"folded": "fol1", "partly disordered": "par1", "disordered": "dis1"}
