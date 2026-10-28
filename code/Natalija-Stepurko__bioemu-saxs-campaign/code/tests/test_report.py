"""Figure conventions of the report stage."""
from bsc import report


def test_every_model_has_a_colour():
    missing = [m for m in report.MODEL_NAME if m not in report.MODEL_COL]
    assert not missing, missing
    assert all(c.startswith("#") and len(c) == 7 for c in report.MODEL_COL.values())
    # the model under study takes its role colour; no figure role is drawn in ink
    assert report.MODEL_COL["bioemu"] == report.ROLE_COL["bioemu"]
    assert report.INK not in report.ROLE_COL.values()
    assert len(set(report.ROLE_COL.values())) == len(report.ROLE_COL)
