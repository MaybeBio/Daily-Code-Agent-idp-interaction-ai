import sys
from scripts import ingest_daily


def test_daily_mode_ingests_yesterday(monkeypatch, tmp_path):
    monkeypatch.setattr(ingest_daily.common, "ROOT", str(tmp_path))
    calls = {}
    monkeypatch.setattr(ingest_daily.search, "run",
                        lambda cfg, day, since, until: calls.setdefault(
                            "search", (day, since, until)))
    monkeypatch.setattr(ingest_daily.merge, "run",
                        lambda cfg, run_day, **kw: calls.setdefault(
                            "merge", (run_day, kw)))
    monkeypatch.setattr(ingest_daily.score, "run",
                        lambda cfg, day: calls.setdefault("score", day) or [])
    ingest_daily.run({"ingest": {"mode": "daily"}}, "2026-10-05")
    assert calls["search"] == ("2026-10-04", ">=2026-10-04", "2026-10-04")
    assert calls["merge"][0] == "2026-10-05"          # feed 选运行当天的
    assert calls["merge"][1]["window_days"] == 0
    assert calls["merge"][1]["include_run_day"] is True
    assert calls["merge"][1]["search_date"] == "2026-10-04"
    assert calls["merge"][1]["out_date"] == "2026-10-04"
    assert calls["score"] == "2026-10-04"


def test_weekly_mode_is_a_noop(monkeypatch):
    called = []
    monkeypatch.setattr(ingest_daily.common, "load_config",
                        lambda p: {"ingest": {"mode": "weekly"}})
    monkeypatch.setattr(ingest_daily, "run",
                        lambda cfg, day: called.append(day) or [])
    monkeypatch.setattr(sys, "argv", ["ingest_daily"])
    ingest_daily.main()
    assert called == []
