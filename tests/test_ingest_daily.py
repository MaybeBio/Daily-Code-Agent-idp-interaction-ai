import sys
from scripts import ingest_daily


def test_daily_mode_ingests_yesterday(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(ingest_daily.common, "ROOT", str(tmp_path))
    calls = {}
    def fake_search(cfg, day, since, until):
        calls["search"] = (day, since, until)
        return []
    def fake_merge(cfg, run_day, **kw):
        calls["merge"] = (run_day, kw)
        return [{"full_name": "a/b"}]
    def fake_score(cfg, day):
        calls["score"] = day
        return [{"full_name": "a/b", "score": 7}]
    monkeypatch.setattr(ingest_daily.search, "run", fake_search)
    monkeypatch.setattr(ingest_daily.merge, "run", fake_merge)
    monkeypatch.setattr(ingest_daily.score, "run", fake_score)
    ingest_daily.run({"ingest": {"mode": "daily"}}, "2026-10-05")
    assert calls["search"] == ("2026-10-04", ">=2026-10-04", "2026-10-04")
    assert calls["merge"][0] == "2026-10-05"          # feed 选运行当天的
    assert calls["merge"][1]["window_days"] == 0
    assert calls["merge"][1]["include_run_day"] is True
    assert calls["merge"][1]["search_date"] == "2026-10-04"
    assert calls["merge"][1]["out_date"] == "2026-10-04"
    assert calls["score"] == "2026-10-04"
    out = capsys.readouterr().out                       # 分阶段日志:超时时能看出卡在哪
    assert "[ingest] search 2026-10-04: 0 (" in out
    assert "[ingest] candidates 2026-10-04: 1 (" in out
    assert "[ingest] scored 2026-10-04: 1 (" in out


def test_weekly_mode_is_a_noop(monkeypatch):
    called = []
    monkeypatch.setattr(ingest_daily.common, "load_config",
                        lambda p: {"ingest": {"mode": "weekly"}})
    monkeypatch.setattr(ingest_daily, "run",
                        lambda cfg, day: called.append(day) or [])
    monkeypatch.setattr(sys, "argv", ["ingest_daily"])
    ingest_daily.main()
    assert called == []
