from scripts import aggregate


def _scored(tmp_path, day, rows):
    aggregate.common.dump_json(str(tmp_path / "data" / "scored" / f"{day}.json"), rows)


def test_aggregates_window_and_dedups_latest(monkeypatch, tmp_path):
    monkeypatch.setattr(aggregate.common, "ROOT", str(tmp_path))
    _scored(tmp_path, "2026-10-01", [
        {"full_name": "a/b", "pushed_at": "2026-10-01", "score": 8, "stars": 1},
        {"full_name": "c/d", "pushed_at": "2026-10-01", "score": 6, "stars": 2},
    ])
    _scored(tmp_path, "2026-10-02", [
        {"full_name": "a/b", "pushed_at": "2026-10-02", "score": 7, "stars": 1},
    ])
    # window_days=3 → [10-01, 10-02, 10-03];10-03 缺失,这里不补跑
    out = aggregate.run({"window_days": 3, "ingest": {"max_backfill_days": 0}},
                        "2026-10-04", backfill_limit=0)
    assert [r["full_name"] for r in out] == ["a/b", "c/d"]    # 按 score 降序
    assert out[0]["pushed_at"] == "2026-10-02"                # 保留最新一次 push 的打分
    saved = aggregate.common.load_json(str(tmp_path / "data" / "scored" / "2026-10-04.json"))
    assert len(saved) == 2


def test_missing_day_triggers_backfill(monkeypatch, tmp_path):
    monkeypatch.setattr(aggregate.common, "ROOT", str(tmp_path))
    _scored(tmp_path, "2026-10-01", [{"full_name": "a/b", "pushed_at": "2026-10-01",
                                      "score": 8, "stars": 0}])
    calls = {}

    def fake_search(cfg, day, since, until):
        calls["search"] = (day, since, until)

    def fake_merge(cfg, run_day, **kw):
        calls["merge"] = (run_day, kw)

    def fake_score(cfg, day):
        calls["score"] = day
        _scored(tmp_path, day, [{"full_name": "e/f", "pushed_at": day,
                                 "score": 9, "stars": 0}])
        return []

    monkeypatch.setattr(aggregate.search, "run", fake_search)
    monkeypatch.setattr(aggregate.merge, "run", fake_merge)
    monkeypatch.setattr(aggregate.score, "run", fake_score)
    # window_days=2 → [10-01, 10-02];10-02 缺失 → 必须补跑,补跑结果参与聚合
    out = aggregate.run({"window_days": 2, "ingest": {"max_backfill_days": 2}}, "2026-10-03")
    assert calls["search"] == ("2026-10-02", ">=2026-10-02", "2026-10-02")
    assert calls["merge"][0] == "2026-10-03"
    assert calls["merge"][1]["window_days"] == 0 and calls["merge"][1]["include_run_day"]
    assert calls["merge"][1]["search_date"] == "2026-10-02"
    assert calls["merge"][1]["out_date"] == "2026-10-02"
    assert calls["score"] == "2026-10-02"
    assert {r["full_name"] for r in out} == {"a/b", "e/f"}


def test_failed_backfill_does_not_abort_aggregation(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(aggregate.common, "ROOT", str(tmp_path))
    _scored(tmp_path, "2026-10-01", [{"full_name": "a/b", "pushed_at": "2026-10-01",
                                      "score": 8, "stars": 0}])
    monkeypatch.setattr(aggregate.search, "run", lambda *a, **k: None)
    monkeypatch.setattr(aggregate.merge, "run", lambda *a, **k: None)

    def boom(*a, **k):
        raise RuntimeError("llm down")

    monkeypatch.setattr(aggregate.score, "run", boom)
    out = aggregate.run({"window_days": 2, "ingest": {"max_backfill_days": 1}}, "2026-10-03")
    assert [r["full_name"] for r in out] == ["a/b"]        # 已有的那天仍被聚合
    err = capsys.readouterr().err
    assert "backfill 2026-10-02 failed: llm down" in err
    assert "still missing: 2026-10-02" in err


def test_backfill_is_capped_and_warns(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(aggregate.common, "ROOT", str(tmp_path))
    calls = []
    monkeypatch.setattr(aggregate.search, "run", lambda *a, **k: calls.append(a[1]))
    monkeypatch.setattr(aggregate.merge, "run", lambda *a, **k: None)

    def fake_score(cfg, day):
        _scored(tmp_path, day, [])      # 补跑会产出当天 scored 文件(即便是空)
        return []

    monkeypatch.setattr(aggregate.score, "run", fake_score)
    # window_days=3,全缺;上限 1 → 只补最早一天,其余告警
    out = aggregate.run({"window_days": 3, "ingest": {"max_backfill_days": 1}}, "2026-10-04")
    assert calls == ["2026-10-01"]
    assert out == []
    assert "still missing: 2026-10-02, 2026-10-03" in capsys.readouterr().err
