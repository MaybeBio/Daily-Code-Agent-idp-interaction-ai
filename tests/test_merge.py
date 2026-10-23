from scripts import merge

def test_run_unions_dedups_and_sorts(monkeypatch, tmp_path):
    monkeypatch.setattr(merge.common, "ROOT", str(tmp_path))
    monkeypatch.setattr(merge.common, "today", lambda: "2026-10-04")
    (tmp_path / "data" / "daily").mkdir(parents=True)
    merge.common.dump_json(str(tmp_path / "data" / "daily" / "2026-10-04.repos.json"),
                           [{"full_name": "a/b", "url": "u", "stars": 1, "language": "",
                             "description": "", "event": "star", "who": ["x"]}])
    merge.common.dump_json(str(tmp_path / "data" / "search" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "stars": 9, "language": "Py",
                             "description": "d", "pushed_at": "2026-10-03"}])
    monkeypatch.setattr(merge.gh, "repo_meta", lambda fn: {"pushed_at": "2026-10-03", "stars": 9,
                                                          "language": "Py", "description": "d", "url": "u"})
    out = merge.run({"window_days": 7}, "2026-10-04")
    assert len(out) == 1                            # 去重为一条
    assert out[0]["full_name"] == "a/b"
    assert out[0]["stars"] == 9                     # search 的星星覆盖 daily 的


def test_daily_run_uses_only_run_day_feed_and_custom_dates(monkeypatch, tmp_path):
    monkeypatch.setattr(merge.common, "ROOT", str(tmp_path))
    for day, name in (("2026-10-04", "prev/prev"), ("2026-10-05", "today/feed")):
        merge.common.dump_json(str(tmp_path / "data" / "daily" / f"{day}.repos.json"),
                               [{"full_name": name, "pushed_at": day}])
    merge.common.dump_json(str(tmp_path / "data" / "search" / "2026-10-04.json"),
                           [{"full_name": "search/hit", "pushed_at": "2026-10-04"}])
    # daily 模式:只取运行当天(10-05)的 feed + 被覆盖日(10-04)的搜索
    out = merge.run({"window_days": 7}, "2026-10-05", window_days=0,
                    include_run_day=True, search_date="2026-10-04",
                    out_date="2026-10-04")
    assert {r["full_name"] for r in out} == {"today/feed", "search/hit"}
    saved = merge.common.load_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"))
    assert {r["full_name"] for r in saved} == {"today/feed", "search/hit"}
