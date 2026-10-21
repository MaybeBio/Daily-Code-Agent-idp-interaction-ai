from scripts import search

def _qcfg(tmp_path, text="item_type: repos\nlimit: 500\n"):
    qdir = tmp_path / "queries"
    qdir.mkdir(exist_ok=True)
    (qdir / "search_repos.yaml").write_text(text)
    return "queries/search_repos.yaml"

def test_run_unions_and_dedups_queries(monkeypatch, tmp_path):
    monkeypatch.setattr(search.common, "ROOT", str(tmp_path))
    cfgpath = _qcfg(tmp_path)
    def fake_search(cfg_path, updated, fields, query=None):
        assert query in ("A in:name", "B in:name")
        base = [{"fullName": "a/b", "url": "u", "language": "Python",
                 "stargazersCount": 7, "description": "d", "pushedAt": "2026-10-03T00:00:00Z"}]
        if query == "B in:name":
            base.append({"fullName": "c/d", "url": "u2", "language": "Go",
                         "stargazersCount": 1, "description": "", "pushedAt": "2026-10-02T00:00:00Z"})
        return base
    monkeypatch.setattr(search.gh, "search_repos", fake_search)
    cfg = {"search": {"config": cfgpath, "queries": ["A in:name", "B in:name"]}}
    rows = search.run(cfg, "2026-10-04", ">=2026-09-27")
    assert sorted(r["full_name"] for r in rows) == ["a/b", "c/d"]   # a/b deduped across queries
    assert (tmp_path / "data" / "search" / "2026-10-04.json").exists()

def test_run_falls_back_to_single_query_in_yaml(monkeypatch, tmp_path):
    monkeypatch.setattr(search.common, "ROOT", str(tmp_path))
    cfgpath = _qcfg(tmp_path, "item_type: repos\nquery: 'X in:name'\n")
    monkeypatch.setattr(search.gh, "search_repos",
                        lambda *a, **k: [{"fullName": "a/b", "url": "u", "language": "",
                                          "stargazersCount": 1, "description": "",
                                          "pushedAt": "2026-10-03T00:00:00Z"}])
    rows = search.run({"search": {"config": cfgpath}}, "2026-10-04", ">=2026-09-27")
    assert rows[0]["full_name"] == "a/b" and rows[0]["stars"] == 1

def test_run_filters_out_stale_pushed(monkeypatch, tmp_path):
    monkeypatch.setattr(search.common, "ROOT", str(tmp_path))
    cfgpath = _qcfg(tmp_path)
    def fake_search(cfg_path, updated, fields, query=None):
        return [
            {"fullName": "a/b", "url": "u", "language": "Py", "stargazersCount": 7,
             "description": "d", "pushedAt": "2026-10-03T00:00:00Z"},   # 窗口内 → 保留
            {"fullName": "old/x", "url": "u2", "language": "Py", "stargazersCount": 9,
             "description": "", "pushedAt": "2025-04-07T00:00:00Z"},     # 旧(被 star 顶 updated) → 丢弃
            {"fullName": "empty/y", "url": "u3", "language": "Py", "stargazersCount": 1,
             "description": "", "pushedAt": ""},                          # 无 pushed_at → 丢弃
        ]
    monkeypatch.setattr(search.gh, "search_repos", fake_search)
    cfg = {"search": {"config": cfgpath, "queries": ["A in:name"]}}
    rows = search.run(cfg, "2026-10-04", ">=2026-09-27")
    assert [r["full_name"] for r in rows] == ["a/b"]

def test_run_skips_failed_query(monkeypatch, tmp_path):
    monkeypatch.setattr(search.common, "ROOT", str(tmp_path))
    cfgpath = _qcfg(tmp_path)
    def fake_search(cfg_path, updated, fields, query=None):
        if query == "bad in:name":
            raise ValueError("unknown flag")
        return [{"fullName": "a/b", "url": "u", "language": "Python",
                 "stargazersCount": 7, "description": "d", "pushedAt": "2026-10-03T00:00:00Z"}]
    monkeypatch.setattr(search.gh, "search_repos", fake_search)
    cfg = {"search": {"config": cfgpath, "queries": ["bad in:name", "good in:name"]}}
    rows = search.run(cfg, "2026-10-04", ">=2026-09-27")
    assert [r["full_name"] for r in rows] == ["a/b"]   # 坏 query 跳过,好 query 仍出结果
