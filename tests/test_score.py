from scripts import score

CFG = {"topic": "demo", "title": "D", "llm": {"concurrency": 2}}
PROMPTS = {"score": {"system": "s", "user": "{readme}"}}

def test_run_scores_and_writes(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "2026-10-03",
                             "stars": 1, "language": "", "description": ""}])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    monkeypatch.setattr(score.agent, "score_repo",
                        lambda *a, **k: {"score": 9, "one_liner": "nice"})
    rows = score.run(CFG, "2026-10-04", client=object(), fetch=lambda fn: "readme text")
    assert rows[0]["score"] == 9 and rows[0]["one_liner"] == "nice"
    saved = score.common.load_json(str(tmp_path / "data" / "scored" / "2026-10-04.json"))
    assert saved[0]["score"] == 9

def test_fetch_failure_without_description_or_tree_is_skipped(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "2026-10-03",
                             "stars": 1, "language": "", "description": ""}])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    called = {"n": 0}
    def should_not_run(*a, **k):
        called["n"] += 1
        return {"score": 9, "one_liner": "x"}
    monkeypatch.setattr(score.agent, "score_repo", should_not_run)
    def boom(fn):
        raise RuntimeError("404")
    rows = score.run(CFG, "2026-10-04", client=object(), fetch=boom, fetch_tree=boom)
    assert rows[0]["score"] == 0
    assert rows[0]["one_liner"] == "无 README/描述/目录树,跳过"
    assert called["n"] == 0          # 没有可读内容,不该花一次 LLM 调用

def test_fetch_failure_falls_back_to_description(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "t",
                             "stars": 1, "language": "Python", "description": "IDP phase separation tool"}])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    seen = {}
    def capture(client, model, prompts, text):
        seen["text"] = text
        return {"score": 8, "one_liner": "n"}
    monkeypatch.setattr(score.agent, "score_repo", capture)
    def boom(fn):
        raise RuntimeError("404")
    rows = score.run(CFG, "2026-10-04", client=object(), fetch=boom, fetch_tree=lambda fn: "")
    assert "IDP phase separation tool" in seen["text"]     # 用 description 兜底喂 LLM
    assert rows[0]["score"] == 8

def test_fetch_failure_falls_back_to_tree(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "t",
                             "stars": 1, "language": "", "description": ""}])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    seen = {}
    def capture(client, model, prompts, text):
        seen["text"] = text
        return {"score": 6, "one_liner": "n"}
    monkeypatch.setattr(score.agent, "score_repo", capture)
    def boom(fn):
        raise RuntimeError("404")
    rows = score.run(CFG, "2026-10-04", client=object(), fetch=boom,
                     fetch_tree=lambda fn: "scripts/calvados_functions.py")
    assert "calvados_functions.py" in seen["text"]         # 无描述也能靠目录树打分
    assert rows[0]["score"] == 6

def test_llm_failure_marked_as_llm_failed(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "t",
                             "stars": 1, "language": "", "description": ""}])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    def boom(*a, **k):
        raise RuntimeError("llm down")
    monkeypatch.setattr(score.agent, "score_repo", boom)
    rows = score.run(CFG, "2026-10-04", client=object(), fetch=lambda fn: "readme")
    assert rows[0]["score"] is None and rows[0]["one_liner"] == "[llm failed]"


def test_run_noagent_passthrough(monkeypatch, tmp_path):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [{"full_name": "a/b", "url": "u", "pushed_at": "t",
                             "stars": 1, "language": "", "description": "IDP tool"},
                            {"full_name": "c/d", "url": "u", "pushed_at": "t",
                             "stars": 2, "language": "", "description": ""}])
    called = {"n": 0}
    def should_not_run(*a, **k):
        called["n"] += 1
        return {"score": 9, "one_liner": "x"}
    monkeypatch.setattr(score.agent, "score_repo", should_not_run)
    rows = score.run({"llm": {"enabled": False}, "topic": "t"}, "2026-10-04",
                     client=object(), fetch=lambda fn: "readme")
    assert called["n"] == 0                         # 不打分,不花 LLM 调用
    assert rows[0]["score"] is None
    assert rows[0]["one_liner"] == "IDP tool"       # description 兜底
    assert rows[1]["one_liner"] == ""
    saved = score.common.load_json(str(tmp_path / "data" / "scored" / "2026-10-04.json"))
    assert saved[0]["score"] is None


def _cand(name, source):
    return {"full_name": name, "url": "u", "pushed_at": "t", "stars": 1,
            "language": "", "description": "", "source": source}


def _run_sources(monkeypatch, tmp_path, score_sources=None):
    monkeypatch.setattr(score.common, "ROOT", str(tmp_path))
    score.common.dump_json(str(tmp_path / "data" / "candidates" / "2026-10-04.json"),
                           [_cand("s/1", ["search"]),
                            _cand("w/2", ["watchlist"]),
                            _cand("b/3", ["watchlist", "search"])])
    monkeypatch.setattr(score, "load_prompts", lambda p: PROMPTS)
    seen = []

    def fake(client, model, prompts, text):
        seen.append(text)
        return {"score": 8, "one_liner": "n"}

    monkeypatch.setattr(score.agent, "score_repo", fake)
    cfg = {**CFG, "llm": dict(CFG["llm"])}
    if score_sources:
        cfg["llm"]["score_sources"] = score_sources
    rows = score.run(cfg, "2026-10-04", client=object(), fetch=lambda fn: fn)
    return {r["full_name"]: r for r in rows}, seen, rows


def test_run_scores_only_search_when_configured(monkeypatch, tmp_path):
    by, seen, rows = _run_sources(monkeypatch, tmp_path, "search")
    assert by["s/1"]["score"] == 8
    assert by["b/3"]["score"] == 8       # 两边都命中,含 search 照样打
    assert by["w/2"]["score"] is None    # 仅 watchlist 被跳过
    assert len(seen) == 2
    assert [r["full_name"] for r in rows] == ["s/1", "w/2", "b/3"]   # 顺序不变


def test_run_default_is_search(monkeypatch, tmp_path):
    by, seen, _ = _run_sources(monkeypatch, tmp_path)
    assert by["s/1"]["score"] == 8
    assert by["b/3"]["score"] == 8       # 含 search
    assert by["w/2"]["score"] is None    # 默认不打量仅 watchlist 的
    assert len(seen) == 2


def test_run_scores_all_when_configured(monkeypatch, tmp_path):
    by, seen, _ = _run_sources(monkeypatch, tmp_path, "all")
    assert all(by[n]["score"] == 8 for n in ("s/1", "w/2", "b/3"))
    assert len(seen) == 3


def test_run_scores_only_watchlist_when_configured(monkeypatch, tmp_path):
    by, seen, _ = _run_sources(monkeypatch, tmp_path, "watchlist")
    assert by["w/2"]["score"] == 8
    assert by["b/3"]["score"] == 8
    assert by["s/1"]["score"] is None
    assert len(seen) == 2
