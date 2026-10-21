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
