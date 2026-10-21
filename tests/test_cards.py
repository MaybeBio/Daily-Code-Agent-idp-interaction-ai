import types
from scripts import cards, agent

PROMPTS = {"code_card": {"system": "s", "user": "{readme}"}}

def _fake_client(content):
    msg = types.SimpleNamespace(content=content)
    resp = types.SimpleNamespace(choices=[types.SimpleNamespace(message=msg)])
    return types.SimpleNamespace(chat=types.SimpleNamespace(
        completions=types.SimpleNamespace(create=lambda **kw: resp)))

def test_card_record_shape():
    row = {"full_name": "a/b", "url": "https://github.com/a/b", "language": "Python",
           "stars": 5, "pushed_at": "t", "score": 8, "one_liner": "x"}
    out = cards.card_record(row, "## 是什么\n...")
    assert out == {"full_name": "a/b", "url": "https://github.com/a/b", "language": "Python",
                   "stars": 5, "pushed_at": "t", "score": 8,
                   "one_liner": "x", "card": "## 是什么\n..."}

def test_run_cards_only_keepers(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    scored = [
        {"full_name": "hi/a", "url": "u", "language": "Python", "stars": 1,
         "score": 9, "one_liner": "x"},
        {"full_name": "lo/b", "url": "u", "language": "Go", "stars": 1,
         "score": 2, "one_liner": "y"},
    ]
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"), scored)
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=lambda fn: "README", fetch_tree=lambda fn: "",
                    card_fn=lambda c, m, p, readme: {"card": "card:" + readme})
    assert [r["full_name"] for r in out] == ["hi/a"]
    saved = cards.common.load_json(str(tmp_path / "data" / "cards" / "2026-10-05.json"))
    assert saved[0]["card"].startswith("card:README")
    assert saved[0]["readme"] == "README"
    assert saved[0]["window"] == {"start": "2026-09-28", "end": "2026-10-04"}   # run-7 → run-1
    assert set(saved[0]) == {"full_name", "url", "language", "stars", "pushed_at",
                             "score", "one_liner", "card", "readme", "window"}

def test_run_card_failure_is_isolated(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    def boom(fn):
        raise RuntimeError("fetch failed")
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=boom, fetch_tree=boom,
                    card_fn=lambda c, m, p, r: {"card": "x"})
    assert out[0]["card"] == ""

def test_run_falls_back_to_description(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "description": "IDP phase separation tool",
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    seen = {}
    def capture(c, m, p, text):
        seen["text"] = text
        return {"card": "C"}
    def boom(fn):
        raise RuntimeError("404")
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=boom, fetch_tree=lambda fn: "", card_fn=capture)
    assert "IDP phase separation tool" in seen["text"]     # README 缺失用 description 兜底
    assert out[0]["card"] == "C"

def test_run_falls_back_to_tree(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    seen = {}
    def capture(c, m, p, text):
        seen["text"] = text
        return {"card": "C"}
    def boom(fn):
        raise RuntimeError("404")
    cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
              client=object(), fetch=boom,
              fetch_tree=lambda fn: "scripts/calvados_functions.py",
              fetch_deps=lambda fn: "", card_fn=capture)
    assert "calvados_functions.py" in seen["text"]         # 无描述也能靠目录树出卡片

def test_run_no_text_skips_llm(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    called = {"n": 0}
    def should_not_run(c, m, p, text):
        called["n"] += 1
        return {"card": "x"}
    def boom(fn):
        raise RuntimeError("404")
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=boom, fetch_tree=boom, card_fn=should_not_run)
    assert out[0]["card"] == "" and called["n"] == 0       # 无内容不花 LLM 调用

def test_run_defaults_min_score_when_absent(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 4, "one_liner": "x"},
                            {"full_name": "c/d", "url": "u", "language": "Py", "stars": 1,
                             "score": 6, "one_liner": "y"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    out = cards.run({"topic": "t"}, "2026-10-05",
                    client=object(), fetch=lambda fn: "r", fetch_tree=lambda fn: "",
                    card_fn=lambda c, m, p, r: {"card": "x"})
    assert [r["full_name"] for r in out] == ["c/d"]

def test_run_missing_scored_yields_empty(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    calls = []
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=lambda fn: calls.append(fn) or "r",
                    card_fn=lambda c, m, p, r: {"card": "x"})
    assert out == [] and calls == []

def test_run_feeds_readme_tree_and_deps(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    seen = {}
    def capture(c, m, p, text):
        seen["text"] = text
        return {"card": "C"}
    cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
              client=object(), fetch=lambda fn: "READMEME",
              fetch_tree=lambda fn: "src/a.py",
              fetch_deps=lambda fn: "### requirements.txt\nnumpy",
              card_fn=capture)
    t = seen["text"]
    assert "READMEME" in t and "src/a.py" in t and "requirements.txt" in t

def test_run_stores_full_readme_not_capped(monkeypatch, tmp_path):
    monkeypatch.setattr(cards.common, "ROOT", str(tmp_path))
    cards.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                           [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "score": 8, "one_liner": "x"}])
    monkeypatch.setattr(cards, "load_prompts", lambda p: PROMPTS)
    big = "x" * (cards._README_CAP + 500)
    seen = {}
    def capture(c, m, p, text):
        seen["text"] = text
        return {"card": "C"}
    out = cards.run({"llm": {"min_score": 5}, "topic": "t"}, "2026-10-05",
                    client=object(), fetch=lambda fn: big,
                    fetch_tree=lambda fn: "", card_fn=capture)
    assert out[0]["readme"] == big                       # 存盘 = 全文
    assert "已截断" in seen["text"]                       # prompt 仍按上限截断

def test_build_code_card_parses_json():
    client = _fake_client('{"card": "## 是什么\\nhello"}')
    out = agent.build_code_card(client, "m", PROMPTS, "README")
    assert out == {"card": "## 是什么\nhello"}
