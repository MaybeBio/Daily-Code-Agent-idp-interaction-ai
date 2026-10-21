import os, json
from scripts import common

def test_load_config_has_required_keys(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("topic: demo\ntitle: Demo\nllm:\n  concurrency: 8\n")
    cfg = common.load_config(str(p))
    assert cfg["topic"] == "demo"
    assert cfg["llm"]["concurrency"] == 8

def test_dump_json_creates_parents(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "ROOT", str(tmp_path))
    path = common.data_dir("daily", "x.json")
    common.dump_json(path, {"k": 1})
    assert json.load(open(path))["k"] == 1
