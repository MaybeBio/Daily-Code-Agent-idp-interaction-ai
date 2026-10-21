import json, subprocess, types
import pytest
from scripts import gh

def test_search_repos_parses_json(monkeypatch):
    payload = json.dumps([{"fullName": "a/b", "stargazersCount": 3}])
    monkeypatch.setattr(gh, "run", lambda cmd, **kw: payload)
    out = gh.search_repos("q.yaml", ">=2026-09-27", ["fullName", "stargazersCount"])
    assert out[0]["fullName"] == "a/b"

def test_search_repos_parses_python_repr(monkeypatch):
    payload = "[{'fullName': 'a/b', 'stargazersCount': 3}]"
    monkeypatch.setattr(gh, "run", lambda cmd, **kw: payload)
    out = gh.search_repos("q.yaml", ">=2026-09-27", ["fullName", "stargazersCount"])
    assert out[0]["fullName"] == "a/b"

def test_safe_repo_meta_falls_back_on_error(monkeypatch):
    def boom(cmd, **kw):
        raise subprocess.CalledProcessError(1, ["gh"])
    monkeypatch.setattr(gh, "run", boom)
    assert gh.safe_repo_meta("a/b") == {"full_name": "a/b", "url": "", "language": "",
                                        "stars": 0, "description": "", "pushed_at": ""}

def test_repo_meta_maps_fields(monkeypatch):
    raw = json.dumps({"full_name": "a/b", "html_url": "u", "language": "Python",
                      "stargazers_count": 5, "description": "d", "pushed_at": "2026-10-03T00:00:00Z"})
    monkeypatch.setattr(gh, "run", lambda cmd, **kw: raw)
    m = gh.repo_meta("a/b")
    assert m == {"full_name": "a/b", "url": "u", "language": "Python",
                 "stars": 5, "description": "d", "pushed_at": "2026-10-03T00:00:00Z"}

def test_fetch_tree_keeps_code_files_and_drops_dirs_and_data(monkeypatch):
    raw = json.dumps({"tree": [
        {"type": "tree", "path": "src"},
        {"type": "blob", "path": "src/a.py"},
        {"type": "blob", "path": "data/big.h5"},
        {"type": "blob", "path": "README.md"},
    ]})
    monkeypatch.setattr(gh, "run", lambda cmd, **kw: raw)
    assert gh.fetch_tree("a/b").splitlines() == ["README.md", "src/a.py"]

def test_fetch_tree_caps_and_marks_truncation(monkeypatch):
    raw = json.dumps({"tree": [{"type": "blob", "path": f"f{i}.py"} for i in range(250)]})
    monkeypatch.setattr(gh, "run", lambda cmd, **kw: raw)
    lines = gh.fetch_tree("a/b", limit=200).splitlines()
    assert len(lines) == 201 and lines[-1] == "..."

def test_fetch_dep_files_reads_common_manifests(monkeypatch):
    def fake_run(cmd, **kw):
        endpoint = cmd[2] if len(cmd) > 2 else ""
        if "git/trees" in endpoint:
            return json.dumps({"tree": [
                {"type": "blob", "path": "requirements.txt"},
                {"type": "blob", "path": "pyproject.toml"},
                {"type": "blob", "path": "README.md"},
            ]})
        if "requirements.txt" in endpoint:
            return "numpy\npandas\n"
        if "pyproject.toml" in endpoint:
            return "[project]\nname = x"
        raise RuntimeError("unexpected")
    monkeypatch.setattr(gh, "run", fake_run)
    out = gh.fetch_dep_files("a/b")
    assert "### requirements.txt" in out and "numpy" in out
    assert "### pyproject.toml" in out and "name = x" in out
    assert "README.md" not in out          # 非依赖文件不抓

def test_fetch_dep_files_none_found(monkeypatch):
    monkeypatch.setattr(gh, "run",
                        lambda cmd, **kw: json.dumps({"tree": [{"type": "blob", "path": "src/a.py"}]}))
    assert gh.fetch_dep_files("a/b") == ""

def test_run_retries_transient_then_succeeds(monkeypatch):
    calls = {"n": 0}
    def fake_run(cmd, **kw):
        calls["n"] += 1
        if calls["n"] < 3:
            raise subprocess.CalledProcessError(1, cmd, output="", stderr="net/http: TLS handshake timeout")
        return types.SimpleNamespace(stdout="ok")
    monkeypatch.setattr(gh.subprocess, "run", fake_run)
    monkeypatch.setattr(gh.time, "sleep", lambda s: None)
    assert gh.run(["gh", "api", "x"]) == "ok"
    assert calls["n"] == 3

def test_run_gives_up_immediately_on_permanent(monkeypatch):
    calls = {"n": 0}
    def fake_run(cmd, **kw):
        calls["n"] += 1
        raise subprocess.CalledProcessError(1, cmd, output="", stderr="gh: Not Found (HTTP 404)")
    monkeypatch.setattr(gh.subprocess, "run", fake_run)
    monkeypatch.setattr(gh.time, "sleep", lambda s: None)
    with pytest.raises(subprocess.CalledProcessError):
        gh.run(["gh", "api", "x"])
    assert calls["n"] == 1
