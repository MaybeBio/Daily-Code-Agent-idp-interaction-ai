import datetime as dt, os, subprocess
from scripts import archive

class _FakeR:
    def __init__(self, code, out):
        self.returncode = code
        self.stdout = out

def test_threshold_filters_and_sorts():
    scored = [
        {"full_name": "a/b", "score": 3},
        {"full_name": "c/d", "score": 9},
        {"full_name": "e/f", "score": None},
        {"full_name": "g/h", "score": 5},
    ]
    out = archive.threshold(scored, 5)
    assert [r["full_name"] for r in out] == ["c/d", "g/h"]

def test_repo_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    assert archive.repo_dir("MaybeBio/GhResearcher") == str(tmp_path / "code" / "MaybeBio__GhResearcher")

def test_run_archives_only_threshold(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    scored = [
        {"full_name": "hi/a", "score": 9, "one_liner": "x", "pushed_at": "t"},
        {"full_name": "lo/b", "score": 2, "one_liner": "y", "pushed_at": "t"},
    ]
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"), scored)
    clone_calls, wiki_calls = [], []
    out = archive.run({"llm": {"min_score": 5}}, "2026-10-05",
                      clone=lambda fn, d: clone_calls.append(fn) or True,
                      wiki=lambda s, fn, d: wiki_calls.append((s, fn)) or True)
    assert [r["full_name"] for r in out] == ["hi/a"]
    assert clone_calls == ["hi/a"]
    assert sorted(s for s, _ in wiki_calls) == ["codewiki", "deepwiki", "zread"]
    saved = archive.common.load_json(str(tmp_path / "data" / "archive" / "2026-10-05.json"))
    assert [r["full_name"] for r in saved] == ["hi/a"]
    assert saved[0]["clone_ok"] is True

def test_run_isolates_clone_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    out = archive.run({"llm": {"min_score": 5}}, "2026-10-05",
                      clone=lambda fn, d: False, wiki=lambda s, fn, d: "ok")
    assert out[0]["clone_ok"] is False
    assert out[0]["wikis"]["deepwiki"] == "ok"   # wiki 仍继续

def test_run_isolates_wiki_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    out = archive.run({"llm": {"min_score": 5}}, "2026-10-05",
                      clone=lambda fn, d: True,
                      wiki=lambda s, fn, d: "ok" if s != "zread" else "failed")
    assert out[0]["wikis"] == {"deepwiki": "ok", "codewiki": "ok", "zread": "failed"}

def test_run_overwrites_existing_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    base = str(tmp_path / "code" / "a__b")
    os.makedirs(os.path.join(base, "code"), exist_ok=True)
    open(os.path.join(base, "code", "stale.txt"), "w").write("old")
    def fake_clone(fn, dest):
        os.makedirs(dest, exist_ok=True)
        open(os.path.join(dest, "fresh.txt"), "w").write("new")
        return True
    def fake_wiki(src, fn, dest):
        os.makedirs(dest, exist_ok=True)
        return True
    archive.run({"llm": {"min_score": 5}}, "2026-10-05", clone=fake_clone, wiki=fake_wiki)
    assert not os.path.exists(os.path.join(base, "code", "stale.txt"))
    assert os.path.exists(os.path.join(base, "code", "fresh.txt"))

def test_run_defaults_min_score_when_absent(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 4, "one_liner": "x"},
                              {"full_name": "c/d", "score": 6, "one_liner": "y"}])
    out = archive.run({}, "2026-10-05", clone=lambda fn, d: True, wiki=lambda s, fn, d: True)
    assert [r["full_name"] for r in out] == ["c/d"]

def test_run_missing_scored_yields_empty(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    calls = []
    out = archive.run({"llm": {"min_score": 5}}, "2026-10-05",
                      clone=lambda fn, d: calls.append(fn) or True,
                      wiki=lambda s, fn, d: calls.append(s) or True)
    assert out == [] and calls == []

def test_clone_repo_command(monkeypatch):
    captured = {}
    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
    monkeypatch.setattr(archive.subprocess, "run", fake_run)
    assert archive.clone_repo("a/b", "/dest") is True
    assert captured["cmd"] == ["degit", "a/b", "/dest", "--force"]

def test_export_wiki_command(monkeypatch):
    captured = []
    monkeypatch.setattr(archive.subprocess, "run",
                        lambda cmd, **kw: captured.append(cmd) or None)
    archive.export_wiki("codewiki", "a/b", "/dest")
    assert captured == [["repowiki-cli", "codewiki", "cp", "a/b", "/dest"]]

def test_wiki_dest_uses_google_code_wiki_for_codewiki(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    dests = {}
    def fake_clone(fn, dest):
        return True
    def fake_wiki(src, fn, dest):
        dests[src] = dest
        return True
    archive.run({"llm": {"min_score": 5}}, "2026-10-05", clone=fake_clone, wiki=fake_wiki)
    assert dests["deepwiki"].endswith(os.path.join("wikis", "2026-10-05", "deepwiki"))
    assert dests["zread"].endswith(os.path.join("wikis", "2026-10-05", "zread"))
    assert dests["codewiki"].endswith(os.path.join("wikis", "2026-10-05", "google_code_wiki"))

def test_archive_one_reuses_stale_wiki(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    base = archive.repo_dir("a/b")
    prev = os.path.join(base, "wikis", "2026-09-28", "deepwiki")
    os.makedirs(prev, exist_ok=True)
    open(os.path.join(prev, "llms.txt"), "w").write("old")
    calls = []
    out = archive.archive_one(
        {"full_name": "a/b", "score": 8, "one_liner": "x"},
        clone=lambda fn, d: True,
        wiki=lambda s, fn, d: calls.append(s) or "ok",
        date="2026-10-05",
        stat_date=lambda s, fn: dt.date(2026, 9, 1),   # 过期 → 复用旧快照
        in_week=lambda d: False)
    assert out["wikis"]["deepwiki"] == "reused"
    assert "deepwiki" not in calls
    assert sorted(calls) == ["codewiki", "zread"]
    assert open(os.path.join(prev, "llms.txt")).read() == "old"
    assert not os.path.isdir(os.path.join(base, "wikis", "2026-10-05", "deepwiki"))

def test_archive_one_refreshes_fresh_wiki(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    base = archive.repo_dir("a/b")
    prev = os.path.join(base, "wikis", "2026-09-28", "deepwiki")
    os.makedirs(prev, exist_ok=True)
    open(os.path.join(prev, "llms.txt"), "w").write("old")
    calls = []
    out = archive.archive_one(
        {"full_name": "a/b", "score": 8, "one_liner": "x"},
        clone=lambda fn, d: True,
        wiki=lambda s, fn, d: calls.append(s) or "ok",
        date="2026-10-05",
        stat_date=lambda s, fn: dt.date(2026, 10, 6),   # 本周 → 重拉新快照
        in_week=lambda d: True)
    assert out["wikis"]["deepwiki"] == "ok"
    assert "deepwiki" in calls
    assert open(os.path.join(prev, "llms.txt")).read() == "old"   # 旧快照仍保留

def test_archive_one_reuses_when_stat_fails_but_dir_exists(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    base = archive.repo_dir("a/b")
    os.makedirs(os.path.join(base, "wikis", "2026-09-28", "zread"), exist_ok=True)
    calls = []
    out = archive.archive_one(
        {"full_name": "a/b", "score": 8, "one_liner": "x"},
        clone=lambda fn, d: True,
        wiki=lambda s, fn, d: calls.append(s) or "ok",
        date="2026-10-05",
        stat_date=lambda s, fn: None,              # stat 失败 → 保守复用旧快照
        in_week=lambda d: False)
    assert out["wikis"]["zread"] == "reused"
    assert "zread" not in calls

def test_archive_one_cp_when_no_old_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    calls = []
    out = archive.archive_one(
        {"full_name": "a/b", "score": 8, "one_liner": "x"},
        clone=lambda fn, d: True,
        wiki=lambda s, fn, d: calls.append(s) or "ok",
        date="2026-10-05",
        stat_date=lambda s, fn: dt.date(2026, 9, 1),   # 过期但无旧快照 → 仍要拉
        in_week=lambda d: False)
    assert out["wikis"]["deepwiki"] == "ok"
    assert sorted(calls) == ["codewiki", "deepwiki", "zread"]

def test_run_appends_snapshot_never_deletes_old(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    for date in ("2026-10-05", "2026-10-12"):
        archive.common.dump_json(str(tmp_path / "data" / "scored" / f"{date}.json"),
                                 [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    def fake_clone(fn, dest):
        os.makedirs(dest, exist_ok=True)
        return True
    def fake_wiki(src, fn, dest):
        os.makedirs(dest, exist_ok=True)
        open(os.path.join(dest, "page.md"), "w").write(src)
        return True
    archive.run({"llm": {"min_score": 5}}, "2026-10-05", clone=fake_clone, wiki=fake_wiki)
    archive.run({"llm": {"min_score": 5}}, "2026-10-12", clone=fake_clone, wiki=fake_wiki)
    base = archive.repo_dir("a/b")
    assert os.path.isfile(os.path.join(base, "wikis", "2026-10-05", "deepwiki", "page.md"))
    assert os.path.isfile(os.path.join(base, "wikis", "2026-10-12", "deepwiki", "page.md"))

def test_in_week():
    assert archive._in_week(dt.date(2026, 10, 5), "2026-10-08") is True   # 本周一
    assert archive._in_week(dt.date(2026, 10, 11), "2026-10-08") is True  # 本周日
    assert archive._in_week(dt.date(2026, 10, 4), "2026-10-08") is False  # 上周日
    assert archive._in_week(None, "2026-10-08") is False
    assert archive._in_week(dt.date(2026, 10, 5), "bad") is False

def test_stat_date_deepwiki(monkeypatch):
    monkeypatch.setattr(archive.subprocess, "run",
                        lambda cmd, **kw: _FakeR(0, '{"command":"stat","repo":"a/b","data":{"last_modified":"2026-10-01T00:00:00Z"}}'))
    assert archive._stat_date("deepwiki", "a/b") == dt.date(2026, 10, 1)

def test_stat_date_zread_epoch(monkeypatch):
    monkeypatch.setattr(archive.subprocess, "run",
                        lambda cmd, **kw: _FakeR(0, '{"repo":"a/b","command":"stat","data":{"updated_at": 1760000000}}'))
    assert archive._stat_date("zread", "a/b") == dt.datetime.fromtimestamp(1760000000, tz=dt.timezone.utc).date()

def test_stat_date_codewiki_resolves_sha(monkeypatch):
    def fake(cmd, **kw):
        if cmd[0] == "repowiki-cli":
            return _FakeR(0, '{"repo":"a/b","command":"stat","commit":"abc123"}')
        return _FakeR(0, "2026-10-02T03:04:05Z")
    monkeypatch.setattr(archive.subprocess, "run", fake)
    assert archive._stat_date("codewiki", "a/b") == dt.date(2026, 10, 2)

def test_clone_repo_falls_back_to_git_then_gives_up(monkeypatch, capsys):
    calls = []
    def boom(cmd, **kw):
        calls.append(cmd)
        raise RuntimeError("net")
    monkeypatch.setattr(archive.subprocess, "run", boom)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.clone_repo("a/b", "/d") is False
    assert len(calls) == 10                             # degit 5 次 + git 回退 5 次
    assert calls[0][0] == "degit" and calls[5][0] == "git"
    assert "[clone failed]" in capsys.readouterr().err  # 不再静默

def test_clone_repo_falls_back_to_git_on_degit_failure(monkeypatch):
    calls = []
    def fail_degit(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "degit":
            raise RuntimeError("degit tarball 404")
    monkeypatch.setattr(archive.subprocess, "run", fail_degit)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.clone_repo("a/b", "/d") is True
    assert calls[-1][:3] == ["git", "clone", "--depth"]

def test_clone_repo_permanent_errors_break_early(monkeypatch):
    calls = []
    def boom(cmd, **kw):
        calls.append(cmd)
        raise subprocess.CalledProcessError(1, cmd, stderr="404 Not Found")
    monkeypatch.setattr(archive.subprocess, "run", boom)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.clone_repo("a/b", "/d") is False
    assert len(calls) == 2                              # degit 1 + git 1,永久失败不再重试

def test_clone_repo_retries_then_succeeds(monkeypatch):
    n = {"i": 0}
    def flaky(cmd, **kw):
        n["i"] += 1
        if n["i"] == 1:
            raise RuntimeError("transient")
    monkeypatch.setattr(archive.subprocess, "run", flaky)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.clone_repo("a/b", "/d") is True
    assert n["i"] == 2

def test_run_archives_all_keepers_concurrently(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    scored = [{"full_name": f"o{i}/r", "score": 8, "one_liner": "x"} for i in range(6)]
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"), scored)
    out = archive.run({"llm": {"min_score": 5}, "archive": {"concurrency": 3}}, "2026-10-05",
                      clone=lambda fn, d: True, wiki=lambda s, fn, d: True)
    assert sorted(r["full_name"] for r in out) == sorted(r["full_name"] for r in scored)

def test_run_isolates_raising_archive(monkeypatch, tmp_path):
    monkeypatch.setattr(archive.common, "ROOT", str(tmp_path))
    archive.common.dump_json(str(tmp_path / "data" / "scored" / "2026-10-05.json"),
                             [{"full_name": "a/b", "score": 8, "one_liner": "x"}])
    def raise_clone(fn, dest):
        raise RuntimeError("boom")
    out = archive.run({"llm": {"min_score": 5}}, "2026-10-05", clone=raise_clone,
                      wiki=lambda s, fn, d: "ok")
    assert out[0]["clone_ok"] is False
    assert out[0]["wikis"] == {"deepwiki": "failed", "codewiki": "failed", "zread": "failed"}

def test_export_wiki_returns_ok(monkeypatch):
    monkeypatch.setattr(archive.subprocess, "run", lambda cmd, **kw: None)
    assert archive.export_wiki("deepwiki", "a/b", "/d") == "ok"

def test_export_wiki_not_indexed(monkeypatch):
    def boom(cmd, **kw):
        raise subprocess.CalledProcessError(1, cmd, stderr="404 not indexed")
    monkeypatch.setattr(archive.subprocess, "run", boom)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.export_wiki("deepwiki", "a/b", "/d") == "not_indexed"

def test_export_wiki_transient_failed(monkeypatch):
    def boom(cmd, **kw):
        raise RuntimeError("net")
    monkeypatch.setattr(archive.subprocess, "run", boom)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    assert archive.export_wiki("deepwiki", "a/b", "/d") == "failed"

def test_export_wiki_zread_submits_when_not_indexed(monkeypatch):
    def fake(cmd, **kw):
        if "cp" in cmd:
            raise subprocess.CalledProcessError(1, cmd, stderr="404 not indexed")
        return None  # submit 成功
    monkeypatch.setattr(archive.subprocess, "run", fake)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    monkeypatch.setenv("ZREAD_TOKEN", "tok")
    assert archive.export_wiki("zread", "a/b", "/d") == "submitted"

def test_export_wiki_zread_submit_skipped_without_token(monkeypatch):
    def fake(cmd, **kw):
        if "cp" in cmd:
            raise subprocess.CalledProcessError(1, cmd, stderr="404 not indexed")
        raise AssertionError("submit 不该被调用")
    monkeypatch.setattr(archive.subprocess, "run", fake)
    monkeypatch.setattr(archive.time, "sleep", lambda s: None)
    monkeypatch.delenv("ZREAD_TOKEN", raising=False)
    assert archive.export_wiki("zread", "a/b", "/d") == "not_indexed"
