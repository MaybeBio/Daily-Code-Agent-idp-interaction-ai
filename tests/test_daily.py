import os, subprocess
from scripts import daily

CFG = {"topic": "demo", "title": "Demo", "daily": {"enable_repo_table": True}}

def test_collect_skips_failing_watchlist(monkeypatch, tmp_path):
    for name in ("users", "users_core", "orgs"):
        d = tmp_path / "monitor" / "lists"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{name}.txt").write_text("# comment only\n", encoding="utf-8")
    monkeypatch.setattr(daily.common, "ROOT", str(tmp_path))

    def boom(cmd, **kw):
        raise subprocess.CalledProcessError(1, ["ghresearcher"])

    monkeypatch.setattr(daily.gh, "run", boom)
    cfg = {"daily": {"watchlists": {"users": "monitor/lists/users.txt",
                                    "users_core": "monitor/lists/users_core.txt",
                                    "orgs": "monitor/lists/orgs.txt"}}}
    assert daily.collect("2026-10-04", cfg) == []

def test_build_daily_issue_table():
    table = [{"full_name": "a/b", "who": ["alice"], "event": "star",
              "stars": 5, "language": "Python", "description": "x", "url": "https://g/a/b"}]
    title, body = daily.build_daily_issue(table, raw_lines=12, cfg=CFG, date="2026-10-04")
    assert "2026-10-04" in title
    assert "a/b" in body and "| repo | stars | language | description |" in body
    assert "| 5 |" in body and "Python" in body
    assert "alice" not in body          # who/event 列已移除

def test_build_daily_issue_disabled_table():
    title, body = daily.build_daily_issue([], raw_lines=0, cfg={"topic": "d", "title": "D",
                                              "daily": {"enable_repo_table": False}}, date="2026-10-04")
    assert "| repo |" not in body

def test_daily_issue_includes_per_list_logs():
    logs = {"users": "Fetching events for target(s): alice\n\n2026-10-04 10:00:00 | ⭐️ alice starred a/b\n",
            "orgs": "", "received": ""}
    _, body = daily.build_daily_issue([], raw_lines=1, cfg=CFG, date="2026-10-04", logs=logs)
    assert "## users 原始动态" in body and "starred a/b" in body
    assert "Fetching events for target(s):" not in body      # 前缀行被去掉
    assert "## orgs 原始动态" not in body                     # 空清单整段省略
    assert "## received 原始动态" not in body

def test_daily_issue_log_lines_use_hard_breaks():
    logs = {"users": "2026-10-04 10:00:00 | ⭐️ alice starred a/b\n"
                     "2026-10-04 11:00:00 | 🚀 alice pushed to a/b\n"}
    _, body = daily.build_daily_issue([], raw_lines=2, cfg=CFG, date="2026-10-04", logs=logs)
    assert "starred a/b  \n" in body                         # 行尾双空格 = 硬换行

def test_daily_issue_size_guard_links():
    big = "2026-10-04 10:00:00 | ⭐️ a starred a/b\n" * 3000   # > 60000 字符
    logs = {"users": big, "orgs": big, "received": big}
    _, body = daily.build_daily_issue([], raw_lines=9000, cfg=CFG, date="2026-10-04",
                                      logs=logs, repo_url="https://github.com/o/r")
    assert "当日日志过大" in body
    assert "https://github.com/o/r/blob/main/monitor/users/2026/10/04.txt" in body
    assert "## users 原始动态" not in body

def test_to_web_normalizes_git_remote():
    assert daily._to_web("git@github.com:o/r.git") == "https://github.com/o/r"
