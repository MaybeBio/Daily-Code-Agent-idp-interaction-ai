import csv
from scripts import issue

SCORED = [
    {"full_name": "a/b", "url": "u", "language": "Py", "stars": 9,
     "pushed_at": "2026-10-03", "score": 9, "one_liner": "hi"},
    {"full_name": "c/d", "url": "u2", "language": "R", "stars": 2,
     "pushed_at": "2026-10-02", "score": 5, "one_liner": "yo"},
]

def test_issue_single_table_sorted_by_score():
    title, body = issue.build_weekly_issue(SCORED, {"title": "T", "topic": "t"}, "2026-09-27", "2026-10-04")
    assert "每周仓库发现" in title and "2 个" in title
    assert "## 候选仓库（2）" in body
    # 按 score 降序:a/b(9) 在 c/d(5) 之前
    assert body.index("a/b") < body.index("c/d")

def test_issue_table_columns_match_header():
    _, body = issue.build_weekly_issue(SCORED, {"topic": "t"}, "2026-09-27", "2026-10-04")
    lines = body.splitlines()
    header = next(ln for ln in lines if ln.startswith("| repo"))
    sep = lines[lines.index(header) + 1]
    row = next(ln for ln in lines if ln.startswith("| ["))
    assert header.count("|") == row.count("|")
    assert sep.count("|") == row.count("|")

def test_csv_columns(tmp_path):
    p = tmp_path / "w.csv"
    issue.write_csv(SCORED, str(p))
    rows = list(csv.DictReader(open(p)))
    assert list(rows[0].keys()) == ["repo", "url", "language", "stars",
                                    "last_commit", "score", "one_liner"]


def test_issue_noagent_drops_ai_columns():
    scored = [
        {"full_name": "a/b", "url": "u", "language": "Py", "stars": 9,
         "pushed_at": "t", "score": None, "one_liner": "d"},
        {"full_name": "c/d", "url": "u2", "language": "R", "stars": 2,
         "pushed_at": "t", "score": None, "one_liner": "d"},
    ]
    title, body = issue.build_weekly_issue(scored, {"topic": "t"}, "2026-09-27",
                                           "2026-10-04", enabled=False)
    assert "每周仓库发现" in title
    header = next(ln for ln in body.splitlines() if ln.startswith("| repo"))
    assert "score" not in header and "one_liner" not in header
    assert "| repo | language | stars | last_commit |" in header
    assert body.index("a/b") < body.index("c/d")     # 按 stars 降序


def test_csv_noagent_drops_ai_columns(tmp_path):
    scored = [
        {"full_name": "a/b", "url": "u", "language": "Py", "stars": 9,
         "pushed_at": "t", "score": None, "one_liner": "d"},
    ]
    p = tmp_path / "w.csv"
    issue.write_csv(scored, str(p), enabled=False)
    rows = list(csv.DictReader(open(p)))
    assert list(rows[0].keys()) == ["repo", "url", "language", "stars", "last_commit"]


_MIXED = [
    {"full_name": "hi/a", "url": "u", "language": "Py", "stars": 1,
     "pushed_at": "t", "score": 9, "one_liner": "x"},
    {"full_name": "edge/b", "url": "u", "language": "Py", "stars": 1,
     "pushed_at": "t", "score": 5, "one_liner": "x"},
    {"full_name": "lo/c", "url": "u", "language": "Py", "stars": 1,
     "pushed_at": "t", "score": 4, "one_liner": "x"},
    {"full_name": "none/d", "url": "u", "language": "Py", "stars": 1,
     "pushed_at": "t", "score": None, "one_liner": "x"},
]


def test_issue_show_high_filters_by_min_score():
    cfg = {"topic": "t", "llm": {"min_score": 5}, "issue": {"show": "high"}}
    title, body = issue.build_weekly_issue(_MIXED, cfg, "2026-09-27", "2026-10-04")
    assert "2 个" in title                        # 9 和 5
    assert "hi/a" in body and "edge/b" in body
    assert "lo/c" not in body and "none/d" not in body


def test_issue_show_all_is_default():
    title, body = issue.build_weekly_issue(_MIXED, {"topic": "t"}, "2026-09-27", "2026-10-04")
    assert "4 个" in title
    for name in ("hi/a", "edge/b", "lo/c", "none/d"):
        assert name in body


def test_issue_show_high_is_noop_when_noagent():
    """noagent 没有分数,high 过滤不能把 Issue 变成空的。"""
    scored = [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 2,
               "pushed_at": "t", "score": None, "one_liner": "d"},
              {"full_name": "c/d", "url": "u", "language": "R", "stars": 1,
               "pushed_at": "t", "score": None, "one_liner": "d"}]
    cfg = {"topic": "t", "issue": {"show": "high"}}
    title, body = issue.build_weekly_issue(scored, cfg, "2026-09-27", "2026-10-04",
                                           enabled=False)
    assert "2 个" in title and "a/b" in body and "c/d" in body
