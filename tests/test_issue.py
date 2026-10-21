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
