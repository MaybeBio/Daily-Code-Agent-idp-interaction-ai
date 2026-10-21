import argparse, csv, os, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common  # noqa: E402

COLS = ["repo", "url", "language", "stars", "last_commit", "score", "one_liner"]

def write_csv(scored, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for r in scored:
            w.writerow([r["full_name"], r["url"], r["language"], r["stars"],
                        r.get("pushed_at", ""), r["score"], r["one_liner"]])

def _row(r):
    return "| [{0}]({1}) | {2} | {3} | {4} | {5} | {6} |".format(
        r["full_name"], r["url"], r["language"], r["stars"], r.get("pushed_at", ""),
        r["score"] if r["score"] is not None else "-", r["one_liner"])

def _section(heading, rows):
    if not rows:
        return []
    out = [f"## {heading}（{len(rows)}）", "",
           "| repo | language | stars | last_commit | score | one_liner |",
           "|---|---|---|---|---|---|"]
    out += [_row(r) for r in rows]
    out.append("")
    return out

def build_weekly_issue(scored, cfg, start, end):
    rows = sorted(scored, key=lambda r: (r["score"] is not None, r["score"] or 0), reverse=True)
    title = f"\U0001F50D 每周仓库发现 {end} · {len(rows)} 个"
    body = [f"# {title}", "", f"- 窗口:{start} → {end}", f"- 课题:{cfg.get('topic','')}", ""]
    body += _section("候选仓库", rows)
    return title, "\n".join(body) + "\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--date", default=None)
    ap.add_argument("--issue-body", required=True)
    ap.add_argument("--issue-title", required=True)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    date = a.date or common.today()
    run_day = _date.fromisoformat(date)
    # 窗口 [run-window_days, run-1],与 search/merge/cards 同源(不含运行当天)。
    start = (run_day - timedelta(days=int(cfg.get("window_days", 7)))).isoformat()
    end = (run_day - timedelta(days=1)).isoformat()
    scored = common.load_json(common.data_dir("scored", f"{date}.json"))
    title, body = build_weekly_issue(scored, cfg, start, end)
    write_csv(scored, common.data_dir("weekly", f"{date}.csv"))
    open(a.issue_body, "w", encoding="utf-8").write(body)
    open(a.issue_title, "w", encoding="utf-8").write(title)

if __name__ == "__main__":
    main()
