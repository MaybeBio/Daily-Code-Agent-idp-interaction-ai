import argparse, csv, os, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common  # noqa: E402

# CSV 有独立 url 列;Markdown 表把 url 折进 repo 链接,故两套列。
CSV_COLS = ["repo", "url", "language", "stars", "last_commit", "score", "one_liner"]
CSV_COLS_NOAI = ["repo", "url", "language", "stars", "last_commit"]
MD_COLS = ["repo", "language", "stars", "last_commit", "score", "one_liner"]
MD_COLS_NOAI = ["repo", "language", "stars", "last_commit"]


def _cols(enabled, md):
    if md:
        return MD_COLS if enabled else MD_COLS_NOAI
    return CSV_COLS if enabled else CSV_COLS_NOAI


def _cell(r, c):
    if c == "repo":
        return r["full_name"]
    if c == "last_commit":
        return r.get("pushed_at", "")
    return r.get(c, "")


def write_csv(scored, path, enabled=True):
    cols = _cols(enabled, md=False)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in scored:
            w.writerow([_cell(r, c) for c in cols])


def _row(r, cols):
    cells = []
    for c in cols:
        if c == "repo":
            cells.append(f"[{r['full_name']}]({r['url']})")
        elif c == "score":
            cells.append(r["score"] if r.get("score") is not None else "-")
        else:
            cells.append(_cell(r, c))
    return "| " + " | ".join(str(x) for x in cells) + " |"


def _section(heading, rows, cols):
    if not rows:
        return []
    out = [f"## {heading}（{len(rows)}）", "",
           "| " + " | ".join(cols) + " |",
           "|" + "---|" * len(cols)]
    out += [_row(r, cols) for r in rows]
    out.append("")
    return out


def build_weekly_issue(scored, cfg, start, end, enabled=True):
    cols = _cols(enabled, md=True)
    if enabled:
        rows = sorted(scored, key=lambda r: (r["score"] is not None, r["score"] or 0), reverse=True)
    else:
        rows = sorted(scored, key=lambda r: r.get("stars") or 0, reverse=True)
    title = f"\U0001F50D 每周仓库发现 {end} · {len(rows)} 个"
    body = [f"# {title}", "", f"- 窗口:{start} → {end}", f"- 课题:{cfg.get('topic','')}", ""]
    body += _section("候选仓库", rows, cols)
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
    enabled = bool(cfg.get("llm", {}).get("enabled", True))
    title, body = build_weekly_issue(scored, cfg, start, end, enabled=enabled)
    write_csv(scored, common.data_dir("weekly", f"{date}.csv"), enabled)
    open(a.issue_body, "w", encoding="utf-8").write(body)
    open(a.issue_title, "w", encoding="utf-8").write(title)


if __name__ == "__main__":
    main()
