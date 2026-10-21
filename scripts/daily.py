import argparse, os, subprocess, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, logs, gh  # noqa: E402

# 每日三个清单:展示名 → monitor 子目录(与 collect() 的 subdir 映射一致)
SECTIONS = (("users", "users"), ("orgs", "orgs"), ("received", "received"))
_MAX_BODY = 60000                     # Issue 正文超此长度则日志改为链接
_PREAMBLE = "Fetching events for target(s):"

def _log_path(date: str, sub: str) -> str:
    return os.path.join(common.ROOT, "monitor", sub,
                        date[:4], date[5:7], f"{date[8:10]}.txt")

def _clean_log(text: str) -> str:
    keep = [ln for ln in text.splitlines()
            if ln.strip() and not ln.startswith(_PREAMBLE)]
    # 行尾双空格 = Markdown 硬换行:GitHub 自带该行为,但本地预览器等不认单换行,
    # 加双空格保证任何渲染器都逐行显示(不会挤成一段)。
    return "  \n".join(keep)

def _to_web(url: str) -> str:
    url = (url or "").strip()
    if url.startswith("git@github.com:"):
        url = "https://github.com/" + url[len("git@github.com:"):]
    if url.endswith(".git"):
        url = url[:-4]
    return url.rstrip("/")

def _repo_web_url(cfg: dict) -> str:
    repo = os.environ.get("GITHUB_REPOSITORY")
    if repo:
        return f"https://github.com/{repo}"
    url = _to_web(cfg.get("repo_url", ""))
    if url:
        return url
    try:
        url = subprocess.run(["git", "remote", "get-url", "origin"],
                             capture_output=True, text=True, check=True).stdout
    except Exception:
        url = ""
    return _to_web(url)

def _log_links(date: str, repo_url: str) -> list[str]:
    out = ["当日日志过大,未内嵌正文,请查看原文件:", ""]
    for name, sub in SECTIONS:
        out.append(f"- [{name}]({repo_url}/blob/main/monitor/{sub}/"
                   f"{date[:4]}/{date[5:7]}/{date[8:10]}.txt)")
    out.append("")
    return out

def collect(date: str, cfg: dict) -> list[dict]:
    root = common.ROOT
    wl = cfg["daily"]["watchlists"]
    y, m = date[:4], date[5:7]
    since = (_date.fromisoformat(date) - timedelta(days=1)).isoformat()  # ghresearcher 需 YYYY-MM-DD
    table: dict = {}
    for key, is_org, received in (("users", False, False), ("orgs", True, False),
                                  ("users_core", False, True)):
        list_path = os.path.join(root, wl[key])
        if not os.path.exists(list_path):
            continue
        subdir = {"users": "users", "orgs": "orgs", "users_core": "received"}[key]
        out_dir = os.path.join(root, "monitor", subdir, y, m)
        os.makedirs(out_dir, exist_ok=True)
        cmd = ["ghresearcher", "monitor", "-f", list_path, "--since", since,
               "--expand-commits"]
        if is_org:
            cmd.append("--org")
        if received:
            cmd += ["-r", "-l", "30"]
        try:
            text = gh.run(cmd)
        except subprocess.CalledProcessError as e:
            print(f"[daily] watchlist {key} skipped: {e}", file=sys.stderr)
            continue
        with open(os.path.join(out_dir, f"{date[8:10]}.txt"), "w", encoding="utf-8") as f:
            f.write(text)
        for full_name, info in logs.repos_from_events(logs.extract_events(text)).items():
            e = table.setdefault(full_name, {"full_name": full_name, "who": [], "event": ""})
            e["who"] = sorted(set(e["who"]) | set(info["who"]))
            e["event"] = info["kinds"][0] if len(info["kinds"]) == 1 else "mixed"
    rows = []
    for full_name, e in table.items():
        meta = gh.safe_repo_meta(full_name)
        rows.append({**e, "stars": meta["stars"], "language": meta["language"],
                     "description": meta["description"], "url": meta["url"]})
    rows.sort(key=lambda r: r["stars"], reverse=True)
    return rows

def build_daily_issue(repo_table: list[dict], raw_lines: int, cfg: dict, date: str,
                      logs: dict | None = None, repo_url: str = "") -> tuple[str, str]:
    logs = logs or {}
    title = f"\U0001F4E1 每日动态 {date} · {cfg.get('title', '')}".strip()
    prefix = [f"# {title}", "",
              f"- 原始动态行数:{raw_lines}",
              f"- 候选仓库:{len(repo_table)}", ""]
    if cfg.get("daily", {}).get("enable_repo_table", True) and repo_table:
        prefix += ["## 今日候选仓库(免费元数据,未做 LLM)", "",
                   "| repo | stars | language | description |",
                   "|---|---|---|---|"]
        for r in repo_table:
            prefix.append("| [{0}]({1}) | {2} | {3} | {4} |".format(
                r["full_name"], r["url"], r["stars"], r["language"],
                (r["description"] or "").replace("|", "\\|")))
        prefix.append("")
    sections = []
    for name, _sub in SECTIONS:
        cleaned = _clean_log(logs.get(name, ""))
        if cleaned:                          # 空清单整段省略(与 AI4Bio 一致)
            sections += [f"## {name} 原始动态", cleaned, ""]
    text = "\n".join(prefix + sections) + "\n"
    if len(text) > _MAX_BODY and repo_url:
        text = "\n".join(prefix + _log_links(date, repo_url)) + "\n"
    return title, text

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--out-dir", default=common.ROOT)
    ap.add_argument("--issue-body", required=True)
    ap.add_argument("--issue-title", required=True)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    date = common.today()
    table = collect(date, cfg)
    logs_by_name, raw_lines = {}, 0
    for name, sub in SECTIONS:
        p = _log_path(date, sub)
        content = ""
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                content = fh.read()
        logs_by_name[name] = content
        raw_lines += sum(1 for ln in content.splitlines() if ln.strip())
    common.dump_json(common.data_dir("daily", f"{date}.repos.json"), table)
    title, body = build_daily_issue(table, raw_lines=raw_lines, cfg=cfg, date=date,
                                    logs=logs_by_name, repo_url=_repo_web_url(cfg))
    open(a.issue_body, "w", encoding="utf-8").write(body)
    open(a.issue_title, "w", encoding="utf-8").write(title)

if __name__ == "__main__":
    main()
