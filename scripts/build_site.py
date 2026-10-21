import argparse, json, os, re, shutil
import datetime as dt
from html import unescape
from urllib.parse import urlparse
import markdown as md
import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATES = os.path.join(os.path.dirname(HERE), "templates")
WIKI_SOURCES = ["deepwiki", "codewiki", "zread"]

def load_config(path):
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise ValueError("config must be a YAML mapping")
    return cfg

def natural_week(iso):
    """由 ISO 日期反推自然周(周一~周日)。仅作 _window_of 的兜底,正常情况下不走这里。"""
    try:
        d = dt.date.fromisoformat((iso or "")[:10])
    except ValueError:
        return "", ""
    monday = d - dt.timedelta(days=d.weekday())
    return monday.isoformat(), (monday + dt.timedelta(days=6)).isoformat()

def _window_of(rec):
    """仓库归入哪个「抓取窗口」,逻辑同 paper 模板的 _window_of。

    优先用 cards.py 写入的 rec["window"](= 本次运行真实抓取的 [run-window_days, run-1]),
    这样同一次运行的仓库整块落进同一个桶;缺 window 的老数据才退回 pushed_at 所在自然周。

    ⚠️ 不能只按 pushed_at 的自然周切:一次运行的仓库会按 pushed_at 落在周日/周一被劈成两周,
    首页(取最新桶)几乎空白、归档也多出半个空周。
    """
    window = rec.get("window") or {}
    start = window.get("start") or ""
    end = window.get("end") or ""
    if start and end:
        return start, end
    return natural_week(rec.get("pushed_at") or rec.get("date"))

def site_base_path(site_base_url):
    url = (site_base_url or "").strip()
    if not url:
        return ""
    return urlparse(url).path.rstrip("/")

def score_tier(score):
    try:
        s = int(score)
    except (TypeError, ValueError):
        return "none"
    if s >= 8:
        return "high"
    if s >= 5:
        return "mid"
    return "low"

def wiki_url(source, full_name):
    if source == "deepwiki":
        return f"https://deepwiki.com/{full_name}"
    if source == "codewiki":
        return f"https://codewiki.google/github.com/{full_name}"
    if source == "zread":
        return f"https://zread.ai/{full_name}"
    return ""

def repo_dir_name(full_name):
    return full_name.replace("/", "__")

def wiki_dates_upto(out_dir, full_name, upto_date):
    """各 wiki 源在 upto_date(含)之前最近一次快照日期(code/{repo}/wikis/{date}/{source}/llms.txt)。

    以卡片自身的运行日期(upto_date=record["date"])为上限,确保周归档页链到「当周冻结」的快照,
    而不是最新的那份。
    """
    wikis_root = os.path.join(out_dir, "code", repo_dir_name(full_name), "wikis")
    result = {}
    if not os.path.isdir(wikis_root):
        return result
    for source in WIKI_SOURCES:
        dates = [d for d in os.listdir(wikis_root)
                 if (not upto_date or d <= upto_date)
                 and os.path.isfile(os.path.join(wikis_root, d, source, "llms.txt"))]
        if dates:
            result[source] = max(dates)
    return result

def wiki_snapshot_url(source, full_name, date, github_repo="", github_branch="main"):
    if not github_repo or not date:
        return ""
    return (f"https://github.com/{github_repo}/blob/{github_branch}/"
            f"code/{repo_dir_name(full_name)}/wikis/{date}/{source}/llms.txt")

def load_cards(out_dir):
    cards_root = os.path.join(out_dir, "data", "cards")
    if not os.path.isdir(cards_root):
        return []
    records = []
    for fname in sorted(os.listdir(cards_root)):
        if not fname.endswith(".json"):
            continue
        with open(os.path.join(cards_root, fname), encoding="utf-8") as f:
            rows = json.load(f)
        for r in rows:
            records.append({**r, "date": fname[:-5]})
    records.sort(key=lambda r: r.get("date") or "", reverse=True)
    return records

def load_issue_url(out_dir):
    path = os.path.join(out_dir, "data", "latest_issue.txt")
    if not os.path.isfile(path):
        return ""
    with open(path, encoding="utf-8") as f:
        return f.read().strip()

def _md_to_html(text):
    return md.markdown(text or "", extensions=["tables", "fenced_code", "sane_lists"])

_TAG_RE = re.compile(r"<[^>]+>")

def _plain_text(md_text):
    rendered = _md_to_html(md_text)
    text = _TAG_RE.sub(" ", rendered).replace("\n", " ")
    return re.sub(r"\s+", " ", unescape(text)).strip()

def _build_search_documents(records, base_path=""):
    """拆成两个索引:head(常驻,仓库名/语言/一句话/描述等轻量字段) + deep(懒加载,Code Card + README 全文)。"""
    head, deep = [], []
    for r in records:
        doc_id = r["full_name"]
        head.append({
            "id": doc_id,
            "title": r["full_name"],
            "subtitle": "",
            "url": f"{base_path}/repos/{repo_dir_name(r['full_name'])}/",
            "meta": [r.get("language") or "",
                     f"★ {r.get('stars')}" if r.get("stars") else "",
                     (r.get("pushed_at") or "")[:10]],
            "tags": [],
            "summary": r.get("one_liner") or "",
            "abstract": r.get("description") or "",
            "date": (r.get("pushed_at") or r.get("date") or "")[:10],
        })
        deep.append({"id": doc_id,
                     "deep": _plain_text((r.get("card") or "") + "\n" + (r.get("readme") or ""))})
    return head, deep

def build_site(out_dir, config=None):
    cfg = config or {}
    base_path = site_base_path(cfg.get("site_base_url") or "")
    site_title = (cfg.get("title") or cfg.get("topic") or "").strip()
    records = load_cards(out_dir)
    env = Environment(loader=FileSystemLoader(TEMPLATES), autoescape=select_autoescape(["html"]))
    env.globals["BASE"] = base_path
    env.globals["SITE_TITLE"] = site_title
    env.globals["SEARCH_PLACEHOLDER"] = (cfg.get("search_placeholder") or "").strip()
    env.globals["score_tier"] = score_tier
    env.globals["wiki_url"] = wiki_url
    env.globals["repo_dir_name"] = repo_dir_name
    env.globals["WIKI_SOURCES"] = WIKI_SOURCES
    github_repo = (cfg.get("github_repo") or "").strip()
    github_branch = (cfg.get("github_branch") or "main").strip()
    env.globals["wiki_snapshot_url"] = (
        lambda s, fn, d: wiki_snapshot_url(s, fn, d, github_repo, github_branch))
    env.filters["md"] = _md_to_html

    site_dir = os.path.join(out_dir, "site")
    data_dir = os.path.join(site_dir, "data")
    assets_dir = os.path.join(site_dir, "assets")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(assets_dir, exist_ok=True)

    # 每个仓库归入其「抓取窗口」(优先用卡片记录的 window,缺则退 pushed_at 自然周);
    # 同一个桶内同一 full_name 去重保留 pushed_at 最新。
    for r in records:
        w_start, w_end = _window_of(r)
        r["window_start"] = w_start
        r["window_end"] = w_end
        r["wiki_dates"] = wiki_dates_upto(out_dir, r["full_name"], r.get("date") or "")

    batches = {}
    for r in records:
        batches.setdefault(r["window_end"] or "unknown", []).append(r)
    for key in batches:
        latest = {}
        for r in batches[key]:
            fn = r["full_name"]
            if fn not in latest or (r.get("pushed_at") or "") > (latest[fn].get("pushed_at") or ""):
                latest[fn] = r
        batches[key] = sorted(latest.values(),
                              key=lambda r: r.get("score") if r.get("score") is not None else -1,
                              reverse=True)

    years = {}
    for key in sorted(batches, reverse=True):
        if key == "unknown":
            continue
        repos = batches[key]
        years.setdefault(key[:4], []).append(
            {"start": repos[0]["window_start"] if repos else "", "end": key, "repos": repos})

    latest_key = max((k for k in batches if k != "unknown"), default="unknown")
    this_week = batches.get(latest_key, [])
    window_start = this_week[0]["window_start"] if this_week else ""
    window_end = latest_key if latest_key != "unknown" else ""
    issue_url = load_issue_url(out_dir)

    latest_by_repo = {}
    for r in records:
        latest_by_repo.setdefault(r["full_name"], r)
    for r in latest_by_repo.values():
        page = env.get_template("code.html").render(repo=r)
        page_dir = os.path.join(site_dir, "repos", repo_dir_name(r["full_name"]))
        os.makedirs(page_dir, exist_ok=True)
        with open(os.path.join(page_dir, "index.html"), "w", encoding="utf-8") as f:
            f.write(page)

    for key, repos in batches.items():
        if key == "unknown":
            continue
        for r in repos:
            page = env.get_template("code.html").render(repo=r)
            page_dir = os.path.join(site_dir, "weeks", key, "repos", repo_dir_name(r["full_name"]))
            os.makedirs(page_dir, exist_ok=True)
            with open(os.path.join(page_dir, "index.html"), "w", encoding="utf-8") as f:
                f.write(page)
        page = env.get_template("week.html").render(
            window_start=repos[0]["window_start"] if repos else "",
            window_end=key, repos=repos)
        week_dir = os.path.join(site_dir, "weeks", key)
        os.makedirs(week_dir, exist_ok=True)
        with open(os.path.join(week_dir, "index.html"), "w", encoding="utf-8") as f:
            f.write(page)

    with open(os.path.join(site_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(env.get_template("index.html").render(
            window_start=window_start, window_end=window_end,
            repos=this_week, issue_url=issue_url))

    with open(os.path.join(site_dir, "archive.html"), "w", encoding="utf-8") as f:
        f.write(env.get_template("archive.html").render(years=years, total=len(latest_by_repo)))

    with open(os.path.join(site_dir, "search.html"), "w", encoding="utf-8") as f:
        f.write(env.get_template("search.html").render())
    head_docs, deep_docs = _build_search_documents(list(latest_by_repo.values()), base_path)
    with open(os.path.join(data_dir, "search.json"), "w", encoding="utf-8") as f:
        json.dump({"version": 1, "documents": head_docs}, f, ensure_ascii=False, indent=2)
    with open(os.path.join(data_dir, "search-deep.json"), "w", encoding="utf-8") as f:
        json.dump({"version": 1, "documents": deep_docs}, f, ensure_ascii=False, indent=2)

    shutil.copy(os.path.join(TEMPLATES, "assets", "style.css"),
                os.path.join(assets_dir, "style.css"))
    shutil.copy(os.path.join(TEMPLATES, "assets", "search.js"),
                os.path.join(assets_dir, "search.js"))

def main():
    parser = argparse.ArgumentParser(description="Build the static site from data/.")
    parser.add_argument("--out-dir", default=".", help="Repo root")
    parser.add_argument("--config", default=None, help="Path to config.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config) if args.config else None
    build_site(args.out_dir, cfg)

if __name__ == "__main__":
    main()
