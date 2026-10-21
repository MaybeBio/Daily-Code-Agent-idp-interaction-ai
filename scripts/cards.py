import argparse, datetime as dt, os, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, agent, gh, archive  # noqa: E402
from scripts.agent import load_prompts  # noqa: E402

def card_record(row, card, readme="", window=None):
    rec = {"full_name": row["full_name"], "url": row.get("url", ""),
            "language": row.get("language", ""), "stars": row.get("stars", 0),
            "pushed_at": row.get("pushed_at", ""),
            "score": row["score"], "one_liner": row.get("one_liner", ""), "card": card}
    if readme:
        rec["readme"] = readme
    if window:
        rec["window"] = window
    return rec

def _window(date, window_days):
    """本次运行抓取的窗口 [run-window_days, run-1],与 paper 同源(不含运行当天)。

    站点(build_site._window_of)靠它把一次运行的仓库整块归入一个桶,而不是按 pushed_at 的自然周。
    """
    try:
        d = dt.date.fromisoformat((date or "")[:10])
    except ValueError:
        return None
    return {"start": (d - dt.timedelta(days=window_days)).isoformat(),
            "end": (d - dt.timedelta(days=1)).isoformat()}

# README 全文喂给 card,但超大 README(awesome-list 类可上万行)设一个宽松上限,防止撑爆上下文/成本。
_README_CAP = 40000

def _card_context(row, fetch, fetch_tree, fetch_deps):
    try:
        readme = fetch(row["full_name"]) or ""
    except Exception:
        readme = ""
    full_readme = readme.strip()
    readme = full_readme
    if len(readme) > _README_CAP:
        readme = readme[:_README_CAP] + "\n...(README 过长,已截断)"
    try:
        tree = fetch_tree(row["full_name"]) or ""
    except Exception:
        tree = ""
    try:
        deps = fetch_deps(row["full_name"]) if tree else ""
    except Exception:
        deps = ""
    parts = []
    if readme:
        parts.append(f"README:\n{readme}")
    elif (row.get("description") or "").strip():
        parts.append(f"描述:{row['description'].strip()}")
    if tree:
        parts.append(f"目录树:\n{tree}")
    if deps:
        parts.append(f"依赖/构建文件:\n{deps}")
    return "\n\n".join(parts), full_readme

def run(cfg, date, client=None, fetch=None, fetch_tree=None, fetch_deps=None, card_fn=None):
    fetch = fetch or gh.fetch_readme
    fetch_tree = fetch_tree or gh.fetch_tree
    fetch_deps = fetch_deps or gh.fetch_dep_files
    client = client if client is not None else agent.make_client()
    card_fn = card_fn or agent.build_code_card
    prompts = load_prompts(os.path.join(common.ROOT, "prompts.yaml"))
    model = agent.model_name()
    min_score = int(cfg.get("llm", {}).get("min_score", 5))
    window = _window(date, int(cfg.get("window_days", 7)))
    scored_path = common.data_dir("scored", f"{date}.json")
    scored = common.load_json(scored_path) if os.path.exists(scored_path) else []
    keepers = archive.threshold(scored, min_score)
    def work(row):
        readme = ""
        try:
            text, readme = _card_context(row, fetch, fetch_tree, fetch_deps)
            if not text.strip():
                return card_record(row, "", readme=readme, window=window)
            res = card_fn(client, model, prompts, text)
            card = res.get("card", "")
        except Exception as e:
            print(f"[card failed] {row['full_name']}: {e}", file=sys.stderr)
            card = ""
        return card_record(row, card, readme=readme, window=window)
    with ThreadPoolExecutor(max_workers=int(cfg.get("llm", {}).get("concurrency", 8))) as ex:
        out = list(ex.map(work, keepers))
    common.dump_json(common.data_dir("cards", f"{date}.json"), out)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--date", default=None)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    date = a.date or common.today()
    out = run(cfg, date)
    print(f"{sum(1 for r in out if r['card'])}/{len(out)} cards")

if __name__ == "__main__":
    main()
