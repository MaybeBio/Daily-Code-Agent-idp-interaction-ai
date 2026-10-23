import argparse, os, sys, time
import yaml
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, gh  # noqa: E402

FIELDS = ["fullName", "url", "language", "stargazersCount", "description", "pushedAt"]

# GitHub 搜索 API 对认证用户限流 30 请求/分钟(gh 内部还会翻页,一次查询可能发多个请求)。
# 查询集有 20+ 条时连着跑会超限,尾部查询直接 HTTP 403。留间隔把整体压到限流以下,
# 代价只是每条查询多等几秒(23 条约 1 分钟),远小于被 403 丢掉一路召回的损失。
QUERY_INTERVAL_SECONDS = 3

def _norm(d: dict) -> dict:
    return {"full_name": d.get("fullName", ""), "url": d.get("url", ""),
            "language": d.get("language") or "", "stars": d.get("stargazersCount", 0),
            "description": d.get("description") or "", "pushed_at": d.get("pushedAt", "")}

def _queries(search_cfg: dict, qcfg_path: str) -> list[str]:
    if search_cfg.get("queries"):
        return list(search_cfg["queries"])
    with open(qcfg_path, encoding="utf-8") as f:      # 兼容旧式单查询
        q = yaml.safe_load(f) or {}
    if q.get("query"):
        return [q["query"]]
    raise ValueError(f"no queries in config search.queries or {qcfg_path}")

def _bound_date(s: str) -> str:
    return (s or "").lstrip("<>=")

def run(cfg: dict, date: str, since: str, until: str | None = None) -> list[dict]:
    qcfg = os.path.join(common.ROOT, cfg["search"]["config"])
    since_date = _bound_date(since)
    until_date = _bound_date(until)
    merged: dict[str, dict] = {}
    for i, query in enumerate(_queries(cfg["search"], qcfg)):
        if i:
            time.sleep(QUERY_INTERVAL_SECONDS)
        try:
            results = gh.search_repos(qcfg, since, FIELDS, query=query)
        except Exception as e:
            print(f"[search failed] {query}: {e}", file=sys.stderr)
            continue
        for d in results:
            n = _norm(d)
            if not n["full_name"]:
                continue
            # gh search 只能按 updated_at 预过滤(--updated 与 in: 可共存),会漏进「被 star/watch
            # 顶掉 updated_at、但代码早已不 push」的旧仓库;这里再按 pushed_at(真·最后一次 push)
            # 客户端过滤,只留窗口内真推过代码的仓库。pushed: 限定符不能和 in: 共存(会被忽略),
            # 且 gh 无上界参数,窗口上界(run-1)同样只能客户端裁。
            p = (n["pushed_at"] or "")[:10]
            if since_date and p < since_date:
                continue
            if until_date and p > until_date:
                continue
            merged[n["full_name"]] = n
    rows = list(merged.values())
    common.dump_json(common.data_dir("search", f"{date}.json"), rows)
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--since", required=True)   # e.g. ">=2026-09-27"
    ap.add_argument("--until", default=None)    # e.g. "2026-10-03"(窗口上界,不含运行当天)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    print(len(run(cfg, common.today(), a.since, a.until)))

if __name__ == "__main__":
    main()
