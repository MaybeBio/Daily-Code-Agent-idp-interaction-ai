import argparse, os, sys, threading
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, agent, gh  # noqa: E402
from scripts.agent import load_prompts  # noqa: E402

def _score_one(row, client, prompts, model, fetch, fetch_tree):
    text = common.repo_text(row, fetch, fetch_tree)
    if not text:
        return {"score": 0, "one_liner": "无 README/描述/目录树,跳过"}
    return agent.score_repo(client, model, prompts, text)

def _wanted(row, cfg):
    """该仓库是否参与打分:llm.score_sources = search(默认) | all | watchlist。

    两者都命中的仓库 source 是 ["watchlist","search"],任一模式都会选中。
    老数据没有 source 字段时保守处理为参与打分,避免静默丢仓库。
    """
    mode = str((cfg.get("llm") or {}).get("score_sources", "search") or "search").lower()
    if mode in ("all", "any"):
        return True
    srcs = common.sources(row)
    return (mode in srcs) if srcs else True

def run(cfg, date, client=None, fetch=None, fetch_tree=None):
    fetch = fetch or gh.fetch_readme
    fetch_tree = fetch_tree or gh.fetch_tree
    cand_path = common.data_dir("candidates", f"{date}.json")
    candidates = common.load_json(cand_path) if os.path.exists(cand_path) else []
    if not bool(cfg.get("llm", {}).get("enabled", True)):
        # noagent 降级:不打分,score=None、one_liner 用描述兜底(供站点搜索摘要)。
        out = [{**row, "score": None,
                "one_liner": (row.get("description") or "").strip()}
               for row in candidates]
        common.dump_json(common.data_dir("scored", f"{date}.json"), out)
        return out
    to_score = [r for r in candidates if _wanted(r, cfg)]
    if not to_score:
        out = [{**row, "score": None,
                "one_liner": (row.get("description") or "").strip()}
               for row in candidates]
        common.dump_json(common.data_dir("scored", f"{date}.json"), out)
        return out
    client = client if client is not None else agent.make_client()
    prompts = load_prompts(os.path.join(common.ROOT, "prompts.yaml"))
    model = agent.model_name()
    done, lock = 0, threading.Lock()
    total = len(to_score)
    def work(row):
        nonlocal done
        try:
            res = _score_one(row, client, prompts, model, fetch, fetch_tree)
        except Exception as e:
            # 只留一行原因:并发下每个失败都打会让日志爆炸,但完全吞掉又没法判断
            # 是限流、超时还是模型名错。
            print(f"[score failed] {row.get('full_name')}: {type(e).__name__}: {str(e)[:120]}",
                  flush=True)
            res = {"score": None, "one_liner": "[llm failed]"}
        with lock:
            done += 1
            if done % 100 == 0 or done == total:
                print(f"[score] {done}/{total}", flush=True)
        return {**row, **res}
    with ThreadPoolExecutor(max_workers=int(cfg.get("llm", {}).get("concurrency", 8))) as ex:
        scored_rows = list(ex.map(work, to_score))
    # 未入选来源的仓库保留在 scored 里(score=None),下游 threshold 会自动排除;
    # 输出顺序仍与 candidates 一致。
    by_name = {r["full_name"]: r for r in scored_rows}
    out = []
    for row in candidates:
        hit = by_name.get(row["full_name"])
        out.append(hit if hit is not None else {
            **row, "score": None,
            "one_liner": (row.get("description") or "").strip()})
    common.dump_json(common.data_dir("scored", f"{date}.json"), out)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    print(len(run(cfg, common.today())))

if __name__ == "__main__":
    main()
