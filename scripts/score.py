import argparse, os, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, agent, gh  # noqa: E402
from scripts.agent import load_prompts  # noqa: E402

def _score_one(row, client, prompts, model, fetch, fetch_tree):
    text = common.repo_text(row, fetch, fetch_tree)
    if not text:
        return {"score": 0, "one_liner": "无 README/描述/目录树,跳过"}
    return agent.score_repo(client, model, prompts, text)

def run(cfg, date, client=None, fetch=None, fetch_tree=None):
    fetch = fetch or gh.fetch_readme
    fetch_tree = fetch_tree or gh.fetch_tree
    client = client if client is not None else agent.make_client()
    prompts = load_prompts(os.path.join(common.ROOT, "prompts.yaml"))
    model = agent.model_name()
    cand_path = common.data_dir("candidates", f"{date}.json")
    candidates = common.load_json(cand_path) if os.path.exists(cand_path) else []
    out = []
    def work(row):
        try:
            res = _score_one(row, client, prompts, model, fetch, fetch_tree)
        except Exception:
            res = {"score": None, "one_liner": "[llm failed]"}
        return {**row, **res}
    with ThreadPoolExecutor(max_workers=int(cfg.get("llm", {}).get("concurrency", 8))) as ex:
        out = list(ex.map(work, candidates))
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
