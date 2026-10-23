import argparse, glob, os, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, gh  # noqa: E402

def _recent_daily(date: str, window_days: int, include_run_day: bool = False) -> list[str]:
    files = sorted(glob.glob(os.path.join(common.ROOT, "data", "daily", "*.repos.json")))
    run_day = _date.fromisoformat(date)
    # 窗口 = [run-window_days, run-1],与 search/cards/paper 同源(不含运行当天),
    # 保证「收集范围」和站点记录的「抓取窗口」完全一致(否则会多出 run 当天这天)。
    # daily 模式通过 include_run_day 把当天刚生成的 feed 纳入:该文件装的是「昨天」的事件,
    # 正好和同一次 daily 搜索覆盖的 pushed_at 对齐。
    cutoff = run_day - timedelta(days=window_days)
    upper = run_day if include_run_day else run_day - timedelta(days=1)
    keep = []
    for f in files:
        try:
            d = _date.fromisoformat(os.path.basename(f)[:10])
        except ValueError:
            continue
        if cutoff <= d <= upper:
            keep.append(f)
    return keep

def run(cfg: dict, date: str, *, window_days: int | None = None,
        include_run_day: bool = False, search_date: str | None = None,
        out_date: str | None = None) -> list[dict]:
    # 池 = 窗口 [run-window_days, run-1] 内的每日 feed ∪ 本周搜索;按 full_name 去重。
    # daily 模式传入 window_days=0/include_run_day=True(只取运行当天那份 feed),并让
    # 搜索文件与输出都按「被覆盖日」(search_date/out_date)命名,weekly 聚合才能直接按
    # [run-7, run-1] 选片。
    wd = int(cfg.get("window_days", 7) if window_days is None else window_days)
    pool: dict[str, dict] = {}
    for f in _recent_daily(date, wd, include_run_day):
        for r in common.load_json(f):
            pool.setdefault(r["full_name"], dict(r))
    search_path = common.data_dir("search", f"{search_date or date}.json")
    if os.path.exists(search_path):
        for r in common.load_json(search_path):
            cur = pool.setdefault(r["full_name"], dict(r))
            # 盲搜带 pushed_at,优先补全
            for k in ("pushed_at", "stars", "language", "description", "url"):
                cur[k] = r.get(k) or cur.get(k, "")
    candidates = list(pool.values())
    for r in candidates:
        if not r.get("pushed_at"):
            meta = gh.safe_repo_meta(r["full_name"])
            r.update({k: meta[k] for k in ("pushed_at", "stars", "language", "description", "url")})
    candidates.sort(key=lambda r: r.get("stars") or 0, reverse=True)
    common.dump_json(common.data_dir("candidates", f"{out_date or date}.json"), candidates)
    return candidates

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    print(f"candidates={len(run(cfg, common.today()))}")

if __name__ == "__main__":
    main()
