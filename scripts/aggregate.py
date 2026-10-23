import argparse, os, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, search, merge, score  # noqa: E402


def _window(date: str, window_days: int) -> list[str]:
    """[run-window_days, ..., run-1],与 merge/cards/issue 的窗口同源。"""
    run_day = _date.fromisoformat(date)
    return [(run_day - timedelta(days=i)).isoformat()
            for i in range(window_days, 0, -1)]


def _scored_path(day: str) -> str:
    return os.path.join(common.ROOT, "data", "scored", f"{day}.json")


def missing_days(days: list[str]) -> list[str]:
    return [d for d in days if not os.path.exists(_scored_path(d))]


def _dedup(rows: list[dict]) -> list[dict]:
    """同一仓库一周内可能多次被 push → 保留 pushed_at 最新的一条;同一天则取高分。"""
    best: dict[str, dict] = {}
    for r in rows:
        name = r.get("full_name")
        if not name:
            continue
        cur = best.get(name)
        key = (r.get("pushed_at") or "", r.get("score") or -1)
        if cur is None or key > (cur.get("pushed_at") or "", cur.get("score") or -1):
            best[name] = dict(r)
    out = list(best.values())
    out.sort(key=lambda r: (r.get("score") or -1, r.get("stars") or 0), reverse=True)
    return out


def _backfill(cfg: dict, day: str) -> None:
    """补跑某个被覆盖日:搜索该日 → 合并(当天 feed + 该日搜索)→ 打分。"""
    run_day = (_date.fromisoformat(day) + timedelta(days=1)).isoformat()
    search.run(cfg, day, f">={day}", day)
    merge.run(cfg, run_day, window_days=0, include_run_day=True,
              search_date=day, out_date=day)
    score.run(cfg, day)


def run(cfg: dict, date: str, window_days: int | None = None,
        backfill_limit: int | None = None) -> list[dict]:
    wd = int(cfg.get("window_days", 7) if window_days is None else window_days)
    days = _window(date, wd)
    missing = missing_days(days)
    if missing:
        limit = int(cfg.get("ingest", {}).get("max_backfill_days", 2)) \
            if backfill_limit is None else backfill_limit
        for day in missing[:limit]:
            print(f"[aggregate] backfill {day}", file=sys.stderr)
            try:
                _backfill(cfg, day)
            except Exception as e:      # 单日补跑失败不拖垮整周,缺片留给下面的告警
                print(f"[aggregate] backfill {day} failed: {e}", file=sys.stderr)
        still = missing_days(days)
        if still:
            print(f"[aggregate] still missing: {', '.join(still)}", file=sys.stderr)
    rows = []
    for day in days:
        p = _scored_path(day)
        if os.path.exists(p):
            rows.extend(common.load_json(p))
    out = _dedup(rows)
    common.dump_json(common.data_dir("scored", f"{date}.json"), out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--date", default=None)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    if cfg.get("ingest", {}).get("mode", "weekly") != "daily":
        print("aggregate only runs in ingest.mode=daily, skip")
        return
    out = run(cfg, a.date or common.today())
    print(f"aggregated={len(out)}")


if __name__ == "__main__":
    main()
