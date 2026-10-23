import argparse, os, sys
from datetime import date as _date, timedelta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common, search, merge, score  # noqa: E402


def run(cfg: dict, run_date: str) -> list[dict]:
    """daily 模式的增量入口:只处理「昨天」,产出 scored/{昨天}.json。

    文件按被覆盖日命名(而不是运行日),weekly 的 aggregate 才能直接按
    [run-window_days, run-1] 选片,且不会和 weekly 自己产出的 scored/{运行日}.json 撞名。
    """
    covered = (_date.fromisoformat(run_date) - timedelta(days=1)).isoformat()
    search.run(cfg, covered, f">={covered}", covered)
    # window_days=0 + include_run_day:只取运行当天那份 feed —— 它装的是「昨天」的事件,
    # 正好和上面这次搜索覆盖的 pushed_at 对齐(不把前一天 feed 重复算一遍)。
    merge.run(cfg, run_date, window_days=0, include_run_day=True,
              search_date=covered, out_date=covered)
    return score.run(cfg, covered)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    if cfg.get("ingest", {}).get("mode", "weekly") != "daily":
        print("ingest mode=weekly, skip")
        return
    print(f"ingested={len(run(cfg, common.today()))}")


if __name__ == "__main__":
    main()
