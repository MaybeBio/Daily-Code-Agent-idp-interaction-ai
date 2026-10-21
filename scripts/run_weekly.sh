#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PY:-python}"
# 窗口天数与 merge.py 同源：统一从 config.yaml 的 window_days 读取
WINDOW_DAYS="$("$PY" -c 'import yaml; print(yaml.safe_load(open("config.yaml")).get("window_days", 7))')"
# 窗口 = [run-window_days, run-1](不含运行当天),与 merge/cards/issue 同源。
SINCE=">=$(date -d "$WINDOW_DAYS days ago" +%Y-%m-%d)"
UNTIL="$(date -d "1 day ago" +%Y-%m-%d)"

echo "[weekly] search";       "$PY" scripts/search.py --config config.yaml --since "$SINCE" --until "$UNTIL"
echo "[weekly] merge";        "$PY" scripts/merge.py  --config config.yaml
echo "[weekly] score";        "$PY" scripts/score.py  --config config.yaml
echo "[weekly] archive";      "$PY" scripts/archive.py --config config.yaml
echo "[weekly] cards";        "$PY" scripts/cards.py   --config config.yaml
echo "[weekly] issue";        "$PY" scripts/issue.py  --config config.yaml \
    --issue-body /tmp/issue.md --issue-title /tmp/issue.title
echo "[weekly] done"
