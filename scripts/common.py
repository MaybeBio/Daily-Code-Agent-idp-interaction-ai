import json, os
from datetime import datetime
from zoneinfo import ZoneInfo
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def _fallback_text(row, tree):
    desc = (row.get("description") or "").strip()
    if not desc and not tree:
        return ""
    parts = [f"仓库:{row['full_name']}"]
    if row.get("language"):
        parts.append(f"语言:{row['language']}")
    if desc:
        parts.append(f"描述:{desc}")
    if tree:
        parts.append(f"目录树:\n{tree}")
    return "\n".join(parts)

def repo_text(row, fetch, fetch_tree) -> str:
    """喂给 LLM 的仓库文本梯子:README 优先,缺失退到 描述+目录树,都无则空串。"""
    try:
        readme = fetch(row["full_name"])
    except Exception:
        readme = ""
    text = (readme or "").strip()
    if text:
        return text
    try:
        tree = fetch_tree(row["full_name"])
    except Exception:
        tree = ""
    return _fallback_text(row, tree)

def data_dir(*parts: str) -> str:
    d = os.path.join(ROOT, "data", *parts[:-1])
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, parts[-1])

def load_json(path: str):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def dump_json(path: str, obj) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)

def today() -> str:
    tz = os.environ.get("TZ", "Asia/Shanghai")
    return datetime.now(ZoneInfo(tz)).strftime("%Y-%m-%d")
