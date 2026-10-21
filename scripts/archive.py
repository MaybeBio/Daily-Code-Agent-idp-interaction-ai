import argparse, datetime as dt, json, os, random, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common  # noqa: E402

WIKI_SOURCES = ["deepwiki", "codewiki", "zread"]

# 仓库在该 wiki 源上尚未被索引(404/not indexed),与瞬时网络故障区分。
_NOT_INDEXED_MARKERS = (
    "404", "not found", "not indexed", "no pages", "no wiki",
    "no documentation", "does not exist", "is empty", "empty repository",
)
# 永久性失败(索引缺失 + 认证/配置错误):重试也救不回,直接放弃。
_PERMANENT_MARKERS = _NOT_INDEXED_MARKERS + (
    "invalid", "authentication failed", "bad credentials", "requires authentication",
)

def threshold(scored, min_score):
    keep = [r for r in scored if r.get("score") is not None and r["score"] >= min_score]
    keep.sort(key=lambda r: r["score"], reverse=True)
    return keep

def repo_dir(full_name):
    return os.path.join(common.ROOT, "code", full_name.replace("/", "__"))

def _wiki_dir(source):
    return "google_code_wiki" if source == "codewiki" else source

def _wiki_dest(base, source, date):
    # wiki 只增不删:按抓取日期落快照,过去周次的快照永不被覆盖。
    return os.path.join(base, "wikis", date, _wiki_dir(source))

def _has_snapshot(base, source):
    wikis_root = os.path.join(base, "wikis")
    if not os.path.isdir(wikis_root):
        return False
    target = _wiki_dir(source)
    return any(os.path.isdir(os.path.join(wikis_root, d, target))
               for d in os.listdir(wikis_root))

def _error_text(err):
    if isinstance(err, subprocess.CalledProcessError):
        return f"{err.stderr or ''} {err.stdout or ''}".lower()
    return str(err).lower()

def _is_not_indexed(err):
    return any(m in _error_text(err) for m in _NOT_INDEXED_MARKERS)

def _is_permanent(err):
    return any(m in _error_text(err) for m in _PERMANENT_MARKERS)

def _run_retry(cmd, tag, attempts=5):
    last = None
    for i in range(attempts):
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            return True, None
        except Exception as e:
            last = e
            if _is_permanent(e):      # 404/无此仓库/认证错 → 别再白试
                break
            if i < attempts - 1:
                time.sleep(2 ** (i + 2) + random.random())
    print(f"[{tag} failed] {' '.join(cmd)}: {last}", file=sys.stderr)
    return False, last

def clone_repo(full_name, dest):
    ok, _ = _run_retry(["degit", full_name, dest, "--force"], "clone")
    if ok:
        return True
    shutil.rmtree(dest, ignore_errors=True)  # degit 可能留下半截目录
    ok, _ = _run_retry(["git", "clone", "--depth", "1",
                        f"https://github.com/{full_name}.git", dest], "clone:git")
    return ok

def _submit_zread(full_name):
    if not os.environ.get("ZREAD_TOKEN"):
        return False, None
    return _run_retry(["repowiki-cli", "zread", "submit", full_name], "wiki:zread:submit")

def export_wiki(source, full_name, dest):
    ok, err = _run_retry(["repowiki-cli", source, "cp", full_name, dest], f"wiki:{source}")
    if ok:
        return "ok"
    if source == "zread" and _is_not_indexed(err):
        sub_ok, _ = _submit_zread(full_name)
        if sub_ok:
            return "submitted"
    return "not_indexed" if _is_not_indexed(err) else "failed"

def _run_json(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed ({r.returncode}): {' '.join(cmd)}")
    return json.loads(r.stdout)

def _run_text(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed ({r.returncode}): {' '.join(cmd)}")
    return r.stdout.strip()

def _parse_iso(v):
    if not isinstance(v, str):
        return None
    try:
        return dt.date.fromisoformat(v.strip()[:10])
    except ValueError:
        return None

def _parse_epoch(v):
    try:
        return dt.datetime.fromtimestamp(int(v), tz=dt.timezone.utc).date()
    except (TypeError, ValueError, OSError):
        return None

def _stat_date(source, full_name):
    """返回该 wiki 源的「内容截止日期」,拿不到返回 None。

    deepwiki 用 last_modified、zread 用 updated_at(epoch 秒)、codewiki 只有
    commit sha,需经 gh api 反查该 commit 的提交日期。
    """
    try:
        j = _run_json(["repowiki-cli", source, "stat", full_name, "--json"])
    except Exception:
        return None
    if not isinstance(j, dict):
        return None
    data = j.get("data") if isinstance(j.get("data"), dict) else {}
    if source == "deepwiki":
        return _parse_iso(data.get("last_modified"))
    if source == "zread":
        return _parse_epoch(data.get("updated_at"))
    if source == "codewiki":
        sha = (j.get("commit") or "").strip()
        if not sha:
            return None
        owner, _, name = full_name.partition("/")
        try:
            out = _run_text(["gh", "api", f"repos/{owner}/{name}/commits/{sha}",
                             "--jq", ".commit.committer.date"])
        except Exception:
            return None
        return _parse_iso(out)
    return None

def _in_week(d, run_date):
    """d 是否落在 run_date 所在的自然周(周一~周日)内。"""
    if d is None:
        return False
    try:
        r = dt.date.fromisoformat((run_date or "")[:10])
    except ValueError:
        return False
    monday = r - dt.timedelta(days=r.weekday())
    return monday <= d <= monday + dt.timedelta(days=6)

def archive_one(row, clone, wiki, date, stat_date=None, in_week=None):
    full_name = row["full_name"]
    base = repo_dir(full_name)
    os.makedirs(base, exist_ok=True)
    code_dir = os.path.join(base, "code")
    shutil.rmtree(code_dir, ignore_errors=True)
    clone_ok = clone(full_name, code_dir)
    wikis = {}
    for s in WIKI_SOURCES:
        dest = _wiki_dest(base, s, date)
        if stat_date is not None and in_week is not None:
            d = stat_date(s, full_name)
            fresh = d is not None and in_week(d)
            if not fresh and _has_snapshot(base, s):   # 没更新(或 stat 失败)且已留过快照 → 复用,不再新打
                wikis[s] = "reused"
                continue
        shutil.rmtree(dest, ignore_errors=True)
        wikis[s] = wiki(s, full_name, dest)
    return {"full_name": full_name, "score": row.get("score"),
            "one_liner": row.get("one_liner", ""), "clone_ok": clone_ok, "wikis": wikis}

def _safe_archive_one(row, clone, wiki, date, stat_date=None, in_week=None):
    try:
        return archive_one(row, clone, wiki, date, stat_date, in_week)
    except Exception as e:
        print(f"[archive failed] {row['full_name']}: {e}", file=sys.stderr)
        return {"full_name": row["full_name"], "score": row.get("score"),
                "one_liner": row.get("one_liner", ""), "clone_ok": False,
                "wikis": {s: "failed" for s in WIKI_SOURCES}}

def run(cfg, date, clone=None, wiki=None, stat_date=None, in_week=None):
    clone = clone or clone_repo
    wiki = wiki or export_wiki
    min_score = int(cfg.get("llm", {}).get("min_score", 5))
    scored_path = common.data_dir("scored", f"{date}.json")
    scored = common.load_json(scored_path) if os.path.exists(scored_path) else []
    keep = threshold(scored, min_score)
    # 低并发:网络工具(degit/repowiki-cli)串行太慢,但太高会触发 repowiki 限流。
    workers = max(1, int(cfg.get("archive", {}).get("concurrency", 4)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        archived = list(ex.map(
            lambda r: _safe_archive_one(r, clone, wiki, date, stat_date, in_week), keep))
    common.dump_json(common.data_dir("archive", f"{date}.json"), archived)
    return archived

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(common.ROOT, "config.yaml"))
    ap.add_argument("--date", default=None)
    a = ap.parse_args()
    cfg = common.load_config(a.config)
    date = a.date or common.today()
    in_week = lambda d: _in_week(d, date)
    out = run(cfg, date, stat_date=_stat_date, in_week=in_week)
    ok_clone = sum(1 for r in out if r["clone_ok"])
    parts = [f"{len(out)} archived; clone ok {ok_clone}/{len(out)}"]
    for s in WIKI_SOURCES:
        counts = {st: sum(1 for r in out if r["wikis"].get(s) == st)
                  for st in ("ok", "submitted", "not_indexed", "failed", "reused")}
        nonzero = " ".join(f"{st}:{n}" for st, n in counts.items() if n)
        parts.append(f"{s} {nonzero or '-'}")
    print("; ".join(parts))

if __name__ == "__main__":
    main()
