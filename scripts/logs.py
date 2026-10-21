import re

_TS = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) \| (.*)$")
_EMOJI_KIND = {"\U0001F680": "push", "⭐": "star", "\U0001F500": "pr",
               "\U0001F41B": "issue", "\U0001F4AC": "issue", "\U0001F195": "branch"}
_KIND_HINTS = [("starred", "star"), ("pushed to", "push"), ("forked", "fork"),
               ("opened PR", "pr"), ("reopened PR", "pr"), ("closed PR", "pr"),
               ("created branch", "branch"), ("opened issue", "issue"),
               ("created issue", "issue")]
# 仓库 = 时间戳行(headline)末尾的 owner/repo。commit 详情在缩进续行上,经
# unwrap_records 并进同一 record,故先按 " - [hash]" 切掉详情再在 headline 上取
# 最后一个 owner/repo:既避开 commit message 里的 "a/b",又让分支名 'x/y' 在前、
# 真仓库在后时取到真仓库。
_DETAIL = re.compile(r"\s+-\s*\[[0-9a-fA-F]{4,}\]")
# owner 以字母数字开头(GitHub 登录名规则);repo 段允许下划线开头,如
# KULL-Centre/_2024_Cao_CALVADOSCOM
_REPO = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*/[A-Za-z0-9._-]+")

def unwrap_records(text: str) -> list[str]:
    out, buf = [], None
    for raw in text.splitlines():
        if _TS.match(raw):
            if buf:
                out.append(buf)
            buf = raw.rstrip()
        elif buf is not None:
            buf += " " + raw.strip()
    if buf:
        out.append(buf)
    return out

def _kind(emoji_and_text: str) -> str:
    for emoji, kind in _EMOJI_KIND.items():
        if emoji in emoji_and_text:
            return kind
    for hint, kind in _KIND_HINTS:
        if hint in emoji_and_text:
            return kind
    return "other"

def parse_event(record: str) -> dict | None:
    m = _TS.match(record)
    if not m:
        return None
    ts, body = m.group(1), m.group(2)
    kind = _kind(body)
    # 去掉 emoji 与首个人名/actor token
    stripped = re.sub(r"^[^\w]+\s*", "", body)          # 去 emoji
    actor = stripped.split()[0] if stripped.split() else ""
    # 仓库:切掉 commit 详情后,取 headline 上最后一个 owner/repo
    headline = _DETAIL.split(body, maxsplit=1)[0]
    matches = _REPO.findall(headline)
    repo = matches[-1] if matches else None
    return {"ts": ts, "kind": kind, "actor": actor, "repo": repo, "text": body}

def extract_events(text: str) -> list[dict]:
    evs = []
    for rec in unwrap_records(text):
        ev = parse_event(rec)
        if ev and ev["repo"]:
            evs.append(ev)
    return evs

def repos_from_events(events: list[dict]) -> dict:
    out: dict = {}
    for e in events:
        r = out.setdefault(e["repo"], {"who": [], "kinds": []})
        if e["actor"] and e["actor"] not in r["who"]:
            r["who"].append(e["actor"])
        if e["kind"] not in r["kinds"]:
            r["kinds"].append(e["kind"])
    return out
