"""build_runnable.py: makes runnable/ from original/ and helpers/. The ONLY change is to path prefixes: absolute paths of the machines
the work was done on become $BSS_REPO (this repository) or $BSS_DATA (the data folder). Python files go through the tokenizer, so nothing
but string literals and comments can change; every other token is carried over as it is. Standard library only.
    python3 -I tools/build_runnable.py   ->   runnable/, runnable/env.sh, audit/runnable_changes.csv
"""
import csv
import io
import re
import shutil
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "runnable"
NU, ZB = "/Volumes/User Homes/Nu Project", "/Volumes/User Homes"
HOMES = ("/Users/taseefrahman", "/Users/taseef", "~")
# longest prefix first: (old prefix, variable, new subfolder)
MAP = [(NU + "/BENDER-DYN/redo_evidence", "BSS_REPO", "redo_evidence"), (NU + "/BENDER-DYN/methods_paper/figures", "BSS_REPO", "methods_paper/figures")]
for h in HOMES:
    MAP += [(h + "/Downloads/Scripts/bender_reanalysis", "BSS_REPO", "bender_reanalysis"), (h + "/Downloads/Scripts", "BSS_DATA", "Scripts"),
            (h + "/Downloads", "BSS_DATA", "Downloads")]
MAP += [(NU, "BSS_DATA", ""), (ZB, "BSS_DATA", "")]
MAP.sort(key=lambda m: -len(m[0]))
LIT = re.compile(r"""^([rRfFbBuU]{0,2})(['"]{1,3})(.*)$""", re.S)
ENV_SH = """# source this file before running anything:  export BSS_DATA=/path/to/data; source runnable/env.sh
export BSS_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export BSS_DATA="${BSS_DATA:?set BSS_DATA to the folder that holds the unpacked data (see README)}"
export BENDER_REANALYSIS_DIR="$BSS_REPO/bender_reanalysis"
export PYTHONPATH="$BSS_REPO/redo_evidence:$BSS_REPO/bender_reanalysis:$BSS_REPO/methods_paper/figures:$PYTHONPATH"
"""


def match(text):
    for old, var, sub in MAP:
        if text == old or text.startswith(old + "/"):
            return var, (("/" + sub) if sub else "") + text[len(old):]
    return None


def as_text(s):
    """Comments and docstrings: the same paths, written as shell variables."""
    for old, var, sub in MAP:
        s = s.replace(old, "$" + var + (("/" + sub) if sub else ""))
    return s


def rewrite_py(src, log, name):
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    orig = list(toks)
    out = []
    prev_sig = None
    for i, t in enumerate(toks):
        tt, s = t.type, t.string
        new = s
        if tt == tokenize.COMMENT:
            new = as_text(s)
        elif tt == tokenize.STRING:
            docstring = prev_sig in (None, ":", "NEWLINE") and toks[i + 1].type in (tokenize.NEWLINE, tokenize.ENDMARKER, tokenize.COMMENT)
            m = LIT.match(s)
            body = m.group(3) if m else ""
            hit = match(body[:-len(m.group(2))]) if m else None
            if docstring:
                new = as_text(s)
            elif hit:
                var, rest = hit
                q, pre = m.group(2), m.group(1)
                new = f"(__import__('os').environ['{var}'] + {pre}{q}{rest}{q})"
            elif as_text(body) != body:
                # a path inside a longer string, such as a shell command printed for the user: left as text with the variable
                new = pre_text = as_text(s)
                log.append((name, t.start[0], "inside a longer string", s[:90]))
                new = pre_text
        elif tt == getattr(tokenize, "FSTRING_START", -1):
            nxt = toks[i + 1]
            hit = match(nxt.string) if nxt.type == getattr(tokenize, "FSTRING_MIDDLE", -2) else None
            if hit:
                var, rest = hit
                new = f"(__import__('os').environ['{var}'] + {s}"
                toks[i + 1] = nxt._replace(string=rest)
                j = next(k for k in range(i + 1, len(toks)) if toks[k].type == tokenize.FSTRING_END and toks[k].start >= nxt.start)
                toks[j] = toks[j]._replace(string=toks[j].string + ")")
        elif tt == getattr(tokenize, "FSTRING_MIDDLE", -2):
            new = as_text(s)
        if new != s:
            log.append((name, t.start[0], "path", s[:90]))
        out.append((t, new))
        if tt in (tokenize.NEWLINE,):
            prev_sig = "NEWLINE"
        elif tt not in (tokenize.NL, tokenize.COMMENT, tokenize.INDENT, tokenize.DEDENT):
            prev_sig = s
    # apply only the edits to the source text, so everything else is copied exactly as it was
    lines = src.splitlines(keepends=True)
    starts = [0]
    for ln in lines:
        starts.append(starts[-1] + len(ln))
    res = src
    for t0, (t, new) in sorted(zip(orig, out), key=lambda z: z[0].start, reverse=True):
        if new == t0.string:
            continue
        (r0, c0), (r1, c1) = t0.start, t0.end
        res = res[:starts[r0 - 1] + c0] + new + res[starts[r1 - 1] + c1:]
    return res


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    shutil.copytree(ROOT / "original", OUT, ignore=shutil.ignore_patterns("SHA256SUMS"))
    for h in (ROOT / "helpers").glob("*.py"):
        shutil.copy2(h, OUT / "redo_evidence" / h.name)
    log = []
    for p in sorted(OUT.rglob("*")):
        if p.suffix == ".py":
            src = p.read_text()
            p.write_text(rewrite_py(src, log, str(p.relative_to(OUT))))
        elif p.suffix == ".sh":
            src = p.read_text()
            new = as_text(src)
            if new != src:
                log.append((str(p.relative_to(OUT)), 0, "shell", ""))
            p.write_text(new)
    (OUT / "env.sh").write_text(ENV_SH)
    with open(ROOT / "audit/runnable_changes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "line", "kind", "was"])
        w.writerows(log)
    print(f"runnable/ written; {len(log)} changes in {len({l[0] for l in log})} files")


if __name__ == "__main__":
    main()
