import ast, json, os, random, subprocess, time

_CODE_EXT = frozenset({
    ".py", ".pyx", ".pxd", ".pyi", ".pyw",
    ".js", ".jsx", ".mjs", ".cjs", ".ts", ".tsx",
    ".html", ".htm", ".css", ".scss", ".sass", ".less",
    ".java", ".kt", ".kts", ".scala", ".groovy", ".gradle",
    ".c", ".cpp", ".cc", ".cxx", ".h", ".hpp", ".hh", ".hxx",
    ".rs", ".go", ".zig",
    ".rb", ".rake", ".php", ".phtml", ".pl", ".pm", ".lua", ".tcl",
    ".sh", ".bash", ".zsh", ".fish", ".ps1", ".psm1", ".psd1",
    ".swift", ".m", ".mm",
    ".cs", ".fs", ".fsx", ".vb",
    ".r", ".jl", ".ipynb",
    ".hs", ".lhs", ".elm", ".ex", ".exs", ".erl", ".hrl",
    ".clj", ".cljs", ".cljc", ".edn",
    ".dart", ".vue", ".svelte",
    ".xsl", ".xslt", ".svg",
    ".proto", ".avsc", ".thrift",
    ".tf", ".tfvars", ".hcl",
    ".cmake", ".mk", ".meson",
    ".sql", ".psql",
    ".graphql", ".gql",
    ".md", ".mdx", ".markdown",
    ".tex", ".sty", ".cls", ".bib",
})

# 永久性失败:重试也救不回(404/空仓库/认证错),与瞬时网络故障(TLS 超时/连接重置/5xx)区分。
_PERMANENT_MARKERS = (
    "404", "not found", "no such", "is empty", "empty repository",
    "bad credentials", "authentication failed", "requires authentication",
    "not authenticated", "unknown flag", "no such option",
)

def _error_text(err):
    if isinstance(err, subprocess.CalledProcessError):
        return f"{err.stderr or ''} {err.stdout or ''}".lower()
    return str(err).lower()

def _is_permanent(err):
    return any(m in _error_text(err) for m in _PERMANENT_MARKERS)

def run(cmd: list[str], **kw) -> str:
    kw.setdefault("capture_output", True)
    kw.setdefault("text", True)
    kw.setdefault("check", True)
    attempts = kw.pop("attempts", 5)
    for i in range(attempts):
        try:
            return subprocess.run(cmd, **kw).stdout
        except Exception as e:
            if _is_permanent(e) or i == attempts - 1:
                raise
            time.sleep(2 ** (i + 2) + random.random())

def search_repos(config_path: str, updated: str, fields: list[str],
                 query: str | None = None) -> list[dict]:
    cmd = ["ghresearcher", "search"]
    if query:
        cmd += ["repos", query]
    cmd += ["--config", config_path, "--updated", updated, "--json", ",".join(fields)]
    out = run(cmd)
    start, end = out.find("["), out.rfind("]")
    if start == -1 or end <= start:
        raise ValueError(f"no result array in ghresearcher stdout: {out[:500]!r}")
    return ast.literal_eval(out[start:end + 1])

def repo_meta(full_name: str) -> dict:
    out = run(["gh", "api", f"repos/{full_name}"])
    d = json.loads(out)
    return {"full_name": d["full_name"], "url": d.get("html_url", ""),
            "language": d.get("language") or "", "stars": d.get("stargazers_count", 0),
            "description": d.get("description") or "", "pushed_at": d.get("pushed_at", "")}

def fetch_readme(full_name: str) -> str:
    return run(["gh", "api", f"repos/{full_name}/readme",
                "-H", "Accept: application/vnd.github.raw"])

def fetch_tree(full_name: str, limit: int = 200) -> str:
    # HEAD + recursive=1 直接取默认分支整棵树,免去先查 default_branch/tree_sha 的两次调用
    out = run(["gh", "api", f"repos/{full_name}/git/trees/HEAD?recursive=1"])
    tree = json.loads(out).get("tree", [])
    paths = sorted(e["path"] for e in tree
                   if e.get("type") == "blob"
                   and os.path.splitext(e.get("path", ""))[1].lower() in _CODE_EXT)
    if len(paths) > limit:
        paths = paths[:limit] + ["..."]
    return "\n".join(paths)

# 常见依赖/构建清单文件名,喂给 card 的「依赖」「复现成本」小节。
_DEP_FILES = frozenset({
    "requirements.txt", "pyproject.toml", "setup.py", "setup.cfg",
    "environment.yml", "environment.yaml", "Pipfile", "tox.ini",
    "package.json", "go.mod", "Cargo.toml", "pom.xml", "build.gradle",
    "build.gradle.kts", "Dockerfile", "Makefile", "CMakeLists.txt",
})

def fetch_dep_files(full_name: str, limit: int = 3000) -> str:
    """抓常见依赖/构建清单文件内容。独立查整棵树:这些文件扩展名(.toml/.txt/.cfg)多不在
    _CODE_EXT 白名单里,复用 fetch_tree 的结果会把它们过滤掉。"""
    out = run(["gh", "api", f"repos/{full_name}/git/trees/HEAD?recursive=1"])
    tree = json.loads(out).get("tree", [])
    found = sorted({e["path"] for e in tree
                    if e.get("type") == "blob"
                    and e["path"].rsplit("/", 1)[-1] in _DEP_FILES})
    blocks = []
    for path in found:
        try:
            content = run(["gh", "api", f"repos/{full_name}/contents/{path}",
                           "-H", "Accept: application/vnd.github.raw"])
        except Exception:
            continue
        content = (content or "").strip()
        if not content:
            continue
        if len(content) > limit:
            content = content[:limit] + "\n...(截断)"
        blocks.append(f"### {path}\n{content}")
    return "\n\n".join(blocks)

def safe_repo_meta(full_name: str) -> dict:
    try:
        return repo_meta(full_name)
    except Exception:
        return {"full_name": full_name, "url": "", "language": "",
                "stars": 0, "description": "", "pushed_at": ""}
