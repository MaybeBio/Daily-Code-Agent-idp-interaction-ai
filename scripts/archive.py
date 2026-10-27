import argparse, base64, datetime as dt, json, os, posixpath, random, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import quote, unquote
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts import common  # noqa: E402

WIKI_SOURCES = ["deepwiki", "codewiki", "zread"]

# ── 代码摘录策略(repo_copy=code) ─────────────────────────────────────
# 归档的目标是攒一份能被 RAG / skill 直接取用的代码语料,不是备份仓库:
# 源码、构建文件、文档保留;数据、二进制、图片、权重一律不留,并在
# MANIFEST.md 里写明省略了什么、去哪取回。
# 后缀名单取自 ghresearcher.scraper._PROGRAMMING_EXTENSIONS(该表只用于目录树
# 渲染,没有大小逻辑),去掉 .svg(多为成组图标/示意图),补上科研常用后缀。
SOURCE_EXTS = frozenset("""
.py .pyx .pxd .pyi .pyw .py3
.js .jsx .mjs .cjs .ts .tsx .vue .svelte
.html .htm .css .scss .sass .less
.java .kt .kts .scala .groovy .gradle
.c .h .cc .cpp .cxx .hpp .hxx .cu .cuh .ino
.rs .go .zig
.rb .php .pl .pm .lua .tcl
.sh .bash .zsh .fish .ps1 .csh .tcsh .awk
.pbs .sbatch .slurm
.swift .m .mm .ml .mli
.cs .fs .fsx .vb
.r .jl .hs .elm .ex .exs .erl .clj .cljs .cljc .dart
.f .f77 .f90 .f95 .f03 .for
.nf .smk .stan
.plt .proto .thrift .graphql .gql
.sql .cmake .mk .meson
.xsl .xslt
.ipynb
""".split())

# 文档格式:人读的说明/论文/讲义。通用文本(.txt/.log/.csv)不在此列。
DOC_EXTS = frozenset({".md", ".mdx", ".markdown", ".rst", ".rmd", ".qmd",
                      ".tex", ".sty", ".cls", ".bib", ".adoc", ".org"})

# 配置/声明式格式:环境与项目配置。大 JSON/XML 数据不在此列(仅按名保留清单)。
BUILD_EXTS = frozenset({".yaml", ".yml", ".toml", ".ini", ".cfg", ".conf", ".lock"})

# 无后缀或特殊名的构建/依赖/环境清单(按文件名精确匹配,大小写不敏感)。
BUILD_NAMES = frozenset({
    "makefile", "gnumakefile", "dockerfile", "containerfile",
    "cmakelists.txt", "snakefile", "justfile", "rakefile",
    "gemfile", "procfile", "vagrantfile", "configure", "build",
    "workspace", "module.bazel", "go.mod", "go.sum",
    "pyproject.toml", "setup.py", "setup.cfg",
    "package.json", "package-lock.json", "yarn.lock",
    "cargo.toml", "cargo.lock",
    "pipfile", "pipfile.lock", "poetry.lock",
    "pom.xml", "build.xml", "settings.xml",
    "tsconfig.json", "jsconfig.json", "manifest.json",
    "composer.json", "composer.lock",
    "build.gradle", "build.gradle.kts", "settings.gradle",
    "meson.build", "makefile.am",
})

# 依赖/环境清单的常见前缀(requirements*.txt / environment*.yml 等)。
BUILD_PREFIXES = ("requirements", "constraints", "environment")

# 按文件名识别的文档(后缀可能是 .txt 或无后缀,如 LICENSE / CHANGELOG.txt)。
DOC_NAMES = frozenset({
    "license", "licence", "copying", "copying.lesser",
    "notice", "changelog", "changes", "history", "news",
    "contributing", "authors", "thanks",
    "acknowledgements", "acknowledgments",
    "security", "code_of_conduct", "codeowners",
    "citation", "citations", "about", "support",
    "maintainers", "governance",
})

SOURCE_CAP = 2 * 1024 * 1024   # 普通源码保险丝(实测最大 0.46 MB,2 MB 只挡 vendored 大包)
DOCS_CAP = 1 * 1024 * 1024     # 文档上限(实测最大 0.27 MB)
BUILD_CAP = 1 * 1024 * 1024    # 构建/配置上限(实测最大 79 KB)
# README 与 .ipynb 不设上限:前者是核心材料,后者是「代码+结果」的完整产物,都该整留。

# README 可带的后缀;语言/区域变体(README_zh.md 等)由 _is_readme 前缀判断。
_PLAIN_EXTS = frozenset({"", ".md", ".markdown", ".rst", ".txt",
                         ".adoc", ".org", ".html", ".htm"})

_CLASS_LABEL = {"source": "源码", "build": "构建/配置", "docs": "文档",
                "readme": "README", "tree": "目录树"}

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

def _run_retry(cmd, tag, attempts=5, env=None):
    last = None
    for i in range(attempts):
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True, env=env)
            return True, None
        except Exception as e:
            last = e
            if _is_permanent(e):      # 404/无此仓库/认证错 → 别再白试
                break
            if i < attempts - 1:
                time.sleep(2 ** (i + 2) + random.random())
    print(f"[{tag} failed] {' '.join(cmd)}: {last}", file=sys.stderr)
    return False, last

def _strip_lfs_attrs(dest):
    """删掉快照里声明 filter=lfs 的 .gitattributes。

    degit 下载的是 tarball,LFS 大文件只剩 ~130B 指针文本;但仓库自带的
    .gitattributes 仍把这些路径标成 filter=lfs。提交时 git 会把指针当 LFS
    对象去上传,而本地根本没有对象,push 直接被拒(整个 weekly 提交丢失)。
    归档只需要源码文本,删掉这些声明即可,指针文件随后按普通文本入库。
    """
    for root, dirs, files in os.walk(dest):
        dirs[:] = [d for d in dirs if d != ".git"]
        if ".gitattributes" not in files:
            continue
        path = os.path.join(root, ".gitattributes")
        try:
            with open(path, encoding="utf-8", errors="ignore") as f:
                text = f.read()
        except OSError:
            continue
        if "filter=lfs" in text.replace(" ", "").replace("\t", ""):
            os.remove(path)

def _git_env():
    # 关掉 LFS smudge:归档只要文本,大权重没必要下载(还会拖慢甚至撞超时)。
    return {**os.environ, "GIT_LFS_SKIP_SMUDGE": "1"}

def clone_repo(full_name, dest):
    """完整浅克隆(degit tarball,失败退回 git)。落盘后删掉 .git。"""
    ok, _ = _run_retry(["degit", full_name, dest, "--force"], "clone")
    if not ok:
        shutil.rmtree(dest, ignore_errors=True)  # degit 可能留下半截目录
        ok, _ = _run_retry(["git", "clone", "--depth", "1",
                            f"https://github.com/{full_name}.git", dest],
                           "clone:git", env=_git_env())
    if ok:
        _strip_lfs_attrs(dest)
        # 语料不需要 .git(历史、对象都在 GitHub 上),删掉省体积也避免嵌套仓库。
        shutil.rmtree(os.path.join(dest, ".git"), ignore_errors=True)
    return ok

# ── 代码摘录:选文件 / 取文件 / 写清单 ──────────────────────────────

def _ext(path):
    return os.path.splitext(posixpath.basename(path))[1].lower()

def _fmt_size(n):
    n = float(n or 0)
    if n < 1024:
        return f"{n:.0f} B"
    for unit in ("KB", "MB", "GB"):
        n /= 1024
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
    return f"{n:.1f} GB"

def _is_readme(path):
    """README 及其变体:README.md / readme.txt / README_zh.md / README-zh-CN.md 等。"""
    name = posixpath.basename(path).lower()
    return name.startswith("readme") and _ext(path) in _PLAIN_EXTS

def _classify(path):
    """→ 类别名;None 表示不在白名单。"""
    name = posixpath.basename(path).lower()
    ext = _ext(path)
    if _is_readme(path):
        return "readme"
    if name in BUILD_NAMES:
        return "build"
    if name.startswith(BUILD_PREFIXES) and \
            ext in (".txt", ".yml", ".yaml", ".toml", ".ini", ".cfg", ".conf", ".lock", ""):
        return "build"
    if ext in BUILD_EXTS:
        return "build"
    stem = name.split(".", 1)[0]
    if (name in DOC_NAMES or stem in DOC_NAMES) and ext in _PLAIN_EXTS:
        return "docs"
    if ext in DOC_EXTS:
        return "docs"
    if ext in SOURCE_EXTS:
        return "source"
    return None

def _cap(cls, path):
    """单文件上限(字节);None 表示不限。"""
    if cls == "readme":
        return None
    if cls == "source" and _ext(path) == ".ipynb":
        return None
    if cls == "docs":
        return DOCS_CAP
    if cls == "build":
        return BUILD_CAP
    return SOURCE_CAP

def select_entries(entries, mode="code"):
    """按策略筛文件。entries: [{"path","size","type"?}]。

    返回 (keep, omitted),元素都带 path/size;keep 另带 class,omitted 带 reason。
    """
    keep, omitted = [], []
    for e in entries:
        path, size = e["path"], int(e.get("size") or 0)
        if e.get("type") == "commit":
            omitted.append({"path": path, "size": size, "reason": "子模块未拉取"})
            continue
        if mode == "readme":
            if _is_readme(path):
                keep.append({"path": path, "size": size, "class": "readme"})
            else:
                omitted.append({"path": path, "size": size, "reason": "仅保留 README"})
            continue
        cls = _classify(path)
        if cls is None:
            omitted.append({"path": path, "size": size, "reason": "后缀不在白名单"})
            continue
        cap = _cap(cls, path)
        if cap is not None and size > cap:
            omitted.append({"path": path, "size": size, "reason": "超过大小上限"})
        else:
            keep.append({"path": path, "size": size, "class": cls})
    return keep, omitted

def _tree_api(full_name):
    """GitHub tree API:一次拿全量路径 + 每个 blob 的 size(不下载内容)。"""
    j = _run_json(["gh", "api", f"repos/{full_name}/git/trees/HEAD?recursive=1"])
    entries = [{"path": e["path"], "size": int(e.get("size") or 0), "type": e["type"]}
               for e in j.get("tree", []) if e.get("type") in ("blob", "commit")]
    return entries, bool(j.get("truncated"))

def _walk_entries(dest):
    entries = []
    for root, dirs, files in os.walk(dest):
        dirs[:] = [d for d in dirs if d != ".git"]
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, dest).replace(os.sep, "/")
            try:
                size = os.path.getsize(p)
            except OSError:
                size = 0
            entries.append({"path": rel, "size": size, "type": "blob"})
    return entries

def _sparse_pattern(path):
    # gitignore 语法:首个 "/" 锚定根目录;通配符要转义,否则文件名会被当模式。
    p = path.replace("\\", "\\\\")
    for ch in "*?[]":
        p = p.replace(ch, "\\" + ch)
    return "/" + p

def _checkout(dest, paths):
    """只把 paths 检出到工作区,其余文件连 blob 都不下载。"""
    if not paths:
        return
    info = os.path.join(dest, ".git", "info")
    os.makedirs(info, exist_ok=True)
    with open(os.path.join(info, "sparse-checkout"), "w", encoding="utf-8") as f:
        f.write("\n".join(_sparse_pattern(p) for p in paths) + "\n")
    subprocess.run(["git", "-C", dest, "config", "core.sparseCheckout", "true"],
                   capture_output=True, text=True, check=True, env=_git_env())
    subprocess.run(["git", "-C", dest, "checkout", "HEAD"],
                   capture_output=True, text=True, check=True, env=_git_env())

def _head_sha(dest):
    try:
        return _run_text(["git", "-C", dest, "rev-parse", "HEAD"]).strip()
    except Exception:
        return ""

# ── README / 文档里的相对链接改写为 GitHub 绝对地址 ─────────────────

_IMG_EXTS = {".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp", ".bmp",
             ".ico", ".tif", ".tiff", ".avif"}
_SCHEME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.\-]*:")
_REF_DEF_RE = re.compile(r"^(\s{0,3}\[[^\]]+\]:\s*)(\S+)(.*)$")
_HTML_IMG_RE = re.compile(r"(<img\b[^>]*?\bsrc\s*=\s*)([\"'])([^\"']+)\2", re.I)
_HTML_A_RE = re.compile(r"(<a\b[^>]*?\bhref\s*=\s*)([\"'])([^\"']+)\2", re.I)
_MD_IMG_RE = re.compile(r"(!\[[^\]]*\]\()([^)\s]+)((?:\s+(?:\"[^\"]*\"|'[^']*'))?\))")
_MD_LINK_RE = re.compile(r"((?<!!)\[[^\]]*\]\()([^)\s]+)((?:\s+(?:\"[^\"]*\"|'[^']*'))?\))")
_FENCE_RE = re.compile(r"^\s*(```|~~~)")

def _rewrite_url(raw, base_dir, full_name, kind):
    """相对路径 → github.com/{repo}/{raw|blob}/HEAD/{path};其余原样返回。"""
    u = raw.strip()
    # "#anchor" / "//cdn" / "https:" / "mailto:" 原样;
    # "@ref" 是 Documenter.jl 的交叉引用,不是路径,也不能改写。
    if not u or u.startswith(("#", "@", "//")) or _SCHEME_RE.match(u):
        return raw
    m = re.match(r"^([^?#]*)([?#].*)?$", u)
    path, suffix = m.group(1), m.group(2) or ""
    if not path:
        return raw
    if path.startswith("/"):                 # GitHub 里 "/x" 指仓库根
        path = path.lstrip("/")
    else:
        path = posixpath.normpath(posixpath.join(base_dir, path))
    if path.startswith("..") or path.startswith("/"):
        return raw                           # 跑出仓库的链接无从改写
    return (f"https://github.com/{full_name}/{kind}/HEAD/"
            f"{quote(unquote(path), safe='/')}{suffix}")

def _rewrite_line(line, base_dir, full_name):
    def ref(m):
        url = m.group(2)
        ext = _ext(url.split("#")[0].split("?")[0])
        kind = "raw" if ext in _IMG_EXTS else "blob"
        return m.group(1) + _rewrite_url(url, base_dir, full_name, kind) + m.group(3)

    def attr(kind):
        def repl(m):
            return (m.group(1) + m.group(2) +
                    _rewrite_url(m.group(3), base_dir, full_name, kind) + m.group(2))
        return repl

    def md(kind):
        def repl(m):
            return (m.group(1) +
                    _rewrite_url(m.group(2), base_dir, full_name, kind) + m.group(3))
        return repl

    line = _REF_DEF_RE.sub(ref, line)
    line = _HTML_IMG_RE.sub(attr("raw"), line)
    line = _HTML_A_RE.sub(attr("blob"), line)
    line = _MD_IMG_RE.sub(md("raw"), line)
    line = _MD_LINK_RE.sub(md("blob"), line)
    return line

def rewrite_md(text, md_path, full_name):
    """改写 Markdown 里的相对图片/链接;围栏代码块内不动。"""
    base_dir = posixpath.dirname(md_path)
    out, fenced = [], False
    for line in text.splitlines(keepends=True):
        if _FENCE_RE.match(line):
            fenced = not fenced
            out.append(line)
        elif fenced:
            out.append(line)
        else:
            out.append(_rewrite_line(line, base_dir, full_name))
    return "".join(out)

def _rewrite_file_links(dest, rel, full_name):
    if not rel.lower().endswith((".md", ".markdown")):
        return
    p = os.path.join(dest, *rel.split("/"))
    try:
        with open(p, encoding="utf-8") as f:
            text = f.read()
    except (OSError, UnicodeDecodeError):
        return
    new = rewrite_md(text, rel, full_name)
    if new != text:
        with open(p, "w", encoding="utf-8") as f:
            f.write(new)

def _rewrite_all_md(dest, keep, full_name):
    for k in keep:
        _rewrite_file_links(dest, k["path"], full_name)

# ── MANIFEST.md ─────────────────────────────────────────────────────

def _group(items, key):
    out = {}
    for it in items:
        k = key(it)
        n, b = out.get(k, (0, 0))
        out[k] = (n + 1, b + it["size"])
    return out

def _ext_key(item):
    name = posixpath.basename(item["path"])
    if name.lower() in BUILD_NAMES:
        return name
    return _ext(name) or "无后缀"

def _cell(path):
    return path.replace("|", "\\|")

def _write_manifest(base, full_name, mode, keep, omitted, commit, date=""):
    policy = {
        "code": ("`repo_copy=code` — **代码摘录,不是完整仓库**。保留源码、构建/配置"
                 "文件、README、文档;数据、二进制、图片、权重不保留。"),
        "readme": "`repo_copy=readme` — 只保留 README 与目录树(TREE.md),其余文件均未下载。",
        "full": "`repo_copy=full` — 完整浅克隆(已删除 `.git`)。",
    }.get(mode, "")
    kept_n = len(keep)
    kept_b = sum(k["size"] for k in keep)
    lines = ["# 代码快照清单", "",
             f"- 上游仓库: https://github.com/{full_name}",
             f"- 快照日期: {date or '未知'}"]
    if commit:
        lines.append(f"- 快照 commit: [`{commit[:12]}`]"
                     f"(https://github.com/{full_name}/commit/{commit})")
    lines += [f"- 归档策略: {policy}", "",
              f"## 保留 ({kept_n} 个文件, {_fmt_size(kept_b)})", ""]
    if keep:
        lines += ["| 类别 | 文件数 | 大小 |", "|---|---:|---:|"]
        for cls, (n, b) in sorted(_group(keep, lambda k: k["class"]).items(),
                                  key=lambda kv: -kv[1][1]):
            lines.append(f"| {_CLASS_LABEL.get(cls, cls)} | {n} | {_fmt_size(b)} |")
        by_ext = sorted(_group(keep, _ext_key).items(), key=lambda kv: -kv[1][1])
        lines += ["", "按后缀: " + "、".join(
            f"`{ext}` {n} ({_fmt_size(b)})" for ext, (n, b) in by_ext)]
        lines += ["", "> README / 文档中的相对图片与链接已改写为 GitHub 绝对地址。"]
    else:
        lines.append("(无)")
    lines += ["", f"## 省略 ({len(omitted)} 个文件, "
              f"{_fmt_size(sum(o['size'] for o in omitted))})", ""]
    if omitted:
        lines += ["| 原因 | 文件数 | 大小 |", "|---|---:|---:|"]
        for reason, (n, b) in sorted(_group(omitted, lambda o: o["reason"]).items(),
                                     key=lambda kv: -kv[1][1]):
            lines.append(f"| {reason} | {n} | {_fmt_size(b)} |")
        lines += ["", "体积最大的省略文件:", "",
                  "| 路径 | 大小 |", "|---|---:|"]
        for o in sorted(omitted, key=lambda o: -o["size"])[:20]:
            url = (f"https://github.com/{full_name}/raw/HEAD/"
                   f"{quote(o['path'], safe='/')}")
            lines.append(f"| [`{_cell(o['path'])}`]({url}) | {_fmt_size(o['size'])} |")
    else:
        lines.append("(无)")
    if mode == "code":
        lines += ["", "> 省略规则:普通源码 >2 MB、文档/构建配置 >1 MB(README 与 .ipynb "
                  "不设限),或后缀不在白名单(图片/权重/数据表/二进制/通用 "
                  ".txt/.log/.csv/.json/.xml 等)、以及子模块。"
                  "需要全文时按上表链接直接取回,或用 git 克隆上游仓库。"]
    lines += ["", "> 本清单由 `scripts/archive.py` 每次归档时自动生成并覆盖。", ""]
    with open(os.path.join(base, "MANIFEST.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

# ── 三种落盘方式 ────────────────────────────────────────────────────

def copy_corpus(full_name, dest, mode="code", date=""):
    """按 repo_copy 策略落成代码语料。

    先用 tree API 拿到全量路径和每个 blob 的大小(不下载内容),筛出保留清单;
    再用 --filter=blob:none 克隆 + sparse-checkout 只拉这些文件。Droplet 那种
    655 MB 的仓库实测只下载 44 MB。tree API 不可用/被截断时退回完整克隆后本地裁剪。
    """
    if mode == "readme":
        return _copy_readme(full_name, dest, date)
    shutil.rmtree(dest, ignore_errors=True)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        entries, truncated = _tree_api(full_name)
    except Exception as e:
        print(f"[tree failed] {full_name}: {e}", file=sys.stderr)
        entries, truncated = None, True
    if truncated or entries is None:
        return _copy_full_fallback(full_name, dest, mode, date)
    try:
        ok, _ = _run_retry(["git", "clone", "--depth", "1", "--filter=blob:none",
                            "--no-checkout",
                            f"https://github.com/{full_name}.git", dest],
                           "clone:partial", env=_git_env())
        if not ok:
            shutil.rmtree(dest, ignore_errors=True)
            return False
        keep, omitted = select_entries(entries, mode)
        _checkout(dest, [k["path"] for k in keep])
        commit = _head_sha(dest)
        _rewrite_all_md(dest, keep, full_name)
        shutil.rmtree(os.path.join(dest, ".git"), ignore_errors=True)
        _write_manifest(os.path.dirname(dest), full_name, mode, keep, omitted,
                        commit, date)
        return True
    except Exception as e:
        print(f"[copy failed] {full_name}: {e}", file=sys.stderr)
        shutil.rmtree(dest, ignore_errors=True)
        return False

def _copy_full_fallback(full_name, dest, mode, date):
    print(f"[copy fallback] {full_name}: 目录树不可用/被截断,"
          f"退回完整克隆后本地裁剪", file=sys.stderr)
    if not clone_repo(full_name, dest):
        return False
    try:
        entries = _walk_entries(dest)
        keep, omitted = select_entries(entries, mode)
        kept = {k["path"] for k in keep}
        for e in entries:
            if e["path"] not in kept:
                try:
                    os.remove(os.path.join(dest, *e["path"].split("/")))
                except OSError:
                    pass
        _prune_empty_dirs(dest)
        _rewrite_all_md(dest, keep, full_name)
        _write_manifest(os.path.dirname(dest), full_name, mode, keep, omitted,
                        "", date)
        return True
    except Exception as e:
        print(f"[copy failed] {full_name}: {e}", file=sys.stderr)
        return False

def _prune_empty_dirs(dest):
    for root, _, _ in os.walk(dest, topdown=False):
        if root == dest:
            continue
        try:
            if not os.listdir(root):
                os.rmdir(root)
        except OSError:
            pass

def _copy_readme(full_name, dest, date=""):
    shutil.rmtree(dest, ignore_errors=True)
    try:
        j = _run_json(["gh", "api", f"repos/{full_name}/readme"])
        name = j.get("name") or "README.md"
        content = base64.b64decode(j.get("content") or "").decode("utf-8", "replace")
    except Exception as e:
        print(f"[copy failed] {full_name}: {e}", file=sys.stderr)
        return False
    if name.lower().endswith((".md", ".markdown")):
        content = rewrite_md(content, name, full_name)
    os.makedirs(dest, exist_ok=True)
    with open(os.path.join(dest, name), "w", encoding="utf-8") as f:
        f.write(content)
    keep = [{"path": name, "size": int(j.get("size") or len(content)),
             "class": "readme"}]
    tree_text = ""
    try:
        entries, truncated = _tree_api(full_name)
        tree_text = _render_tree([e["path"] for e in entries
                                  if e.get("type") == "blob"], truncated)
        with open(os.path.join(dest, "TREE.md"), "w", encoding="utf-8") as f:
            f.write(tree_text)
        keep.append({"path": "TREE.md", "size": len(tree_text.encode()),
                     "class": "tree"})
    except Exception as e:
        print(f"[tree failed] {full_name}: {e}", file=sys.stderr)
    commit = ""
    try:
        commit = _run_text(["gh", "api", f"repos/{full_name}/commits/HEAD",
                            "--jq", ".sha"])
    except Exception:
        pass
    _write_manifest(os.path.dirname(dest), full_name, "readme", keep, [], commit, date)
    return True

def _render_tree(paths, truncated=False, limit=3000):
    """把路径列表渲染成缩进目录树。"""
    paths = sorted(set(paths))
    if not paths:
        return "# 目录树\n\n(空仓库)\n"
    lines = ["# 目录树", ""]
    total = len(paths)
    if truncated:
        lines.append("> 目录过大,已被 tree API 截断,以下不完整。")
    if total > limit:
        lines.append(f"> 共 {total} 个文件,只显示前 {limit} 个。")
        paths = paths[:limit]
    root = {"dirs": {}, "files": []}
    for p in paths:
        parts = p.split("/")
        node = root
        for d in parts[:-1]:
            node = node["dirs"].setdefault(d, {"dirs": {}, "files": []})
        node["files"].append(parts[-1])

    def walk(node, prefix):
        dirs = sorted(node["dirs"])
        files = sorted(node["files"])
        entries = [(d, True) for d in dirs] + [(f, False) for f in files]
        for i, (name, is_dir) in enumerate(entries):
            last = i == len(entries) - 1
            lines.append(prefix + ("└── " if last else "├── ") + name + ("/" if is_dir else ""))
            if is_dir:
                walk(node["dirs"][name], prefix + ("    " if last else "│   "))

    walk(root, "")
    return "\n".join(lines) + "\n"

def _make_copy(mode, date):
    if mode == "full":
        return clone_repo
    return lambda full_name, dest: copy_corpus(full_name, dest, mode, date)

# ── 归档主流程 ──────────────────────────────────────────────────────

def archive_one(row, copy, wiki, date, stat_date=None, in_week=None):
    full_name = row["full_name"]
    base = repo_dir(full_name)
    os.makedirs(base, exist_ok=True)
    code_dir = os.path.join(base, "code")
    shutil.rmtree(code_dir, ignore_errors=True)   # 语料每次重建,不叠加旧文件
    copy_ok = copy(full_name, code_dir)
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
            "one_liner": row.get("one_liner", ""), "copy_ok": copy_ok, "wikis": wikis}

def _safe_archive_one(row, copy, wiki, date, stat_date=None, in_week=None):
    try:
        return archive_one(row, copy, wiki, date, stat_date, in_week)
    except Exception as e:
        print(f"[archive failed] {row['full_name']}: {e}", file=sys.stderr)
        return {"full_name": row["full_name"], "score": row.get("score"),
                "one_liner": row.get("one_liner", ""), "copy_ok": False,
                "wikis": {s: "failed" for s in WIKI_SOURCES}}

def run(cfg, date, copy=None, wiki=None, stat_date=None, in_week=None):
    wiki = wiki or export_wiki
    mode = str(cfg.get("archive", {}).get("repo_copy", "code"))
    if mode not in ("code", "full", "readme"):
        mode = "code"
    copy = copy or _make_copy(mode, date)
    min_score = int(cfg.get("llm", {}).get("min_score", 5))
    scored_path = common.data_dir("scored", f"{date}.json")
    scored = common.load_json(scored_path) if os.path.exists(scored_path) else []
    keep = threshold(scored, min_score)
    # 低并发:网络工具(tree API/degit/repowiki-cli)串行太慢,但太高会触发 repowiki 限流。
    workers = max(1, int(cfg.get("archive", {}).get("concurrency", 4)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        archived = list(ex.map(
            lambda r: _safe_archive_one(r, copy, wiki, date, stat_date, in_week), keep))
    for r in archived:
        r["copy"] = mode
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
    ok_copy = sum(1 for r in out if r["copy_ok"])
    parts = [f"{len(out)} archived; copy ok {ok_copy}/{len(out)}"]
    for s in WIKI_SOURCES:
        counts = {st: sum(1 for r in out if r["wikis"].get(s) == st)
                  for st in ("ok", "submitted", "not_indexed", "failed", "reused")}
        nonzero = " ".join(f"{st}:{n}" for st, n in counts.items() if n)
        parts.append(f"{s} {nonzero or '-'}")
    print("; ".join(parts))

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

def export_wiki(source, full_name, dest):
    ok, err = _run_retry(["repowiki-cli", source, "cp", full_name, dest], f"wiki:{source}")
    if ok:
        return "ok"
    if source == "zread" and _is_not_indexed(err):
        sub_ok, _ = _submit_zread(full_name)
        if sub_ok:
            return "submitted"
    return "not_indexed" if _is_not_indexed(err) else "failed"

def _submit_zread(full_name):
    if not os.environ.get("ZREAD_TOKEN"):
        return False, None
    return _run_retry(["repowiki-cli", "zread", "submit", full_name], "wiki:zread:submit")

if __name__ == "__main__":
    main()
