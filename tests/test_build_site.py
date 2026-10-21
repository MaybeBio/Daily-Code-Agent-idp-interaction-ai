import json, os
from scripts import build_site

def test_score_tier():
    assert build_site.score_tier(8) == "high"
    assert build_site.score_tier(5) == "mid"
    assert build_site.score_tier(3) == "low"
    assert build_site.score_tier(None) == "none"

def test_wiki_url():
    assert build_site.wiki_url("deepwiki", "a/b") == "https://deepwiki.com/a/b"
    assert build_site.wiki_url("codewiki", "a/b") == "https://codewiki.google/github.com/a/b"
    assert build_site.wiki_url("zread", "a/b") == "https://zread.ai/a/b"
    assert build_site.wiki_url("nope", "a/b") == ""

def test_load_cards_empty(tmp_path):
    assert build_site.load_cards(str(tmp_path)) == []

def test_load_issue_url_missing(tmp_path):
    assert build_site.load_issue_url(str(tmp_path)) == ""

def test_site_base_path_and_base_prefix():
    assert build_site.site_base_path("") == ""
    assert build_site.site_base_path("https://example.com/repos/foo") == "/repos/foo"

def test_natural_week():
    assert build_site.natural_week("2026-10-01") == ("2026-09-28", "2026-10-04")
    assert build_site.natural_week("2026-10-05") == ("2026-10-05", "2026-10-11")
    assert build_site.natural_week("2026-10-01T00:00:00Z") == ("2026-09-28", "2026-10-04")
    assert build_site.natural_week("bad") == ("", "")
    assert build_site.natural_week("") == ("", "")

def _write_cards(tmp_path, rows):
    d = str(tmp_path / "data" / "cards")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "2026-10-05.json"), "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False)

def test_build_site_renders(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "https://github.com/a/b",
                             "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z",
                             "score": 8, "one_liner": "one",
                             "card": "## 是什么\nhello"}])
    open(str(tmp_path / "data" / "latest_issue.txt"), "w").write("https://github.com/x/y/issues/1")
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    idx = open(str(tmp_path / "site" / "index.html"), encoding="utf-8").read()
    assert "当周 Issue" in idx and "https://github.com/x/y/issues/1" in idx
    arch = open(str(tmp_path / "site" / "archive.html"), encoding="utf-8").read()
    assert "2026-09-28" in arch and "2026-10-04" in arch          # 自然周归档链接
    week = open(str(tmp_path / "site" / "weeks" / "2026-10-04" / "index.html"), encoding="utf-8").read()
    assert "a/b" in week
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "https://deepwiki.com/a/b" in detail and "hello" in detail

def test_build_site_empty_no_crash(tmp_path):
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    idx = open(str(tmp_path / "site" / "index.html"), encoding="utf-8").read()
    assert "本周推荐" in idx

def test_build_site_no_issue_link(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "pushed_at": "t", "score": 8,
                             "one_liner": "x", "card": "## 是什么\nhi"}])
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    idx = open(str(tmp_path / "site" / "index.html"), encoding="utf-8").read()
    assert "当周 Issue" not in idx

def test_build_site_empty_card(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "pushed_at": "t", "score": 8,
                             "one_liner": "x", "card": ""}])
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "Code Card" not in detail

def test_build_site_weeks_dedup_within_week_keep_across(tmp_path):
    d = str(tmp_path / "data" / "cards")
    os.makedirs(d, exist_ok=True)
    def rec(pushed, score):
        return {"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                "pushed_at": pushed, "score": score, "one_liner": "x", "card": ""}
    # 周A(09-28~10-04) 一条;周B(10-05~10-11) 两条(同周去重保留 pushed_at 最新)
    with open(os.path.join(d, "2026-10-03.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-01T00:00:00Z", 7)], f, ensure_ascii=False)
    with open(os.path.join(d, "2026-10-08.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-05T00:00:00Z", 9)], f, ensure_ascii=False)
    with open(os.path.join(d, "2026-10-09.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-06T00:00:00Z", 8)], f, ensure_ascii=False)
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    week_a = open(str(tmp_path / "site" / "weeks" / "2026-10-04" / "index.html"), encoding="utf-8").read()
    week_b = open(str(tmp_path / "site" / "weeks" / "2026-10-11" / "index.html"), encoding="utf-8").read()
    assert "a/b" in week_a and "a/b" in week_b           # 跨自然周都保留
    assert week_b.count('class="repo-card"') == 1        # 同自然周内去重,只剩一条

def test_build_site_week_detail_preserves_old_card(tmp_path):
    d = str(tmp_path / "data" / "cards")
    os.makedirs(d, exist_ok=True)
    def rec(pushed, one, card):
        return {"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                "pushed_at": pushed, "score": 8, "one_liner": one, "card": card}
    with open(os.path.join(d, "2026-09-30.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-01T00:00:00Z", "old", "## 是什么\nOLD")], f, ensure_ascii=False)
    with open(os.path.join(d, "2026-10-06.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-06T00:00:00Z", "new", "## 是什么\nNEW")], f, ensure_ascii=False)
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    week_a = open(str(tmp_path / "site" / "weeks" / "2026-10-04" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    week_b = open(str(tmp_path / "site" / "weeks" / "2026-10-11" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    latest = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "OLD" in week_a and "NEW" not in week_a          # 旧 card 沉淀在当周详情页
    assert "NEW" in week_b and "OLD" not in week_b          # 新周是新的 card
    assert "NEW" in latest and "OLD" not in latest          # 本周/详情页永远是最新
    week_a_index = open(str(tmp_path / "site" / "weeks" / "2026-10-04" / "index.html"), encoding="utf-8").read()
    assert "weeks/2026-10-04/repos/a__b/" in week_a_index   # 周列表链接到当周详情,而非最新

def test_build_site_latest_card_wins(tmp_path):
    d = str(tmp_path / "data" / "cards")
    os.makedirs(d, exist_ok=True)
    older = {"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
             "pushed_at": "t", "score": 7, "one_liner": "old",
             "card": "## 是什么\nOLD"}
    newer = {"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
             "pushed_at": "t", "score": 9, "one_liner": "new",
             "card": "## 是什么\nNEW"}
    with open(os.path.join(d, "2026-09-28.json"), "w", encoding="utf-8") as f:
        json.dump([older], f, ensure_ascii=False)
    with open(os.path.join(d, "2026-10-05.json"), "w", encoding="utf-8") as f:
        json.dump([newer], f, ensure_ascii=False)
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "NEW" in detail and "OLD" not in detail

def test_plain_text_strips_markdown():
    assert build_site._plain_text("**bold** and [link](http://x)") == "bold and link"

def test_build_search_documents_includes_readme(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z", "score": 8,
                             "one_liner": "one", "card": "## 是什么\nhello",
                             "readme": "## Usage\nlatent-space denoising diffusion"}])
    records = build_site.load_cards(str(tmp_path))
    head, deep = build_site._build_search_documents(records)
    assert "latent-space denoising diffusion" in deep[0]["deep"]

def test_build_search_documents_includes_full_readme(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z", "score": 8,
                             "one_liner": "one", "card": "## 是什么\nhello",
                             "readme": "head " + "x" * 50000 + " TAIL"}])
    records = build_site.load_cards(str(tmp_path))
    head, deep = build_site._build_search_documents(records)
    assert "TAIL" in deep[0]["deep"]              # 全文进索引,不截断
    assert len(deep[0]["deep"]) > 50000

def test_build_site_readme_section(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "pushed_at": "t", "score": 8, "one_liner": "x",
                             "card": "## 是什么\nhi", "readme": "## Install\npip install x"}])
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "README 原文" in detail and "pip install x" in detail

def test_build_site_no_readme_section(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "pushed_at": "t", "score": 8, "one_liner": "x",
                             "card": "## 是什么\nhi"}])
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "README 原文" not in detail

def test_build_search_documents(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z", "score": 8,
                             "one_liner": "one", "card": "## 是什么\nhello"}])
    records = build_site.load_cards(str(tmp_path))
    head, deep = build_site._build_search_documents(records)
    assert len(head) == 1 and len(deep) == 1
    d = head[0]
    assert d["id"] == "a/b" and d["title"] == "a/b"
    assert d["url"] == "/repos/a__b/"
    assert "Python" in d["meta"] and "★ 5" in d["meta"]
    assert d["summary"] == "one"
    assert d["date"] == "2026-10-01"
    assert deep[0]["id"] == "a/b"
    assert "hello" in deep[0]["deep"]

def test_build_site_writes_search_files(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z", "score": 8,
                             "one_liner": "one", "card": "## 是什么\nhello"}])
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    site = str(tmp_path / "site")
    assert os.path.isfile(os.path.join(site, "search.html"))
    with open(os.path.join(site, "data", "search.json"), encoding="utf-8") as f:
        head = json.load(f)
    assert head["version"] == 1 and len(head["documents"]) == 1
    with open(os.path.join(site, "data", "search-deep.json"), encoding="utf-8") as f:
        deep = json.load(f)
    assert deep["version"] == 1 and len(deep["documents"]) == 1
    assert os.path.isfile(os.path.join(site, "assets", "search.js"))
    html = open(os.path.join(site, "search.html"), encoding="utf-8").read()
    assert "data-search-root" in html

def _write_wiki(tmp_path, full_name, source, date):
    d = str(tmp_path / "code" / build_site.repo_dir_name(full_name) / "wikis" / date / source)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "llms.txt"), "w", encoding="utf-8").write(f"# {source} {date}")

def test_wiki_dates_upto(tmp_path):
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-10-05")
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-10-12")
    _write_wiki(tmp_path, "a/b", "zread", "2026-10-05")
    assert build_site.wiki_dates_upto(str(tmp_path), "a/b", "") == {
        "deepwiki": "2026-10-12", "zread": "2026-10-05"}
    # 以 upto_date 为上界:冻结在 10-05,看不到 10-12 的快照
    assert build_site.wiki_dates_upto(str(tmp_path), "a/b", "2026-10-05") == {
        "deepwiki": "2026-10-05", "zread": "2026-10-05"}

def test_wiki_dates_upto_missing(tmp_path):
    assert build_site.wiki_dates_upto(str(tmp_path), "a/b", "") == {}

def test_wiki_snapshot_url():
    assert build_site.wiki_snapshot_url(
        "deepwiki", "a/b", "2026-10-12", "MaybeBio/Daily-Code-Agent-Template", "main") == \
        "https://github.com/MaybeBio/Daily-Code-Agent-Template/blob/main/code/a__b/wikis/2026-10-12/deepwiki/llms.txt"
    assert build_site.wiki_snapshot_url("deepwiki", "a/b", "2026-10-12", "") == ""

def test_build_site_wiki_snapshot_link(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "https://github.com/a/b",
                             "language": "Python", "stars": 5,
                             "pushed_at": "2026-10-01T00:00:00Z",
                             "score": 8, "one_liner": "one",
                             "card": "## 是什么\nhello"}])
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-10-05")
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": "",
                                          "github_repo": "MaybeBio/Daily-Code-Agent-Template",
                                          "github_branch": "main"})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "MaybeBio/Daily-Code-Agent-Template/blob/main/code/a__b/wikis/2026-10-05/deepwiki/llms.txt" in detail
    assert "deepwiki 快照" in detail

def test_build_site_no_wiki_snapshot_link_without_github_repo(tmp_path):
    _write_cards(tmp_path, [{"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                             "pushed_at": "t", "score": 8, "one_liner": "x", "card": "## 是什么\nhi"}])
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-10-05")
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": ""})
    detail = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    assert "快照" not in detail

def test_build_site_week_detail_links_frozen_wiki(tmp_path):
    d = str(tmp_path / "data" / "cards")
    os.makedirs(d, exist_ok=True)
    def rec(pushed):
        return {"full_name": "a/b", "url": "u", "language": "Py", "stars": 1,
                "pushed_at": pushed, "score": 8, "one_liner": "x", "card": ""}
    with open(os.path.join(d, "2026-09-30.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-01T00:00:00Z")], f, ensure_ascii=False)
    with open(os.path.join(d, "2026-10-06.json"), "w", encoding="utf-8") as f:
        json.dump([rec("2026-10-06T00:00:00Z")], f, ensure_ascii=False)
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-09-30")
    _write_wiki(tmp_path, "a/b", "deepwiki", "2026-10-06")
    build_site.build_site(str(tmp_path), {"title": "T", "site_base_url": "",
                                          "github_repo": "MaybeBio/Daily-Code-Agent-Template",
                                          "github_branch": "main"})
    week_a = open(str(tmp_path / "site" / "weeks" / "2026-10-04" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    week_b = open(str(tmp_path / "site" / "weeks" / "2026-10-11" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    latest = open(str(tmp_path / "site" / "repos" / "a__b" / "index.html"), encoding="utf-8").read()
    # 周归档页链到「当周冻结」的快照,而非最新那份
    assert "wikis/2026-09-30/deepwiki/llms.txt" in week_a
    assert "wikis/2026-10-06/deepwiki/llms.txt" not in week_a
    assert "wikis/2026-10-06/deepwiki/llms.txt" in week_b
    assert "wikis/2026-09-30/deepwiki/llms.txt" not in week_b
    assert "wikis/2026-10-06/deepwiki/llms.txt" in latest
