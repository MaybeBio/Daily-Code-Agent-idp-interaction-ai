from scripts import logs

WRAPPED = (
    "Fetching events for target(s): a, b...\n\n"
    "2026-08-14 07:52:25 | \U0001F680 thematrixmaster pushed to \n"
    "thematrixmaster/thematrixmaster.github.io\n"
    "    - [9177129] (expanded) Fix author self-highlighting\n"
    "2026-08-13 23:32:28 | ⭐️ QizhiPei starred deepseek-ai/deepseek-harness\n"
)

def test_unwrap_rejoins_wrapped_record():
    recs = logs.unwrap_records(WRAPPED)
    joined = [r for r in recs if r.startswith("2026-08-14")]
    assert len(joined) == 1
    assert "thematrixmaster.github.io" in joined[0]

def test_parse_push_and_star():
    recs = logs.unwrap_records(WRAPPED)
    evs = [logs.parse_event(r) for r in recs]
    evs = [e for e in evs if e]
    push = [e for e in evs if e["kind"] == "push"][0]
    assert push["repo"] == "thematrixmaster/thematrixmaster.github.io"
    assert push["actor"] == "thematrixmaster"
    star = [e for e in evs if e["kind"] == "star"][0]
    assert star["repo"] == "deepseek-ai/deepseek-harness"

def test_branch_name_not_mistaken_for_repo():
    # 分支名含 '/'，不应被当作仓库；真实仓库在 'at' 之后
    rec = "2026-08-13 22:02:24 | \U0001F195 Shenggan created branch 'feat/task-timing-slot' at Shenggan/pypto"
    ev = logs.parse_event(rec)
    assert ev["kind"] == "branch"
    assert ev["repo"] == "Shenggan/pypto"

def test_repos_from_events_dedups_and_groups():
    recs = ["2026-08-13 23:32:28 | ⭐️ A starred x/y",
            "2026-08-13 23:40:00 | \U0001F680 A pushed to x/y"]
    out = logs.repos_from_events(logs.extract_events("\n".join(recs)))
    assert out["x/y"]["who"] == ["A"]
    assert set(out["x/y"]["kinds"]) == {"star", "push"}

def test_push_repo_not_confused_by_slash_in_commit():
    rec = ("2026-08-13 13:56:16 | \U0001F680 jnwei pushed to aqlaboratory/openfold-3 "
           "- [fb027a4] (expanded) Merge pull request #358 from CesarPuentes/feat/user-default-runner-yaml")
    ev = logs.parse_event(rec)
    assert ev["repo"] == "aqlaboratory/openfold-3"

def test_repo_at_end_when_issue_title_has_colon():
    # 仓库在 'in' 之后、标题之前,取 headline 上最后一个 owner/repo
    rec = "2026-10-05 20:04:09 | \U0001F4AC milot-mirdita created issue 'Riboseek web server jobs failing' in steineggerlab/foldseek"
    assert logs.parse_event(rec)["repo"] == "steineggerlab/foldseek"

def test_review_event_on_repo_parsed():
    # 'on' 句式(旧锚点表未覆盖,现应解析到仓库)
    rec = "2026-10-05 22:52:29 | \U0001F539 jnwei performed PullRequestReviewEvent on aqlaboratory/openfold-3"
    assert logs.parse_event(rec)["repo"] == "aqlaboratory/openfold-3"
