import pytest

from scripts import agent


@pytest.fixture(autouse=True)
def _disable_llm_rate_limit(monkeypatch):
    """测试里不要真的按 7.5s 排队。

    限速逻辑由 test_agent 的 test_rate_limit_spaces_requests 单独验证;
    其余测试若共享模块级 _next_slot,会被前一个测试推后的时间槽拖住几十秒。
    """
    monkeypatch.setenv("LLM_MIN_INTERVAL", "0")
    monkeypatch.setattr(agent, "_next_slot", 0.0)
