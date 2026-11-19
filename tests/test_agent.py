import types
import httpx
import openai
import pytest
from scripts import agent

def _resp(content):
    msg = types.SimpleNamespace(content=content)
    return types.SimpleNamespace(choices=[types.SimpleNamespace(message=msg)])

def _client(seq):
    """seq:依次返回一个 content 字符串,或一个要抛出的异常实例。"""
    calls = []
    def create(**kw):
        item = seq[len(calls)]
        calls.append(item)
        if isinstance(item, Exception):
            raise item
        return _resp(item)
    client = types.SimpleNamespace(chat=types.SimpleNamespace(
        completions=types.SimpleNamespace(create=create)))
    return client, calls

def _status(code, headers=None):
    return openai.APIStatusError("e", body=None,
                                 response=httpx.Response(
                                     code, headers=headers or {},
                                     request=httpx.Request("POST", "http://x")))

def test_retry_json_retries_transient_then_succeeds(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([RuntimeError("net"), '{"ok": 1}'])
    assert agent._retry_json(client, "m", [], lambda o: o) == {"ok": 1}
    assert len(calls) == 2

def test_retry_json_gives_up_after_attempts(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([RuntimeError("net")] * 8)
    with pytest.raises(RuntimeError):
        agent._retry_json(client, "m", [], lambda o: o)
    assert len(calls) == 8                                  # 默认 8 次尝试

def test_retry_json_short_circuits_permanent(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([_status(401), '{"ok": 1}'])
    with pytest.raises(openai.APIStatusError):
        agent._retry_json(client, "m", [], lambda o: o)
    assert len(calls) == 1                                  # key 错 → 不重试

def test_retry_json_retries_429(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([_status(429), '{"ok": 1}'])
    assert agent._retry_json(client, "m", [], lambda o: o) == {"ok": 1}
    assert len(calls) == 2                                  # 限流 → 重试

@pytest.mark.parametrize("code,perm", [(400, True), (401, True), (403, True), (404, True),
                                       (408, False), (429, False), (500, False)])
def test_is_permanent_by_status(code, perm):
    assert agent._is_permanent(_status(code)) is perm

def test_is_permanent_ignores_non_openai_error():
    assert agent._is_permanent(RuntimeError("x")) is False

def test_make_client_sets_request_timeout_and_disables_sdk_retries(monkeypatch):
    seen = {}
    def fake_openai(**kw):
        seen.update(kw)
        return object()
    monkeypatch.setattr(agent, "OpenAI", fake_openai)
    monkeypatch.setenv("LLM_BASE_URL", "http://x")
    monkeypatch.setenv("LLM_API_KEY", "k")
    agent.make_client()
    # 不设的话 SDK 默认 600s/请求 + 2 次自动重试,卡住的请求会拖垮整步
    assert seen["timeout"] == agent.REQUEST_TIMEOUT_SECONDS
    assert seen["max_retries"] == 0


def test_retry_json_honors_retry_after(monkeypatch):
    slept = []
    monkeypatch.setattr(agent, "_rate_limit", lambda: None)
    monkeypatch.setattr(agent.time, "sleep", lambda s: slept.append(round(s, 3)))
    client, calls = _client([_status(429, {"retry-after": "17"}), '{"ok": 1}'])
    assert agent._retry_json(client, "m", [], lambda o: o) == {"ok": 1}
    assert slept == [17.0]                                  # 按服务端给的间隔等


def test_retry_json_retries_connection_then_gives_up(monkeypatch):
    monkeypatch.setattr(agent, "_rate_limit", lambda: None)
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    conn = openai.APIConnectionError(request=httpx.Request("POST", "http://x"))
    client, calls = _client([conn] * 5)
    with pytest.raises(openai.APIConnectionError):
        agent._retry_json(client, "m", [], lambda o: o)
    assert len(calls) == 3                                  # 连接错误最多试 3 次


def test_rate_limit_spaces_requests(monkeypatch):
    monkeypatch.setattr(agent, "_next_slot", 0.0)
    monkeypatch.setenv("LLM_MIN_INTERVAL", "7.5")
    clock = {"t": 100.0}
    slept = []

    def fake_sleep(s):
        slept.append(round(s, 3))
        clock["t"] += s

    monkeypatch.setattr(agent.time, "monotonic", lambda: clock["t"])
    monkeypatch.setattr(agent.time, "sleep", fake_sleep)
    agent._rate_limit()            # 第一个:立即放行,占 100→107.5
    clock["t"] = 101.0
    agent._rate_limit()            # 第二个:等到 107.5,睡 6.5
    assert slept == [6.5]
