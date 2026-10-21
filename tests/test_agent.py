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

def _status(code):
    return openai.APIStatusError("e", body=None,
                                 response=httpx.Response(code, request=httpx.Request("POST", "http://x")))

def test_retry_json_retries_transient_then_succeeds(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([RuntimeError("net"), '{"ok": 1}'])
    assert agent._retry_json(client, "m", [], lambda o: o) == {"ok": 1}
    assert len(calls) == 2

def test_retry_json_gives_up_after_attempts(monkeypatch):
    monkeypatch.setattr(agent.time, "sleep", lambda s: None)
    client, calls = _client([RuntimeError("net")] * 5)
    with pytest.raises(RuntimeError):
        agent._retry_json(client, "m", [], lambda o: o)
    assert len(calls) == 5                                  # 5 次尝试

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
