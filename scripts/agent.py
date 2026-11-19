import json, os, random, threading, time
import yaml
from openai import OpenAI, APIStatusError, APIConnectionError

DEFAULT_MODEL = "deepseek-chat"

# 单次请求超时。不传 timeout 时 openai SDK 默认 600s 且自带 2 次自动重试:
# 一个卡住的请求最坏能占 30 分钟,并发里几个一起卡,整步就会静默跑满 timeout-minutes。
# 这里压到 90s 并关掉 SDK 重试,退避统一交给 _retry_json。
REQUEST_TIMEOUT_SECONDS = 90

# 全局限速:所有线程共享同一个最小间隔。服务商实测限额约 10 次/分钟,
# 取 7.5s/次(8 次/分钟)留余量。
# 注意:真正的吞吐由这里封顶,并发数只决定「同时等待的线程数」。
# 因此不要把 concurrency 调很低(请求通常 5~15s,2 个线程喂不满 8 次/分钟);
# 用 LLM_MIN_INTERVAL 可以覆盖间隔(秒),设为 0 关闭限速。
MIN_REQUEST_INTERVAL_SECONDS = 7.5
_rate_lock = threading.Lock()
_next_slot = 0.0

def load_prompts(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def make_client() -> OpenAI:
    return OpenAI(base_url=os.environ["LLM_BASE_URL"], api_key=os.environ["LLM_API_KEY"],
                  timeout=REQUEST_TIMEOUT_SECONDS, max_retries=0)

def model_name() -> str:
    return os.environ.get("LLM_MODEL", DEFAULT_MODEL).strip() or DEFAULT_MODEL

def _rate_limit():
    """全局限速:按最小间隔给每次请求排一个时间槽,并发再高也不会同时打满限额。"""
    interval = float(os.environ.get("LLM_MIN_INTERVAL", MIN_REQUEST_INTERVAL_SECONDS))
    if interval <= 0:
        return
    global _next_slot
    with _rate_lock:
        now = time.monotonic()
        wait = max(0.0, _next_slot - now)
        _next_slot = max(now, _next_slot) + interval
    if wait > 0:
        time.sleep(wait)

def _retry_after(err):
    """429 响应里的 Retry-After(秒);没有就返回 None。"""
    headers = getattr(getattr(err, "response", None), "headers", None)
    raw = headers.get("retry-after") if headers else None
    try:
        return max(0.0, float(raw)) if raw is not None else None
    except (TypeError, ValueError):
        return None

# 4xx 是客户端错误(key 错/model 名错/请求格式错),再试也一样,立即放弃;
# 429 限流、408 超时、5xx 服务端错都是瞬时,留给重试。
def _is_permanent(err) -> bool:
    return (isinstance(err, APIStatusError)
            and 400 <= err.status_code < 500
            and err.status_code not in (408, 429))

def _retry_json(client, model, msgs, parse, attempts=8):
    """JSON 模式调 LLM:全局限速 + 指数退避,429 优先遵循 Retry-After。

    连接层错误(APIConnectionError,含超时)最多试 3 次:端点不可达时再试只是白等;
    429/5xx 按 attempts(默认 8)退避重试。
    """
    conn_attempts = min(3, attempts)
    conn_fail = 0
    for attempt in range(attempts):
        try:
            _rate_limit()
            r = client.chat.completions.create(model=model, messages=msgs, temperature=0.0,
                                               response_format={"type": "json_object"})
            return parse(json.loads(r.choices[0].message.content))
        except Exception as e:      # 网络/限流/解析
            if _is_permanent(e):
                raise
            if isinstance(e, APIConnectionError):
                conn_fail += 1
                if conn_fail >= conn_attempts or attempt == attempts - 1:
                    raise
                time.sleep(min(2 ** conn_fail, 30) + random.random())
                continue
            if attempt == attempts - 1:
                raise
            delay = _retry_after(e)
            if delay is None:
                delay = min(2 ** (attempt + 2), 60) + random.random()
            time.sleep(delay)

def score_repo(client, model, prompts, readme) -> dict:
    p = prompts["score"]
    msgs = [{"role": "system", "content": p["system"]},
            {"role": "user", "content": p["user"].format(readme=readme[:10000])}]
    return _retry_json(client, model, msgs,
                       lambda obj: {"score": int(obj["score"]),
                                    "one_liner": str(obj.get("one_liner", ""))})

def build_code_card(client, model, prompts, readme) -> dict:
    p = prompts["code_card"]
    msgs = [{"role": "system", "content": p["system"]},
            {"role": "user", "content": p["user"].format(readme=readme)}]
    return _retry_json(client, model, msgs, lambda obj: {"card": str(obj.get("card", ""))})
