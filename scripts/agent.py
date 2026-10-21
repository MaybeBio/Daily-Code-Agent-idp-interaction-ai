import json, os, random, time
import yaml
from openai import OpenAI, APIStatusError

DEFAULT_MODEL = "deepseek-chat"

def load_prompts(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def make_client() -> OpenAI:
    return OpenAI(base_url=os.environ["LLM_BASE_URL"], api_key=os.environ["LLM_API_KEY"])

def model_name() -> str:
    return os.environ.get("LLM_MODEL", DEFAULT_MODEL).strip() or DEFAULT_MODEL

# 4xx 是客户端错误(key 错/model 名错/请求格式错),再试也一样,立即放弃;
# 429 限流、408 超时、5xx 服务端错都是瞬时,留给重试。
def _is_permanent(err) -> bool:
    return (isinstance(err, APIStatusError)
            and 400 <= err.status_code < 500
            and err.status_code not in (408, 429))

def _retry_json(client, model, msgs, parse, attempts=5):
    """JSON 模式调 LLM:5 次指数退避重试,叠加 jitter 打散高并发下的惊群。"""
    for attempt in range(attempts):
        try:
            r = client.chat.completions.create(model=model, messages=msgs, temperature=0.0,
                                               response_format={"type": "json_object"})
            return parse(json.loads(r.choices[0].message.content))
        except Exception as e:      # 网络/限流/解析
            if _is_permanent(e) or attempt == attempts - 1:
                raise
            time.sleep(2 ** (attempt + 2) + random.random())

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
