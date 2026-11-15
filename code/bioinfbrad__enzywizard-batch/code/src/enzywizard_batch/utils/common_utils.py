from __future__ import annotations
from typing import List, Dict, Tuple, Set, Any
import numpy as np
import json
from dataclasses import dataclass
import re

def vec_zeros(n: int) -> List[int]:
    return [0] * n

def one_hot_vec_generator(length:int,base:str,index:Dict[str,int])-> List[int]:
    v = vec_zeros(length)
    if base in index:
        v[index[base]] += 1
    return v

def multi_hot_vec_generator(length: int, base: str,base_to_indices: Dict[str, List[int]]) -> List[int]:
    v: List[int] = vec_zeros(length)

    base = base.upper()
    for idx in base_to_indices.get(base, []):
        v[idx] = 1

    return v

def convert_to_json_serializable(obj):
    if isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}

    if isinstance(obj, (list, tuple, set)):
        return [convert_to_json_serializable(v) for v in obj]

    if isinstance(obj, np.ndarray):
        return convert_to_json_serializable(obj.tolist())

    if isinstance(obj, np.generic):
        return obj.item()

    return obj

@dataclass(frozen=True)
class RawJSON:
    """Wrap a pre-rendered JSON snippet that should be embedded as-is."""
    raw: str

class InlineJSONEncoder(json.JSONEncoder):
    _RAW_PREFIX = "__RAWJSON__:"
    _RAW_SUFFIX = ":__ENDRAW__"

    def default(self, obj):
        if isinstance(obj, RawJSON):
            return self._RAW_PREFIX + obj.raw + self._RAW_SUFFIX
        return super().default(obj)

    def _raw_pattern(self) -> re.Pattern:
        # match the *JSON-encoded string token* which includes surrounding quotes
        # e.g. "__RAWJSON__:[1,2,3]:__ENDRAW__"  (with escapes if any)
        return re.compile(
            r"\"%s.*?%s\"" % (re.escape(self._RAW_PREFIX), re.escape(self._RAW_SUFFIX))
        )

    def _unpack_rawjson_token(self, token_with_quotes: str) -> str:
        """
        token_with_quotes is a valid JSON string token (including the surrounding quotes).
        We json.loads it to unescape, then strip prefix/suffix, and return raw JSON.
        """
        val = json.loads(token_with_quotes)  # -> "__RAWJSON__:[...]:__ENDRAW__"
        raw = val[len(self._RAW_PREFIX):-len(self._RAW_SUFFIX)]
        return raw

    def encode(self, obj):
        # Keep a non-stream fallback; safe for small objects
        s = super().encode(obj)
        pat = self._raw_pattern()

        def repl(m: re.Match) -> str:
            return self._unpack_rawjson_token(m.group(0))

        return pat.sub(repl, s)

    def iterencode(self, obj, _one_shot=False):
        """
        Stream-safe: do NOT materialize all chunks.
        We keep a small buffer to handle matches that cross chunk boundaries.
        """
        pat = self._raw_pattern()
        start_seq = '"' + self._RAW_PREFIX  # potential match start

        buf = ""

        for chunk in super().iterencode(obj, _one_shot=_one_shot):
            buf += chunk

            # consume complete matches from the buffer
            while True:
                m = pat.search(buf)
                if not m:
                    break
                # yield text before the token
                if m.start() > 0:
                    yield buf[:m.start()]
                # yield unpacked raw JSON
                yield self._unpack_rawjson_token(m.group(0))
                # drop consumed part
                buf = buf[m.end():]

            # If no match, flush safe prefix of buffer, keep only possible tail
            # that could be the beginning of a RawJSON token.
            idx = buf.rfind(start_seq)
            if idx == -1:
                if buf:
                    yield buf
                    buf = ""
            else:
                # flush everything before the last possible start
                if idx > 0:
                    yield buf[:idx]
                    buf = buf[idx:]

        # flush remaining buffer
        if buf:
            # if something incomplete remains, just output it (shouldn't happen normally)
            # or you can be strict and raise
            yield buf

def wrap_leaf_lists_as_rawjson(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: wrap_leaf_lists_as_rawjson(v) for k, v in obj.items()}

    if isinstance(obj, tuple):
        obj = list(obj)

    if isinstance(obj, list):
        is_leaf_list = all(not isinstance(x, (dict, list, tuple)) for x in obj)
        if is_leaf_list:
            raw = json.dumps(obj, ensure_ascii=False, separators=(", ", ": "))
            return RawJSON(raw)
        return [wrap_leaf_lists_as_rawjson(v) for v in obj]

    return obj

def get_clean_filename(name: str) -> str:
    name = re.sub(r"[^\w\-\.]", "_", name)
    return name[:150]


def get_optimized_filename(name: str) -> str:
    if not isinstance(name, str):
        return ""

    name = name.strip()

    name = re.sub(r"[,;:=+\s]+", "_", name)

    name = re.sub(r"[^\w\-\.]", "_", name)

    name = re.sub(r"_+", "_", name)

    name = name.strip("_.")

    return name[:150]