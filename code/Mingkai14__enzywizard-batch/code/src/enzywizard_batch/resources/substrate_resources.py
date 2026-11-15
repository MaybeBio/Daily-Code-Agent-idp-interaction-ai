import re

BRACKET_RE = re.compile(r"[\(\[\{]([^)\]\}]+)[\)\]\}]")
SEP_RE = re.compile(r"[\u2010\u2011\u2012\u2013\u2014\u2015\-\_\/]+")
MID_PUNCT_RE = re.compile(r"[,:;]+")
TRAILING_PUNCT_RE = re.compile(r"[，,;；。\.:\s]+$")

PREFIX_PATTERNS = [
    r"^l[\-\s]+",
    r"^d[\-\s]+",
    r"^dl[\-\s]+",
]

SUFFIX_PHRASES = [
    r"monohydrate",
    r"dihydrate",
    r"trihydrate",
    r"tetrahydrate",
    r"hydrate",
    r"sodium\s+salt",
    r"potassium\s+salt",
    r"lithium\s+salt",
    r"ammonium\s+salt",
    r"hydrochloride",
    r"chloride",
]