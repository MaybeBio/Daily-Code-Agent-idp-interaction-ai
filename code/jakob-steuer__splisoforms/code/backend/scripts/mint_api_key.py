#!/usr/bin/env python3
"""
mint_api_key.py — issue a public API key for /api/v1.

Generates a random key, stores only its SHA-256 hash, and prints the plaintext
ONCE (it can never be recovered afterwards). Hand the printed key to the user.

Usage:
    python scripts/mint_api_key.py --label "Dr. Smith lab" --owner smith@uni.edu
    python scripts/mint_api_key.py --label "high volume" --rate-per-min 120 --rate-per-day 20000
    python scripts/mint_api_key.py --revoke spl_live_ab12   # revoke by prefix
"""
from __future__ import annotations

import argparse
import secrets
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from app.database import SessionLocal              # noqa: E402
from app.models.bio import ApiKey                  # noqa: E402
from app.api.security import hash_key, KEY_PREFIX  # noqa: E402


def mint(label: str, owner: str | None, tier: str, per_min: int, per_day: int) -> None:
    raw = KEY_PREFIX + secrets.token_urlsafe(32)
    db = SessionLocal()
    try:
        key = ApiKey(
            key_hash=hash_key(raw),
            key_prefix=raw[: len(KEY_PREFIX) + 4],
            label=label,
            owner_email=owner,
            tier=tier,
            rate_per_min=per_min,
            rate_per_day=per_day,
        )
        db.add(key)
        db.commit()
        print("\n  API key created — copy it now, it will NOT be shown again:\n")
        print(f"    {raw}\n")
        print(f"  prefix : {key.key_prefix}")
        print(f"  label  : {label}")
        print(f"  limits : {per_min}/min, {per_day}/day")
        print("\n  Usage:  curl -H 'Authorization: Bearer <key>' https://<host>/api/v1/genes/TP53\n")
    finally:
        db.close()


def revoke(prefix: str) -> None:
    db = SessionLocal()
    try:
        n = db.query(ApiKey).filter(ApiKey.key_prefix == prefix).update({ApiKey.revoked: True})
        db.commit()
        print(f"revoked {n} key(s) with prefix {prefix}")
    finally:
        db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", help="Human-readable label for the key.")
    ap.add_argument("--owner", help="Owner email.")
    ap.add_argument("--tier", default="standard")
    ap.add_argument("--rate-per-min", type=int, default=600)
    ap.add_argument("--rate-per-day", type=int, default=100000)
    ap.add_argument("--revoke", metavar="PREFIX", help="Revoke key(s) by prefix instead of minting.")
    args = ap.parse_args()

    if args.revoke:
        revoke(args.revoke)
    else:
        if not args.label:
            ap.error("--label is required when minting")
        mint(args.label, args.owner, args.tier, args.rate_per_min, args.rate_per_day)


if __name__ == "__main__":
    main()
