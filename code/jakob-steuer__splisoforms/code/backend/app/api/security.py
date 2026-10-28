"""
Optional API-key authentication + per-key rate limiting for the public /api/v1 surface.

The public API is OPEN: anonymous requests are allowed and throttled by the nginx
per-IP edge limit. A key is the higher-volume tier — supply `Authorization: Bearer
<key>` to lift the per-IP cap (keyed requests bypass that nginx zone) up to the key's
own per-minute/per-day quota. Keys are issued by scripts/mint_api_key.py; only the
SHA-256 hash is stored. Use `optional_api_key` for open endpoints, `require_api_key`
where a key must be mandatory.
"""

from __future__ import annotations

import hashlib
from datetime import datetime

from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sqlalchemy.orm import Session

from app.database import get_db
from app.models.bio import ApiKey

KEY_PREFIX = "spl_live_"          # issued keys look like spl_live_<random>
_bearer = HTTPBearer(auto_error=False, description="Public API key as a Bearer token")


def hash_key(raw: str) -> str:
    """SHA-256 hex of a raw key. The raw key is never stored."""
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _authenticate(db: Session, raw_key: str) -> ApiKey:
    """
    Validate a Bearer key, enforce its per-key quota, and bump usage counters.
    Raises 401 (invalid), 403 (revoked), or 429 (quota exceeded).
    """
    key_hash = hash_key(raw_key.strip())

    # Row lock so concurrent requests on the same key increment safely.
    key: ApiKey | None = (
        db.query(ApiKey)
        .filter(ApiKey.key_hash == key_hash)
        .with_for_update()
        .first()
    )
    if key is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid API key.",
            headers={"WWW-Authenticate": "Bearer"},
        )
    if key.revoked:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="API key has been revoked.")

    now = datetime.utcnow()
    minute = now.replace(second=0, microsecond=0)
    today = now.date()

    # Lazily roll the windows over.
    if key.minute_window != minute:
        key.minute_window = minute
        key.minute_count = 0
    if key.day_window != today:
        key.day_window = today
        key.day_count = 0

    if key.minute_count >= key.rate_per_min:
        db.commit()  # persist any window reset
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=f"Rate limit exceeded ({key.rate_per_min}/min).",
            headers={"Retry-After": "60"},
        )
    if key.day_count >= key.rate_per_day:
        db.commit()
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=f"Daily quota exceeded ({key.rate_per_day}/day).",
            headers={"Retry-After": "3600"},
        )

    key.minute_count += 1
    key.day_count += 1
    key.request_count += 1
    key.last_used_at = now
    db.commit()
    return key


def optional_api_key(
    creds: HTTPAuthorizationCredentials | None = Depends(_bearer),
    db: Session = Depends(get_db),
) -> ApiKey | None:
    """
    Public read access is OPEN: anonymous requests are allowed and throttled by the
    nginx per-IP edge limit (zone=public_api). If a Bearer key IS supplied it is
    validated and its per-key quota is applied (the higher-volume tier); a
    supplied-but-invalid key is rejected rather than silently downgraded to anonymous.
    """
    if creds is None or not creds.credentials:
        return None
    return _authenticate(db, creds.credentials)


def require_api_key(
    creds: HTTPAuthorizationCredentials | None = Depends(_bearer),
    db: Session = Depends(get_db),
) -> ApiKey:
    """Strict variant: a valid Bearer key is mandatory (401 if missing)."""
    if creds is None or not creds.credentials:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key. Send 'Authorization: Bearer <key>'.",
            headers={"WWW-Authenticate": "Bearer"},
        )
    return _authenticate(db, creds.credentials)
