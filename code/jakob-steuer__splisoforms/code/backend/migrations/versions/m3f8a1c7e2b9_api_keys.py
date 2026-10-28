"""public API access — api_keys table (auth + per-key rate-limit counters)

Revision ID: m3f8a1c7e2b9
Revises: k7e4b3a2c9d1
Create Date: 2026-05-28 16:00:00.000000

Backs the public REST API (/api/v1).  Each consumer gets a key whose SHA-256
hash is stored here (the plaintext is shown once at issuance and never
persisted).  Per-key minute/day counters provide an application-level quota on
top of nginx edge rate-limiting; the counters reset lazily when their window
bucket rolls over.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "m3f8a1c7e2b9"
down_revision: Union[str, Sequence[str], None] = "k7e4b3a2c9d1"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "api_keys",
        sa.Column("id", sa.Integer(), primary_key=True),
        # sha256 hex of the raw key; the raw key is never stored.
        sa.Column("key_hash", sa.String(64), nullable=False, unique=True, index=True),
        # Short non-secret prefix for identification in UIs/logs (e.g. "spl_live_ab12").
        sa.Column("key_prefix", sa.String(20), nullable=False, index=True),
        sa.Column("label", sa.String(255), nullable=True),
        sa.Column("owner_email", sa.String(255), nullable=True),
        sa.Column("tier", sa.String(20), nullable=False, server_default="standard"),
        # Quotas (per key).
        sa.Column("rate_per_min", sa.Integer(), nullable=False, server_default="60"),
        sa.Column("rate_per_day", sa.Integer(), nullable=False, server_default="5000"),
        # Lazily-reset sliding-bucket counters.
        sa.Column("minute_window", sa.DateTime(), nullable=True),
        sa.Column("minute_count", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("day_window", sa.Date(), nullable=True),
        sa.Column("day_count", sa.Integer(), nullable=False, server_default="0"),
        # Lifetime accounting.
        sa.Column("request_count", sa.BigInteger(), nullable=False, server_default="0"),
        sa.Column("revoked", sa.Boolean(), nullable=False, server_default=sa.false(), index=True),
        sa.Column("created_at", sa.DateTime(), server_default=sa.func.now(), nullable=False),
        sa.Column("last_used_at", sa.DateTime(), nullable=True),
    )


def downgrade() -> None:
    op.drop_table("api_keys")
