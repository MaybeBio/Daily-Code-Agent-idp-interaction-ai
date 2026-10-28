"""pfam_domains — shift every stored coordinate by −1 (review DOM-2)

Revision ID: a4c2e7d9b1f3
Revises: f1a2d6c93e07
Create Date: 2026-10-19

``app/analysis/pfam_domains.py`` used to add +1 to pyhmmer's alignment
coordinates (``hmm_from``/``hmm_to``/``target_from``/``target_to``), which are
already 1-based, inclusive (they mirror HMMER's ``P7_ALIDISPLAY`` and equal the
``hmmscan --domtblout`` "hmm from/to" and "ali from/to" columns).  Every row
written by that code is therefore shifted by +1 on all four columns.

DB signature of the bug (2026-10-19, 436,318 rows): 192,478 rows with
``hmm_end = hmm_length + 1``; zero rows with ``seq_start = 1`` (min = 2).
Example: isoform 1211 Pkinase — ``hmmscan`` says ali 8–324 / hmm 1–255,
the table held 9–325 / 2–256.

This is a one-off DATA migration.  It is guarded so it only fires on a table
that still carries the +1 signature (no ``seq_start = 1`` row AND at least one
``hmm_end = hmm_length + 1`` row); a table scanned with the corrected code, or
already migrated, is left untouched.  That makes it safe both on a DB upgraded
in place and on a prod DB that is later overwritten by a ``pg_restore`` of a
local dump (whichever state the dump is in, ``alembic upgrade head`` does the
right thing).

``domain_coverage`` is a difference of HMM coordinates and is unaffected.
``alteration_type`` (truncated_n / truncated_c) is set at scan time from
``hmm_start``/``hmm_end`` vs 10 %/90 % of ``hmm_length`` and is NOT recomputed
here; the −1 shift can in principle move a handful of ``truncated`` rows that
sat exactly on those thresholds, so recompute the label (or re-scan) if an
exact truncated_n/_c count matters.
"""
from alembic import op
import sqlalchemy as sa

revision = "a4c2e7d9b1f3"
down_revision = "f1a2d6c93e07"
branch_labels = None
depends_on = None


def _has_plus_one_signature(conn) -> bool:
    row = conn.execute(sa.text(
        "SELECT "
        "  count(*) FILTER (WHERE seq_start = 1)               AS starts_at_one, "
        "  count(*) FILTER (WHERE hmm_end = hmm_length + 1)    AS over_by_one "
        "FROM pfam_domains"
    )).one()
    return row.starts_at_one == 0 and row.over_by_one > 0


def upgrade():
    conn = op.get_bind()
    if not _has_plus_one_signature(conn):
        print("[a4c2e7d9b1f3] pfam_domains has no +1 signature — nothing to shift")
        return
    res = conn.execute(sa.text(
        "UPDATE pfam_domains SET "
        "  seq_start = seq_start - 1, "
        "  seq_end   = seq_end   - 1, "
        "  hmm_start = hmm_start - 1, "
        "  hmm_end   = hmm_end   - 1"
    ))
    print(f"[a4c2e7d9b1f3] shifted {res.rowcount} pfam_domains rows by -1")


def downgrade():
    conn = op.get_bind()
    # Only re-apply the +1 if the table is in the corrected state (otherwise
    # we would double-shift a never-migrated table).
    if _has_plus_one_signature(conn):
        print("[a4c2e7d9b1f3] pfam_domains already carries the +1 signature — not shifting")
        return
    res = conn.execute(sa.text(
        "UPDATE pfam_domains SET "
        "  seq_start = seq_start + 1, "
        "  seq_end   = seq_end   + 1, "
        "  hmm_start = hmm_start + 1, "
        "  hmm_end   = hmm_end   + 1"
    ))
    print(f"[a4c2e7d9b1f3] re-shifted {res.rowcount} pfam_domains rows by +1")
