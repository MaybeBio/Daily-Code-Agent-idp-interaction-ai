"""neoantigen_predictions — re-assert the atlas partial indexes; covering stats index

Revision ID: r1b6e3d8c5a2
Revises: b7e2c4d8f1a0
Create Date: 2026-10-27

The partial indexes from c3f8a1b0d2e5 (ix_neo_representative, ix_neo_repr_peptide,
ix_neo_repr_isoform) existed only as migration SQL. A database built by reset_db.py
(metadata.create_all + `alembic stamp head`) never gets them — that is how the 2026-10
v2 rebuild, and prod after the dump swap, lost all three while alembic reported head.
/neoantigens/top fell from ~0.85 s to 2–4 s: the default sort became a full sort of
the 2.56 M representative rows and the stats block's COUNT(DISTINCT peptide_sequence)
an on-disk merge sort.

The indexes are now also declared on the model (NeoantigenPrediction.__table_args__),
so future rebuilds create them. This migration re-creates them idempotently and
CONCURRENTLY (no lock on the live atlas; ~4 s each on 7.9 M rows), and replaces
ix_neo_repr_peptide with the covering ix_neo_repr_stats: same leading key, INCLUDEs
the stats-block columns, so FILTER counts and COUNT(DISTINCT peptide_sequence) are
index-only scans for every filter on an included column (incl. flag=proteome_unique).

If a CONCURRENT build is interrupted it leaves an INVALID index that a later
IF NOT EXISTS would skip, so the upgrade fails loudly on any invalid index of ours;
`DROP INDEX <name>` and re-run `alembic upgrade head` in that case.
"""
from alembic import op
import sqlalchemy as sa

revision = "r1b6e3d8c5a2"
down_revision = "b7e2c4d8f1a0"
branch_labels = None
depends_on = None

_ATLAS_INDEXES = ("ix_neo_representative", "ix_neo_repr_isoform", "ix_neo_repr_stats")

_STATS_INCLUDE = ("priority_tier, presentation_tier, is_immunogenic, is_expressed, "
                  "is_ba_rescue, is_proteome_unique, netmhcpan_rank_el, ptm_impact_status")


def _fail_on_invalid(names):
    rows = op.get_bind().execute(sa.text(
        "SELECT c.relname FROM pg_index i JOIN pg_class c ON c.oid = i.indexrelid "
        "WHERE NOT i.indisvalid AND c.relname = ANY(:names)"
    ), {"names": list(names)}).scalars().all()
    if rows:
        raise RuntimeError(
            f"invalid (interrupted CONCURRENTLY) index(es) {rows}: DROP INDEX them and re-run"
        )


def upgrade():
    # CREATE/DROP INDEX CONCURRENTLY cannot run inside a transaction block.
    with op.get_context().autocommit_block():
        op.execute(
            "CREATE INDEX CONCURRENTLY IF NOT EXISTS ix_neo_representative "
            "ON neoantigen_predictions (priority_tier, consensus_rank) WHERE is_representative"
        )
        op.execute(
            "CREATE INDEX CONCURRENTLY IF NOT EXISTS ix_neo_repr_isoform "
            "ON neoantigen_predictions (isoform_id) WHERE is_representative"
        )
        op.execute(
            "CREATE INDEX CONCURRENTLY IF NOT EXISTS ix_neo_repr_stats "
            f"ON neoantigen_predictions (peptide_sequence) INCLUDE ({_STATS_INCLUDE}) "
            "WHERE is_representative"
        )
        _fail_on_invalid(_ATLAS_INDEXES)
        # Subsumed by ix_neo_repr_stats (same leading key, same predicate).
        op.execute("DROP INDEX CONCURRENTLY IF EXISTS ix_neo_repr_peptide")


def downgrade():
    with op.get_context().autocommit_block():
        op.execute(
            "CREATE INDEX CONCURRENTLY IF NOT EXISTS ix_neo_repr_peptide "
            "ON neoantigen_predictions (peptide_sequence) WHERE is_representative"
        )
        _fail_on_invalid(("ix_neo_repr_peptide",))
        op.execute("DROP INDEX CONCURRENTLY IF EXISTS ix_neo_repr_stats")
        # ix_neo_representative / ix_neo_repr_isoform stay: c3f8a1b0d2e5 owns them.
