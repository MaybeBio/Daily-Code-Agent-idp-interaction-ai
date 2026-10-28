"""add functional_class + alteration_type on pfam_domains; backfill from existing data

Revision ID: c3e9a4b7d2f1
Revises: a7c4f9d2e1b3
Create Date: 2026-05-19 11:00:00.000000

Two new orthogonal annotation axes on pfam_domains:

  functional_class    catalytic_kinase | catalytic_phosphatase |
                      catalytic_protease | catalytic_other |
                      ub_ligase | dna_binding | rna_binding |
                      signaling_adapter | oligomerization | other

  alteration_type     intact | truncated_n | truncated_c | truncated_both |
                      internal_disruption | lost | gained | unknown |
                      truncated_unknown

Both fields are derivable from existing pfam_domains columns; this migration
backfills the entire table in-place — no HMMER re-scan is needed.

The legacy `is_kinase` boolean is left intact (set wherever
functional_class == 'catalytic_kinase' AND the legacy criteria matched).
"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "c3e9a4b7d2f1"
down_revision: Union[str, Sequence[str], None] = "a7c4f9d2e1b3"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("pfam_domains", sa.Column("functional_class", sa.String(40), nullable=True))
    op.add_column("pfam_domains", sa.Column("alteration_type",  sa.String(30), nullable=True))
    op.create_index(
        "ix_pfam_domains_functional_class", "pfam_domains", ["functional_class"],
    )
    op.create_index(
        "ix_pfam_domains_alteration_type",  "pfam_domains", ["alteration_type"],
    )

    # ── alteration_type — pure SQL, derived from existing columns ────────────
    # Logic mirrors functional_classes.alteration_type():
    #   lost      → 'lost'
    #   gained    → 'gained'
    #   present   → 'intact'
    #   truncated → n_cut + c_cut + coverage decide between
    #               truncated_n / truncated_c / truncated_both /
    #               internal_disruption / truncated_unknown
    op.execute("""
        UPDATE pfam_domains
        SET alteration_type = CASE
            WHEN domain_status = 'lost'    THEN 'lost'
            WHEN domain_status = 'gained'  THEN 'gained'
            WHEN domain_status = 'present' THEN 'intact'
            WHEN domain_status = 'truncated' THEN CASE
                WHEN hmm_length IS NULL OR hmm_length <= 0 THEN 'truncated_unknown'
                WHEN hmm_start > hmm_length * 0.10
                 AND hmm_end   < hmm_length * 0.90 THEN
                    CASE WHEN domain_coverage < 0.50
                         THEN 'internal_disruption'
                         ELSE 'truncated_both'
                    END
                WHEN hmm_start > hmm_length * 0.10 THEN 'truncated_n'
                WHEN hmm_end   < hmm_length * 0.90 THEN 'truncated_c'
                ELSE 'internal_disruption'
            END
            ELSE 'unknown'
        END;
    """)

    # ── functional_class — load the curated mapping in Python, push it via a
    # parametrised UPDATE per (accession, clan) combination.  This keeps the
    # migration deterministic and easy to audit.
    bind = op.get_bind()

    from app.analysis.functional_classes import (
        ACCESSION_TO_CLASS, CLAN_TO_CLASS,
    )

    # 1. Apply accession-level overrides (most specific)
    for acc, cls in ACCESSION_TO_CLASS.items():
        bind.execute(
            sa.text("UPDATE pfam_domains SET functional_class = :cls WHERE pfam_accession = :acc"),
            {"cls": cls, "acc": acc},
        )

    # 2. Apply clan-level mapping where accession didn't already set it
    for clan, cls in CLAN_TO_CLASS.items():
        bind.execute(
            sa.text("""UPDATE pfam_domains SET functional_class = :cls
                       WHERE pfam_clan = :clan AND functional_class IS NULL"""),
            {"cls": cls, "clan": clan},
        )

    # 3. Default everything else to 'other'
    bind.execute(sa.text(
        "UPDATE pfam_domains SET functional_class = 'other' WHERE functional_class IS NULL"
    ))


def downgrade() -> None:
    op.drop_index("ix_pfam_domains_alteration_type", table_name="pfam_domains")
    op.drop_index("ix_pfam_domains_functional_class", table_name="pfam_domains")
    op.drop_column("pfam_domains", "alteration_type")
    op.drop_column("pfam_domains", "functional_class")
