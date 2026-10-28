"""add presentation_score and processing_score to neoantigen_predictions

Revision ID: a1b2c3d4e5f6
Revises: f1a2b3c4d5e6
Branch labels: None
Depends on: None

MHCflurry 2.0 upgrade: stores the combined antigen-presentation probability
(presentation_score) and the proteasomal-cleavage component (processing_score)
in addition to the existing binding-affinity fields.

These columns are nullable so existing rows (predicted under the old affinity-only
pipeline) remain valid — they will simply show NULL for the two new columns.
"""
from alembic import op
import sqlalchemy as sa

revision     = "a1b2c3d4e5f6"
down_revision = "f1a2b3c4d5e6"
branch_labels = None
depends_on    = None


def upgrade() -> None:
    op.add_column(
        "neoantigen_predictions",
        sa.Column("presentation_score", sa.Float(), nullable=True,
                  comment="MHCflurry 2.0 combined presentation probability (0–1). "
                          "Primary ranking metric; accounts for both MHC binding affinity "
                          "and antigen-processing (proteasomal cleavage). "
                          "NULL for rows predicted without the presentation model."),
    )
    op.add_column(
        "neoantigen_predictions",
        sa.Column("processing_score", sa.Float(), nullable=True,
                  comment="MHCflurry 2.0 antigen-processing / proteasomal cleavage score (0–1). "
                          "Derived from flanking sequence context around the epitope. "
                          "NULL for rows predicted without flanking-sequence data."),
    )
    # Index on presentation_score for efficient 'top-N presented neoantigens' queries
    op.create_index(
        "ix_neoantigen_predictions_presentation_score",
        "neoantigen_predictions",
        ["presentation_score"],
    )


def downgrade() -> None:
    op.drop_index("ix_neoantigen_predictions_presentation_score",
                  table_name="neoantigen_predictions")
    op.drop_column("neoantigen_predictions", "processing_score")
    op.drop_column("neoantigen_predictions", "presentation_score")
