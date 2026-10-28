"""Add multi-tool neoantigen ensemble scores (NetMHCpan-4.2 + BigMHC)

Revision ID: b9c3d7e2f1a0
Revises: a1b2c3d4e5f6
Branch labels: None
Depends on: None

Adds 7 columns to neoantigen_predictions for the v5 multi-tool ensemble:
  - netmhcpan_affinity_nm  : NetMHCpan-4.2 BA IC50 in nM
  - netmhcpan_rank_el      : NetMHCpan-4.2 EL %rank
  - netmhcpan_rank_ba      : NetMHCpan-4.2 BA %rank (rescue signal for cysteine/atypical binders)
  - bigmhc_el_score        : BigMHC EL presentation score (0–1)
  - bigmhc_im_score        : BigMHC IM immunogenicity score (0–1); T-cell reactivity layer
  - consensus_rank         : Median %rank across available tools
  - is_ba_rescue           : Flag for peptides rescued by BA-mode but missed by EL tools

All columns are nullable — existing rows remain valid with NULL values.
Run backfill_multi_tool_neoantigens.py to populate.
"""
from alembic import op
import sqlalchemy as sa

revision     = "b9c3d7e2f1a0"
down_revision = "b2c3d4e5f6a1"
branch_labels = None
depends_on    = None


def upgrade() -> None:
    op.add_column("neoantigen_predictions",
        sa.Column("netmhcpan_affinity_nm", sa.Float(), nullable=True,
                  comment="NetMHCpan-4.2 binding affinity IC50 in nM (BA mode)"))
    op.add_column("neoantigen_predictions",
        sa.Column("netmhcpan_rank_el", sa.Float(), nullable=True,
                  comment="NetMHCpan-4.2 EL %rank — eluted ligand percentile (lower = better presented)"))
    op.add_column("neoantigen_predictions",
        sa.Column("netmhcpan_rank_ba", sa.Float(), nullable=True,
                  comment="NetMHCpan-4.2 BA %rank — binding affinity percentile; orthogonal to EL, rescues cysteine-rich peptides"))
    op.add_column("neoantigen_predictions",
        sa.Column("bigmhc_el_score", sa.Float(), nullable=True,
                  comment="BigMHC EL score (0–1): deep-learning presentation probability"))
    op.add_column("neoantigen_predictions",
        sa.Column("bigmhc_im_score", sa.Float(), nullable=True,
                  comment="BigMHC IM score (0–1): T-cell immunogenicity probability via transfer learning"))
    op.add_column("neoantigen_predictions",
        sa.Column("consensus_rank", sa.Float(), nullable=True,
                  comment="Median %rank across MHCflurry, NetMHCpan-EL, BigMHC-EL; primary ranking metric post-v5"))
    op.add_column("neoantigen_predictions",
        sa.Column("is_ba_rescue", sa.Boolean(), nullable=True, server_default="false",
                  comment="True when NetMHCpan-BA predicts binding but EL tools miss; flags potential MS-bias-missed binders"))

    op.create_index("ix_neoantigen_predictions_consensus_rank",
                    "neoantigen_predictions", ["consensus_rank"])
    op.create_index("ix_neoantigen_predictions_bigmhc_im_score",
                    "neoantigen_predictions", ["bigmhc_im_score"])


def downgrade() -> None:
    op.drop_index("ix_neoantigen_predictions_bigmhc_im_score",
                  table_name="neoantigen_predictions")
    op.drop_index("ix_neoantigen_predictions_consensus_rank",
                  table_name="neoantigen_predictions")
    for col in ["is_ba_rescue", "consensus_rank", "bigmhc_im_score",
                "bigmhc_el_score", "netmhcpan_rank_ba", "netmhcpan_rank_el",
                "netmhcpan_affinity_nm"]:
        op.drop_column("neoantigen_predictions", col)
