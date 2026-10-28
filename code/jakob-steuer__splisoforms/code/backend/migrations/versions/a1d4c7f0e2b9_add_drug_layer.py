"""Add drug–isoform layer (drugs, drug_complexes, drug_pocket_fate)

Revision ID: a1d4c7f0e2b9
Revises: 933fa3182e51
Create Date: 2026-07-08

Drugs overlay existing genes/isoforms (not a new isoform Source). Anchored to the
DB by protein_md5. See backend/scripts/build_drug_metadata.py / ingest_drugs.py.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = 'a1d4c7f0e2b9'
down_revision: Union[str, Sequence[str], None] = '933fa3182e51'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        'drugs',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('pubchem_cid', sa.String(length=20), nullable=True),
        sa.Column('broad_id', sa.String(length=50), nullable=True),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.Column('slug', sa.String(length=255), nullable=False),
        sa.Column('chembl_id', sa.String(length=20), nullable=True),
        sa.Column('moa', sa.Text(), nullable=True),
        sa.Column('smiles', sa.Text(), nullable=True),
        sa.Column('inchikey', sa.String(length=27), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
    )
    op.create_index(op.f('ix_drugs_id'), 'drugs', ['id'], unique=False)
    op.create_index(op.f('ix_drugs_pubchem_cid'), 'drugs', ['pubchem_cid'], unique=True)
    op.create_index(op.f('ix_drugs_slug'), 'drugs', ['slug'], unique=False)

    op.create_table(
        'drug_complexes',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('drug_id', sa.Integer(), nullable=False),
        sa.Column('gene_id', sa.Integer(), nullable=True),
        sa.Column('target_isoform_id', sa.Integer(), nullable=True),
        sa.Column('target_uniprot_id', sa.String(length=100), nullable=True),
        sa.Column('target_enst', sa.String(length=30), nullable=True),
        sa.Column('seq_md5', sa.String(length=32), nullable=True),
        sa.Column('job_name', sa.String(length=255), nullable=False),
        sa.Column('structure_dir', sa.String(length=500), nullable=True),
        sa.Column('target_length', sa.Integer(), nullable=True),
        sa.Column('has_structure', sa.Boolean(), nullable=False),
        sa.Column('iptm', sa.Float(), nullable=True),
        sa.Column('ptm', sa.Float(), nullable=True),
        sa.Column('interface_pae_min', sa.Float(), nullable=True),
        sa.Column('ligand_plddt', sa.Float(), nullable=True),
        sa.Column('pocket_plddt', sa.Float(), nullable=True),
        sa.Column('holo_grade', sa.String(length=1), nullable=True),
        sa.Column('n_contact_residues', sa.Integer(), nullable=True),
        sa.Column('contact_residues_json', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['drug_id'], ['drugs.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['gene_id'], ['genes.id']),
        sa.ForeignKeyConstraint(['target_isoform_id'], ['isoforms.id']),
        sa.PrimaryKeyConstraint('id'),
    )
    op.create_index(op.f('ix_drug_complexes_id'), 'drug_complexes', ['id'], unique=False)
    op.create_index(op.f('ix_drug_complexes_drug_id'), 'drug_complexes', ['drug_id'], unique=False)
    op.create_index(op.f('ix_drug_complexes_gene_id'), 'drug_complexes', ['gene_id'], unique=False)
    op.create_index(op.f('ix_drug_complexes_target_isoform_id'), 'drug_complexes', ['target_isoform_id'], unique=False)
    op.create_index(op.f('ix_drug_complexes_target_uniprot_id'), 'drug_complexes', ['target_uniprot_id'], unique=False)
    op.create_index(op.f('ix_drug_complexes_seq_md5'), 'drug_complexes', ['seq_md5'], unique=False)
    op.create_index(op.f('ix_drug_complexes_job_name'), 'drug_complexes', ['job_name'], unique=True)
    op.create_index(op.f('ix_drug_complexes_has_structure'), 'drug_complexes', ['has_structure'], unique=False)
    op.create_index(op.f('ix_drug_complexes_holo_grade'), 'drug_complexes', ['holo_grade'], unique=False)

    op.create_table(
        'drug_pocket_fate',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('drug_complex_id', sa.Integer(), nullable=False),
        sa.Column('isoform_id', sa.Integer(), nullable=False),
        sa.Column('n_pocket_residues', sa.Integer(), nullable=True),
        sa.Column('n_retained', sa.Integer(), nullable=True),
        sa.Column('frac_retained', sa.Float(), nullable=True),
        sa.Column('fate', sa.String(length=20), nullable=True),
        sa.Column('mechanism_class', sa.String(length=30), nullable=True),
        sa.Column('decidable', sa.Boolean(), nullable=False),
        sa.Column('pocket_plddt', sa.Float(), nullable=True),
        sa.Column('retained_residues_json', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['drug_complex_id'], ['drug_complexes.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['isoform_id'], ['isoforms.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
    )
    op.create_index(op.f('ix_drug_pocket_fate_id'), 'drug_pocket_fate', ['id'], unique=False)
    op.create_index(op.f('ix_drug_pocket_fate_drug_complex_id'), 'drug_pocket_fate', ['drug_complex_id'], unique=False)
    op.create_index(op.f('ix_drug_pocket_fate_isoform_id'), 'drug_pocket_fate', ['isoform_id'], unique=False)
    op.create_index(op.f('ix_drug_pocket_fate_fate'), 'drug_pocket_fate', ['fate'], unique=False)


def downgrade() -> None:
    op.drop_table('drug_pocket_fate')
    op.drop_table('drug_complexes')
    op.drop_table('drugs')
