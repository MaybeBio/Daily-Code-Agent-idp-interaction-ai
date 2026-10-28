"""
Cell-type specificity lookup via the Human Protein Atlas (HPA).

Fetches per-gene single cell type expression data from the HPA JSON API.
The data is gene-level (not isoform-level), so expression profiles are
linked to the gene and inherited by all isoforms of that gene.

HPA endpoint: https://www.proteinatlas.org/{ENSEMBL_GENE_ID}.json
Key fields parsed:
  - "RNA single cell type specific nCPM": dict of cell_type → expression value
  - "RNA single cell type specificity": classification string
"""
import logging
import re
from typing import Dict, List, Optional, Tuple

import httpx
from sqlalchemy.orm import Session

from app.models.bio import Gene, Isoform, TranscriptMetadata, CellTypeExpression

logger = logging.getLogger(__name__)

HPA_BASE_URL = "https://www.proteinatlas.org"
HPA_TIMEOUT = 10.0  # seconds


def resolve_ensembl_gene_id(db: Session, gene: Gene) -> Optional[str]:
    """
    Resolve an Ensembl gene ID (ENSG...) for a Gene record.

    Strategies:
    1. Gene.ensembl_gene_id if already set
    2. Direct UniProt ID lookup for Ensembl cross-reference
    3. UniProt gene symbol search for the reviewed SwissProt entry
       (handles TrEMBL IDs that lack Ensembl cross-refs)
    """
    if gene.ensembl_gene_id:
        return gene.ensembl_gene_id

    ensg = None

    # Try direct UniProt lookup first (fastest, no search needed)
    if gene.uniprot_id:
        ensg = _lookup_ensembl_via_uniprot(gene.uniprot_id)

    # Fallback: search by gene symbol (works for TrEMBL IDs and genes without UniProt)
    if not ensg and gene.symbol and not gene.symbol.startswith("PB."):
        ensg = _search_ensembl_by_gene_symbol(gene.symbol)

    if ensg:
        gene.ensembl_gene_id = ensg
        db.commit()
        return ensg

    return None


def _lookup_ensembl_via_uniprot(uniprot_id: str) -> Optional[str]:
    """
    Look up Ensembl gene ID from UniProt ID via the UniProt REST API.
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        with httpx.Client(timeout=HPA_TIMEOUT) as client:
            resp = client.get(url)
            if resp.status_code != 200:
                return None

            data = resp.json()
            # Look through cross-references for Ensembl
            for xref in data.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "Ensembl":
                    # Properties contain the gene ID
                    for prop in xref.get("properties", []):
                        if prop.get("key") == "GeneId":
                            gene_id = prop.get("value", "")
                            if gene_id.startswith("ENSG"):
                                # Strip version suffix (ENSG00000074800.17 → ENSG00000074800)
                                return gene_id.split(".")[0]
    except Exception as e:
        logger.debug(f"UniProt lookup failed for {uniprot_id}: {e}")

    return None


def _search_ensembl_by_gene_symbol(gene_symbol: str) -> Optional[str]:
    """
    Search UniProt for the reviewed human entry of a gene symbol
    and extract the Ensembl gene ID from its cross-references.

    Handles the common case where Gene.uniprot_id is a TrEMBL accession
    that lacks Ensembl cross-references, but the reviewed SwissProt entry
    for the same gene does have them.
    """
    try:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
            "format": "json",
            "fields": "accession,xref_ensembl",
            "size": 3,
        }
        with httpx.Client(timeout=HPA_TIMEOUT) as client:
            resp = client.get(url, params=params)
            if resp.status_code != 200:
                return None

            data = resp.json()
            for entry in data.get("results", []):
                for xref in entry.get("uniProtKBCrossReferences", []):
                    if xref.get("database") == "Ensembl":
                        for prop in xref.get("properties", []):
                            if prop.get("key") == "GeneId":
                                gene_id = prop.get("value", "")
                                if gene_id.startswith("ENSG"):
                                    return gene_id.split(".")[0]
    except Exception as e:
        logger.debug(f"Ensembl search by gene symbol failed for {gene_symbol}: {e}")

    return None


def fetch_hpa_data(ensembl_gene_id: str) -> Tuple[Optional[Dict[str, float]], Optional[str]]:
    """
    Fetch single cell type expression data from the HPA JSON API.

    Args:
        ensembl_gene_id: Ensembl gene ID (e.g., "ENSG00000169083")

    Returns:
        Tuple of (cell_type_expression_dict, specificity_category).
        Expression dict maps cell type names to nCPM values.
        Returns (None, None) if data unavailable.
    """
    url = f"{HPA_BASE_URL}/{ensembl_gene_id}.json"

    try:
        with httpx.Client(timeout=HPA_TIMEOUT) as client:
            resp = client.get(url)
            if resp.status_code != 200:
                logger.warning(f"HPA returned {resp.status_code} for {ensembl_gene_id}")
                return None, None

            data = resp.json()

    except Exception as e:
        logger.error(f"HPA fetch failed for {ensembl_gene_id}: {e}")
        return None, None

    # Parse single cell type expression
    expression_raw = data.get("RNA single cell type specific nCPM")
    specificity = data.get("RNA single cell type specificity")

    if not expression_raw or not isinstance(expression_raw, dict):
        logger.debug(f"No single cell type data for {ensembl_gene_id}")
        return None, specificity

    # Convert string values to floats
    expression: Dict[str, float] = {}
    for cell_type, value in expression_raw.items():
        try:
            expression[cell_type] = float(value)
        except (ValueError, TypeError):
            continue

    # Normalize specificity category
    if specificity:
        specificity = specificity.strip().lower().replace(" ", "_")

    return expression, specificity


def analyze_gene(
    db: Session,
    gene: Gene,
    force: bool = False,
) -> List[CellTypeExpression]:
    """
    Fetch and store cell-type expression data for a gene from HPA.

    Args:
        db: SQLAlchemy session.
        gene: The gene to look up.
        force: If True, delete existing results and re-fetch.

    Returns:
        List of CellTypeExpression ORM objects.
    """
    # Check existing
    existing = db.query(CellTypeExpression).filter(
        CellTypeExpression.gene_id == gene.id
    ).first()
    if existing and not force:
        logger.info(f"Cell-type expression already exists for {gene.symbol}, skipping")
        return db.query(CellTypeExpression).filter(
            CellTypeExpression.gene_id == gene.id
        ).all()

    if force:
        db.query(CellTypeExpression).filter(
            CellTypeExpression.gene_id == gene.id
        ).delete()
        db.flush()

    # Resolve Ensembl ID
    ensg = resolve_ensembl_gene_id(db, gene)
    if not ensg:
        logger.warning(f"Could not resolve Ensembl gene ID for {gene.symbol}")
        return []

    # Fetch from HPA
    expression, specificity = fetch_hpa_data(ensg)
    if not expression:
        logger.info(f"No HPA single cell data for {gene.symbol} ({ensg})")
        return []

    # Store results
    results: List[CellTypeExpression] = []
    for cell_type, ncpm in expression.items():
        record = CellTypeExpression(
            gene_id=gene.id,
            cell_type=cell_type,
            expression_ncpm=ncpm,
            source="HPA_v25",
            specificity_category=specificity,
        )
        results.append(record)

    db.bulk_save_objects(results)
    db.commit()

    logger.info(
        f"Cell-type expression for {gene.symbol}: "
        f"{len(results)} cell types, specificity={specificity}"
    )

    return results


def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
) -> List[CellTypeExpression]:
    """
    Convenience wrapper that runs cell-type analysis for an isoform's gene.

    This is gene-level data, so all isoforms of the same gene share
    the same cell-type expression profile.
    """
    if not isoform.gene_id:
        return []

    gene = db.query(Gene).filter(Gene.id == isoform.gene_id).first()
    if not gene:
        return []

    return analyze_gene(db, gene, force=force)
