"""
SPLISOFORMS public REST API — v1.

A small, curated, read-only surface for programmatic access to genes, isoforms,
summary metrics, and structures. Every endpoint requires an API key
(`Authorization: Bearer <key>`) and is rate-limited per key. Data is read
through a SELECT-only Postgres role.

Deliberately NOT exposed: bulk/unscoped dumps, raw per-residue mass data,
write/compute endpoints, internal workflow tables (prediction_jobs), and the
proprietary per-tool neoantigen ensemble internals.
"""

from __future__ import annotations

from typing import List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import RedirectResponse
from pydantic import BaseModel
from sqlalchemy import func, or_
from sqlalchemy.orm import Session

from app.api.security import optional_api_key
from app.database import get_public_db
from app.models import bio as m

# Auth + rate limit applied to every endpoint in this router.
router = APIRouter(dependencies=[Depends(optional_api_key)])

MAX_PAGE = 100
MAX_SEARCH = 25
MAX_NEOANTIGENS = 50


# ── Response models (field-whitelisted; no internal/IP fields) ─────────────────

class IsoformBrief(BaseModel):
    transcript_id: str
    ensembl_transcript_id: Optional[str] = None
    ensembl_protein_id: Optional[str] = None
    length: Optional[int] = None
    delta_length: Optional[int] = None
    is_canonical: Optional[bool] = None
    is_mane_select: Optional[bool] = None
    has_structure: bool = False


class GeneSummary(BaseModel):
    symbol: str
    ensembl_gene_id: Optional[str] = None
    chromosome: Optional[str] = None
    strand: Optional[str] = None
    n_isoforms: int
    isoforms: List[IsoformBrief]


class IsoformMetrics(BaseModel):
    nmd_status: Optional[str] = None
    pfam_domains_present: int = 0
    pfam_domains_lost: int = 0
    ted_domains: int = 0
    neoantigens: int = 0
    ptm_sites: int = 0
    mean_plddt: Optional[float] = None


class IsoformDetail(BaseModel):
    transcript_id: str
    gene_symbol: Optional[str] = None
    ensembl_transcript_id: Optional[str] = None
    ensembl_protein_id: Optional[str] = None
    length: Optional[int] = None
    delta_length: Optional[int] = None
    is_canonical: Optional[bool] = None
    is_mane_select: Optional[bool] = None
    structure_provenance: Optional[str] = None
    gencode_version: Optional[str] = None
    has_structure: bool = False
    sequence: Optional[str] = None
    metrics: IsoformMetrics


class DomainOut(BaseModel):
    source: str            # "pfam" | "ted"
    name: Optional[str] = None
    accession: Optional[str] = None
    seq_start: int
    seq_end: int
    status: Optional[str] = None
    alteration_type: Optional[str] = None


class NeoantigenOut(BaseModel):
    peptide: str
    hla_allele: str
    binding_affinity_nm: float
    binding_category: str
    presentation_score: Optional[float] = None
    percentile_rank: Optional[float] = None
    agreement_tier: Optional[str] = None


class SearchGene(BaseModel):
    symbol: str
    ensembl_gene_id: Optional[str] = None


class SearchIsoform(BaseModel):
    transcript_id: str
    gene_symbol: Optional[str] = None


class SearchResult(BaseModel):
    genes: List[SearchGene]
    isoforms: List[SearchIsoform]


# ── Helpers ────────────────────────────────────────────────────────────────────

def _brief(db: Session, iso: m.Isoform) -> IsoformBrief:
    has_struct = db.query(m.Structure.id).filter(m.Structure.isoform_id == iso.id).first() is not None
    return IsoformBrief(
        transcript_id=iso.transcript_id,
        ensembl_transcript_id=iso.ensembl_transcript_id,
        ensembl_protein_id=iso.ensembl_protein_id,
        length=len(iso.sequence) if iso.sequence else None,
        delta_length=iso.delta_length,
        is_canonical=iso.is_canonical,
        is_mane_select=iso.is_mane_select,
        has_structure=has_struct,
    )


def _get_gene_or_404(db: Session, symbol: str) -> m.Gene:
    gene = db.query(m.Gene).filter(m.Gene.symbol.ilike(symbol)).first()
    if not gene:
        raise HTTPException(status_code=404, detail=f"Gene '{symbol}' not found.")
    return gene


def _get_isoform_or_404(db: Session, transcript_id: str) -> m.Isoform:
    # Accept either the stable transcript_id (ENST…/PB…) or the display transcript_name
    # (e.g. F2-201, TP53-201-NOV-…) — the latter is what users actually know.
    iso = (
        db.query(m.Isoform)
        .filter(or_(m.Isoform.transcript_id == transcript_id,
                    m.Isoform.transcript_name == transcript_id))
        .first()
    )
    if not iso:
        raise HTTPException(status_code=404, detail=f"Isoform '{transcript_id}' not found.")
    return iso


# ── Endpoints ──────────────────────────────────────────────────────────────────

@router.get("/genes/{symbol}", response_model=GeneSummary, summary="Gene summary + its isoforms")
def get_gene(symbol: str, db: Session = Depends(get_public_db)):
    gene = _get_gene_or_404(db, symbol)
    isos = (
        db.query(m.Isoform)
        .filter(m.Isoform.gene_id == gene.id)
        .order_by(m.Isoform.is_canonical.desc(), m.Isoform.transcript_id)
        .limit(MAX_PAGE)
        .all()
    )
    return GeneSummary(
        symbol=gene.symbol,
        ensembl_gene_id=gene.ensembl_gene_id,
        chromosome=gene.chromosome,
        strand=gene.strand,
        n_isoforms=db.query(func.count(m.Isoform.id)).filter(m.Isoform.gene_id == gene.id).scalar() or 0,
        isoforms=[_brief(db, i) for i in isos],
    )


@router.get("/genes/{symbol}/isoforms", response_model=List[IsoformBrief], summary="Paginated isoforms of a gene")
def list_gene_isoforms(
    symbol: str,
    limit: int = Query(50, ge=1, le=MAX_PAGE),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_public_db),
):
    gene = _get_gene_or_404(db, symbol)
    isos = (
        db.query(m.Isoform)
        .filter(m.Isoform.gene_id == gene.id)
        .order_by(m.Isoform.is_canonical.desc(), m.Isoform.transcript_id)
        .offset(offset).limit(limit).all()
    )
    return [_brief(db, i) for i in isos]


@router.get("/isoforms/{transcript_id}", response_model=IsoformDetail, summary="Isoform core record + summary metrics")
def get_isoform(transcript_id: str, db: Session = Depends(get_public_db)):
    iso = _get_isoform_or_404(db, transcript_id)
    iid = iso.id

    def count(model, *extra):
        q = db.query(func.count(model.id)).filter(model.isoform_id == iid)
        for e in extra:
            q = q.filter(e)
        return q.scalar() or 0

    metrics = IsoformMetrics(
        nmd_status=db.query(m.NmdPrediction.nmd_status).filter(m.NmdPrediction.isoform_id == iid).scalar(),
        pfam_domains_present=count(m.PfamDomain, m.PfamDomain.domain_status == "present"),
        pfam_domains_lost=count(m.PfamDomain, m.PfamDomain.domain_status == "lost"),
        ted_domains=count(m.TedDomain),
        neoantigens=count(m.NeoantigenPrediction),
        ptm_sites=count(m.PtmAccessibility),
        mean_plddt=db.query(m.StructuralMetrics.plddt_mean).filter(m.StructuralMetrics.isoform_id == iid).scalar(),
    )
    has_struct = db.query(m.Structure.id).filter(m.Structure.isoform_id == iid).first() is not None
    return IsoformDetail(
        transcript_id=iso.transcript_id,
        gene_symbol=iso.gene.symbol if iso.gene else None,
        ensembl_transcript_id=iso.ensembl_transcript_id,
        ensembl_protein_id=iso.ensembl_protein_id,
        length=len(iso.sequence) if iso.sequence else None,
        delta_length=iso.delta_length,
        is_canonical=iso.is_canonical,
        is_mane_select=iso.is_mane_select,
        structure_provenance=iso.structure_provenance,
        gencode_version=iso.gencode_version,
        has_structure=has_struct,
        sequence=iso.sequence,
        metrics=metrics,
    )


@router.get("/isoforms/{transcript_id}/structure", summary="Redirect to the isoform's predicted structure (mmCIF)")
def get_isoform_structure(transcript_id: str, request: Request, db: Session = Depends(get_public_db)):
    iso = _get_isoform_or_404(db, transcript_id)
    struct = db.query(m.Structure).filter(m.Structure.isoform_id == iso.id).first()
    if not struct or not struct.af3_cif_url:
        raise HTTPException(status_code=404, detail="No structure available for this isoform.")
    url = struct.af3_cif_url
    if not url.startswith("http"):
        url = str(request.base_url).rstrip("/") + url
    return RedirectResponse(url=url, status_code=307)


@router.get("/isoforms/{transcript_id}/domains", response_model=List[DomainOut], summary="Pfam + TED domains")
def get_isoform_domains(transcript_id: str, db: Session = Depends(get_public_db)):
    iso = _get_isoform_or_404(db, transcript_id)
    out: List[DomainOut] = []
    for d in db.query(m.PfamDomain).filter(m.PfamDomain.isoform_id == iso.id).order_by(m.PfamDomain.seq_start).limit(MAX_PAGE).all():
        out.append(DomainOut(source="pfam", name=d.pfam_name, accession=d.pfam_accession,
                             seq_start=d.seq_start, seq_end=d.seq_end,
                             status=d.domain_status, alteration_type=d.alteration_type))
    for t in db.query(m.TedDomain).filter(m.TedDomain.isoform_id == iso.id).order_by(m.TedDomain.seq_start).limit(MAX_PAGE).all():
        out.append(DomainOut(source="ted", name=t.cath_label, accession=t.ted_id,
                             seq_start=t.seq_start, seq_end=t.seq_end))
    return out


@router.get("/isoforms/{transcript_id}/neoantigens", response_model=List[NeoantigenOut], summary="Top predicted neoantigens (capped)")
def get_isoform_neoantigens(
    transcript_id: str,
    limit: int = Query(25, ge=1, le=MAX_NEOANTIGENS),
    db: Session = Depends(get_public_db),
):
    iso = _get_isoform_or_404(db, transcript_id)
    rows = (
        db.query(m.NeoantigenPrediction)
        .filter(m.NeoantigenPrediction.isoform_id == iso.id)
        .order_by(m.NeoantigenPrediction.presentation_score.desc().nullslast())
        .limit(limit).all()
    )
    return [
        NeoantigenOut(
            peptide=r.peptide_sequence,
            hla_allele=r.hla_allele,
            binding_affinity_nm=r.binding_affinity_nm,
            binding_category=r.binding_category,
            presentation_score=r.presentation_score,
            percentile_rank=r.percentile_rank,
            agreement_tier=r.agreement_tier,
        )
        for r in rows
    ]


@router.get("/search", response_model=SearchResult, summary="Search genes + isoforms")
def search(
    query: str = Query(..., min_length=2, max_length=64),
    limit: int = Query(10, ge=1, le=MAX_SEARCH),
    db: Session = Depends(get_public_db),
):
    term = f"%{query}%"
    genes = (
        db.query(m.Gene)
        .filter(or_(m.Gene.symbol.ilike(term), m.Gene.hgnc_symbol.ilike(term), m.Gene.uniprot_id.ilike(term)))
        .limit(limit).all()
    )
    isos = (
        db.query(m.Isoform)
        .filter(or_(m.Isoform.transcript_id.ilike(term), m.Isoform.ensembl_transcript_id.ilike(term)))
        .limit(limit).all()
    )
    return SearchResult(
        genes=[SearchGene(symbol=g.symbol, ensembl_gene_id=g.ensembl_gene_id) for g in genes],
        isoforms=[SearchIsoform(transcript_id=i.transcript_id, gene_symbol=i.gene.symbol if i.gene else None) for i in isos],
    )
