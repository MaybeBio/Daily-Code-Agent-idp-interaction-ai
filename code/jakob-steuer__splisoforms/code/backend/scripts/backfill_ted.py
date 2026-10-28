import os
import sys
import csv
import logging
from sqlalchemy.orm import sessionmaker

# Ensure backend imports work
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app.database import engine
from app.models.bio import Base, Isoform, TedDomain, TranscriptMetadata, Gene

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

def backfill_ted_domains(tsv_path: str, limit: int = None):
    """
    Parses a Biorxiv TED / TEDLH structural domain mapping tab-separated file.
    Maps these pre-computed consensus boundaries to Ensembl canonical isoforms using UniProt IDs.
    Ingests the chopping boundaries, average pLDDT, and CATH labels.
    """
    if not os.path.exists(tsv_path):
        # Hard-fail (not a silent return) so the backfill chain records FAIL rather
        # than a misleading PASS when the TED source file is missing/misplaced.
        sys.exit(f"[!] Critical Error: TED source file not found at {tsv_path}")

    print(f"[*] Extracting and matching TED Domains from {tsv_path}...")
    db = SessionLocal()
    matched_count = 0
    missed_count = 0
    length_mismatch_count = 0
    failed = False

    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for index, row in enumerate(reader):
                if limit and matched_count >= limit:
                    break
                    
                if not row or len(row) < 14:
                    continue
                    
                domain_id = row[0]
                if not domain_id.startswith("AF-"):
                    continue
                
                parts = domain_id.split('-')
                if len(parts) < 2:
                    continue
                
                uniprot_id = parts[1]
                chopping = row[3]
                try:
                    plddt = float(row[6]) if row[6] and row[6] != '-' else None
                except ValueError:
                    plddt = None
                
                cath_label = row[13] if row[13] and row[13] != '-' else None
                domain_type = row[2].upper() if row[2] else "TED"
                
                isoforms = []
                gene = db.query(Gene).filter(Gene.uniprot_id == uniprot_id).first()
                if gene:
                    canonical_iso_id = gene.canonical_transcript_id
                    if canonical_iso_id:
                        iso = db.query(Isoform).filter(Isoform.id == canonical_iso_id).first()
                        if iso:
                            isoforms.append(iso)
                    if not isoforms:
                        isos = db.query(Isoform).filter(Isoform.gene_id == gene.id, Isoform.is_canonical == True).all()
                        isoforms.extend(isos)

                if not isoforms:
                    meta = db.query(TranscriptMetadata).filter(TranscriptMetadata.uniprot_id == uniprot_id).all()
                    for m in meta:
                        iso = db.query(Isoform).filter(Isoform.id == m.isoform_id).first()
                        if iso and iso not in isoforms:
                            # Prefer Canonical isoform if possible
                            if iso.is_canonical:
                                isoforms.insert(0, iso)
                            else:
                                isoforms.append(iso)
                            break

                if not isoforms:
                    missed_count += 1
                    continue
                
                for isoform in isoforms:
                    existing = db.query(TedDomain).filter(
                        TedDomain.isoform_id == isoform.id,
                        TedDomain.ted_id == domain_id
                    ).first()

                    if existing:
                        continue

                    seq_start = 0
                    seq_end = 0
                    if chopping and chopping != '-':
                        segments = chopping.split('_')
                        if segments:
                            first_seg = segments[0]
                            last_seg = segments[-1]
                            try:
                                seq_start = int(first_seg.split('-')[0])
                                seq_end = int(last_seg.split('-')[1])
                            except Exception:
                                pass

                    # The TED domain is defined on the AFDB model of the UniProt
                    # canonical sequence; if our designated canonical isoform is a
                    # different (shorter) protein the boundaries do not transfer.
                    if isoform.sequence and seq_end > len(isoform.sequence):
                        length_mismatch_count += 1
                        continue

                    ted = TedDomain(
                        isoform_id=isoform.id,
                        ted_id=domain_id,
                        domain_type=domain_type,
                        seq_start=seq_start,
                        seq_end=seq_end,
                        score=0.0,
                        chopping=chopping,
                        plddt=plddt,
                        cath_label=cath_label
                    )
                    db.add(ted)
                
                matched_count += 1

                if matched_count % 500 == 0:
                    db.commit()
                    print(f"  ... Cached {matched_count} TED domains")

        db.commit()
        print(f"\n✅ Success! Loaded {matched_count} TED domains. (Missed {missed_count} uniprot "
              f"alignments; {length_mismatch_count} domains skipped because they extend past "
              f"the designated canonical isoform).")

    except Exception as e:
        print(f"[!] Extraction Failure: {e}")
        db.rollback()
        failed = True
    finally:
        db.close()
    # A swallowed failure used to exit 0 and read as PASS in the backfill chain
    # (review ENG-2 / DOM-5). Fail loudly when nothing was loaded.
    if failed or matched_count == 0:
        sys.exit(1)

def _default_ted_path():
    # The 29 MB TED source lives in db/ted_domains/ (a sibling of db/splisoforms/),
    # NOT under the repo's data/. Probe both so a future move fails loudly.
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(here))   # db/splisoforms
    db_dir = os.path.dirname(repo)                   # db/
    cands = [os.path.join(db_dir, "ted_domains", "human_ted_domains.tsv"),
             os.path.join(repo, "data", "ted_domains", "human_ted_domains.tsv")]
    return next((p for p in cands if os.path.exists(p)), cands[0])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Populate TED Structural Domains")
    parser.add_argument("--tsv-path", type=str, default=_default_ted_path(), help="Path to structural mapping TSV")
    parser.add_argument("--limit", type=int, help="Limit number of domains to process")
    args = parser.parse_args()

    print("=== SPLISOFORMS TED Domain Extraction Engine ===")
    backfill_ted_domains(args.tsv_path, args.limit)
