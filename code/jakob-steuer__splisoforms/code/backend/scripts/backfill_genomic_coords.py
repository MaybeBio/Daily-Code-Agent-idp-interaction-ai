"""
Backfill transcript_metadata.chrom and strand from two sources:
  1. SQANTI3 classification files (per-isoform chrom/strand)
  2. Gene table fallback (gene.chromosome / gene.strand)

Also updates transcript_metadata.length if missing.
"""
import os, sys, csv, glob, logging
from sqlalchemy.orm import sessionmaker

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app.database import engine
from app.models.bio import TranscriptMetadata, Isoform, Gene

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
Session = sessionmaker(bind=engine)

SQANTI_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "ccrcc_metadata", "SQANTI3"
)

def load_sqanti_coords():
    """Return {transcript_id: (chrom, strand, length)} from all SQANTI classification files."""
    coords = {}
    pattern = os.path.join(SQANTI_DIR, "**", "*_classification.txt")
    files = glob.glob(pattern, recursive=True)
    log.info(f"Found {len(files)} SQANTI classification files")
    for f in files:
        with open(f, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                tid = row.get("isoform", "").strip()
                chrom = row.get("chrom", "").strip()
                strand = row.get("strand", "").strip()
                try:
                    length = int(row.get("length", 0))
                except ValueError:
                    length = None
                if tid and chrom and chrom not in ("", "NA", "unknown"):
                    coords[tid] = (chrom, strand, length)
    log.info(f"Loaded {len(coords)} transcript coords from SQANTI")
    return coords

def run():
    db = Session()
    sqanti = load_sqanti_coords()

    # Load gene chrom/strand fallback
    gene_map = {}  # gene_id → (chromosome, strand)
    for g in db.query(Gene).all():
        if g.chromosome:
            gene_map[g.id] = (g.chromosome, g.strand)

    rows = db.query(TranscriptMetadata).filter(
        (TranscriptMetadata.chrom == None) |
        (TranscriptMetadata.chrom == "unknown") |
        (TranscriptMetadata.chrom == "")
    ).all()
    log.info(f"Rows to update: {len(rows)}")

    updated_sqanti = 0
    updated_gene = 0
    skipped = 0

    for meta in rows:
        # Get the isoform's gene_id for fallback
        iso = db.query(Isoform).filter(Isoform.id == meta.isoform_id).first()
        gene_id = iso.gene_id if iso else None

        if meta.transcript_id in sqanti:
            chrom, strand, length = sqanti[meta.transcript_id]
            meta.chrom = chrom
            meta.strand = strand
            if length and not meta.length:
                meta.length = length
            updated_sqanti += 1
        elif gene_id and gene_id in gene_map:
            chrom, strand = gene_map[gene_id]
            meta.chrom = f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom)
            meta.strand = strand
            updated_gene += 1
        else:
            skipped += 1

        if (updated_sqanti + updated_gene) % 1000 == 0 and (updated_sqanti + updated_gene) > 0:
            db.commit()
            log.info(f"  Progress: {updated_sqanti} from SQANTI, {updated_gene} from gene, {skipped} skipped")

    db.commit()
    db.close()
    log.info(f"Done. SQANTI: {updated_sqanti}, gene fallback: {updated_gene}, skipped: {skipped}")

if __name__ == "__main__":
    run()
