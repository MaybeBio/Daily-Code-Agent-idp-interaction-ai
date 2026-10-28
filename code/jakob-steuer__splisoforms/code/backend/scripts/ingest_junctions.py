import os
import sys

# Ensure backend imports work regardless of cwd (must precede `app.*` imports).
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
from sqlalchemy.orm import Session
from sqlalchemy import or_
from app.database import SessionLocal
from app.models import bio as bio_models

import argparse

def ingest_junctions(data_dir: str, limit: int = None):
    db = SessionLocal()
    try:
        sqanti_dir = os.path.join(data_dir, "sqanti3")
        if not os.path.isdir(sqanti_dir):
            print(f"Error: sqanti3 directory not found at {sqanti_dir}")
            return

        # Pre-fetch ALL isoforms for fast PB.xxx mapping
        print("Caching isoforms for fast mapping...")
        isoforms = db.query(bio_models.Isoform).filter(
            bio_models.Isoform.transcript_id.like("PB.%")
        ).all()

        # Map PB IDs -> Isoform object (handle versioning: PB.581.114 -> PB.581)
        iso_map = {}
        for iso in isoforms:
            # Store both full ID and stem (PB.xxx)
            iso_map[iso.transcript_id] = iso
            stem = ".".join(iso.transcript_id.split(".")[:2])
            iso_map[stem] = iso
        print(f"Cached {len(isoforms)} PacBio isoforms ({len(iso_map)} lookup keys).")

        # Discover sample groups from SQANTI3 subdirectories
        groups = sorted([
            d for d in os.listdir(sqanti_dir)
            if os.path.isdir(os.path.join(sqanti_dir, d))
        ])
        print(f"Found {len(groups)} sample groups: {groups}")
        processed_count = 0

        for group in groups:
            print(f"Processing group: {group}")
            group_dir = os.path.join(sqanti_dir, group)
            class_file = os.path.join(group_dir, "scisoseq_transcripts.sorted.filtered_lite_classification.txt")
            
            if not os.path.exists(class_file):
                print(f"  Warning: Classification file not found for {group}: {class_file}")
                continue

            # Load classification
            df_class = pd.read_csv(class_file, sep="\t")
            
            for idx, row in df_class.iterrows():
                if limit and processed_count >= limit:
                    break
                    
                pb_id = str(row['isoform'])
                isoform = iso_map.get(pb_id)
                if not isoform:
                    continue

                # Clear existing junctions for this isoform to prevent duplicates
                db.query(bio_models.IsoformJunction).filter(bio_models.IsoformJunction.isoform_id == isoform.id).delete()

                # Derived Logic
                try:
                    exons = int(row['exons']) if not pd.isna(row['exons']) else 0
                    ref_exons = int(row['ref_exons']) if not pd.isna(row['ref_exons']) else 0
                except (ValueError, TypeError):
                    exons = 0
                    ref_exons = 0
                
                struct_cat = str(row['structural_category']) if not pd.isna(row['structural_category']) else "unknown"
                sub_cat = str(row['subcategory']) if not pd.isna(row['structural_category']) else "none"
                
                junction_type = str(row.get('structural_category', 'unknown'))
                is_novel = "novel" in struct_cat.lower() or "atypical" in sub_cat.lower()

                # Create Junction Event
                junc_event = bio_models.IsoformJunction(
                    isoform_id=isoform.id,
                    chrom=str(row['chrom']),
                    strand=str(row['strand']),
                    junction_type=junction_type,
                    subcategory=sub_cat,
                    exon_count_novel=exons,
                    exon_count_canonical=ref_exons,
                    is_novel=is_novel
                )
                db.add(junc_event)
                processed_count += 1

            db.commit()
            print(f"  Processed junctions for {group}")
            if limit and processed_count >= limit:
                break

        print(f"✅ Ingestion Complete. Total events created: {processed_count}")

    except Exception as e:
        db.rollback()
        print(f"❌ Error during ingestion: {e}")
    finally:
        db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ingest SQANTI3 Splice Junctions")
    parser.add_argument("--data-dir", type=str, default="data/sources/ccrcc", help="Base directory containing sqanti3/<sample>/ results")
    parser.add_argument("--limit", type=int, help="Limit number of isoforms to process for testing")
    args = parser.parse_args()

    ingest_junctions(args.data_dir, args.limit)
