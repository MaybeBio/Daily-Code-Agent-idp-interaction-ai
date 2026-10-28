import os
import sys
import glob
import subprocess
import csv
import logging
from sqlalchemy.orm import sessionmaker

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app.database import engine
from app.models.bio import Base, Isoform, TedDomain, Structure
from app.analysis.ptm_accessibility import analyze_isoform

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

def generate_novel_ted_domains(af3_results_dir: str, chainsaw_script_dir: str):
    print("[*] Monitoring novel AF3 structures...")
    
    if not os.path.isdir(af3_results_dir):
        print(f"[!] Directory not found: {af3_results_dir}")
        return

    db = SessionLocal()
    
    # 1. Prepare flat directory for Chainsaw
    flat_dir = os.path.join(af3_results_dir, "flat_cifs")
    os.makedirs(flat_dir, exist_ok=True)
    
    print("[*] Finding and symlinking .cif files...")
    structures_dir = os.path.join(af3_results_dir, "structures")
    cif_files = glob.glob(os.path.join(structures_dir, "**", "*.cif"), recursive=True)
    cif_count = 0
    for cif in cif_files:
        filename = os.path.basename(cif)
        symlink_path = os.path.join(flat_dir, filename)
        if not os.path.exists(symlink_path):
            os.symlink(cif, symlink_path)
        cif_count += 1
        
    print(f"[*] Symlinked {cif_count} .cif files.")
    if cif_count == 0:
        print("[!] No .cif files found. Exiting.")
        return

    # 2. Run Chainsaw
    results_tsv = os.path.join(af3_results_dir, "chainsaw_results.tsv")
    
    if not os.path.exists(results_tsv):
        print("[*] Running Chainsaw segmentation...")
        cmd = [
            sys.executable,
            os.path.join(chainsaw_script_dir, "get_predictions.py"),
            "--structure_directory", flat_dir,
            "-o", results_tsv
        ]
        
        env = os.environ.copy()
        env["PYTHONPATH"] = chainsaw_script_dir
        
        try:
            subprocess.run(cmd, env=env, check=True, cwd=chainsaw_script_dir)
            print("[*] Chainsaw prediction complete.")
        except subprocess.CalledProcessError as e:
            print(f"[!] Chainsaw prediction failed: {e}")
            return
    else:
        print("[*] Found existing chainsaw_results.tsv, skipping segmentation...")
        
    print("[*] Parsing Chainsaw outputs and storing novel domains...")
    if not os.path.exists(results_tsv):
        print("[!] No novel domains were segmented.")
        return
        
    # Map Chainsaw chain_id ("<md5>_model") → isoform_id in O(1). The af3_cif_url
    # basename is "<md5>_model.cif"; its stem equals the chain_id the segmenter
    # emits, so we key on that directly (the old code did an O(structures) substring
    # scan per row — billions of comparisons at full scale).
    struct_map = {}
    for url, iso_id in (
        db.query(Structure.af3_cif_url, Structure.isoform_id)
          .filter(Structure.af3_cif_url.isnot(None)).all()
    ):
        stem = os.path.basename(url)
        if stem.endswith(".cif"):
            stem = stem[:-4]
        struct_map[stem] = iso_id

    print(f"[*] Fetched {len(struct_map)} existing structures mapping.")
        
    # 3. Parse outputs and map to database
    with open(results_tsv, 'r', encoding='utf-8') as f:
        # columns: chain_id, sequence_md5, nres, ndom, chopping, uncertainty
        reader = csv.DictReader(f, delimiter='\t')
        row_count = 0
        novel_domain_objects = []
        
        for row in reader:
            chain_id = row['chain_id']
            chopping = row['chopping']
            uncertainty = row['uncertainty']
            
            target_isoform_id = struct_map.get(chain_id)
            if not target_isoform_id:
                # No structure matches this chain_id — skip.
                continue
                
            if chopping == 'NULL' or not chopping:
                continue
                
            # Handle possible existing domains
            existing = db.query(TedDomain).filter(
                TedDomain.isoform_id == target_isoform_id,
                TedDomain.domain_type.in_(['TED_NOVEL', 'CHAINSAW'])
            ).all()
            for ex in existing:
                db.delete(ex)
            db.commit()

            domains_choppings = chopping.split(',')
            for idx, d_chop in enumerate(domains_choppings):
                # Chainsaw choppings are 0-based, inclusive residue indices
                # (adjacent domains read "0-85,86-171"; first residue is 0).
                # Everything else in the DB (TED consensus, Pfam, PTM positions)
                # is 1-based, so shift by +1 on both ends (review DOM-6). The
                # stored `chopping` string is rewritten in 1-based form too.
                seq_start = 0
                seq_end = 0
                segments_1based = []
                try:
                    for seg in d_chop.split('_'):
                        a, b = seg.split('-')
                        segments_1based.append((int(a) + 1, int(b) + 1))
                    seq_start = segments_1based[0][0]
                    seq_end = segments_1based[-1][1]
                except Exception:
                    segments_1based = []
                chop_1based = "_".join(f"{a}-{b}" for a, b in segments_1based) or d_chop

                ted = TedDomain(
                    isoform_id=target_isoform_id,
                    ted_id=f"NOVEL-{chain_id}-{idx+1}",
                    domain_type="CHAINSAW",
                    seq_start=seq_start,
                    seq_end=seq_end,
                    # Chainsaw reports an UNCERTAINTY (lower = better), stored in
                    # the `score` column for lack of a better one — label it so.
                    score=float(uncertainty) if uncertainty else 0.0,
                    chopping=chop_1based,
                    plddt=None,
                    cath_label=None  # No CATH classification from single-tool prediction
                )
                novel_domain_objects.append(ted)
            
            row_count += 1
            if len(novel_domain_objects) >= 5000:
                db.add_all(novel_domain_objects)
                db.commit()
                print(f"[*] Committed 5000 novel domains. (Rows processed: {row_count})")
                novel_domain_objects = []

    if novel_domain_objects:
        db.add_all(novel_domain_objects)
        db.commit()
        print(f"[*] Committed remaining {len(novel_domain_objects)} domains. (Total rows: {row_count})")
    db.close()
    print("[*] Novel domain mapping complete!")


if __name__ == "__main__":
    app_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sources_dir = os.path.join(app_root, "data", "sources")
    chainsaw_dir = os.path.join(app_root, "tools", "ted_domains", "ted-tools-main", "ted_consensus_1.0", "programs", "chainsaw")
    
    if os.path.exists(sources_dir):
        for folder_name in os.listdir(sources_dir):
            af3_dir = os.path.join(sources_dir, folder_name)
            if os.path.isdir(af3_dir):
                print(f"\\n{'='*50}\\n[*] Processing source directory: {folder_name}\\n{'='*50}")
                generate_novel_ted_domains(af3_dir, chainsaw_dir)
    else:
        print(f"[!] Sources directory not found at {sources_dir}")
