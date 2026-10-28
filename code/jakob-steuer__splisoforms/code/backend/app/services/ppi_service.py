import os
import json
import requests
import Bio.PDB
from Bio.PDB.NeighborSearch import NeighborSearch
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import text
from app.models import bio as bio_models
from app.database import SessionLocal
import logging

logger = logging.getLogger(__name__)

# Constants
# Use the AFDB uniprot/summary endpoint — it returns ALL structures (monomer + complex)
# including oligomeric_state ("HOMODIMER" / "MONOMER") and the correct CIF URLs.
# The legacy /api/prediction/{id} endpoint only returns monomers when queried by UniProt ID;
# complexes (e.g. NVIDIA AlphaFold-Multimer batch, AF-000000XXXXXXXXX) are invisible there.
AFDB_SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/uniprot/summary/{uniprot_id}.json"
STRING_API_URL = "https://string-db.org/api/json/network"
STRUCTURE_CACHE_DIR = "backend/data/structures/complexes"

class PPIService:
    def __init__(self, db: Session):
        self.db = db
        os.makedirs(STRUCTURE_CACHE_DIR, exist_ok=True)

    def analyze_isoform_ppi(self, isoform_id: int):
        """Main entry point for PPI analysis for a given isoform."""
        isoform = self.db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
        if not isoform or not isoform.gene:
            return

        uniprot_id = isoform.gene.uniprot_id
        if not uniprot_id:
            return

        # 1. AFDB Homodimer Analysis
        self._run_homodimer_analysis(isoform, uniprot_id)

        # 2. STRING Interaction Analysis
        self._run_string_analysis(isoform, uniprot_id)

    def _find_complex_entries(self, uniprot_id: str):
        """
        Query the AFDB uniprot/summary endpoint and return a list of all
        HOMODIMER and HETERODIMER structures.
        """
        url = AFDB_SUMMARY_URL.format(uniprot_id=uniprot_id)
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            return []
        resp.raise_for_status()
        data = resp.json()

        complexes = []
        seen_ids = set()
        
        for struct in data.get("structures", []):
            summary = struct.get("summary", {})
            state = summary.get("oligomeric_state", "")
            model_id = summary.get("model_identifier", "")
            
            if state in ("HOMODIMER", "HETERODIMER") and model_id and model_id not in seen_ids:
                seen_ids.add(model_id)
                
                target_chains = []
                partner_chains = []
                partner_uniprot_id = None
                # Must be initialised here: for a HOMODIMER every POLYMER entity is
                # the target itself, so the `else` branch below never runs and
                # `partner_name` would otherwise be unbound at the append
                # (review PPI-1 — 1,164 UnboundLocalError per backfill_complex run).
                partner_name = None

                for entity in summary.get("entities", []):
                    if entity.get("entity_type") != "POLYMER":
                        continue
                    if entity.get("identifier") == uniprot_id:
                        target_chains.extend(entity.get("chain_ids", []))
                    else:
                        partner_chains.extend(entity.get("chain_ids", []))
                        if not partner_uniprot_id:
                            partner_uniprot_id = entity.get("identifier")
                            partner_name = entity.get("description")
                
                if not target_chains:
                    continue
                
                target_chain = target_chains[0]
                
                if state == "HOMODIMER":
                    # For homodimers, the partner chain is just the second chain of the same protein
                    partner_chain = target_chains[1] if len(target_chains) > 1 else (partner_chains[0] if partner_chains else None)
                else:
                    partner_chain = partner_chains[0] if partner_chains else None
                
                if not partner_chain:
                    continue
                    
                complexes.append({
                    "model_identifier": model_id,
                    "model_url": summary.get("model_url"),
                    "oligomeric_state": state,
                    "target_chain": target_chain,
                    "partner_chain": partner_chain,
                    "partner_uniprot_id": partner_uniprot_id,
                    "partner_name": partner_name
                })
                
        return complexes

    def _run_complex_analysis(self, isoform, uniprot_id: str):
        """Check for homodimer/heterodimers in AFDB and analyze interface disruption."""
        # Skip if already checked
        if isoform.gene.afdb_complex_checked:
            return

        try:
            complexes = self._find_complex_entries(uniprot_id)

            if not complexes:
                isoform.gene.afdb_complex_checked = True
                self.db.commit()
                return

            canonical_seq = isoform.gene.canonical_sequence or ""
            all_isoforms = (
                self.db.query(bio_models.Isoform)
                .options(selectinload(bio_models.Isoform.ptm_accessibility))
                .filter(bio_models.Isoform.gene_id == isoform.gene_id)
                .all()
            )

            for cpx in complexes:
                accession = cpx["model_identifier"]    # e.g. AF-0000000065778530
                cif_url   = cpx["model_url"]           # MMCIF URL for the complex
                state     = cpx["oligomeric_state"]
                target_ch = cpx["target_chain"]
                partner_ch = cpx["partner_chain"]
                
                # Download and cache the complex CIF
                cif_path = os.path.join(STRUCTURE_CACHE_DIR, f"{accession}.cif")
                if not os.path.exists(cif_path):
                    cif_res = requests.get(cif_url, timeout=120)
                    cif_res.raise_for_status()
                    with open(cif_path, "wb") as f:
                        f.write(cif_res.content)

                # Extract canonical interface residue positions (1-based, from target chain)
                interface_residues = self._extract_interface_residues(cif_path, target_ch, partner_ch)
                if not interface_residues:
                    logger.warning(f"No interface residues extracted from {cif_path} ({target_ch} vs {partner_ch})")
                    continue

                for iso in all_isoforms:
                    if not iso.sequence:
                        continue
                    self._analyze_complex_disruption(iso, cpx, interface_residues, canonical_seq)

            isoform.gene.afdb_complex_checked = True
            self.db.commit()

        except Exception as e:
            logger.error(f"Error in complex analysis for {uniprot_id}: {e}")

    def _extract_interface_residues(self, cif_path: str, target_ch_id: str, partner_ch_id: str):
        """Extract residues on target chain that are within 5Å of partner chain."""
        parser = Bio.PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure("complex", cif_path)

        model = structure[0]
        try:
            chain_a = model[target_ch_id]
            chain_b = model[partner_ch_id]
        except KeyError:
            return []

        atoms_a = [atom for atom in chain_a.get_atoms() if atom.element != 'H']
        atoms_b = [atom for atom in chain_b.get_atoms() if atom.element != 'H']

        ns = NeighborSearch(atoms_b)
        interface_residues = set()
        for atom_a in atoms_a:
            if ns.search(atom_a.coord, 5.0):
                res = atom_a.get_parent()
                interface_residues.add(res.get_id()[1])

        return sorted(list(interface_residues))

    def _map_canonical_positions(self, canonical_seq: str, isoform_seq: str,
                                  interface_positions: list) -> dict:
        """
        For each canonical interface position (1-based), determine its status in
        the isoform using SequenceMatcher (handles exon-skip deletions well).

        Returns {pos: 'present' | 'absent' | 'perturbed'}
          - present   : same amino acid at the aligned position
          - perturbed : different amino acid at the aligned position (substitution)
          - absent    : position falls in a gap (deleted region)
        """
        from difflib import SequenceMatcher

        if canonical_seq == isoform_seq:
            return {pos: "present" for pos in interface_positions}

        # Build mapping: canonical_0based_idx -> isoform_0based_idx
        matcher = SequenceMatcher(None, canonical_seq, isoform_seq, autojunk=False)
        canon_to_iso: dict[int, int] = {}
        for a, b, size in matcher.get_matching_blocks():
            for offset in range(size):
                canon_to_iso[a + offset] = b + offset

        result = {}
        for pos in interface_positions:
            idx = pos - 1  # 0-based
            if idx not in canon_to_iso:
                result[pos] = "absent"
            else:
                iso_idx = canon_to_iso[idx]
                c_aa = canonical_seq[idx] if idx < len(canonical_seq) else None
                i_aa = isoform_seq[iso_idx] if iso_idx < len(isoform_seq) else None
                result[pos] = "present" if c_aa == i_aa else "perturbed"

        return result

    def _analyze_complex_disruption(self, isoform, cpx: dict,
                                       interface_residues: list, canonical_seq: str):
        """
        Compare each canonical complex interface position against the isoform
        sequence via alignment.  SASA perturbation is used as a secondary signal
        when PTM accessibility data is available.
        """
        total = len(interface_residues)
        if total == 0:
            return

        iso_seq = isoform.sequence or ""

        # --- Primary signal: sequence alignment ---------------------------------
        position_status = self._map_canonical_positions(
            canonical_seq, iso_seq, interface_residues
        )

        # --- Secondary signal: SASA change from PTM accessibility data ----------
        ptm_data = {ptm.residue_position: ptm for ptm in (isoform.ptm_accessibility or [])}

        residue_status = []
        absent_count = 0
        perturbed_count = 0

        for res_num in interface_residues:
            status = position_status.get(res_num, "present")
            sasa_c = sasa_n = None
            ptm = ptm_data.get(res_num)

            if ptm:
                sasa_c = ptm.canonical_sasa
                sasa_n = ptm.novel_sasa

            if status == "present":
                # Upgrade to "perturbed" if SASA changed significantly
                if sasa_c is not None and sasa_n is not None:
                    if abs((sasa_n - sasa_c) / (sasa_c + 0.1)) > 0.2:
                        status = "perturbed"
                        perturbed_count += 1
            elif status == "absent":
                absent_count += 1
            elif status == "perturbed":
                perturbed_count += 1

            residue_status.append({
                "pos": res_num,
                "status": status,
                "sasa_canonical": sasa_c,
                "sasa_novel": sasa_n,
            })

        disruption_score = (absent_count + perturbed_count) / total
        accession = cpx["model_identifier"]
        is_homodimer = cpx["oligomeric_state"] == "HOMODIMER"

        if is_homodimer:
            existing = self.db.query(bio_models.HomodimerAnalysis).filter(
                bio_models.HomodimerAnalysis.isoform_id == isoform.id
            ).first()
            if existing:
                self.db.delete(existing)
                self.db.flush()

            analysis = bio_models.HomodimerAnalysis(
                isoform_id=isoform.id,
                afdb_accession=accession,
                total_interface_residues=total,
                count_absent=absent_count,
                count_sasa_changed=perturbed_count,
                disruption_score=disruption_score,
                interface_residues_json=json.dumps(residue_status),
            )
            self.db.add(analysis)
        else:
            existing = self.db.query(bio_models.HeterodimerAnalysis).filter(
                bio_models.HeterodimerAnalysis.isoform_id == isoform.id,
                bio_models.HeterodimerAnalysis.afdb_accession == accession
            ).first()
            if existing:
                self.db.delete(existing)
                self.db.flush()

            analysis = bio_models.HeterodimerAnalysis(
                isoform_id=isoform.id,
                afdb_accession=accession,
                partner_uniprot_id=cpx["partner_uniprot_id"],
                partner_name=cpx.get("partner_name"),
                total_interface_residues=total,
                count_absent=absent_count,
                count_sasa_changed=perturbed_count,
                disruption_score=disruption_score,
                interface_residues_json=json.dumps(residue_status),
            )
            self.db.add(analysis)

    def _run_string_analysis(self, isoform, uniprot_id: str):
        """Fetch STRING interactions and check for domain loss."""
        try:
            # Clear existing interactions to prevent infinite accumulation
            self.db.query(bio_models.StringInteraction).filter(
                bio_models.StringInteraction.isoform_id == isoform.id
            ).delete()
            self.db.commit()

            params = {
                "identifiers": uniprot_id,
                "species": 9606,
                "required_score": 700
            }
            # 30 s timeout — STRING occasionally hangs and the rest of this
            # service has timeouts on every other HTTP call (audit-5 fix).
            response = requests.get(STRING_API_URL, params=params, timeout=30)
            response.raise_for_status()
            partners = response.json()
            
            # Get Pfam domains for domain loss detection
            lost_domains = {d.pfam_name for d in isoform.pfam_domains if d.domain_status == 'lost'}
            gene_symbol = isoform.gene.symbol

            seen_partners = set()

            for p in partners:
                # The STRING API network call returns interactions for the identifier provided.
                # We identify which one is the partner.
                p_a = p.get('preferredName_A')
                p_b = p.get('preferredName_B')
                
                if p_a == gene_symbol and p_b == gene_symbol:
                    continue # skip self
                
                partner_name = p_b if p_a == gene_symbol else p_a
                
                # Deduplicate: Only keep the highest scoring interaction per partner
                if partner_name in seen_partners:
                    continue
                seen_partners.add(partner_name)

                # STRING's stringId_A / stringId_B fields are "9606.ENSP*";
                # split on '.' yields the Ensembl protein ID (ENSP*), not a
                # UniProt accession.  Stored in partner_ensp_id (renamed in
                # migration i6f1a3d8e5b2).
                partner_ensp = p.get('stringId_B' if p_a == gene_symbol else 'stringId_A', '').split('.')[-1]

                interaction = bio_models.StringInteraction(
                    isoform_id=isoform.id,
                    partner_gene_symbol=partner_name,
                    partner_ensp_id=partner_ensp,
                    string_score=int(float(p['score']) * 1000) if 'score' in p else 700,
                    evidence_json=json.dumps(p),
                    domain_loss_possible=len(lost_domains) > 0
                )
                self.db.add(interaction)

            self.db.commit()

        except Exception as e:
            logger.error(f"Error in STRING analysis for {uniprot_id}: {e}")

if __name__ == "__main__":
    # Test runner
    db = SessionLocal()
    service = PPIService(db)
    # Test with isoform 36
    service.analyze_isoform_ppi(36)
    db.close()
