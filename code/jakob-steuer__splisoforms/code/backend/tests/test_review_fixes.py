"""
Regression checks for the four code bugs fixed after the 2026-10-19 methods
review (October 2026 methods review): PPI-1, STR-1, DOM-2, STR-3.

No database is needed — sessions are faked.  Run from ``backend/``:

    .venv/bin/python -m unittest tests.test_review_fixes -v

The two data-dependent checks (DOM-2 against ``hmmscan --cut_ga``, STR-3 on a
real AF3 CIF) are skipped when Pfam-A / a local structure is not present.
"""
import importlib.util
import math
import os
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

BACKEND = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BACKEND))
os.environ.setdefault("DATABASE_URL", "postgresql://x:y@localhost:1/never_connected")

from app.services import ppi_service  # noqa: E402


# ── PPI-1 ─────────────────────────────────────────────────────────────────────
class TestPPI1HomodimerFirst(unittest.TestCase):
    """A gene whose FIRST AFDB dimer model is a homodimer must not raise."""

    FIXTURE = {"structures": [
        {"summary": {"oligomeric_state": "HOMODIMER", "model_identifier": "AF-HOMO-1",
                     "model_url": "https://x/homo.cif",
                     "entities": [{"entity_type": "POLYMER", "identifier": "P12345",
                                   "chain_ids": ["A", "B"], "description": "Target"}]}},
        {"summary": {"oligomeric_state": "HETERODIMER", "model_identifier": "AF-HET-1",
                     "model_url": "https://x/het.cif",
                     "entities": [{"entity_type": "POLYMER", "identifier": "P12345", "chain_ids": ["A"]},
                                  {"entity_type": "POLYMER", "identifier": "Q99999", "chain_ids": ["B"],
                                   "description": "Partner protein"}]}},
    ]}

    def test_homodimer_first_parses(self):
        resp = MagicMock(status_code=200)
        resp.json.return_value = self.FIXTURE
        svc = ppi_service.PPIService.__new__(ppi_service.PPIService)  # no db / mkdir
        with patch.object(ppi_service.requests, "get", return_value=resp):
            out = svc._find_complex_entries("P12345")

        self.assertEqual([c["oligomeric_state"] for c in out], ["HOMODIMER", "HETERODIMER"])
        self.assertIsNone(out[0]["partner_name"])
        self.assertIsNone(out[0]["partner_uniprot_id"])
        self.assertEqual((out[0]["target_chain"], out[0]["partner_chain"]), ("A", "B"))
        self.assertEqual(out[1]["partner_name"], "Partner protein")
        self.assertEqual(out[1]["partner_uniprot_id"], "Q99999")


# ── STR-1 ─────────────────────────────────────────────────────────────────────
class TestSTR1MetricsAttachedBeforeAnalyze(unittest.TestCase):
    """backfill_structural_metrics_v4 must attach the new StructuralMetrics to
    ``isoform.structural_metrics`` before structural_impact.analyze_isoform runs,
    otherwise analyze_isoform inserts a second row (UniqueViolation)."""

    def _load_script(self):
        path = BACKEND / "scripts" / "backfill_structural_metrics_v4.py"
        spec = importlib.util.spec_from_file_location("v4_under_test", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_attach_then_analyze(self):
        v4 = self._load_script()

        class FakeIso:
            id = 1
            transcript_id = "ENST_TEST"
            structural_metrics = None
        iso = FakeIso()

        class FakeQuery:
            def __init__(self, model): self.model = model
            def join(self, *a, **k): return self
            def outerjoin(self, *a, **k): return self
            def filter(self, *a, **k): return self
            def limit(self, n): return self
            def all(self): return [iso]
            def get(self, i): return iso
            def first(self): return None  # no Structure → PAE step skipped

        class FakeDB:
            def query(self, model): return FakeQuery(model)
            def add(self, o): pass
            def flush(self): pass
            def commit(self): pass
            def rollback(self): pass
            def close(self): pass

        seen = {}

        def fake_analyze(db, isoform):
            seen["attached"] = isoform.structural_metrics is not None
            if isoform.structural_metrics is None:
                raise RuntimeError("second structural_metrics row → UniqueViolation")
            isoform.structural_metrics.plddt_series = [90.0]

        with patch.object(v4, "SessionLocal", FakeDB), \
             patch.object(v4.structural_impact, "analyze_isoform", fake_analyze), \
             patch.object(v4, "PAEService"):
            processed, errors = v4.backfill_v4_metrics(limit=1)

        self.assertTrue(seen.get("attached"))
        self.assertEqual((processed, errors), (1, 0))


# ── DOM-2 ─────────────────────────────────────────────────────────────────────
class TestDOM2PfamCoordinates(unittest.TestCase):
    """pyhmmer alignment coordinates are already 1-based: scan_sequence must
    report exactly what ``hmmscan --cut_ga`` reports (isoform 1211, Pkinase)."""

    # ENST00000619836 (DB isoform 1211); hmmscan --cut_ga Pfam-A (2026-10-19):
    #   Pkinase        hmm 1-255  ali 8-324   (len 262)
    #   PK_Tyr_Ser-Thr hmm 3-86   ali 10-88   (len 258)
    #   PK_Tyr_Ser-Thr hmm 75-250 ali 114-321 (len 258)
    SEQ = ("MVSSQPKYDLIREVGRGSYGVVYEAVIRKTSARVAVKKIRCHAPENVELALREFWALSSIKSQHPNVIHLEECILQKDGM"
           "VQKMSHGSNSSLYLQLVETSLKGEIAFDPRSAYYLWFVMDFCDGGDMNEYLLSRKPNRKTNTSFMLQLSSALAFLHKNQI"
           "IHRDLKPDNILISQTRLDTSDLEPTLKVADFGLSKVCSASGQNPEEPVSVNKCFLSTACGTDFYMAPEVWEGHYTAKADI"
           "FALGIIIWAMLERITFIDTETKKELLGSYVKQGTEIVPVGEALLENPKMELLIPVKKKSMNGRMKQLIKEMLAANPQDRP"
           "DAFELELRLVQIAFKDSSWET")
    EXPECTED = {("PF00069", 1, 255, 8, 324), ("PF07714", 3, 86, 10, 88), ("PF07714", 75, 250, 114, 321)}

    def test_matches_hmmscan(self):
        from app.analysis.pfam_db import get_pfam_hmm_path
        try:
            get_pfam_hmm_path()
        except Exception:
            self.skipTest("Pfam-A.hmm not available (set PFAM_HMM_PATH)")
        from app.analysis.pfam_domains import scan_sequence
        got = {(h.pfam_accession, h.hmm_start, h.hmm_end, h.seq_start, h.seq_end)
               for h in scan_sequence(self.SEQ, "iso1211")}
        self.assertEqual(got, self.EXPECTED)


# ── STR-3 ─────────────────────────────────────────────────────────────────────
class TestSTR3PhiPsi(unittest.TestCase):
    """phi/psi must be populated (degrees) for every non-terminal residue and
    agree with an independent dihedral computation."""

    CIF = Path(os.environ.get(
        "SPLISOFORMS_TEST_CIF",
        BACKEND / "static" / "data" / "structures" / "gencode_v47"
        / "2fecdc87f4b369f8464c94d42fc8aad5" / "2fecdc87f4b369f8464c94d42fc8aad5_model.cif"))

    def test_phi_psi_populated_in_degrees(self):
        if not self.CIF.exists():
            self.skipTest("local AF3 structure not present")
        import gemmi
        from Bio.PDB.vectors import Vector, calc_dihedral
        from app.analysis.structural_impact import extract_residue_features

        feats = extract_residue_features(self.CIF)
        self.assertGreater(len(feats), 10)
        self.assertIsNone(feats[0].phi)
        self.assertIsNone(feats[-1].psi)
        self.assertEqual(sum(f.phi is None for f in feats), 1)
        self.assertEqual(sum(f.psi is None for f in feats), 1)
        self.assertTrue(all(-180 <= f.phi <= 180 for f in feats if f.phi is not None))

        poly = gemmi.read_structure(str(self.CIF))[0][0].get_polymer()

        def P(res, name):
            a = res.find_atom(name, "\0")
            return Vector(a.pos.x, a.pos.y, a.pos.z)

        for i in (1, 10, len(poly) // 2, len(poly) - 2):
            phi = math.degrees(calc_dihedral(P(poly[i - 1], "C"), P(poly[i], "N"), P(poly[i], "CA"), P(poly[i], "C")))
            psi = math.degrees(calc_dihedral(P(poly[i], "N"), P(poly[i], "CA"), P(poly[i], "C"), P(poly[i + 1], "N")))
            self.assertAlmostEqual(feats[i].phi, phi, delta=0.01)
            self.assertAlmostEqual(feats[i].psi, psi, delta=0.01)


# ── Atlas latency (2026-10-27) ────────────────────────────────────────────────
class TestAtlasPartialIndexesDeclared(unittest.TestCase):
    """The /neoantigens/top partial indexes must live on the model, not only in a
    migration: reset_db.py (create_all + stamp head) drops migration-only SQL at every
    rebuild — the v2 rebuild lost them and the atlas went 0.85 s → 2–4 s."""

    def test_partial_indexes_on_model(self):
        from app.models.bio import NeoantigenPrediction
        by_name = {ix.name: ix for ix in NeoantigenPrediction.__table__.indexes}
        for name, cols in {
            "ix_neo_representative": ["priority_tier", "consensus_rank"],
            "ix_neo_repr_isoform": ["isoform_id"],
            "ix_neo_repr_stats": ["peptide_sequence"],
        }.items():
            self.assertIn(name, by_name, name)
            ix = by_name[name]
            self.assertEqual([c.name for c in ix.columns], cols, name)
            self.assertEqual(str(ix.dialect_options["postgresql"]["where"]), "is_representative", name)
        # every column the stats block counts is covered, so it stays index-only
        include = set(by_name["ix_neo_repr_stats"].dialect_options["postgresql"]["include"])
        self.assertTrue({"priority_tier", "presentation_tier", "is_immunogenic", "is_expressed",
                         "is_ba_rescue", "is_proteome_unique", "netmhcpan_rank_el",
                         "ptm_impact_status"} <= include, include)

    def test_atlas_stats_sql_keeps_distinct_separate(self):
        """COUNT(DISTINCT peptide) folded into the FILTER-count statement forces an
        on-disk sort of the whole representative set (~3.5 s); it must be its own
        statement so it reads the peptide-leading index presorted."""
        from app.api import bio as bio_api
        db = MagicMock()
        db.execute.return_value.mappings.return_value.first.return_value = {"total": 7, "p1": 1}
        db.execute.return_value.scalar.return_value = 3
        total, stats = bio_api._atlas_count_and_stats(db, "n.is_representative", "", {})
        sqls = [str(c.args[0]) for c in db.execute.call_args_list]
        self.assertEqual(len(sqls), 3)
        self.assertIn("FILTER", sqls[0]); self.assertNotIn("DISTINCT", sqls[0])
        self.assertIn("COUNT(DISTINCT n.peptide_sequence)", sqls[1])
        self.assertEqual((total, stats["unique_peptides"], stats["unique_genes"]), (7, 3, 3))

    def test_atlas_stats_cache_ttl(self):
        from app.api import bio as bio_api
        bio_api._atlas_stats_cache.clear()
        key = ("n.is_representative", ())
        with patch.object(bio_api, "_ATLAS_STATS_TTL_S", 600), \
             patch.object(bio_api.time, "monotonic", side_effect=[1000.0, 1001.0, 1601.0]):
            bio_api._atlas_stats_put(key, 5, {"total": 5})
            self.assertEqual(bio_api._atlas_stats_get(key), (5, {"total": 5}))
            self.assertIsNone(bio_api._atlas_stats_get(key))   # expired
        with patch.object(bio_api, "_ATLAS_STATS_TTL_S", 0):
            bio_api._atlas_stats_put(key, 5, {"total": 5})
            self.assertIsNone(bio_api._atlas_stats_get(key))   # disabled


if __name__ == "__main__":
    unittest.main()
