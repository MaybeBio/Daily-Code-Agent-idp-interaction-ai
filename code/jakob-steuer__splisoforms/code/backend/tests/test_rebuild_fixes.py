"""
Regression checks for the ingest / analysis fixes made for the 2026-10 rebuild
(October 2026 methods review): ING-1/ING-2, NMD-2, DOM-3, EXP-1, ENG-2.

No database is needed. Run from ``backend/``:

    .venv/bin/python -m unittest tests.test_rebuild_fixes -v
"""
import importlib.util
import os
import sys
import unittest
from pathlib import Path

BACKEND = Path(__file__).resolve().parent.parent
REPO = BACKEND.parent
sys.path.insert(0, str(BACKEND))
sys.path.insert(0, str(BACKEND / "scripts"))
os.environ.setdefault("DATABASE_URL", "postgresql://x:y@localhost:1/never_connected")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── ING-1 / ING-2: ingest priority order ─────────────────────────────────────
class TestIngestPriority(unittest.TestCase):
    def setUp(self):
        self.ing = _load(BACKEND / "scripts" / "ingest_from_metadata.py", "ingest_fm")

    def _row(self, **kw):
        base = {"is_mane_select": "False", "is_canonical": "False", "has_structure": "True",
                "gene_symbol": "GENE", "gencode_tags": ""}
        base.update(kw)
        return base

    def test_mane_beats_file_order(self):
        rows = [self._row(gene_symbol="ZNF268"),                       # non-MANE, first in file
                self._row(gene_symbol="ZNF268", is_mane_select="True", is_canonical="True")]
        order = sorted(range(2), key=lambda i: self.ing.row_priority(rows[i], i))
        self.assertEqual(order, [1, 0])

    def test_readthrough_and_bare_ensg_lose(self):
        rows = [self._row(gene_symbol="ENSG00000289694", gencode_tags="readthrough_transcript,basic",
                          is_canonical="True"),                     # the readthrough gene row
                self._row(gene_symbol="CDC42", is_canonical="True")]
        order = sorted(range(2), key=lambda i: self.ing.row_priority(rows[i], i))
        self.assertEqual(order, [1, 0], "CDC42 must survive over the readthrough gene")

    def test_file_order_is_the_last_tiebreak(self):
        rows = [self._row(), self._row()]
        order = sorted(range(2), key=lambda i: self.ing.row_priority(rows[i], i))
        self.assertEqual(order, [0, 1])

    def test_structure_beats_no_structure(self):
        rows = [self._row(has_structure="False"), self._row(has_structure="True")]
        order = sorted(range(2), key=lambda i: self.ing.row_priority(rows[i], i))
        self.assertEqual(order, [1, 0])


# ── EXP-1: ccRCC sample provenance ───────────────────────────────────────────
class TestCcrccSampleProvenance(unittest.TestCase):
    def setUp(self):
        self.bsm = _load(REPO / "scripts" / "build_source_metadata.py", "build_source_metadata")

    def test_pep_file_maps_to_sample(self):
        f = self.bsm.ccrcc_sample_from_source_file
        self.assertEqual(f("p28443_o31598_1_threecells_combined_best_ORFsequences_x.fa"), "Normal")
        self.assertEqual(f("p28443_o31598_2_threecells_combined_best_ORFsequences_x.fa"), "ccRCC2")
        self.assertEqual(f("p28443_o31598_5_best_ORFsequences_Novel_Complete_forAF_NoCutoff_pep.fa"), "ccRCC5")
        self.assertEqual(f(""), "")
        self.assertEqual(f("something_else.fa"), "")

    def test_schema_has_optional_columns(self):
        self.assertEqual(self.bsm.COLUMNS[-2:], ["sample_source", "gencode_tags"])
        self.assertEqual(len(self.bsm.COLUMNS), 31)


# ── NMD-2: CDS beyond the transcript end must be unknown, not escape ─────────
class TestNmdCdsBeyondEnd(unittest.TestCase):
    def setUp(self):
        from app.analysis import nmd_prediction as nmd
        self.nmd = nmd

    def test_cds_end_past_transcript_is_unknown(self):
        # two exons of 100 nt each → 200-nt transcript; CDS claims to end at 230
        status, has_ptc, _ptc, dist, conf, evidence = self.nmd.classify_nmd_50nt_rule(
            "[[1000,1099],[2000,2099]]", "+", 1, 230)
        self.assertEqual(status, "unknown")
        self.assertFalse(has_ptc)
        self.assertIsNone(dist)
        self.assertIn("beyond the transcript end", evidence)

    def test_regular_escape_still_escape(self):
        # stop in the last exon → escape (CDS 1..180 of a 200-nt transcript)
        # 5'-proximal exception must not fire: start → stop distance 177 nt > 150
        status, *_ = self.nmd.classify_nmd_50nt_rule("[[1000,1099],[2000,2099]]", "+", 1, 180)
        self.assertEqual(status, "nmd_escape")

    def test_regular_target_still_target(self):
        # three exons of 300 nt; CDS 1..600 → stop at nt 597 (0-based), last EJC at 600,
        # distance = 600-597 = 3 → escape; CDS 1..300 → stop 297, last EJC 600 → distance 303 → target
        status, has_ptc, _p, dist, *_ = self.nmd.classify_nmd_50nt_rule(
            "[[1,300],[1000,1299],[2000,2299]]", "+", 1, 300)
        self.assertEqual(status, "nmd_target")
        self.assertTrue(has_ptc)
        self.assertGreater(dist, 50)


# ── DOM-3: Pfam clan-overlap resolution ──────────────────────────────────────
class TestPfamClanOverlap(unittest.TestCase):
    def setUp(self):
        from app.analysis import pfam_domains as pf
        self.pf = pf

    def _hit(self, acc, clan, s, e, ev, bits):
        return self.pf.DomainHit(pfam_accession=acc, pfam_name=acc, pfam_description="",
                                 pfam_clan=clan, seq_start=s, seq_end=e, hmm_start=1,
                                 hmm_end=23, hmm_length=23, domain_coverage=1.0,
                                 e_value=ev, bit_score=bits)

    def test_same_clan_overlap_keeps_best_evalue(self):
        hits = [self._hit("PF13894", "CL0361", 361, 383, 8.9e-4, 20.4),   # zf-C2H2_4
                self._hit("PF00096", "CL0361", 361, 383, 6.8e-6, 26.8),   # zf-C2H2 (best)
                self._hit("PF13912", "CL0361", 361, 384, 5.7e-2, 14.0)]   # zf-C2H2_6
        kept = self.pf.resolve_clan_overlaps(hits)
        self.assertEqual([h.pfam_accession for h in kept], ["PF00096"])

    def test_non_overlapping_same_clan_both_kept(self):
        hits = [self._hit("PF00096", "CL0361", 307, 327, 1e-5, 25.0),
                self._hit("PF00096", "CL0361", 333, 355, 3e-3, 18.3)]
        self.assertEqual(len(self.pf.resolve_clan_overlaps(hits)), 2)

    def test_different_clan_or_no_clan_untouched(self):
        hits = [self._hit("PF00069", "CL0016", 10, 260, 1e-49, 168.1),
                self._hit("PF04382", None, 100, 150, 1e-21, 77.0),       # no clan → kept
                self._hit("PF01352", "CL0999", 120, 160, 2e-23, 82.0)]   # other clan → kept
        self.assertEqual(len(self.pf.resolve_clan_overlaps(hits)), 3)

    def test_output_sorted_by_start(self):
        hits = [self._hit("PF00096", "CL0361", 333, 355, 3e-3, 18.3),
                self._hit("PF00096", "CL0361", 307, 327, 1e-5, 25.0)]
        kept = self.pf.resolve_clan_overlaps(hits)
        self.assertEqual([h.seq_start for h in sorted(kept, key=lambda h: h.seq_start)],
                         [307, 333])


# ── ENG-2: exit status for the backfill chain ────────────────────────────────
class TestExitStatus(unittest.TestCase):
    def test_exit_status(self):
        from _backfill_io import exit_status
        self.assertEqual(exit_status(0, 0), 0)
        self.assertEqual(exit_status(100, 0), 0)
        self.assertEqual(exit_status(100, 10), 0)
        self.assertEqual(exit_status(100, 51), 1)
        self.assertEqual(exit_status(100, 100), 1)
        self.assertEqual(exit_status(1, 1), 1)


if __name__ == "__main__":
    unittest.main()
