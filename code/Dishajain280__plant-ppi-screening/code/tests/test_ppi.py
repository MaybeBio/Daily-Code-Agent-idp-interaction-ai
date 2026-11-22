"""Tests for the screening toolkit. Run with: python -m unittest discover -s tests"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cli  # noqa: E402  (also exercises the import path)
from ppi import batch, curate, metrics, report, results, uniprot  # noqa: E402


class TestMetrics(unittest.TestCase):
    def test_perfect_ranking(self):
        self.assertAlmostEqual(metrics.auc([0.9, 0.8, 0.2, 0.1], [1, 1, 0, 0]), 1.0)

    def test_inverted_ranking(self):
        self.assertAlmostEqual(metrics.auc([0.1, 0.2, 0.8, 0.9], [1, 1, 0, 0]), 0.0)

    def test_ties_count_as_half(self):
        self.assertAlmostEqual(metrics.auc([0.5, 0.5], [1, 0]), 0.5)

    def test_single_class_is_nan(self):
        value = metrics.auc([0.9, 0.8], [1, 1])
        self.assertNotEqual(value, value)  # NaN

    def test_recall_in_top_fraction(self):
        scores = [0.9, 0.4, 0.8, 0.3]
        labels = [1, 0, 1, 0]
        recall, k = metrics.recall_in_top_fraction(scores, labels, fraction=0.5)
        self.assertEqual(k, 2)
        self.assertAlmostEqual(recall, 1.0)

    def test_precision_at_threshold(self):
        precision, selected = metrics.precision_at([0.9, 0.7, 0.5], [1, 0, 0], 0.6)
        self.assertEqual(selected, 2)
        self.assertAlmostEqual(precision, 0.5)

    def test_bootstrap_bounds_are_ordered(self):
        low, high = metrics.bootstrap_auc([0.9, 0.8, 0.3, 0.2, 0.7, 0.1], [1, 1, 0, 0, 1, 0],
                                          n=200)
        self.assertLessEqual(low, high)


class TestSequences(unittest.TestCase):
    def test_full_sequence_passthrough(self):
        self.assertEqual(uniprot.slice_sequence("ABCDE", None, None), "ABCDE")

    def test_range_is_one_based_inclusive(self):
        self.assertEqual(uniprot.slice_sequence("ABCDE", 2, 4), "BCD")

    def test_out_of_bounds_range_rejected(self):
        with self.assertRaises(ValueError):
            uniprot.slice_sequence("ABCDE", 3, 9)


class TestResolve(unittest.TestCase):
    """Resolver behaviour, exercised without touching the network."""

    @staticmethod
    def _hit(accession: str, name: str = "protein", length: int = 100) -> dict:
        return {"primaryAccession": accession, "sequence": {"length": length},
                "proteinDescription": {"recommendedName": {"fullName": {"value": name}}}}

    def test_synonym_fallback_finds_entry_missed_by_exact_name(self):
        import io
        from unittest.mock import patch

        calls = []

        def fake_get(url, timeout=30):
            calls.append(url)
            hits = [] if "gene_exact" in url else [self._hit("P49177", "G beta")]
            return json.dumps({"results": hits}).encode()

        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(uniprot, "CACHE_DIR", Path(tmp)), \
                patch.object(uniprot, "_get", fake_get):
            entry = uniprot.resolve("AGB1")
        self.assertEqual(entry.accession, "P49177")
        self.assertEqual(len(calls), 2)  # exact pass, then synonym pass

    def test_two_hits_are_reported_as_ambiguous(self):
        from unittest.mock import patch

        def fake_get(url, timeout=30):
            hits = [self._hit("Q9MA75", "VIP1 transcription factor"),
                    self._hit("F4J8C6", "VIP1 kinase")]
            return json.dumps({"results": hits}).encode()

        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(uniprot, "CACHE_DIR", Path(tmp)), \
                patch.object(uniprot, "_get", fake_get):
            with self.assertRaises(uniprot.ResolutionError) as ctx:
                uniprot.resolve("VIP1")
        self.assertIn("ambiguous", str(ctx.exception))
        self.assertIn("Q9MA75", str(ctx.exception))

    def test_second_call_is_served_from_cache(self):
        from unittest.mock import patch

        calls = []

        def fake_get(url, timeout=30):
            calls.append(url)
            return json.dumps({"results": [self._hit("O22476", "BRI1")]}).encode()

        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(uniprot, "CACHE_DIR", Path(tmp)), \
                patch.object(uniprot, "_get", fake_get):
            uniprot.resolve("BRI1")
            first = len(calls)
            uniprot.resolve("BRI1")
            self.assertEqual(len(calls), first)  # no extra requests

    def test_missing_symbol_raises_with_guidance(self):
        from unittest.mock import patch

        def fake_get(url, timeout=30):
            return json.dumps({"results": []}).encode()

        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(uniprot, "CACHE_DIR", Path(tmp)), \
                patch.object(uniprot, "_get", fake_get):
            with self.assertRaises(uniprot.ResolutionError) as ctx:
                uniprot.resolve("NOTAGENE")
        self.assertIn("candidates.csv", str(ctx.exception))


class TestBatch(unittest.TestCase):
    def test_plan_slices_and_counts_tokens(self):
        sequences = {"P1": "A" * 100, "P2": "C" * 50}
        rows = [{"pair_id": "x", "acc_a": "P1", "acc_b": "P2",
                 "range_a": "1-20", "range_b": ""}]
        plan_obj = batch.plan(rows, sequences)
        self.assertEqual(len(plan_obj.jobs), 1)
        self.assertEqual(plan_obj.jobs[0].tokens, 70)
        self.assertEqual(plan_obj.jobs[0].chain_a, "A" * 20)

    def test_oversize_job_is_skipped_with_reason(self):
        sequences = {"P1": "A" * 4000, "P2": "C" * 2000}
        plan_obj = batch.plan([{"pair_id": "big", "acc_a": "P1", "acc_b": "P2",
                                "range_a": "", "range_b": ""}], sequences)
        self.assertEqual(plan_obj.jobs, [])
        self.assertIn("5000", plan_obj.skipped[0][1])

    def test_missing_sequence_is_reported(self):
        plan_obj = batch.plan([{"pair_id": "m", "acc_a": "P1", "acc_b": "P9",
                                "range_a": "", "range_b": ""}], {"P1": "AAA"})
        self.assertEqual(plan_obj.jobs, [])
        self.assertIn("P9", plan_obj.skipped[0][1])

    def test_batches_are_chunked_at_limit(self):
        sequences = {"P1": "A" * 10, "P2": "C" * 10}
        rows = [{"pair_id": f"p{i}", "acc_a": "P1", "acc_b": "P2",
                 "range_a": "", "range_b": ""} for i in range(205)]
        plan_obj = batch.plan(rows, sequences)
        self.assertEqual(len(plan_obj.files(max_jobs=100)), 3)

    def test_written_batch_matches_server_schema(self):
        sequences = {"P1": "MSTN", "P2": "GGGG"}
        plan_obj = batch.plan([{"pair_id": "j1", "acc_a": "P1", "acc_b": "P2",
                                "range_a": "", "range_b": ""}], sequences)
        with tempfile.TemporaryDirectory() as tmp:
            paths = batch.write_batches(plan_obj, Path(tmp))
            payload = json.loads(paths[0].read_text(encoding="utf-8"))
        self.assertEqual(payload[0]["name"], "j1")
        self.assertEqual(payload[0]["dialect"], "alphafoldserver")
        self.assertEqual(payload[0]["modelSeeds"], [])  # empty list per the server docs
        self.assertEqual(payload[0]["sequences"][0]["proteinChain"]["sequence"], "MSTN")

    def test_quota_tracker_counts_today(self):
        with tempfile.TemporaryDirectory() as tmp:
            tracker = batch.QuotaTracker(Path(tmp) / "q.json", daily_limit=5)
            self.assertEqual(tracker.remaining(), 5)
            tracker.record(2)
            self.assertEqual(tracker.remaining(), 3)


class TestResults(unittest.TestCase):
    def test_parse_summary_reads_scores(self):
        payload = {"iptm": 0.83, "ptm": 0.71, "ranking_score": 0.80,
                   "chain_pair_pae_min": [[0.0, 4.2], [5.1, 0.0]],
                   "fraction_disordered": 0.12, "has_clash": False}
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "VIP1_BZIP29"
            folder.mkdir()
            path = folder / "VIP1_BZIP29_summary_confidences.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            score = results.parse_summary(path)
        self.assertEqual(score.pair_id, "VIP1_BZIP29")
        self.assertAlmostEqual(score.iptm, 0.83)
        self.assertAlmostEqual(score.min_chain_pair_pae, 4.2)

    def test_server_fold_prefix_is_stripped(self):
        payload = {"iptm": 0.72, "ptm": 0.6}
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "fold_VIP1_BZIP29"
            folder.mkdir()
            path = folder / "fold_VIP1_BZIP29_summary_confidences.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            score = results.parse_summary(path)
        self.assertEqual(score.pair_id, "VIP1_BZIP29")

    def test_collect_writes_one_row_per_pair(self):
        payload = {"iptm": 0.5, "ptm": 0.4}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("pairA", "pairB"):
                folder = root / name
                folder.mkdir()
                (folder / f"{name}_summary_confidences.json").write_text(
                    json.dumps(payload), encoding="utf-8")
            out = root / "scores.csv"
            scores = results.collect(root, out)
            with out.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(scores), 2)
        self.assertEqual(len(rows), 2)


class TestReport(unittest.TestCase):
    def test_svg_marks_positives_and_negatives(self):
        rows = [{"pair_id": "a", "label": 1, "score": 0.9, "min_chain_pair_pae": 3.0},
                {"pair_id": "b", "label": 0, "score": 0.2, "min_chain_pair_pae": 20.0}]
        svg = report.svg_chart(rows, threshold=0.6)
        self.assertIn("<svg", svg)
        self.assertIn("#1a7f37", svg)  # filled positive
        self.assertIn("#8250df", svg)  # hollow negative

    def test_markdown_table_sorted_by_score(self):
        rows = [{"pair_id": "low", "label": 0, "score": 0.2, "min_chain_pair_pae": 20.0},
                {"pair_id": "high", "label": 1, "score": 0.9, "min_chain_pair_pae": 3.0}]
        table = report.markdown_table(rows)
        self.assertLess(table.index("high"), table.index("low"))


class TestCurate(unittest.TestCase):
    @staticmethod
    def _row(**overrides) -> dict:
        row = {"pair_id": "p1", "gene_a": "A", "gene_b": "B", "acc_a": "", "acc_b": "",
               "label": "1", "label_type": "positive_curated",
               "source": "Author et al. (2020) Plant Cell, Y2H; PMID 12345678",
               "verified": "yes"}
        row.update(overrides)
        return row

    def _errors(self, rows, **kwargs):
        return [f for f in curate.validate(rows, **kwargs) if f.level == "error"]

    def _warnings(self, rows, **kwargs):
        return [f for f in curate.validate(rows, **kwargs) if f.level == "warning"]

    def test_clean_rows_pass(self):
        rows = [self._row(pair_id="p1", gene_b="B"),
                self._row(pair_id="p2", gene_a="C", gene_b="D", label="0",
                          label_type="decoy_localization")]
        self.assertEqual(self._errors(rows), [])

    def test_reversed_pair_is_rejected(self):
        rows = [self._row(pair_id="p1", gene_a="VIP1", gene_b="BZIP29"),
                self._row(pair_id="p2", gene_a="BZIP29", gene_b="VIP1")]
        messages = " ".join(f.message for f in self._errors(rows))
        self.assertIn("already listed", messages)

    def test_duplicate_pair_id_is_rejected(self):
        rows = [self._row(pair_id="same"), self._row(pair_id="same", gene_a="C", gene_b="D")]
        self.assertTrue(any("duplicate pair_id" in f.message for f in self._errors(rows)))

    def test_placeholder_source_is_flagged(self):
        rows = [self._row(source="TODO")]
        self.assertTrue(any("no source recorded" in f.message for f in self._warnings(rows)))

    def test_label_type_mismatch_is_an_error(self):
        rows = [self._row(label="0", label_type="positive_curated")]
        self.assertTrue(any("curated positive" in f.message for f in self._errors(rows)))

    def test_decoy_heavy_file_is_flagged(self):
        rows = [self._row(pair_id=f"d{i}", gene_a=f"X{i}", gene_b=f"Y{i}", label="0",
                          label_type="decoy_localization") for i in range(5)]
        rows.append(self._row(pair_id="pos", gene_a="P", gene_b="Q"))
        self.assertTrue(any("decoys" in f.message for f in self._warnings(rows)))

    def test_small_file_warns_about_interval_width(self):
        rows = [self._row(pair_id="a"), self._row(pair_id="b", gene_a="C", gene_b="D")]
        self.assertTrue(any("confidence interval" in f.message for f in self._warnings(rows)))

    def test_cli_reports_errors(self):
        fields = ["pair_id", "gene_a", "gene_b", "acc_a", "acc_b", "range_a", "range_b",
                  "label", "label_type", "source", "verified"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cands.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(self._row(label="7"))
            code = cli.main(["check-labels", "--candidates", str(path)])
        self.assertEqual(code, 1)


class TestEvaluateGuard(unittest.TestCase):
    def _write(self, path: Path, rows: list[dict], fields: list[str]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_unverified_labels_are_refused(self):
        fields = ["pair_id", "gene_a", "gene_b", "acc_a", "acc_b", "range_a", "range_b",
                  "label", "label_type", "source", "verified"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root / "cands.csv",
                        [{"pair_id": "p1", "label": "1", "verified": "no",
                          "gene_a": "A", "gene_b": "B", "acc_a": "", "acc_b": "",
                          "range_a": "", "range_b": "", "label_type": "", "source": ""}],
                        fields)
            self._write(root / "scores.csv",
                        [{"pair_id": "p1", "iptm": "0.8", "min_chain_pair_pae": "4.0"}],
                        ["pair_id", "iptm", "min_chain_pair_pae"])
            code = cli.main(["evaluate", "--candidates", str(root / "cands.csv"),
                             "--results", str(root / "scores.csv"),
                             "--out", str(root / "report")])
        self.assertEqual(code, 2)

    def test_calibration_only_reports_and_exits_cleanly(self):
        fields = ["pair_id", "gene_a", "gene_b", "acc_a", "acc_b", "range_a", "range_b",
                  "label", "label_type", "source", "verified"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root / "cands.csv",
                        [{"pair_id": "BRI1_BAK1_full", "label": "1", "verified": "yes",
                          "label_type": "calibration", "gene_a": "BRI1", "gene_b": "BAK1",
                          "acc_a": "", "acc_b": "", "range_a": "", "range_b": "",
                          "source": "benchmark control"},
                         {"pair_id": "BRI1_BAK1_domains", "label": "1", "verified": "yes",
                          "label_type": "calibration", "gene_a": "BRI1", "gene_b": "BAK1",
                          "acc_a": "", "acc_b": "", "range_a": "", "range_b": "",
                          "source": "benchmark control"}],
                        fields)
            self._write(root / "scores.csv",
                        [{"pair_id": "BRI1_BAK1_full", "iptm": "0.38",
                          "min_chain_pair_pae": "18.0"},
                         {"pair_id": "BRI1_BAK1_domains", "iptm": "0.90",
                          "min_chain_pair_pae": "4.0"}],
                        ["pair_id", "iptm", "min_chain_pair_pae"])
            code = cli.main(["evaluate", "--candidates", str(root / "cands.csv"),
                             "--results", str(root / "scores.csv"),
                             "--out", str(root / "report")])
        self.assertEqual(code, 0)

    def test_verified_labels_produce_a_report(self):
        fields = ["pair_id", "gene_a", "gene_b", "acc_a", "acc_b", "range_a", "range_b",
                  "label", "label_type", "source", "verified"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root / "cands.csv",
                        [{"pair_id": "p1", "label": "1", "verified": "yes",
                          "gene_a": "A", "gene_b": "B", "acc_a": "", "acc_b": "",
                          "range_a": "", "range_b": "", "label_type": "", "source": ""},
                         {"pair_id": "p2", "label": "0", "verified": "yes",
                          "gene_a": "A", "gene_b": "C", "acc_a": "", "acc_b": "",
                          "range_a": "", "range_b": "", "label_type": "", "source": ""}],
                        fields)
            self._write(root / "scores.csv",
                        [{"pair_id": "p1", "iptm": "0.8", "min_chain_pair_pae": "4.0"},
                         {"pair_id": "p2", "iptm": "0.2", "min_chain_pair_pae": "21.0"}],
                        ["pair_id", "iptm", "min_chain_pair_pae"])
            code = cli.main(["evaluate", "--candidates", str(root / "cands.csv"),
                             "--results", str(root / "scores.csv"),
                             "--out", str(root / "report")])
            self.assertEqual(code, 0)
            self.assertTrue((root / "report" / "report.md").exists())
            self.assertTrue((root / "report" / "figure.svg").exists())


if __name__ == "__main__":
    unittest.main()
