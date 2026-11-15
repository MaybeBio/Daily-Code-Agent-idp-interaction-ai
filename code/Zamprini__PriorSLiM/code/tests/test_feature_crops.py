"""CPU-only crop-planning checks; does not load ESM2 or MotifAE."""

from pathlib import Path
import unittest

import pandas as pd

from src.features.build_contextual_esm2 import (
    MAX_RESIDUES, build_contextual_plan, make_idr_scan_windows,
)
from src.features.build_motifae_sparse_windows import validate_context_metadata


class ContextualESM2PlanTests(unittest.TestCase):
    def make_plan(self, protein_length, intervals, crop_mode="idr_coverage"):
        sequence = "A" * protein_length
        windows = make_idr_scan_windows({"P": sequence}, {"P": intervals})
        window_meta, chunks = build_contextual_plan(
            windows, {"P": sequence}, {"P": "P"}, {"P": 0}, Path("/tmp/contextual-plan-test"),
            crop_mode=crop_mode, idrs_by_base={"P": intervals}, context_margin=128,
        )
        return windows, window_meta, chunks

    def test_short_protein_has_one_full_length_chunk(self):
        windows, mapped, chunks = self.make_plan(800, [(1, 800)])
        self.assertEqual((len(windows), len(mapped), len(chunks)), (800, 800, 1))
        self.assertEqual((int(chunks.iloc[0].crop_start), int(chunks.iloc[0].crop_end)), (1, 800))
        self.assertEqual(int(mapped.iloc[0].left_padding), 14)
        self.assertEqual(int(mapped.iloc[-1].right_padding), 14)

    def test_long_idr_uses_738_center_cores_and_complete_windows(self):
        windows, mapped, chunks = self.make_plan(2200, [(200, 1900)])
        self.assertEqual(len(windows), 1701)
        self.assertEqual(len(chunks), 3)
        self.assertEqual(chunks[["core_start", "core_end"]].values.tolist(),
                         [[200, 937], [938, 1675], [1676, 1900]])
        self.assertTrue((chunks.crop_length == MAX_RESIDUES).all())
        self.assertTrue((mapped.window_center.between(200, 1900)).all())
        self.assertTrue((mapped.valid_protein_start >= mapped.crop_start).all())
        self.assertTrue((mapped.valid_protein_end <= mapped.crop_end).all())
        self.assertTrue((mapped.left_padding + mapped.valid_length + mapped.right_padding == 29).all())
        self.assertTrue((mapped.groupby("window_id").size() == 1).all())
        left_context = mapped.valid_protein_start - mapped.crop_start
        right_context = mapped.crop_end - mapped.valid_protein_end
        self.assertTrue((left_context >= (mapped.valid_protein_start - 1).clip(upper=128)).all())
        self.assertTrue((right_context >= (mapped.protein_length - mapped.valid_protein_end).clip(upper=128)).all())
        validate_context_metadata(mapped, chunks)

    def test_sub_1022_idr_still_splits_when_margin_requires_it(self):
        _windows, _mapped, chunks = self.make_plan(1800, [(400, 1299)])
        self.assertEqual(len(chunks), 2)
        self.assertEqual(chunks[["core_start", "core_end"]].values.tolist(),
                         [[400, 1137], [1138, 1299]])

    def test_window_mode_keeps_legacy_grouping(self):
        rows = pd.DataFrame([
            {"uniprot": "P", "window_start": start, "window_end": start + 28}
            for start in (200, 300, 700)
        ])
        mapped, chunks = build_contextual_plan(
            rows, {"P": "A" * 2200}, {"P": "P"}, {"P": 0},
            Path("/tmp/contextual-plan-test"), crop_mode="window",
        )
        self.assertEqual((len(mapped), len(chunks)), (3, 1))
        self.assertNotIn("context_margin", chunks.columns)

    def test_region_file_mismatch_is_rejected(self):
        windows = make_idr_scan_windows({"P": "A" * 800}, {"P": [(100, 200)]})
        with self.assertRaisesRegex(ValueError, "disagrees with region file"):
            build_contextual_plan(
                windows, {"P": "A" * 800}, {"P": "P"}, {"P": 0},
                Path("/tmp/contextual-plan-test"), crop_mode="idr_coverage",
                idrs_by_base={"P": [(99, 200)]},
            )


if __name__ == "__main__":
    unittest.main()
