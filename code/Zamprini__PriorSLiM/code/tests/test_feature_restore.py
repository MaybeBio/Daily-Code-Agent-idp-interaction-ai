"""Lossless sparse restoration, partial reruns and integrity failures."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import restore_features as restore


class FeatureRestoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.store = self.root / 'store'
        self.store.mkdir()
        self.ids = np.array([7, 40959, 2], dtype='<i4')
        self.values = np.array([0.125, 9.5, 0.33333334], dtype='<f4')
        self.offsets = np.array([0, 2, 2, 3], dtype='<i8')
        self.references = np.array([0, 1, 2, 0, 2, 1], dtype='<i4')
        self.archive = self.store / 'motifae_residue_dictionary.npz'
        self.save_dictionary()
        counts = np.diff(self.offsets)[self.references]
        self.ptr = np.r_[np.int64(0), counts.cumsum()]
        self.ptr_path = self.store / 'motifae_residue_indptr.int64.npy'
        np.save(self.ptr_path, self.ptr)
        self.expected = []
        self.outputs = []
        for name, source in zip(restore.OUTPUT_NAMES, [self.ids, self.values]):
            data = np.concatenate([source[self.offsets[j]:self.offsets[j + 1]] for j in self.references]).tobytes()
            self.expected.append(data)
            self.outputs.append(dict(path='store/' + name, bytes=len(data),
                                     sha256=hashlib.sha256(data).hexdigest()))
        self.record = dict(**self.file_record(self.archive), format=restore.FORMAT,
                           residue_indptr=self.file_record(self.ptr_path), restores=self.outputs)

    def file_record(self, path):
        return dict(path=path.relative_to(self.root).as_posix(), bytes=path.stat().st_size,
                    sha256=restore.file_sha256(path))

    def save_dictionary(self):
        np.savez_compressed(self.archive, feature_indices=self.ids, feature_values=self.values,
                            unique_indptr=self.offsets, residue_rows=self.references)

    def test_exact_roundtrip(self):
        self.assertEqual(restore.restore_store(self.root, self.record, block_rows=2), 'restored')
        for output, expected in zip(self.outputs, self.expected):
            self.assertEqual((self.root / output['path']).read_bytes(), expected)

    def test_repeated_run_verifies_existing_files(self):
        restore.restore_store(self.root, self.record)
        times = [(self.root / r['path']).stat().st_mtime_ns for r in self.outputs]
        self.assertEqual(restore.restore_store(self.root, self.record), 'verified')
        self.assertEqual(times, [(self.root / r['path']).stat().st_mtime_ns for r in self.outputs])

    def test_partial_run_restores_missing_file_only(self):
        existing = self.root / self.outputs[0]['path']
        existing.write_bytes(self.expected[0])
        timestamp = existing.stat().st_mtime_ns
        restore.restore_store(self.root, self.record)
        self.assertEqual(existing.stat().st_mtime_ns, timestamp)
        self.assertEqual((self.root / self.outputs[1]['path']).read_bytes(), self.expected[1])

    def test_invalid_existing_file_is_not_overwritten(self):
        target = self.root / self.outputs[0]['path']
        target.write_bytes(b'incorrect')
        with self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            restore.restore_store(self.root, self.record)
        self.assertEqual(target.read_bytes(), b'incorrect')

    def test_corrupt_archive_is_rejected(self):
        self.archive.write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            restore.restore_store(self.root, self.record)
        self.assertFalse((self.root / self.outputs[0]['path']).exists())

    def test_invalid_row_reference_is_rejected(self):
        self.references[0] = 99
        self.save_dictionary()
        self.record.update(self.file_record(self.archive))
        with self.assertRaisesRegex(ValueError, 'outside the dictionary'):
            restore.restore_store(self.root, self.record)

    def test_wrong_original_offsets_are_rejected(self):
        self.ptr[1] = 1
        np.save(self.ptr_path, self.ptr)
        self.record['residue_indptr'] = self.file_record(self.ptr_path)
        with self.assertRaisesRegex(ValueError, 'row sizes disagree'):
            restore.restore_store(self.root, self.record)
        self.assertFalse(list(self.store.glob('*.part')))

    def test_wrong_output_hash_leaves_no_partial_outputs(self):
        self.outputs[1]['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'Restored checksum mismatch'):
            restore.restore_store(self.root, self.record)
        self.assertFalse((self.root / self.outputs[0]['path']).exists())
        self.assertFalse(list(self.store.glob('*.part')))

    def test_verify_only_requires_restored_files(self):
        with self.assertRaises(FileNotFoundError):
            restore.restore_store(self.root, self.record, verify_only=True)
        restore.restore_store(self.root, self.record)
        self.assertEqual(restore.restore_store(self.root, self.record, verify_only=True), 'verified')

    def test_paths_cannot_escape_data_root(self):
        bad = json.loads(json.dumps(self.record))
        bad['restores'][0]['path'] = '../outside.bin'
        with self.assertRaisesRegex(ValueError, 'outside the data directory'):
            restore.restore_store(self.root, bad)


if __name__ == '__main__':
    unittest.main()
