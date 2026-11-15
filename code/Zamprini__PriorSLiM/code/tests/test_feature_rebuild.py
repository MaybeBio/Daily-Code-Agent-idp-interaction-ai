import gzip
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SPEC = importlib.util.spec_from_file_location('rebuild_features', Path(__file__).parents[1] / 'src/rebuild_features.py')
rebuild = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rebuild)


class FeatureRebuildTests(unittest.TestCase):
    def test_sequence_and_crop_checks(self):
        sequence = 'ACDEFGHIKLMNPQRSTVWY' * 3
        windows = pd.DataFrame([dict(window_id='w1', uniprot='P1', window_start=2,
                                     window_end=30, sequence=sequence[1:30])])
        context = pd.DataFrame([dict(window_id='w1')])
        chunks = pd.DataFrame([dict(protein_id='P1', chunk_id=0, protein_length=60,
                                    crop_start=2, crop_end=50,
                                    sequence_sha256=hashlib.sha256(sequence[1:50].encode()).hexdigest())])
        rebuild.validate_sequences(windows, context, chunks, {'P1': sequence})
        with self.assertRaisesRegex(ValueError, 'Window/sequence'):
            rebuild.validate_sequences(windows, context, chunks, {'P1': 'A' * 60})
        chunks.loc[0, 'sequence_sha256'] = 'changed'
        with self.assertRaisesRegex(ValueError, 'Crop/sequence'):
            rebuild.validate_sequences(windows, context, chunks, {'P1': sequence})

    def test_terminal_padding(self):
        seq = 'ACDEFGHIKLMNPQRSTVWY' * 2
        window = pd.DataFrame([dict(window_id='w', uniprot='P', window_start=-13,
                                    window_end=15, sequence='X' * 14 + seq[:15])])
        chunks = pd.DataFrame([dict(protein_id='P', chunk_id=0, protein_length=40,
                                    crop_start=1, crop_end=40,
                                    sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())])
        rebuild.validate_sequences(window, window[['window_id']], chunks, {'P': seq})

    def test_paths_stay_inside_data_root(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.assertEqual(rebuild.under(root, 'feature_inputs/file'), root.resolve() / 'feature_inputs/file')
            with self.assertRaises(ValueError):
                rebuild.under(root, '../elsewhere')

    def test_fasta_duplicates(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'sequences.fasta.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write('>P1\nACD\nEF\n>P2\nHHH\n')
            self.assertEqual(rebuild.read_sequences(path), {'P1': 'ACDEF', 'P2': 'HHH'})
            with gzip.open(path, 'wt') as handle:
                handle.write('>P1\nACD\n>P1\nEF\n')
            with self.assertRaises(ValueError):
                rebuild.read_sequences(path)


if __name__ == '__main__':
    unittest.main()
