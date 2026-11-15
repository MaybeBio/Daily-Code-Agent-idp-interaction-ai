"""Model loading, numerical replay and CV selection regression tests."""
import os
from pathlib import Path
import sys
import unittest
import tempfile
import hashlib
import io
import json
from contextlib import redirect_stdout
from unittest.mock import patch
import numpy as np
import torch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import model
import predict
import train
import download

TEST_CHECKPOINT = Path(os.environ.get('PRIORSLIM_TEST_CHECKPOINT', str(predict.DEFAULT_CHECKPOINT)))


class CoreTests(unittest.TestCase):
    def test_model_dimensions(self):
        self.assertEqual(model.WINDOW_LENGTH, 29)
        self.assertEqual(model.LATENT_DIM, 40960)
        self.assertEqual(set(model.MOTIF_INDICES) | set(model.CHEMISTRY_INDICES), set(range(58)))
        self.assertFalse(set(model.MOTIF_INDICES) & set(model.CHEMISTRY_INDICES))

    def test_masked_summaries(self):
        values = torch.zeros(2, 29, 7)
        mask = torch.ones(2, 29, dtype=torch.bool)
        summary = model.compute_pssm_center_summaries(values, mask)
        self.assertEqual(tuple(summary.shape), (2, 22))
        self.assertTrue(torch.isfinite(summary).all())

    def test_resource_checksum(self):
        content = b'PriorSLiM checksum test'
        entry = {'sha256': hashlib.sha256(content).hexdigest(), 'bytes': len(content), 'filename': 'model.pt', 'url': None}
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / 'model.pt'
            target.write_bytes(content)
            self.assertEqual(download.fetch(entry, target), target)
            target.write_bytes(b'changed')
            with self.assertRaises(ValueError):
                download.fetch(entry, target)

    def test_application_model_manifest(self):
        path = Path(download.__file__).parent / 'configs/resources.json'
        registry = json.loads(path.read_text())
        application = registry['assets']['model']
        source = registry['experiment_models'][application['source_asset']]
        self.assertEqual(application['filename'], 'model.pt')
        self.assertEqual(application['source_asset'], 'fold_08_K619')
        for key in ('sha256', 'bytes'):
            self.assertEqual(application[key], source[key])

    def test_experiment_models_are_archive_members(self):
        path = Path(download.__file__).parent / 'configs/resources.json'
        registry = json.loads(path.read_text())
        self.assertEqual(set(registry['assets']), {'model', 'models'})
        self.assertEqual(registry['assets']['models']['filename'], 'models.tar.gz')
        self.assertEqual(len(registry['experiment_models']), 17)
        for name, entry in registry['experiment_models'].items():
            fold, k = name.removeprefix('fold_').split('_K')
            self.assertEqual(entry['archive_asset'], 'models')
            self.assertEqual(entry['member'],
                             f'train/outputs_jwarm_cv10/seed42/fold_{fold}/K{k}/jwarm_bundle.pt')

    def test_experiment_download_redirects_to_archive(self):
        with patch.object(sys, 'argv', ['download.py', '--asset', 'fold_06_K781']), \
                patch.object(download, 'fetch') as fetch, patch('sys.stderr', new_callable=io.StringIO) as error:
            with self.assertRaises(SystemExit) as exit_status:
                download.main()
        self.assertEqual(exit_status.exception.code, 2)
        self.assertIn('--asset models', error.getvalue())
        fetch.assert_not_called()

    def test_models_archive_download(self):
        with patch.object(sys, 'argv', ['download.py', '--asset', 'models', '--output-dir', 'data/archives']), \
                patch.object(download, 'fetch') as fetch, redirect_stdout(io.StringIO()):
            download.main()
        entry, destination, base_url = fetch.call_args.args
        self.assertEqual(entry['filename'], 'models.tar.gz')
        self.assertEqual(destination, Path('data/archives/models.tar.gz'))

    def test_prediction_default_and_override(self):
        self.assertEqual(predict.parse_args(['--example']).checkpoint, predict.DEFAULT_CHECKPOINT)
        overridden = predict.parse_args(['--checkpoint', 'another.pt', '--example'])
        self.assertEqual(overridden.checkpoint, Path('another.pt'))

    def test_download_defaults_to_application_model(self):
        with patch.object(sys, 'argv', ['download.py']), patch.object(download, 'fetch') as fetch:
            with redirect_stdout(io.StringIO()):
                download.main()
        entry, destination, base_url = fetch.call_args.args
        self.assertEqual(entry['filename'], 'model.pt')
        self.assertEqual(destination, predict.DEFAULT_CHECKPOINT)
        self.assertIsNone(base_url)

    def test_summary_gradients(self):
        branch = model.PSSMSummaryBranch().eval()
        inputs = torch.zeros(2, 29, 7)
        embedding, _ = branch(inputs)
        embedding.sum().backward()
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in branch.parameters() if p.grad is not None))

    def test_validation_only_k_selection(self):
        points = [dict(K=k, eligible=True, converged=True, lambda_ratio=0.1, val={'ap': ap})
                  for k, ap in [(100, 0.50), (200, 0.504), (300, 0.506)]]
        points.append(dict(K=1, eligible=True, converged=False, lambda_ratio=0.2, val={'ap': 0.99}))
        chosen = train.choose_automatic_points(points)
        self.assertEqual(chosen['performance']['K'], 300)
        self.assertEqual(chosen['efficiency']['K'], 200)

    def test_cv_defaults(self):
        opts = train.args_from_cli(['check', '--device', 'cpu'])
        self.assertEqual(opts.seed, 42)
        self.assertEqual(opts.folds, list(range(1, 11)))
        self.assertFalse(opts.evaluate_test)

    @unittest.skipUnless(TEST_CHECKPOINT.is_file(), 'Download model.pt to enable prediction replay')
    def test_frozen_prediction(self):
        fixture = Path(__file__).resolve().parents[1] / 'examples/synthetic_windows.npz'
        logits = predict.predict_example(TEST_CHECKPOINT, fixture)
        self.assertEqual(logits.shape, (32,))
        self.assertTrue(np.isfinite(logits).all())


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
