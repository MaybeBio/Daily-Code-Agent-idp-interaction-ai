"""Resumable external ESM2 downloads and cache placement."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import download_esm2


class ESM2DownloadTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.cache = Path(self.folder.name)
        self.content = b'complete reference checkpoint'
        self.entry = dict(filename='weights.pt', bytes=len(self.content),
                          sha256=hashlib.sha256(self.content).hexdigest(),
                          url='https://dl.fbaipublicfiles.com/fair-esm/models/weights.pt')
        self.target = self.cache / self.entry['filename']
        self.partial = self.cache / (self.entry['filename'] + '.part')

    def test_cache_environment(self):
        with patch.dict(os.environ, {'TORCH_HOME': '/tmp/torch-home'}, clear=True):
            self.assertEqual(download_esm2.cache_directory(), Path('/tmp/torch-home/hub/checkpoints'))
        with patch.dict(os.environ, {'XDG_CACHE_HOME': '/tmp/cache-home'}, clear=True):
            self.assertEqual(download_esm2.cache_directory(), Path('/tmp/cache-home/torch/hub/checkpoints'))
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(download_esm2.cache_directory(), Path('~/.cache/torch/hub/checkpoints').expanduser())

    def test_official_resources(self):
        registry = json.loads((Path(download_esm2.__file__).parent / 'configs/resources.json').read_text())
        resource = registry['external_resources']['esm2']
        self.assertEqual(resource['model_name'], 'esm2_t33_650M_UR50D')
        self.assertEqual(resource['representation_layer'], 33)
        self.assertEqual(resource['embedding_dim'], 1280)
        files = resource['files']
        self.assertEqual({e['filename'] for e in files},
                         {'esm2_t33_650M_UR50D.pt', 'esm2_t33_650M_UR50D-contact-regression.pt'})
        for entry in files:
            self.assertTrue(entry['url'].startswith('https://dl.fbaipublicfiles.com/fair-esm/'))
            self.assertEqual(len(entry['sha256']), 64)
            int(entry['sha256'], 16)
            self.assertGreater(entry['bytes'], 0)

    def test_existing_valid_checkpoint_is_reused(self):
        self.target.write_bytes(self.content)
        with patch.object(subprocess, 'run') as run:
            self.assertEqual(download_esm2.fetch(self.entry, self.cache), self.target)
        run.assert_not_called()

    def test_existing_invalid_checkpoint_is_not_overwritten(self):
        self.target.write_bytes(b'incorrect')
        with patch.object(subprocess, 'run') as run, self.assertRaises(ValueError):
            download_esm2.fetch(self.entry, self.cache)
        run.assert_not_called()
        self.assertEqual(self.target.read_bytes(), b'incorrect')

    def test_resume_and_install_after_verification(self):
        self.partial.write_bytes(self.content[:5])

        def finish(command, check):
            self.assertTrue(check)
            self.assertEqual(command[command.index('--continue-at') + 1], '-')
            self.assertEqual(Path(command[command.index('--output') + 1]), self.partial)
            self.assertFalse(self.target.exists())
            self.assertEqual(self.partial.read_bytes(), self.content[:5])
            self.partial.write_bytes(self.content)

        with patch.object(download_esm2.shutil, 'which', return_value='/usr/bin/curl'), \
                patch.object(subprocess, 'run', side_effect=finish):
            download_esm2.fetch(self.entry, self.cache)
        self.assertEqual(self.target.read_bytes(), self.content)
        self.assertFalse(self.partial.exists())

    def test_interruption_preserves_partial_file(self):
        self.partial.write_bytes(self.content[:5])
        with patch.object(download_esm2.shutil, 'which', return_value='/usr/bin/curl'), \
                patch.object(subprocess, 'run', side_effect=subprocess.CalledProcessError(18, 'curl')), \
                self.assertRaises(subprocess.CalledProcessError):
            download_esm2.fetch(self.entry, self.cache)
        self.assertEqual(self.partial.read_bytes(), self.content[:5])
        self.assertFalse(self.target.exists())

    def test_invalid_download_is_not_installed(self):
        def finish(*args, **kwargs):
            self.partial.write_bytes(b'x' * len(self.content))

        with patch.object(download_esm2.shutil, 'which', return_value='/usr/bin/curl'), \
                patch.object(subprocess, 'run', side_effect=finish), self.assertRaises(ValueError):
            download_esm2.fetch(self.entry, self.cache)
        self.assertFalse(self.target.exists())

    def test_completed_partial_file_needs_no_download(self):
        self.partial.write_bytes(self.content)
        with patch.object(subprocess, 'run') as run:
            download_esm2.fetch(self.entry, self.cache)
        run.assert_not_called()
        self.assertEqual(self.target.read_bytes(), self.content)


if __name__ == '__main__':
    unittest.main()
