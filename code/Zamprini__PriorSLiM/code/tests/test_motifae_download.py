"""Version-pinned MotifAE downloads and offline verification."""
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
import download_motifae
import download_esm2


class MotifAEDownloadTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.directory = Path(self.folder.name)
        self.content = b'reference MotifAE checkpoint'
        self.entry = dict(filename='MotifAE_k3_step_80000.pt', bytes=len(self.content),
                          sha256=hashlib.sha256(self.content).hexdigest(),
                          url='https://zenodo.org/records/20547228/files/MotifAE_k3_step_80000.pt?download=1')
        self.registry = {'external_resources': {'motifae': {'files': [self.entry]}}}

    def test_official_resource_is_pinned(self):
        registry = json.loads((Path(download_motifae.__file__).parent / 'configs/resources.json').read_text())
        resource = registry['external_resources']['motifae']
        self.assertEqual(resource['input_dim'], 1280)
        self.assertEqual(resource['latent_dim'], 40960)
        entry, = resource['files']
        self.assertEqual(entry['filename'], 'MotifAE_k3_step_80000.pt')
        self.assertEqual(entry['bytes'], 419601512)
        self.assertEqual(entry['sha256'], 'e322c859a593c6c3dba461a2fae6c648ccd064a46422221f05f25c6eac03f28d')
        self.assertEqual(entry['url'], self.entry['url'])

    def test_default_directory(self):
        with patch.dict(os.environ, {}, clear=True):
            expected = Path(download_motifae.__file__).resolve().parents[1] / 'data/motifae'
            self.assertEqual(download_motifae.default_directory(), expected)
        with patch.dict(os.environ, {'PRIORSLIM_DATA_ROOT': '/tmp/prior-data'}, clear=True):
            self.assertEqual(download_motifae.default_directory(), Path('/tmp/prior-data/motifae'))

    def test_existing_checkpoint_is_verified_without_network(self):
        path = self.directory / self.entry['filename']
        path.write_bytes(self.content)
        with patch.object(subprocess, 'run') as run:
            self.assertEqual(download_motifae.fetch(self.entry, self.directory), path)
        run.assert_not_called()

    def test_resume_keeps_original_weights(self):
        target = self.directory / self.entry['filename']
        partial = target.with_name(target.name + '.part')
        partial.write_bytes(self.content[:5])

        def finish(command, check):
            self.assertTrue(check)
            self.assertEqual(command[-1], self.entry['url'])
            self.assertEqual(command[command.index('--continue-at') + 1], '-')
            self.assertEqual(partial.read_bytes(), self.content[:5])
            partial.write_bytes(self.content)

        with patch.object(download_esm2.shutil, 'which', return_value='/usr/bin/curl'), \
                patch.object(subprocess, 'run', side_effect=finish):
            self.assertEqual(download_motifae.fetch(self.entry, self.directory), target)
        self.assertEqual(target.read_bytes(), self.content)
        self.assertFalse(partial.exists())

    def test_check_only_valid(self):
        (self.directory / self.entry['filename']).write_bytes(self.content)
        with patch.object(download_motifae.json, 'loads', return_value=self.registry), \
                patch.object(download_motifae, 'fetch') as fetch:
            download_motifae.main(['--output-dir', str(self.directory), '--check-only'])
        fetch.assert_not_called()

    def test_check_only_missing_does_not_download(self):
        with patch.object(download_motifae.json, 'loads', return_value=self.registry), \
                patch.object(download_motifae, 'fetch') as fetch, self.assertRaises(SystemExit) as error:
            download_motifae.main(['--output-dir', str(self.directory), '--check-only'])
        self.assertEqual(error.exception.code, 1)
        fetch.assert_not_called()

    def test_download_failure_reports_resume(self):
        with patch.object(download_motifae.json, 'loads', return_value=self.registry), \
                patch.object(download_motifae, 'fetch', side_effect=subprocess.CalledProcessError(18, 'curl')), \
                self.assertRaises(SystemExit) as error:
            download_motifae.main(['--output-dir', str(self.directory)])
        self.assertEqual(error.exception.code, 1)


if __name__ == '__main__':
    unittest.main()
