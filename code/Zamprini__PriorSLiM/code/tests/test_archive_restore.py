import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import download
import restore_archives as restore


class ArchiveRestoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.content = b'first chunk\x00second chunk\xffthird'
        self.entry = {'filename': 'models.tar.gz', 'bytes': len(self.content),
                      'sha256': hashlib.sha256(self.content).hexdigest(), 'parts': []}
        self.parts = {}
        for number, offset in enumerate(range(0, len(self.content), 10), 1):
            name = f'models.tar.gz.part{number:04d}'
            data = self.content[offset:offset + 10]
            (self.root / name).write_bytes(data)
            self.parts[name] = data
            self.entry['parts'].append({'filename': name, 'bytes': len(data),
                                        'sha256': hashlib.sha256(data).hexdigest()})
        self.document = {'format': restore.FORMAT, 'archives': [self.entry]}
        (self.root / restore.MANIFEST).write_text(json.dumps(self.document))

    def test_exact_restore(self):
        path = restore.restore_archive(self.root, self.entry)
        self.assertEqual(path.read_bytes(), self.content)
        self.assertFalse(list(self.root.glob('*.assembling')))

    def test_existing_valid_archive_is_reused(self):
        path = restore.restore_archive(self.root, self.entry)
        before = path.stat().st_mtime_ns
        restore.restore_archive(self.root, self.entry)
        self.assertEqual(before, path.stat().st_mtime_ns)

    def test_existing_invalid_archive_is_not_overwritten(self):
        target = self.root / self.entry['filename']
        target.write_bytes(b'invalid')
        with self.assertRaises(ValueError):
            restore.restore_archive(self.root, self.entry)
        self.assertEqual(target.read_bytes(), b'invalid')

    def test_corrupt_part_is_rejected(self):
        (self.root / self.entry['parts'][0]['filename']).write_bytes(b'invalid')
        with self.assertRaises(ValueError):
            restore.restore_archive(self.root, self.entry)
        self.assertFalse((self.root / self.entry['filename']).exists())
        self.assertFalse(list(self.root.glob('*.assembling')))

    def test_missing_part_is_reported(self):
        (self.root / self.entry['parts'][0]['filename']).unlink()
        with self.assertRaises(FileNotFoundError):
            restore.restore_archive(self.root, self.entry)

    def test_wrong_final_hash_is_rejected(self):
        entry = dict(self.entry, sha256='0' * 64)
        with self.assertRaises(ValueError):
            restore.restore_archive(self.root, entry)
        self.assertFalse(list(self.root.glob('*.assembling')))

    def test_part_order_is_checked(self):
        entry = dict(self.entry, parts=list(reversed(self.entry['parts'])))
        with self.assertRaises(ValueError):
            restore.restore_archive(self.root, entry)

    def test_path_escape_is_rejected(self):
        with self.assertRaises(ValueError):
            restore.restore_archive(self.root, dict(self.entry, filename='../outside'))

    def test_downloaded_parts_are_verified(self):
        for name in self.parts:
            (self.root / name).unlink()
        with patch.object(restore, 'urlopen', side_effect=lambda url, **kw: io.BytesIO(self.parts[url.rsplit('/', 1)[1]])):
            path = restore.restore_archive(self.root, self.entry, 'https://example.org/files')
        self.assertEqual(path.read_bytes(), self.content)

    def test_registered_asset_identity_is_required(self):
        with self.assertRaises(ValueError):
            restore.fetch_split(dict(self.entry, sha256='0' * 64), self.root / 'models.tar.gz', 'https://example.org/files')

    def test_downloader_falls_back_on_404(self):
        url = 'https://example.org/files/models.tar.gz'
        error = HTTPError(url, 404, 'Not Found', {}, None)
        with patch.object(download, 'urlopen', side_effect=error):
            path = download.fetch(self.entry, self.root / 'models.tar.gz', 'https://example.org/files')
        self.assertEqual(path.read_bytes(), self.content)

    def test_permission_errors_do_not_trigger_fallback(self):
        error = HTTPError('https://example.org/files/models.tar.gz', 403, 'Forbidden', {}, None)
        with patch.object(download, 'urlopen', side_effect=error), patch.object(download, 'fetch_split') as fallback:
            with self.assertRaises(HTTPError):
                download.fetch(self.entry, self.root / 'models.tar.gz', 'https://example.org/files')
        fallback.assert_not_called()


if __name__ == '__main__':
    unittest.main()
