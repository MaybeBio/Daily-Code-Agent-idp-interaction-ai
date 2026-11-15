"""Reassemble split release archives and verify every part and the complete file."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.parse import quote, urlsplit
from urllib.request import urlopen

MANIFEST = 'archive_parts.json'
FORMAT = 'priorslim_split_archives_v1'


def local_file(directory, name):
    directory = Path(directory).resolve()
    if not name or Path(name).name != name or name in {'.', '..'} or '\\' in name:
        raise ValueError('Invalid release filename')
    path = directory / name
    if path.resolve().parent != directory:
        raise ValueError('Release path leaves the selected directory')
    return path


def sha256(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            result.update(block)
    return result.hexdigest()


def verify(path, entry):
    if path.stat().st_size != entry['bytes'] or sha256(path) != entry['sha256']:
        raise ValueError(f'Checksum or size mismatch: {path}')


def asset_url(base_url, filename):
    parsed = urlsplit(base_url)
    if parsed.scheme != 'https' or not parsed.netloc or parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ValueError('Supply an HTTPS release directory without credentials or query parameters')
    return base_url.rstrip('/') + '/' + quote(filename, safe='')


def load_manifest(directory, base_url=None):
    path = local_file(directory, MANIFEST)
    if path.exists():
        document = json.loads(path.read_text())
    elif base_url:
        with urlopen(asset_url(base_url, MANIFEST), timeout=60) as response:
            content = response.read(4_000_001)
        if len(content) > 4_000_000:
            raise ValueError('Archive manifest is too large')
        document = json.loads(content)
    else:
        raise FileNotFoundError(f'Download {MANIFEST} into {directory}')
    if document.get('format') != FORMAT:
        raise ValueError('Unsupported split-archive manifest')
    return document


def obtain_part(directory, entry, base_url):
    path = local_file(directory, entry['filename'])
    if path.exists():
        verify(path, entry)
        return path
    if not base_url:
        raise FileNotFoundError(f'Missing archive part: {path.name}')
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=directory, suffix='.part', delete=False) as handle:
            temporary = Path(handle.name)
            with urlopen(asset_url(base_url, path.name), timeout=60) as response:
                for block in iter(lambda: response.read(1024**2), b''):
                    handle.write(block)
        verify(temporary, entry)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def restore_archive(directory, entry, base_url=None):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    target = local_file(directory, entry['filename'])
    if target.exists():
        verify(target, entry)
        return target
    parts = entry['parts']
    if not parts or sum(p['bytes'] for p in parts) != entry['bytes']:
        raise ValueError('Archive part sizes do not sum to the complete size')
    for number, part in enumerate(parts, 1):
        if part['filename'] != f'{target.name}.part{number:04d}' or part['bytes'] <= 0:
            raise ValueError('Invalid archive part order or size')
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=directory, suffix='.assembling', delete=False) as output:
            temporary = Path(output.name)
            for number, part in enumerate(parts, 1):
                source = obtain_part(directory, part, base_url)
                with source.open('rb') as handle:
                    for block in iter(lambda: handle.read(8 * 1024**2), b''):
                        output.write(block)
                print(f'{target.name}: verified part {number}/{len(parts)}', flush=True)
        verify(temporary, entry)
        temporary.replace(target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return target


def fetch_split(entry, destination, base_url):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    document = load_manifest(destination.parent, base_url)
    candidates = [a for a in document['archives'] if a['filename'] == entry['filename']]
    if len(candidates) != 1:
        raise ValueError('Requested asset is absent or duplicated in the split manifest')
    archive = candidates[0]
    if any(archive[k] != entry[k] for k in ['bytes', 'sha256']):
        raise ValueError('Split archive differs from the registered asset')
    result = restore_archive(destination.parent, archive, base_url)
    if result != destination.resolve():
        result.replace(destination)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('.'))
    parser.add_argument('--archive', action='append', help='Restore only this archive; may be repeated')
    parser.add_argument('--base-url', help='Optional HTTPS release directory for downloading missing parts')
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    document = load_manifest(args.directory, args.base_url)
    entries = document['archives']
    if args.archive:
        entries = [a for a in entries if a['filename'] in args.archive]
        if {a['filename'] for a in entries} != set(args.archive):
            parser.error('Unknown archive name')
    for entry in entries:
        print(f'Verified archive: {restore_archive(args.directory, entry, args.base_url)}')
    print('Archive parts are retained; they may be removed after successful verification.')


if __name__ == '__main__':
    main()
