"""Download a versioned asset and verify its SHA-256 before installation."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.parse import urlparse
from urllib.error import HTTPError
from urllib.request import urlopen

from restore_archives import fetch_split


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(entry, destination, base_url=None):
    destination = Path(destination)
    if destination.exists():
        if sha256(destination) != entry['sha256']:
            raise ValueError(f'Existing asset checksum mismatch: {destination}')
        return destination
    url = entry.get('url')
    if base_url:
        url = base_url.rstrip('/') + '/' + entry['filename']
    if not url:
        raise ValueError('This asset has no published URL. Supply the verified release base URL with --base-url.')
    if urlparse(url).scheme != 'https':
        raise ValueError('Asset URLs must use HTTPS')
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as output:
            temporary = Path(output.name)
            with urlopen(url, timeout=60) as response:
                for chunk in iter(lambda: response.read(1024 * 1024), b''):
                    output.write(chunk)
        if sha256(temporary) != entry['sha256']:
            raise ValueError('Downloaded asset checksum mismatch')
        if temporary.stat().st_size != entry['bytes']:
            raise ValueError('Downloaded asset size mismatch')
        temporary.replace(destination)
    except HTTPError as error:
        if error.code != 404 or not entry['filename'].endswith('.tar.gz'):
            raise
        return fetch_split(entry, destination, url.rsplit('/', 1)[0])
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', help='Default: the version in the resource manifest')
    parser.add_argument('--asset', default='model',
                        help='model: prediction checkpoint (default); models: experiment and GussIDR models; '
                             'feature_inputs: experiment feature-reconstruction inputs')
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'data/models')
    parser.add_argument('--base-url', help='Optional HTTPS mirror; defaults to the published Zenodo URLs')
    args = parser.parse_args()
    manifest = json.loads((Path(__file__).parent / 'configs/resources.json').read_text())
    if args.version is not None and args.version != manifest['version']:
        parser.error('Unknown asset version: ' + args.version)
    if args.asset not in manifest['assets']:
        if args.asset in manifest.get('experiment_models', {}):
            parser.error('Experiment checkpoints are in models.tar.gz. Use --asset models and follow reproduce/README.md.')
        parser.error('Unknown asset; choose: ' + ', '.join(manifest['assets']))
    entry = manifest['assets'][args.asset]
    print(fetch(entry, args.output_dir / entry['filename'], args.base_url))


if __name__ == '__main__':
    main()
