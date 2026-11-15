"""Download official ESM2 weights into the Torch Hub cache, with resume and SHA-256 checks."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

if __package__:
    from .download import sha256
else:
    from download import sha256


def cache_directory():
    torch_home = os.environ.get('TORCH_HOME')
    if torch_home is None:
        torch_home = os.path.join(os.environ.get('XDG_CACHE_HOME', '~/.cache'), 'torch')
    return Path(torch_home).expanduser() / 'hub' / 'checkpoints'


def verify(path, entry):
    return path.stat().st_size == entry['bytes'] and sha256(path) == entry['sha256']


def fetch(entry, cache_dir):
    destination = Path(cache_dir) / entry['filename']
    if destination.exists():
        if not verify(destination, entry):
            raise ValueError(f'Checksum mismatch: {destination}. Move this file aside and rerun.')
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    if partial.exists() and partial.stat().st_size >= entry['bytes']:
        if not verify(partial, entry):
            raise ValueError(f'Invalid completed download: {partial}. Remove it and rerun.')
    else:
        if shutil.which('curl') is None:
            raise RuntimeError('curl is required for resumable downloads. Install curl and rerun.')
        subprocess.run([
            'curl', '--fail', '--location', '--retry', '5', '--connect-timeout', '30',
            '--speed-limit', '1024', '--speed-time', '120', '--continue-at', '-',
            '--output', str(partial), entry['url'],
        ], check=True)
        if not verify(partial, entry):
            raise ValueError(f'Download checksum mismatch: {partial}. Remove it and rerun.')
    partial.replace(destination)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    registry = json.loads((Path(__file__).parent / 'configs/resources.json').read_text())
    cache_dir = cache_directory()
    print(f'ESM2 cache: {cache_dir}', flush=True)
    try:
        for entry in registry['external_resources']['esm2']['files']:
            path = fetch(entry, cache_dir)
            print(f'Verified: {path}', flush=True)
    except subprocess.CalledProcessError:
        parser.exit(1, 'Download interrupted. Rerun this command to resume the .part file.\n')
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()
