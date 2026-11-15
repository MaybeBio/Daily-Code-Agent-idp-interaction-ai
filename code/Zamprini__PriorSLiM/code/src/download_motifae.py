"""Download the official MotifAE k3 checkpoint with resume and SHA-256 checks."""
import argparse
import json
import os
from pathlib import Path
import subprocess

if __package__:
    from .download_esm2 import fetch, verify
else:
    from download_esm2 import fetch, verify


def default_directory():
    root = os.environ.get('PRIORSLIM_DATA_ROOT', Path(__file__).resolve().parents[1] / 'data')
    return Path(root).expanduser() / 'motifae'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=default_directory())
    parser.add_argument('--check-only', action='store_true', help='Verify an existing checkpoint without downloading')
    args = parser.parse_args(argv)
    registry = json.loads((Path(__file__).parent / 'configs/resources.json').read_text())
    entry = registry['external_resources']['motifae']['files'][0]
    directory = args.output_dir.expanduser().resolve()
    try:
        if args.check_only:
            path = directory / entry['filename']
            if not path.is_file() or not verify(path, entry):
                raise ValueError(f'Missing checkpoint or checksum mismatch: {path}')
        else:
            path = fetch(entry, directory)
        print(f'Verified MotifAE k3 checkpoint: {path}', flush=True)
    except subprocess.CalledProcessError:
        parser.exit(1, 'Download interrupted. Rerun this command to resume the .part file.\n')
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()
