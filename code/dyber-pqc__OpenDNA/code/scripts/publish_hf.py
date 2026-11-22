#!/usr/bin/env python3
"""Upload the trained heads and dataset cards to Hugging Face, and optionally
make the repositories public.

Run it yourself. It reads whatever login `hf auth login` stored and never
asks for, prints or writes a token.

Nothing happens without a flag. By default it checks and prints what it would
do, which is the mode to run first:

    python scripts/publish_hf.py                  # check and show the plan
    python scripts/publish_hf.py --upload         # send the files
    python scripts/publish_hf.py --public         # make the repositories public

The order matters. The cards in those repositories are whatever was uploaded
last; if they predate the current ones, making the repositories public
publishes superseded documents. Upload first, then flip.

Before uploading anything it verifies each head's weights against the SHA-256
recorded in its own metadata, which is the same check `opendna models` makes.
A file that does not match is not a file worth publishing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ORG = "dyber-pqc"

# The four trained heads. Each repository holds the weights, the metadata, the
# metrics and the card, and the card has to land as README.md at the root for
# Hugging Face to render it.
HEADS = ["disorder", "burial", "stability", "solubility"]

# The datasets, including `structures`, which trains no head of its own but is
# where burial's splits come from.
DATASETS = ["disorder", "burial", "stability", "solubility", "structures"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def head_files(task: str) -> list[tuple[Path, str]]:
    """Local path paired with the name it takes in the repository."""
    models = ROOT / "models"
    return [
        (models / f"{task}.README.md", "README.md"),
        (models / f"{task}.safetensors", f"{task}.safetensors"),
        (models / f"{task}.metadata.json", f"{task}.metadata.json"),
        (models / f"{task}.metrics.json", f"{task}.metrics.json"),
    ]


def dataset_files(task: str, with_data: bool) -> list[tuple[Path, str]]:
    cards = ROOT / "datasets" / "cards" / task
    files = [(cards / "README.md", "README.md"), (cards / "manifest.json", "manifest.json")]
    if with_data:
        # Parquet is not in git; it is in the builder's output directory when
        # the dataset has been built on this machine.
        out = ROOT / "datasets" / "out" / task / "data"
        for split in ("train", "validation", "test"):
            parquet = out / f"{split}.parquet"
            if parquet.is_file():
                files.append((parquet, f"data/{split}.parquet"))
    return files


def check_head(task: str) -> list[str]:
    """Problems with one head, as sentences. Empty means it is fit to publish."""
    problems = []
    for path, _ in head_files(task):
        if not path.is_file():
            problems.append(f"missing {path.relative_to(ROOT)}")
    if problems:
        return problems

    metadata = json.loads((ROOT / "models" / f"{task}.metadata.json").read_text(encoding="utf-8"))
    recorded = metadata.get("weights_sha256")
    if not recorded:
        problems.append("metadata records no weights_sha256, so the weights cannot be verified")
    else:
        actual = sha256(ROOT / "models" / f"{task}.safetensors")
        if actual != recorded:
            problems.append(
                f"weights hash {actual[:16]} but the metadata records {recorded[:16]}; "
                "do not publish this"
            )

    # The card has to open with front matter or the Hub renders none of its
    # metadata, licence included.
    card = (ROOT / "models" / f"{task}.README.md").read_text(encoding="utf-8")
    if not card.startswith("---\n"):
        problems.append("the card has no YAML front matter")
    elif "\n---\n" not in card[4:]:
        problems.append("the card's front matter is not closed")
    return problems


def describe(task: str, kind: str, files: list[tuple[Path, str]], problems: list[str]) -> None:
    repo = f"{ORG}/opendna-{task}"
    print(f"\n{kind} {repo}")
    for path, name in files:
        size = path.stat().st_size if path.is_file() else 0
        mark = " " if path.is_file() else "!"
        print(f"  {mark} {size / 1e6:8.2f} MB  {path.relative_to(ROOT)}  ->  {name}")
    for problem in problems:
        print(f"  ! {problem}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--upload", action="store_true", help="send the files")
    parser.add_argument("--public", action="store_true", help="make the repositories public")
    parser.add_argument(
        "--with-data",
        action="store_true",
        help="also send each dataset's Parquet files, if they have been built here",
    )
    parser.add_argument("--only", nargs="*", help="limit to these task names")
    args = parser.parse_args()

    wanted = set(args.only) if args.only else None
    heads = [t for t in HEADS if not wanted or t in wanted]
    datasets = [t for t in DATASETS if not wanted or t in wanted]

    print(f"Repository root: {ROOT}")
    blocked = False
    plan: list[tuple[str, str, list[tuple[Path, str]]]] = []

    for task in heads:
        problems = check_head(task)
        files = head_files(task)
        describe(task, "model  ", files, problems)
        if problems:
            blocked = True
        else:
            plan.append(("model", task, files))

    for task in datasets:
        files = dataset_files(task, args.with_data)
        problems = [f"missing {p.relative_to(ROOT)}" for p, _ in files if not p.is_file()]
        describe(task, "dataset", files, problems)
        if problems:
            blocked = True
        else:
            plan.append(("dataset", task, files))
            if not args.with_data:
                print("    note: the Parquet files are not in git; --with-data sends them if built here")

    if blocked:
        print("\nSomething above is missing or does not verify. Nothing has been sent.")
        return 1
    if not (args.upload or args.public):
        print("\nChecks pass. Re-run with --upload to send these, then --public to open them.")
        return 0

    from huggingface_hub import HfApi  # imported here so the checks run without it

    api = HfApi()
    try:
        who = api.whoami()
    except Exception as error:  # noqa: BLE001 - the message is what matters
        print(f"\nNot logged in ({type(error).__name__}). Run `hf auth login` first.")
        return 1
    print(f"\nLogged in as {who.get('name', '?')}.")

    if args.upload:
        for kind, task, files in plan:
            repo = f"{ORG}/opendna-{task}"
            for path, name in files:
                print(f"  {repo}  <-  {name}")
                api.upload_file(
                    path_or_fileobj=str(path),
                    path_in_repo=name,
                    repo_id=repo,
                    repo_type="dataset" if kind == "dataset" else "model",
                    commit_message=f"Update {name} from the OpenDNA repository",
                )
        print("Uploaded.")

    if args.public:
        print("\nMaking these public. This cannot be quietly undone: anything fetched or")
        print("indexed while they are open stays out there.")
        print("opendna-solubility carries CC BY-SA 2.1 JP, so its share-alike obligation")
        print("starts travelling with the weights from this point.")
        if input("Type 'public' to continue: ").strip() != "public":
            print("Left private.")
            return 0
        for kind, task, _ in plan:
            repo = f"{ORG}/opendna-{task}"
            api.update_repo_settings(
                repo_id=repo,
                private=False,
                repo_type="dataset" if kind == "dataset" else "model",
            )
            print(f"  public: {repo}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
