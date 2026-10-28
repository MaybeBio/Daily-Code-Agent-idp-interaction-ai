from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

from idrfeat.annotations import (
    load_dbptm,
    load_drllps,
    load_elm,
    load_phasepro,
)
from idrfeat.complexity import seg_masks, segmasker_available
from idrfeat.features import Annotations, build_feature_table
from idrfeat.io import load_config, read_fasta, write_table

STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true&format=fasta&query=%28proteome%3A{proteome}%29"
)


def _download(proteome: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = STREAM_URL.format(proteome=proteome)
    req = urllib.request.Request(url, headers={"User-Agent": "idrfeat-build"})
    with urllib.request.urlopen(req, timeout=600) as resp, open(dest, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    return dest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="build_proteome")
    parser.add_argument("--fasta", help="pre-downloaded FASTA (plain or .gz); skips the download")
    parser.add_argument("--proteome", default="UP000005640", help="UniProt proteome id")
    parser.add_argument("--out", default="data/processed/proteome_idr_features.parquet")
    parser.add_argument("--config")
    parser.add_argument("--backend", choices=["metapredict", "aiupred"])
    parser.add_argument("--limit", type=int, help="cap number of proteins for a smoke run")
    parser.add_argument("--elm", help="ELM instances TSV")
    parser.add_argument("--dbptm", nargs="+", help="one or more dbPTM files")
    parser.add_argument("--phasepro", help="PhaSePro full JSON")
    parser.add_argument("--drllps", help="DrLLPS LLPS protein table")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    if args.fasta:
        fasta = Path(args.fasta)
    else:
        fasta = Path("data/raw") / f"{args.proteome}.fasta.gz"
        if not fasta.exists():
            print(f"downloading {args.proteome} to {fasta}")
            _download(args.proteome, fasta)

    seqs = read_fasta(fasta)
    if args.limit:
        seqs = {k: seqs[k] for k in list(seqs)[: args.limit]}
    if not seqs:
        print(f"no sequences in {fasta}", file=sys.stderr)
        return 1

    cfg = load_config(args.config)
    if args.backend:
        cfg["disorder"]["primary"] = args.backend
    annotations = Annotations(
        elm=load_elm(args.elm) if args.elm else None,
        ptm=load_dbptm(args.dbptm) if args.dbptm else None,
        phasepro=load_phasepro(args.phasepro) if args.phasepro else None,
        drllps=load_drllps(args.drllps) if args.drllps else None,
        lowcomplexity=seg_masks(seqs) if segmasker_available() else None,
    )
    df = build_feature_table(seqs, cfg, annotations=annotations)
    parquet_path, csv_path = write_table(df, args.out)
    print(f"{len(df)} IDR segments from {len(seqs)} proteins")
    print(f"wrote {parquet_path} and {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
