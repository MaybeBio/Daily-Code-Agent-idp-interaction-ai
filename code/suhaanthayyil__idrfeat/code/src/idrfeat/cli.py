"""Command line interface for idrfeat."""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .annotations import (
    load_dbptm,
    load_drllps,
    load_elm,
    load_phasepro,
    load_plddt_dir,
)
from .complexity import seg_masks, segmasker_available
from .features import Annotations, build_feature_table
from .io import load_config, read_fasta, write_table


def _build_annotations(args: argparse.Namespace, seqs: dict[str, str]) -> Annotations:
    accessions = list(seqs)
    return Annotations(
        elm=load_elm(args.elm) if args.elm else None,
        ptm=load_dbptm(args.dbptm) if args.dbptm else None,
        phasepro=load_phasepro(args.phasepro) if args.phasepro else None,
        drllps=load_drllps(args.drllps) if args.drllps else None,
        plddt=load_plddt_dir(args.plddt_dir, accessions) if args.plddt_dir else None,
        lowcomplexity=seg_masks(seqs) if segmasker_available() else None,
    )


def _cmd_features(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    if args.backend:
        cfg["disorder"]["primary"] = args.backend
    seqs = read_fasta(args.fasta)
    if not seqs:
        print(f"no sequences found in {args.fasta}", file=sys.stderr)
        return 1
    annotations = _build_annotations(args, seqs)
    df = build_feature_table(seqs, cfg, annotations=annotations)
    parquet_path, csv_path = write_table(df, args.out)
    active = [
        name
        for name, value in [
            ("elm", annotations.elm),
            ("dbptm", annotations.ptm),
            ("phasepro", annotations.phasepro),
            ("drllps", annotations.drllps),
            ("plddt", annotations.plddt),
            ("lowcomplexity(SEG)", annotations.lowcomplexity),
        ]
        if value
    ]
    backend = df["disorder_backend"].iloc[0] if len(df) else "n/a"
    print(
        f"{len(df)} IDR segments from {len(seqs)} proteins "
        f"(disorder backend: {backend}; annotations: {', '.join(active) or 'none'})"
    )
    print(f"wrote {parquet_path} and {csv_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="idrfeat", description=__doc__)
    parser.add_argument("--version", action="version", version=f"idrfeat {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    feats = sub.add_parser("features", help="annotate IDRs and compute per-segment features")
    feats.add_argument("--fasta", required=True, help="protein FASTA (plain or .gz)")
    feats.add_argument("--out", required=True, help="output path; parquet plus CSV alongside")
    feats.add_argument("--config", help="optional YAML config overriding the defaults")
    feats.add_argument(
        "--backend",
        choices=["metapredict", "aiupred"],
        help="disorder backend that defines IDR segments (default: metapredict)",
    )
    feats.add_argument("--elm", help="ELM instances TSV")
    feats.add_argument("--dbptm", nargs="+", help="one or more dbPTM files")
    feats.add_argument("--phasepro", help="PhaSePro full JSON")
    feats.add_argument("--drllps", help="DrLLPS LLPS protein table")
    feats.add_argument("--plddt-dir", dest="plddt_dir", help="directory of {accession}.pdb models")
    feats.set_defaults(func=_cmd_features)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv if argv is not None else sys.argv[1:])
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
