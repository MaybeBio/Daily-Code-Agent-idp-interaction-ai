#!/usr/bin/env python3
"""Command line for the plant PPI structural screening pipeline.

Typical use:
    python cli.py resolve
    python cli.py batch
    # upload batches/*.json to the AlphaFold Server, download the outputs
    python cli.py ingest --results results/
    python cli.py evaluate --results results/scores.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from ppi import batch, curate, metrics, report, results, uniprot

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
CANDIDATES = DATA / "candidates.csv"
SEQUENCES = DATA / "sequences"
STATE = ROOT / ".state" / "quota.json"

FIELDS = ["pair_id", "gene_a", "gene_b", "acc_a", "acc_b", "range_a", "range_b",
          "label", "label_type", "source", "verified"]


def _read_csv(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict], fields: list[str] = FIELDS) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def cmd_resolve(args: argparse.Namespace) -> int:
    rows = _read_csv(Path(args.candidates))
    failures = 0
    for row in rows:
        for side in ("a", "b"):
            acc_key, gene_key = f"acc_{side}", f"gene_{side}"
            if row.get(acc_key):
                uniprot.fetch_fasta(row[acc_key], SEQUENCES)
                continue
            try:
                entry = uniprot.resolve(row[gene_key])
            except uniprot.ResolutionError as exc:
                print(f"  ! {row['pair_id']}: {exc}", file=sys.stderr)
                failures += 1
                continue
            row[acc_key] = entry.accession
            uniprot.fetch_fasta(entry.accession, SEQUENCES)
            print(f"  {row[gene_key]:>10} -> {entry.accession} ({entry.length} aa)")
    out = Path(args.out)
    _write_csv(out, rows)
    print(f"resolved rows written to {out}")
    return 1 if failures else 0


def cmd_batch(args: argparse.Namespace) -> int:
    rows = _read_csv(Path(args.candidates))
    usable = [r for r in rows if r.get("acc_a") and r.get("acc_b")]
    if len(usable) != len(rows):
        print(f"  note: {len(rows) - len(usable)} rows have unresolved accessions and are skipped")

    sequences = {}
    for row in usable:
        for side in ("a", "b"):
            acc = row[f"acc_{side}"]
            path = SEQUENCES / f"{acc}.fasta"
            if path.exists():
                sequences[acc] = uniprot.read_fasta(path)
    plan_obj = batch.plan(usable, sequences)

    print(f"planned {len(plan_obj.jobs)} jobs, {plan_obj.total_tokens} tokens total")
    for pair_id, reason in plan_obj.skipped:
        print(f"  ! skipped {pair_id}: {reason}", file=sys.stderr)

    remaining = batch.QuotaTracker(STATE).remaining()
    if len(plan_obj.jobs) > remaining:
        print(f"  note: {len(plan_obj.jobs)} jobs planned, {remaining} left in today's "
              f"quota of {batch.DAILY_JOB_QUOTA} - submit over several days")

    written = batch.write_batches(plan_obj, ROOT / "batches")
    for path in written:
        print(f"  wrote {path.relative_to(ROOT)}")
    if not written:
        print("nothing to write")
        return 1
    print("upload these files at https://alphafoldserver.com and download the outputs")
    return 0


def cmd_ingest(args: argparse.Namespace) -> int:
    scores = results.collect(Path(args.results), Path(args.out))
    print(f"parsed {len(scores)} summary files -> {args.out}")
    missing = [s.pair_id for s in scores if s.iptm != s.iptm]
    if missing:
        print(f"  ! no ipTM found for: {', '.join(missing)}", file=sys.stderr)
    return 0


def cmd_evaluate(args: argparse.Namespace) -> int:
    rows = _read_csv(Path(args.candidates))
    if not Path(args.results).exists():
        print(f"no scores file at {args.results}", file=sys.stderr)
        print("run `python cli.py ingest --results results/` after downloading the jobs, "
              "or `python cli.py demo` for an offline run on synthetic scores", file=sys.stderr)
        return 1
    scores = {r["pair_id"]: r for r in results.read_scores(Path(args.results))}

    joined, unverified, calibration = [], [], []
    for row in rows:
        scored = scores.get(row["pair_id"])
        if not scored:
            continue
        try:
            score = float(scored["iptm"])
        except (TypeError, ValueError):
            continue
        if score != score:
            continue
        entry = {
            "pair_id": row["pair_id"],
            "label": int(float(row["label"])),
            "score": score,
            "min_chain_pair_pae": _float_or_none(scored.get("min_chain_pair_pae")),
        }
        # Calibration pairs validate the pipeline; they must never enter the metrics.
        if (row.get("label_type") or "").strip() == curate.CALIBRATION:
            calibration.append(entry)
            continue
        if row.get("verified", "").strip().lower() != "yes":
            unverified.append(row["pair_id"])
        joined.append(entry)

    if unverified and not args.allow_unverified:
        print("refusing to score: these pairs are marked unverified in candidates.csv:", file=sys.stderr)
        for pair_id in unverified:
            print(f"  - {pair_id}", file=sys.stderr)
        print("check each label against its source, set verified=yes, or pass "
              "--allow-unverified for a dry run", file=sys.stderr)
        return 2

    if not joined and calibration:
        _report_calibration(calibration)
        print("\ncalibration only: no verified benchmark pairs yet, so no metrics were computed.")
        print("add labelled pairs in candidates.csv, or verify the existing ones, then rerun.")
        return 0

    if not joined:
        print("no overlapping pairs between candidates.csv and the results file", file=sys.stderr)
        return 1

    labels = [row["label"] for row in joined]
    values = [row["score"] for row in joined]
    summary = metrics.summarise(values, labels, fraction=args.fraction, threshold=args.threshold)
    note = args.note
    if args.allow_unverified:
        note = (note + " " if note else "") + \
            "Includes unverified labels; treat as a dry run, not a result."
    paths = report.write_report(joined, summary, Path(args.out), args.threshold, note)

    if calibration:
        _report_calibration(calibration)

    print(f"{len(joined)} benchmark pairs ({sum(labels)} positive)")
    for name in ("auc", "recall_top", "precision_at_threshold"):
        if name in summary:
            print(f"  {name}: {summary[name]}")
    print(f"  wrote {paths[0]} and {paths[1]}")
    return 0


def cmd_check_labels(args: argparse.Namespace) -> int:
    """Validate the candidate list before any predictions are scored."""
    path = Path(args.candidates)
    if not path.exists():
        fallback = CANDIDATES
        print(f"  note: {path} not found, checking {fallback.name} instead")
        path = fallback
    rows = _read_csv(path)
    counts = curate.summary(rows)
    print(f"{path.name}: {counts['pairs']} pairs | {counts['positives']} positive | "
          f"{counts['negatives']} negative ({counts['decoys']} decoys, "
          f"{counts['curated_negatives']} curated) | {counts['verified']} verified | "
          f"{counts['calibration']} calibration")

    findings = curate.validate(rows, min_pairs=args.min_pairs)
    errors = [f for f in findings if f.level == "error"]
    warnings = [f for f in findings if f.level == "warning"]
    for finding in errors + warnings:
        print(f"  {finding}")

    if errors:
        print("fix the errors above before running predictions")
        return 1
    if warnings:
        print("warnings only: scoring will run, but read them before trusting the numbers")
    else:
        print("label file looks consistent")
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    """End-to-end run on synthetic scores so the pipeline can be checked offline."""
    demo_dir = ROOT / "demo"
    pairs = [
        ("BRI1_BAK1_domains", 1, 0.91), ("VIP1_BZIP29", 1, 0.78),
        ("BZIP29_BZIP30", 1, 0.74), ("VIP1_AGB1", 1, 0.69),
        ("BZIP30_AGB1", 1, 0.66), ("VIP1_BZIP30", 1, 0.63),
        ("BRI1_VIP1", 0, 0.41), ("BAK1_VIP1", 0, 0.33),
        ("BRI1_BZIP29", 0, 0.28), ("AGB1_GAPDH", 0, 0.22),
        ("BZIP29_NPR1", 0, 0.19), ("VIP1_RCN1", 0, 0.12),
    ]
    rows = [{"pair_id": p, "gene_a": p.split("_")[0], "gene_b": p.split("_")[1],
             "acc_a": "", "acc_b": "", "range_a": "", "range_b": "",
             "label": str(label), "label_type": "demo", "source": "synthetic",
             "verified": "no"} for p, label, _ in pairs]
    _write_csv(demo_dir / "candidates.csv", rows)

    score_rows = [{"pair_id": p, "iptm": f"{s}", "ptm": f"{s - 0.05:.2f}",
                   "min_chain_pair_pae": f"{30 - 25 * s:.2f}", "ranking_score": f"{s:.2f}",
                   "fraction_disordered": "0.1", "has_clash": "False", "source_file": "demo"}
                  for p, _, s in pairs]
    _write_csv(demo_dir / "scores.csv", score_rows, list(score_rows[0].keys()))

    namespace = argparse.Namespace(candidates=str(demo_dir / "candidates.csv"),
                                   results=str(demo_dir / "scores.csv"),
                                   out=str(demo_dir / "report"), threshold=0.6,
                                   fraction=0.25, allow_unverified=True,
                                   note="Synthetic demonstration data, not real predictions.")
    return cmd_evaluate(namespace)


def _report_calibration(calibration: list[dict]) -> None:
    """Show calibration scores and check the expected full-length versus domain gap."""
    print("calibration pairs (excluded from metrics):")
    for entry in sorted(calibration, key=lambda e: -e["score"]):
        print(f"  {entry['pair_id']}: ipTM {entry['score']:.3f}")

    if len(calibration) < 2:
        print("  one calibration pair is enough to confirm the pipeline runs; two are "
              "needed to check the full-length versus domain gap")
        return

    full = next((e for e in calibration if "full" in e["pair_id"].lower()), None)
    truncated = next((e for e in calibration if "domain" in e["pair_id"].lower()), None)
    if full and truncated:
        gap = truncated["score"] - full["score"]
        print(f"  gap (domains - full-length): {gap:+.3f}")
        if gap >= 0.2:
            print("  behaves as expected: the truncated construct scores clearly higher")
        elif gap > -0.05:
            print("  inconclusive: the expected gap did not appear; check that range_a/range_b "
                  "cover the interacting domains and not the disordered flanks")
        else:
            print("  unexpected: full-length scored higher than the truncated construct; "
                  "verify the domain boundaries before trusting anything downstream")


def _float_or_none(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("resolve", help="map gene symbols to UniProt accessions")
    p.add_argument("--candidates", default=str(CANDIDATES))
    p.add_argument("--out", default=str(DATA / "resolved.csv"))
    p.set_defaults(func=cmd_resolve)

    p = sub.add_parser("batch", help="build AlphaFold Server batch files")
    p.add_argument("--candidates", default=str(DATA / "resolved.csv"))
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("ingest", help="parse downloaded AlphaFold Server outputs")
    p.add_argument("--results", default=str(ROOT / "results"))
    p.add_argument("--out", default=str(ROOT / "results" / "scores.csv"))
    p.set_defaults(func=cmd_ingest)

    p = sub.add_parser("evaluate", help="score the screening against known labels")
    p.add_argument("--candidates", default=str(DATA / "resolved.csv"))
    p.add_argument("--results", default=str(ROOT / "results" / "scores.csv"))
    p.add_argument("--out", default=str(ROOT / "report"))
    p.add_argument("--threshold", type=float, default=0.6)
    p.add_argument("--fraction", type=float, default=0.1)
    p.add_argument("--note", default=None)
    p.add_argument("--allow-unverified", action="store_true")
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("check-labels", help="validate the candidate list before scoring")
    p.add_argument("--candidates", default=str(CANDIDATES))
    p.add_argument("--min-pairs", type=int, default=15)
    p.set_defaults(func=cmd_check_labels)

    p = sub.add_parser("demo", help="offline end-to-end run on synthetic scores")
    p.set_defaults(func=cmd_demo)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
