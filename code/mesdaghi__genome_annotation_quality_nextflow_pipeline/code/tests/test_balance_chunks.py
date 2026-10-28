import subprocess
from pathlib import Path

BALANCE_CHUNKS = (
    Path(__file__).resolve().parent.parent
    / "genome_annotation_quality_nextflow_pipeline"
    / "bin"
    / "balance_chunks.sh"
)


def run_balance(fasta, chunks, outdir):
    return subprocess.run(
        ["bash", str(BALANCE_CHUNKS), "--input", str(fasta), "--chunks", str(chunks), "--outdir", str(outdir)],
        capture_output=True,
        text=True,
    )


def test_splits_multi_record_file_preserving_headers_and_sequences(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seq1 desc\nACGT\nACGT\n>seq2\nAC\n")
    outdir = tmp_path / "out"

    result = run_balance(fasta, 2, outdir)

    assert result.returncode == 0, result.stderr
    chunk_files = sorted(outdir.glob("chunk_*.fa"))
    assert len(chunk_files) == 2
    combined = "".join(f.read_text() for f in chunk_files)
    assert ">seq1 desc\nACGT\nACGT\n" in combined
    assert ">seq2\nAC\n" in combined


def test_empty_file_errors(tmp_path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")
    outdir = tmp_path / "out"

    result = run_balance(fasta, 2, outdir)

    assert result.returncode != 0
    assert "no sequences found" in result.stderr


def test_balances_chunks_by_residue_count(tmp_path):
    fasta = tmp_path / "in.fa"
    records = [(">long", "A" * 100)] + [(f">short{i}", "A" * 10) for i in range(10)]
    fasta.write_text("".join(f"{h}\n{s}\n" for h, s in records))
    outdir = tmp_path / "out"

    result = run_balance(fasta, 2, outdir)

    assert result.returncode == 0, result.stderr
    chunk_files = sorted(outdir.glob("chunk_*.fa"))
    assert len(chunk_files) == 2

    totals = []
    for chunk_file in chunk_files:
        lines = chunk_file.read_text().splitlines()
        totals.append(sum(len(line.strip()) for line in lines if not line.startswith(">")))

    assert sum(totals) == 200
    assert abs(totals[0] - totals[1]) <= 10


def test_chunks_capped_at_record_count(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seq1\nAC\n>seq2\nAC\n")
    outdir = tmp_path / "out"

    result = run_balance(fasta, 10, outdir)

    assert result.returncode == 0, result.stderr
    chunk_files = sorted(outdir.glob("chunk_*.fa"))
    assert len(chunk_files) == 2
