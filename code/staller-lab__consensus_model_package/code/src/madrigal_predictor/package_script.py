import argparse
from pathlib import Path
from parrot import encode_sequence
from madrigal_predictor.run_query import run
from madrigal_predictor.process_runs import get_runs, get_runs_coords

def analyze_fasta_file(file_path: str | Path, threshold: float = 0.6):
    """
    Reads protein sequences of amino acids from FASTA File.

    Args:
        file_path: path to FASTA file (str or path)
        threshold: the threshold of likeliness for a protein sequence to be an AD

    Returns:
        titles: a list of the titles of the FASTA proteins
        above_runs: a list of lists, with inner lists having protein sequences likely to be ADs
        below_runs: a list of lists, with inner lists having protein sequences unlikely to be ADs
        above_starts: a list of lists, with inner lists having coordinates of starts of likely ADs
        above_ends: a list of lists, with inner lists having coordinates of ends of likely ADs

    Coordinates are 0-based and inclusive, so a likely AD is
    sequence[above_starts[i][j] : above_ends[i][j] + 1]. Add 1 to each start and
    end to report 1-based residue numbers.

    """
    # Read file path
    contents = Path(file_path).read_text(encoding='utf-8')
    return read_sequences(contents, threshold)


def read_sequence(data: str, threshold: float = 0.6):
    """Reads a protein sequence (string) of amino acids and analyzes what is likely an AD

    Args:
        data: protein sequence to be analyzed, represented by amino acids
        threshold: the threshold of a likeliness for a protein sequence to be an AD

    Returns:
        above_runs: list having protein sequences likely to be ADs
        below_runs: list having protein sequences unlikely to be ADs
        above_starts: list having coordinates of starts of likely ADs
        above_ends: list having coordinates of ends of likely ADs

    Coordinates are 0-based and inclusive, so a likely AD is
    sequence[above_starts[j] : above_ends[j] + 1]. Add 1 to each start and end
    to report 1-based residue numbers.

    """
    data = sanitize(data)
    query_input = encode_sequence.one_hot(data).float().unsqueeze(0)
    output = run(query_input).squeeze().tolist()
    sequence = ''.join(data.upper().split())
    above_runs, below_runs, above_starts, below_starts = get_runs_coords(sequence, output, threshold)
    # get_runs_coords returns the starts of the below-threshold stretches as its
    # fourth value. Those do not line up with the ADs: a sequence that begins
    # below threshold produces a below-run before the first AD, so the two lists
    # are offset. Each AD's end is derived from its own start and length instead.
    above_ends = [start + len(run_seq) - 1 for start, run_seq in zip(above_starts, above_runs)]
    return above_runs, below_runs, above_starts, above_ends

def sanitize(sequence):
    return ''.join(sequence.upper().split())

def read_sequences(data, threshold=0.6):
    lines = data.splitlines()

    titles = []
    above_runs = []
    below_runs = []
    above_starts = []
    above_ends = []

    seq = ''
    for line in lines:
        if line.startswith('>'):
            if seq != '':
                seq = sanitize(seq)
                above_run, below_run, above_start, above_end = read_sequence(seq, threshold)
                above_runs.append(above_run)
                below_runs.append(below_run)
                above_starts.append(above_start)
                above_ends.append(above_end)
                seq = ''
            titles.append(line[1:])
        else:
            seq += line

    if seq != '':
        seq = sanitize(seq)
        above_run, below_run, above_start, above_end = read_sequence(seq, threshold)
        above_runs.append(above_run)
        below_runs.append(below_run)
        above_starts.append(above_start)
        above_ends.append(above_end)
        seq = ''


    return titles, above_runs, below_runs, above_starts, above_ends
