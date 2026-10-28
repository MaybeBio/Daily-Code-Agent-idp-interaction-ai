import pandas as pd
from argparse import ArgumentParser
from typing import List, Dict, Tuple

# --- Configure target compositions here ---
TARGET_FRACTIONS = [0.9, 0.75, 0.5, 0.25, 0.1]
GUIDELINE_TOTAL = 36000
MAX_FRACTION = 0.95
# ------------------------------------------

def unique_pairs(df: pd.DataFrame) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(len(df) - 1) for j in range(i + 1, len(df))]

def generate_compositions(L1: int, L2: int) -> List[Dict]:
    n1_ref = round((GUIDELINE_TOTAL * 0.5) / L1)
    n2_ref = round((GUIDELINE_TOTAL * 0.5) / L2)
    target_total = n1_ref * L1 + n2_ref * L2

    compositions, used_n1 = [], set()

    for frac in TARGET_FRACTIONS:
        if not (1 - MAX_FRACTION) <= frac <= MAX_FRACTION:
            continue

        n1_ideal = (frac * target_total) / L1
        best, min_error = None, float('inf')

        for n1 in range(max(1, round(n1_ideal) - 10), round(n1_ideal) + 11):
            if n1 in used_n1:
                continue
            remainder = target_total - n1 * L1
            if remainder > 0 and remainder % L2 == 0:
                n2 = remainder // L2
                actual_frac = (n1 * L1) / (n1 * L1 + n2 * L2)
                if actual_frac <= MAX_FRACTION and (1 - actual_frac) <= MAX_FRACTION:
                    error = abs(n1 - n1_ideal)
                    if error < min_error:
                        min_error, best = error, (n1, n2)

        if best:
            n1, n2 = best
            total_aa = n1 * L1 + n2 * L2
            compositions.append({
                'nuse1': n1, 'nuse2': n2, 'total_aa': total_aa,
                'fraction_seq1': round(n1 * L1 / total_aa, 2),
                'fraction_seq2': round(n2 * L2 / total_aa, 2),
            })
            used_n1.add(n1)

    for i, c in enumerate(compositions):
        c['composition'] = i
    return compositions

def write_fasta(df: pd.DataFrame, fasta_path: str) -> None:
    seen = set()
    with open(fasta_path, 'w') as f:
        for col in ('seq_name1', 'seq_name2'):
            seq_col = 'seq1' if col == 'seq_name1' else 'seq2'
            for name, seq in zip(df[col], df[seq_col]):
                if name not in seen:
                    f.write(f'>{name}\n{seq}\n')
                    seen.add(name)

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('-input',  help='Path to input CSV (e.g. gin_samples.csv)')
    parser.add_argument('-output', help='Path to output CSV (e.g. gin_prep.csv)')
    args = parser.parse_args()

    df_prot = pd.read_csv(args.input)
    pair_list = unique_pairs(df_prot)
    print(f"Found {len(pair_list)} unique pairs")

    results = []
    for idx, (i, j) in enumerate(pair_list):
        row1, row2 = df_prot.iloc[i], df_prot.iloc[j]
        name1, name2 = row1["seq_name"], row2["seq_name"]
        seq1,  seq2  = row1["fasta"],    row2["fasta"]

        comps = generate_compositions(len(seq1), len(seq2))
        for c in comps:
            c.update({'seq_name1': name1, 'seq_name2': name2, 'seq1': seq1, 'seq2': seq2})
        results.extend(comps)

    df_out = pd.DataFrame(results)
    df_out.to_csv(args.output, index=False)

    fasta_path = args.output.replace('.csv', '.fasta')
    write_fasta(df_out, fasta_path)

    print(f"Saved {len(results)} compositions -> {args.output}")
    print(f"Saved sequences -> {fasta_path}")

if __name__ == "__main__":
    main()