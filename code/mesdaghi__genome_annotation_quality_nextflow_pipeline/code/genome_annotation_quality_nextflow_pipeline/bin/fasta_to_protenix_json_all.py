import json
import sys

def extract_accession(header: str) -> str:
    """
    Extract UniProt accession from FASTA header.
    Examples:
      >sp|P12345|PROT_HUMAN ...
      >tr|Q9XYZ1|...
    Returns:
      P12345
    """
    header = header.lstrip(">").strip()
    parts = header.split("|")

    if len(parts) >= 3:
        return parts[1].strip()   # UniProt accession: second field

    # fallback — take first token if no pipe structure
    return header.split()[0]

def fasta_to_protenix_json_all(fasta_file, json_file):
    sequences = []
    header = None
    seq = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous sequence block
                if header and seq:
                    accession = extract_accession(header)
                    sequences.append({
                        "sequences": [
                            {
                                "proteinChain": {
                                    "sequence": "".join(seq),
                                    "count": 1
                                }
                            }
                        ],
                        "name": accession
                    })

                header = line  # Start new header
                seq = []       # Reset sequence buffer

            else:
                seq.append(line)

        # Save the last sequence (EOF case)
        if header and seq:
            accession = extract_accession(header)
            sequences.append({
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": "".join(seq),
                            "count": 1
                        }
                    }
                ],
                "name": accession
            })

    # Write full dataset
    with open(json_file, "w") as out:
        json.dump(sequences, out, indent=4)

    print(f"✅ JSON written: {len(sequences)} sequences → {json_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fasta_to_protenix_json_all.py input.fasta output.json")
        sys.exit(1)

    fasta_to_protenix_json_all(sys.argv[1], sys.argv[2])


