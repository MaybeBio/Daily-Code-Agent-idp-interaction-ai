import json
import sys
from pathlib import Path

if len(sys.argv) != 4:
    print("Usage: python split_json.py <input.json> <chunk_size> <out_dir>")
    sys.exit(1)

input_json = Path(sys.argv[1])
chunk_size = int(sys.argv[2])
out_dir = Path(sys.argv[3])

out_dir.mkdir(parents=True, exist_ok=True)

with input_json.open() as f:
    data = json.load(f)

data = [
    item for item in data
    if len(item["sequences"][0]["proteinChain"]["sequence"]) <= 1000
]

for i in range(0, len(data), chunk_size):
    chunk = data[i:i + chunk_size]
    chunk_id = i // chunk_size + 1
    out_file = out_dir / f"{input_json.stem}_chunk_{chunk_id}.json"
    with out_file.open("w") as f_out:
        json.dump(chunk, f_out, indent=2)

print(f"Created {len(data) // chunk_size + 1} chunks")

