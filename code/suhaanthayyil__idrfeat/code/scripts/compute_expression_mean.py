import csv
import re
import sys

out_path = sys.argv[1]
rdr = csv.reader(sys.stdin)
header = next(rdr)
genes = header[1:]
sums = [0.0] * len(genes)
cnts = [0] * len(genes)

for row in rdr:
    if not row:
        continue
    for i, v in enumerate(row[1:]):
        if v:
            try:
                sums[i] += float(v)
                cnts[i] += 1
            except ValueError:
                continue


def strip(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", name).strip()


with open(out_path, "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["gene", "expr_mean"])
    for g, s, n in zip(genes, sums, cnts):
        if n > 0:
            w.writerow([strip(g), s / n])

print("genes with expression:", sum(1 for n in cnts if n > 0), file=sys.stderr)
