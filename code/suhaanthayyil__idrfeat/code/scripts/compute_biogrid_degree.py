import csv
import io
import zipfile
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "depmap"

zf = zipfile.ZipFile(DATA / "biogrid_organism.tab3.zip")
member = [n for n in zf.namelist() if "Homo_sapiens" in n][0]

adj: dict[str, set[str]] = {}
with zf.open(member) as fh:
    text = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
    header = next(text).rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    iA = idx["Official Symbol Interactor A"]
    iB = idx["Official Symbol Interactor B"]
    iType = idx.get("Experimental System Type")
    iOA = idx.get("Organism ID Interactor A")
    iOB = idx.get("Organism ID Interactor B")
    for line in text:
        f = line.rstrip("\n").split("\t")
        if iType is not None and f[iType] != "physical":
            continue
        if iOA is not None and iOB is not None and (f[iOA] != "9606" or f[iOB] != "9606"):
            continue
        a, b = f[iA], f[iB]
        if a == "-" or b == "-" or a == b:
            continue
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

with open(DATA / "biogrid_degree.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["gene", "hub_degree"])
    for g, s in adj.items():
        w.writerow([g, len(s)])

print("genes with physical degree:", len(adj))
