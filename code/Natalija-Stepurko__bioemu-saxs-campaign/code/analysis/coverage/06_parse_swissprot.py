"""Parse the Swiss-Prot flat file into one row per entry: length, organism group, disordered
residues (REGION /note="Disordered", mostly MobiDB-lite), Gene3D superfamilies and match
counts, PDB and Pfam cross-reference counts."""
import gzip, re
import pandas as pd
from common import *

rows = []
cur = None


def flush(c):
    if c is None:
        return
    L = c["length"]
    mask = bytearray(L)
    for a, b in c["dis"]:
        for i in range(max(a, 1) - 1, min(b, L)):
            mask[i] = 1
    lin = [x.strip() for x in " ".join(c["oc"]).rstrip(".").split(";")]
    rows.append(dict(acc=c["acc"], length=L, taxid=c["taxid"], org_group=org_category(c["taxid"], c["os"], lin),
                     kingdom=lin[0] if lin else "", dis_res=sum(mask), longest_dis=max([b - a + 1 for a, b in c["dis"]] or [0]),
                     n_pdb=c["pdb"], n_pfam=c["pfam"], gene3d=";".join(c["g3d"]), n_gene3d_matches=c["g3dn"]))


pend_region = None
with gzip.open(DL / "uniprot" / "uniprot_sprot.dat.gz", "rt") as fh:
    for line in fh:
        tag = line[:5]
        if tag == "ID   ":
            flush(cur)
            cur = dict(acc=None, length=int(re.search(r"(\d+) AA\.", line).group(1)), os="", oc=[], taxid=None,
                       dis=[], pdb=0, pfam=0, g3d=[], g3dn=0)
            pend_region = None
        elif tag == "AC   " and cur["acc"] is None:
            cur["acc"] = line[5:].split(";")[0].strip()
        elif tag == "OS   ":
            cur["os"] += line[5:].strip() + " "
        elif tag == "OC   ":
            cur["oc"].append(line[5:].strip())
        elif tag == "OX   " and cur["taxid"] is None:
            cur["taxid"] = int(re.search(r"NCBI_TaxID=(\d+)", line).group(1))
        elif tag == "DR   ":
            db = line[5:].split(";")[0]
            if db == "PDB":
                cur["pdb"] += 1
            elif db == "Pfam":
                cur["pfam"] += 1
            elif db == "Gene3D":
                f = [x.strip() for x in line[5:].split(";")]
                cur["g3d"].append(f[1]); cur["g3dn"] += int(f[3].rstrip(".") or 1)
        elif line.startswith("FT   REGION"):
            m = re.search(r"(\d+)\.\.>?(\d+)", line)
            pend_region = (int(m.group(1)), int(m.group(2))) if m else None
        elif line.startswith("FT                   /note=\"Disordered\"") and pend_region:
            cur["dis"].append(pend_region); pend_region = None
        elif line.startswith("FT   ") and not line.startswith("FT    "):
            pend_region = None
    flush(cur)

df = pd.DataFrame(rows)
df.to_csv(TAB / "swissprot_summary.csv.gz", index=False)
rel = (DL / "uniprot" / "reldate.txt").read_text().splitlines()
record_version("UniProtKB/Swiss-Prot flat file", {"release": rel[1].strip(), "n_entries": len(df),
               "source": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"})
print(len(df)); print(df.org_group.value_counts(normalize=True).round(4))
