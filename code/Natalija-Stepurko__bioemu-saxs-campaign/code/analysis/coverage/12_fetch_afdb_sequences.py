"""Sequences and annotations of the AFDB pre-training clusters.

The cluster definitions (script 11) map each of the ~50k clusters to its structure representatives
(AFDB v4 entries with mean pLDDT and its standard deviation). BioEmu-1 feeds the sequence of the member with
the highest mean pLDDT as input for every structure drawn from the cluster (Lewis et al. 2025, Materials and
Methods, "Pre-training on AFDB"); the other members are at least 80% identical to it.

For every member accession this fetches, through the UniProt ID-mapping service: the sequence, length,
organism and lineage, and the CATH-Gene3D cross-references. Accessions deleted from UniProtKB since AFDB v4
(about three quarters, mostly TrEMBL redundancy removal) are read from UniParc: sequence, first listed
organism (lineage from the UniProt taxonomy service) and CATH-Gene3D matches.
Writes tables/afdb_members.csv.gz (one row per member, flag `input` for the sequence used by the model).
"""
import gzip
import json
import time

import pandas as pd
from common import BT, S, TAB, record_version

D = BT / "afdb"
CHUNK = 90000      # ids per ID-mapping job (the service takes at most 100,000)
FIELDS = "accession,length,organism_id,organism_name,lineage,xref_gene3d,reviewed,sequence"


def idmapping(ids, to, out, fields=None):
    """Run one UniProt ID-mapping job (at most 100,000 ids) and cache its TSV result."""
    if out.exists():
        return out
    r = S.post("https://rest.uniprot.org/idmapping/run", data={"from": "UniProtKB_AC-ID", "to": to,
                                                              "ids": ",".join(ids)}, timeout=300)
    r.raise_for_status()
    job = r.json()["jobId"]
    print(f"  {out.name}: job {job}, {len(ids)} ids", flush=True)
    for _ in range(720):
        st = S.get(f"https://rest.uniprot.org/idmapping/status/{job}", timeout=120, allow_redirects=False)
        js = st.json() if st.headers.get("content-type", "").startswith("application/json") else {}
        if st.status_code in (302, 303) or "results" in js or js.get("jobStatus") == "FINISHED":
            break
        if js.get("jobStatus") in ("ERROR", "FAILURE"):
            raise RuntimeError(f"idmapping job {job}: {js}")
        time.sleep(5)
    else:
        raise RuntimeError(f"idmapping job {job} did not finish within an hour")
    db = "uniprotkb" if to == "UniProtKB" else "uniparc"
    params = {"format": "tsv"}
    if fields:
        params["fields"] = fields
    res = S.get(f"https://rest.uniprot.org/idmapping/{db}/results/stream/{job}", params=params, timeout=1800)
    res.raise_for_status()
    out.write_bytes(gzip.compress(res.content, mtime=0))
    if to == "UniProtKB":
        record_version("UniProtKB (AFDB cluster members)", {
            "release": res.headers.get("X-UniProt-Release"), "release_date": res.headers.get("X-UniProt-Release-Date"),
            "source": "https://rest.uniprot.org/idmapping (UniProtKB_AC-ID -> UniProtKB)"})
    return out


def idmapping_all(ids, to, stem, fields):
    """ID mapping in jobs of CHUNK ids. The result stream of a large job can end early without an error, so
    ids missing from a result are submitted again until a round adds nothing."""
    parts = []
    for i in range(0, len(ids), CHUNK):
        todo, k = ids[i:i + CHUNK], 0
        while todo:
            name = f"{stem}_{i // CHUNK}.tsv.gz" if k == 0 else f"{stem}_{i // CHUNK}_retry{k}.tsv.gz"
            df = pd.read_csv(idmapping(todo, to, D / name, fields), sep="\t", compression="gzip")
            parts.append(df)
            got = set(df.From)
            left = [x for x in todo if x not in got]
            if len(left) == len(todo):
                break
            todo, k = left, k + 1
    return pd.concat(parts).drop_duplicates("From")


defs = json.load(gzip.open(D / "afdb_cluster_definitions.json.gz"))
rows = []
for rep, members in defs.items():
    best = max(members, key=lambda m: m[1])[0]
    for name, plddt, plddt_sd in members:
        rows.append({"cluster": rep.split("-")[1], "member": name.split("-")[1], "afdb_entry": name,
                     "plddt_mean": plddt, "plddt_sd": plddt_sd, "input": name == best})
mem = pd.DataFrame(rows)
accs = sorted(set(mem.member))
print(len(defs), "clusters;", len(mem), "members;", len(accs), "accessions")

up = idmapping_all(accs, "UniProtKB", "uniprotkb", FIELDS)
up = up[up.Sequence.notna()]
# Entries deleted from UniProtKB since AFDB v4 (mostly TrEMBL redundancy removal) come back without data;
# their sequences, organisms and CATH-Gene3D matches are read from UniParc.
inactive = sorted(set(accs) - set(up.From))
print(len(up), "active UniProtKB entries;", len(inactive), "inactive, read from UniParc")
record_version("UniProtKB (AFDB cluster members): counts", {"n_ids": len(accs), "n_active": int(len(up))})
upa = idmapping_all(inactive, "UniParc", "uniparc", "upi,organism_id,organism,length,sequence,Gene3D")
upa = upa[upa.Sequence.notna()]
print(len(upa), "of them found in UniParc")
# UniParc lists every organism the sequence occurs in; keep the first, and look its lineage up
upa["Organism (ID)"] = upa["Organism ID"].astype(str).str.split(";").str[0].str.strip()
upa["Organism"] = upa["Organisms"].astype(str).str.split(";").str[0].str.strip()
upa["Gene3D"] = upa["Gene3D"].fillna("").astype(str).str.replace("G3DSA:", "", regex=False)
taxa = sorted({t for t in upa["Organism (ID)"] if t.isdigit()})
lin_file = D / "taxonomy_lineages.tsv.gz"
if not lin_file.exists():
    chunks = []
    for i in range(0, len(taxa), 50):
        qy = " OR ".join(f"tax_id:{t}" for t in taxa[i:i + 50])
        r = S.get("https://rest.uniprot.org/taxonomy/stream", params={"query": qy, "fields": "id,scientific_name,lineage",
                                                                     "format": "tsv"}, timeout=600)
        r.raise_for_status()
        chunks.append(r.text if not chunks else r.text.split("\n", 1)[1])
    lin_file.write_bytes(gzip.compress("".join(chunks).encode(), mtime=0))
lin = pd.read_csv(lin_file, sep="\t", compression="gzip", dtype=str)
lin_map = dict(zip(lin.iloc[:, 0], lin["Lineage"], strict=True))
upa["Taxonomic lineage"] = upa["Organism (ID)"].map(lin_map)
upa = upa.rename(columns={"Entry": "uniparc"})

cols = {"From": "member", "Entry": "uniprot_entry", "Length": "length", "Organism (ID)": "taxid", "Organism": "organism",
        "Taxonomic lineage": "lineage", "Gene3D": "gene3d", "Reviewed": "reviewed", "Sequence": "sequence"}
up = up.rename(columns=cols)[[c for c in cols.values()]]
up["source"] = "UniProtKB"
upa = upa.rename(columns=cols)[["member", "uniparc", "length", "taxid", "organism", "lineage", "gene3d", "sequence"]]
upa["source"] = "UniParc"
seqs = pd.concat([up, upa], ignore_index=True)
out = mem.merge(seqs, on="member", how="left")
out.to_csv(TAB / "afdb_members.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
print(out.source.value_counts(dropna=False).to_dict())
print("input sequences without a sequence:", int((out.input & out.sequence.isna()).sum()))
record_version("UniParc (AFDB cluster members deleted from UniProtKB)", {
    "source": "https://rest.uniprot.org/idmapping (UniProtKB_AC-ID -> UniParc)", "n_ids": len(inactive),
    "n_found": int(len(upa))})
