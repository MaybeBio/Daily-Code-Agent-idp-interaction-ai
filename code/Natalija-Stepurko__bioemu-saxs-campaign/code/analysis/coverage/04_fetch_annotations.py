"""Per accession: InterPro Gene3D (CATH superfamily) and Pfam locations, PDBe SIFTS best
structures, TED (AFDB domain chopping with CATH labels). All cached as raw JSON."""
import json
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from common import *

m = pd.read_csv(TAB / "mapping_step1.csv")
accs = sorted(set(m.uniprot.dropna()))


def paged(url, path):
    if path.exists():
        return json.loads(path.read_text())
    out, nxt = [], url
    while nxt:
        r = get(nxt, ok404=True)
        if r.status_code in (204, 404) or not r.content:
            break
        d = r.json(); out += d["results"]; nxt = d.get("next")
    path.write_text(json.dumps(out))
    return out


def one(acc):
    base = "https://www.ebi.ac.uk/interpro/api/entry"
    paged(f"{base}/cathgene3d/protein/uniprot/{acc}?page_size=200", DL / "interpro" / f"gene3d_{acc}.json")
    paged(f"{base}/pfam/protein/uniprot/{acc}?page_size=200", DL / "interpro" / f"pfam_{acc}.json")
    cached_json(DL / "pdbe" / f"best_{acc}.json",
                f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}", ok404=True)
    cached_json(DL / "ted" / f"ted_{acc}.json",
                f"https://ted.cathdb.info/api/v1/uniprot/summary/{acc}", ok404=True)
    return acc


(DL / "ted").mkdir(exist_ok=True)
with ThreadPoolExecutor(4) as ex:
    done = list(ex.map(one, accs))
print(len(done), "accessions annotated")
r = get("https://www.ebi.ac.uk/interpro/api/").json()["databases"]
for k in ["interpro", "cathgene3d", "pfam"]:
    record_version(f"InterPro:{k}", {"name": r[k]["name"], "version": r[k]["version"],
                                      "release_date": r[k]["releaseDate"],
                                      "source": "https://www.ebi.ac.uk/interpro/api/"})
record_version("PDBe SIFTS best_structures", {"source": "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/<ACC>",
                                             "version": "live API (weekly SIFTS update)"})
record_version("TED", {"source": "https://ted.cathdb.info/api/v1/uniprot/summary/<ACC>",
                       "version": "TED v1 (AFDB v4 models, CATH labels)"})
