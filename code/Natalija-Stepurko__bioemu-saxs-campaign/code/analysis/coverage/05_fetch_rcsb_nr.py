"""Length and source organism for one representative (first listed member) of every RCSB
30%-identity polymer-entity cluster, via the RCSB GraphQL data API."""
import json, time, requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from common import *

lines = (DL / "rcsb" / "clusters-by-entity-30.txt").read_text().split("\n")
# The file also clusters computed structure models (AF_* AlphaFold DB, MA_* ModelArchive);
# keep experimental entities only and use the first experimental member as representative.
all_clusters = [l.split() for l in lines if l.strip()]
clusters = [[m for m in c if not m.startswith(("AF_", "MA_"))] for c in all_clusters]
clusters = [c for c in clusters if c]
reps = [c[0] for c in clusters]
sizes = [len(c) for c in clusters]
print(len(all_clusters), "clusters in file;", len(clusters), "with >=1 experimental entity")

Q = """query($ids:[String!]!){ polymer_entities(entity_ids:$ids){ rcsb_id
  entity_poly{ rcsb_sample_sequence_length rcsb_entity_polymer_type }
  rcsb_entity_source_organism{ ncbi_taxonomy_id ncbi_scientific_name taxonomy_lineage{ id name } } } }"""

out_dir = DL / "rcsb" / "graphql_experimental"; out_dir.mkdir(exist_ok=True)


def batch(i):
    p = out_dir / f"batch_{i:05d}.json"
    if p.exists():
        return json.loads(p.read_text())
    ids = reps[i:i + 400]
    for t in range(5):
        try:
            r = S.post("https://data.rcsb.org/graphql", json={"query": Q, "variables": {"ids": ids}}, timeout=180)
        except requests.RequestException:
            time.sleep(5 * (t + 1)); continue
        if r.status_code == 200 and "data" in r.json():
            d = r.json()["data"]["polymer_entities"] or []
            p.write_text(json.dumps(d)); return d
    raise RuntimeError(f"batch {i} failed: {r.status_code} {r.text[:300]}")


with ThreadPoolExecutor(4) as ex:
    res = [e for b in ex.map(batch, range(0, len(reps), 400)) for e in b]

rows = []
for e in res:
    org = (e.get("rcsb_entity_source_organism") or [{}])[0] or {}
    lin = org.get("taxonomy_lineage") or []
    rows.append(dict(entity=e["rcsb_id"], length=(e.get("entity_poly") or {}).get("rcsb_sample_sequence_length"),
                     poly_type=(e.get("entity_poly") or {}).get("rcsb_entity_polymer_type"),
                     n_source_org=len(e.get("rcsb_entity_source_organism") or []),
                     taxid=org.get("ncbi_taxonomy_id"), sci_name=org.get("ncbi_scientific_name"),
                     lineage_names="|".join(x.get("name") or "" for x in lin),
                     lineage_ids="|".join(str(x.get("id")) for x in lin)))
df = pd.DataFrame(rows)
df["cluster_size"] = df.entity.map(dict(zip(reps, sizes)))
df.to_csv(TAB / "rcsb_nr30_representatives.csv.gz", index=False)
print(len(df), "representatives with metadata;", df.poly_type.value_counts().to_dict())
hdr = (DL / "rcsb" / "clusters30.headers").read_text()
lm = [l for l in hdr.splitlines() if l.lower().startswith("last-modified")]
record_version("RCSB clusters-by-entity-30", {"source": "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt",
                                              "file_last_modified": lm[0].split(":", 1)[1].strip() if lm else None,
                                              "n_clusters_in_file": len(all_clusters), "n_clusters_with_experimental_entity": len(clusters)})
