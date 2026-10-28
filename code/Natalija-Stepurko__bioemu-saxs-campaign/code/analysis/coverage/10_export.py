"""Copy the small outputs the main repository tracks into results/coverage/: the per-entry mapping, the
summary of all comparisons, the database versions, and a slim CATH superfamily rank table (one row per
superfamily: S35 clusters, rank, presence in the SAXS set) for the rank figure."""
import json
import shutil

import pandas as pd
from common import EXPORT, TAB, VERSIONS

EXPORT.mkdir(parents=True, exist_ok=True)
shutil.copyfile(TAB / "per_protein_mapping.csv", EXPORT / "per_protein_mapping.csv")
shutil.copyfile(TAB / "summary.json", EXPORT / "summary.json")
(EXPORT / "versions.json").write_text(json.dumps(json.loads(VERSIONS.read_text()), indent=2) + "\n")
rk = pd.read_csv(TAB / "cath_superfamily_rank.csv")
rk = rk[["sf", "s35_clusters", "rank", "in_ours_gene3d", "in_ours_ted", "name"]]
rk.to_csv(EXPORT / "cath_superfamily_rank.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
print("exported to", EXPORT, sorted(p.name for p in EXPORT.iterdir()))
