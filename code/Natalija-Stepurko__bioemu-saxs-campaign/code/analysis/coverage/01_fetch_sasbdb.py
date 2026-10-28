"""Fetch SASBDB entry summaries for all 439 PeptoneDB-SAXS labels (cached per entry)."""
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from common import *

tab = pd.read_csv(SAXS_TABLE)
labels = tab.label.tolist()

def one(lab):
    return lab, cached_json(DL / "sasbdb" / f"{lab}.json",
                            f"https://www.sasbdb.org/rest-api/entry/summary/{lab}/?format=json", ok404=True)

with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(one, labels))
bad = [k for k, v in res.items() if "__status__" in v]
print(len(res), "fetched;", len(bad), "missing:", bad)
record_version("SASBDB", {"source": "https://www.sasbdb.org/rest-api/entry/summary/<LABEL>/", "n_entries": len(res) - len(bad)})
