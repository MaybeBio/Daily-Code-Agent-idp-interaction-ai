"""One row per PeptoneDB-SAXS entry with the SASBDB record's measurement metadata: curve type, size-exclusion
coupling, instrument and source type, sample conditions (cell temperature, buffer, pH, additives, concentration
range) and the molecule's organism and name. Read from the entry summaries cached by 01_fetch_sasbdb.py (the same
REST endpoint; entries missing from the cache are fetched there first).

Columns
  curve          SASBDB "type_of_curve" (Merged, Single concentration, Extrapolated, SEC-SAXS, ...)
  sec            a size-exclusion column is recorded (experiment.seccolumn) or the curve type names SEC
  source_type    X-ray synchrotron, X-ray in house, ...; instrument and beamline as recorded
  temp           experiment.cell_temperature (deg C)
  buffer, buffer_conc_mM, buffer_additive, buffer_comment, ph   experiment.sample.buffer
  cmin, cmax     experiment.concentration_min / max (mg/ml)
  organisms, mol_names, uniprot   per molecule in the sample, joined with " | "
  n_mol          number of molecule records in the sample
Writes tables/sasbdb_metadata.csv and exports it to results/coverage/sasbdb_metadata.csv.
"""
import json
import subprocess
import sys

import pandas as pd
from common import DL, EXPORT, SAXS_TABLE, TAB

labels = pd.read_csv(SAXS_TABLE).label.tolist()
if any(not (DL / "sasbdb" / f"{lab}.json").exists() for lab in labels):
    subprocess.run([sys.executable, "01_fetch_sasbdb.py"], check=True)


def join(xs):
    return " | ".join("" if x is None else str(x) for x in xs)


rows = []
for lab in labels:
    e = json.loads((DL / "sasbdb" / f"{lab}.json").read_text())
    if "__status__" in e:
        rows.append({"label": lab})
        continue
    x = e.get("experiment") or {}
    ins = x.get("instrument") or {}
    smp = x.get("sample") or {}
    buf = smp.get("buffer") or {}
    mol = smp.get("molecule") or []
    curve = e.get("type_of_curve")
    rows.append({
        "label": lab, "curve": curve,
        "sec": bool(x.get("seccolumn")) or ("SEC" in (curve or "").upper()),
        "source_type": ins.get("type_of_source"), "instrument": ins.get("name"), "beamline": ins.get("beamline_name"),
        "temp": x.get("cell_temperature"), "buffer": buf.get("name"), "buffer_conc_mM": buf.get("concentration")
        if buf.get("concentration_unit") == "mM" else None, "buffer_additive": buf.get("additive"),
        "buffer_comment": buf.get("comment"), "ph": buf.get("ph"),
        "cmin": x.get("concentration_min"), "cmax": x.get("concentration_max"),
        "concentration_unit": x.get("concentration_unit"),
        "n_mol": len(mol), "organisms": join(m.get("organism") for m in mol),
        "mol_names": join(m.get("long_name") for m in mol), "uniprot": join(m.get("uniprot_code") for m in mol),
        "molecule_source": join(m.get("molecule_source") for m in mol),
        "project": (e.get("project") or {}).get("title"),
        "contributors": join(c.get("contributor_surname") for c in e.get("contributor") or x.get("contributor") or []),
    })
df = pd.DataFrame(rows)
df.to_csv(TAB / "sasbdb_metadata.csv", index=False)
EXPORT.mkdir(parents=True, exist_ok=True)
df.to_csv(EXPORT / "sasbdb_metadata.csv", index=False)
print(len(df), "entries;", df.curve.value_counts(dropna=False).to_dict(), df.source_type.value_counts(dropna=False).to_dict())
