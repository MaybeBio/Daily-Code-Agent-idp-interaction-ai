"""Per-entry table: construct region, CATH-Gene3D / TED / Pfam domains within the construct
(>=50% of the domain inside the construct region), PDB coverage of the construct region
(SIFTS), UniProt 'Disordered' region content, organism group."""
import json
import numpy as np
import pandas as pd
from common import *

m = pd.read_csv(TAB / "mapping_step1.csv")
ATOMIC = {"X-ray diffraction", "Electron Microscopy", "Solution NMR", "Solid-state NMR",
          "Electron crystallography", "Neutron Diffraction", "X-ray powder diffraction", "Fiber diffraction"}
MIN_OVL = 0.5


def segs_from_interpro(results, acc):
    out = []
    for r in results:
        ida = r["metadata"]["accession"]
        for p in r["proteins"]:
            for loc in p.get("entry_protein_locations") or []:
                out.append((ida, [(f["start"], f["end"]) for f in loc["fragments"]]))
    return out


def overlap_frac(segs, a, b):
    n = sum(e - s + 1 for s, e in segs)
    ov = sum(max(0, min(e, b) - max(s, a) + 1) for s, e in segs)
    return ov / n if n else 0


def ted_domains(acc):
    d = json.loads((DL / "ted" / f"ted_{acc}.json").read_text())
    out = []
    for x in d.get("data") or []:
        segs = [tuple(int(v) for v in s.split("-")) for s in x["chopping"].split("_")]
        out.append((x.get("cath_label") or "-", x.get("cath_assignment_level"), segs))
    return out


def taxonomy_by_name(name):
    """Organism group for unmapped entries from the SASBDB organism name via UniProt taxonomy."""
    if not isinstance(name, str) or not name.strip():
        return "Other", None
    if "synthetic" in name.lower():
        return "Other", "synthetic construct"
    nm = name.split("/")[0].strip()
    nm = " ".join(nm.split()[:2])
    d = cached_json(DL / "uniprot" / f"taxonomy_{nm.replace(' ', '_')}.json",
                    "https://rest.uniprot.org/taxonomy/search",
                    params={"query": f'scientific:"{nm}"', "fields": "id,scientific_name,lineage", "size": 1})
    res = d.get("results") or []
    if not res:
        return "Other", nm
    r = res[0]
    lin = [x["scientificName"] for x in r.get("lineage", [])]
    return org_category(r["taxonId"], r["scientificName"], lin), r["scientificName"]


rows = []
for _, r in m.iterrows():
    row = r.to_dict()
    acc = r.uniprot
    if isinstance(acc, str) and not np.isnan(r.u_first):
        a, b = int(r.u_first), int(r.u_last)
        e = json.loads((DL / "uniprot" / f"{acc}.json").read_text())
        L = e["sequence"]["length"]
        # CATH-Gene3D
        g3 = segs_from_interpro(json.loads((DL / "interpro" / f"gene3d_{acc}.json").read_text()), acc)
        g3_in = [(ida.replace("G3DSA:", ""), s) for ida, s in g3 if overlap_frac(s, a, b) >= MIN_OVL]
        row["gene3d_superfamilies"] = ";".join(sorted(x for x, _ in g3_in))
        row["n_gene3d_domains"] = len(g3_in)
        row["n_gene3d_domains_fulllength"] = len(g3)
        # Pfam
        pf = segs_from_interpro(json.loads((DL / "interpro" / f"pfam_{acc}.json").read_text()), acc)
        pf_in = sorted({ida for ida, s in pf if overlap_frac(s, a, b) >= MIN_OVL})
        row["pfam"] = ";".join(pf_in); row["n_pfam"] = len(pf_in)
        # TED (AFDB v4 chopping); labels at H level only count as superfamilies
        td = [(lab, lev, s) for lab, lev, s in ted_domains(acc) if overlap_frac(s, a, b) >= MIN_OVL]
        row["ted_superfamilies"] = ";".join(sorted(lab for lab, lev, _ in td if lev == "H"))
        row["ted_labels_any_level"] = ";".join(sorted(lab for lab, lev, _ in td if lab != "-"))
        row["n_ted_domains"] = len(td)
        # PDB coverage of the construct region (SIFTS best_structures, atomic-model methods only)
        sb = json.loads((DL / "pdbe" / f"best_{acc}.json").read_text())
        chains = [c for c in sb.get(acc, []) if c["experimental_method"] in ATOMIC] if "__status__" not in sb else []
        cov = np.zeros(b - a + 1, bool)
        for c in chains:
            s, t = max(c["unp_start"], a), min(c["unp_end"], b)
            if t >= s:
                cov[s - a:t - a + 1] = True
        row["n_pdb_chains_overlapping"] = sum(1 for c in chains if min(c["unp_end"], b) >= max(c["unp_start"], a))
        row["pdb_region_coverage"] = cov.mean()
        row["has_pdb_any_overlap"] = bool(cov.any())
        row["has_pdb_50"] = bool(cov.mean() >= 0.5)
        row["uniprot_has_pdb_xref"] = any(x["database"] == "PDB" for x in e.get("uniProtKBCrossReferences", []))
        # UniProt 'Disordered' regions (MobiDB-lite and curated)
        dmask = np.zeros(L + 1, bool)
        for f in e.get("features", []):
            if f["type"] == "Region" and f.get("description") == "Disordered":
                s, t = f["location"]["start"]["value"], f["location"]["end"]["value"]
                if s and t:
                    dmask[s:t + 1] = True
        row["up_disordered_frac_construct"] = dmask[a:b + 1].mean()
        row["up_disordered_frac_fulllength"] = dmask[1:].mean()
        runs, cur = [], 0
        for v in dmask[a:b + 1]:
            cur = cur + 1 if v else 0; runs.append(cur)
        row["up_longest_disordered_construct"] = max(runs or [0])
        runs, cur = [], 0
        for v in dmask[1:]:
            cur = cur + 1 if v else 0; runs.append(cur)
        row["up_longest_disordered_fulllength"] = max(runs or [0])
        row["org_group"] = org_category(e["organism"]["taxonId"], e["organism"]["scientificName"],
                                        e["organism"].get("lineage", []))
        row["organism"] = e["organism"]["scientificName"]
        row["organism_source"] = "UniProt"
    else:
        g, nm = taxonomy_by_name(r.sasbdb_organism)
        row.update(org_group=g, organism=nm or r.sasbdb_organism, organism_source="SASBDB name",
                   n_gene3d_domains=0, n_ted_domains=np.nan, n_pfam=np.nan, gene3d_superfamilies="",
                   ted_superfamilies="", has_pdb_any_overlap=np.nan, has_pdb_50=np.nan)
    rows.append(row)

df = pd.DataFrame(rows)
for c in ["gene3d_superfamilies", "ted_superfamilies", "pfam", "ted_labels_any_level"]:
    df[c] = df[c].fillna("")
df["has_cath_gene3d"] = df.n_gene3d_domains > 0
df["has_cath_ted"] = df.ted_superfamilies != ""
df["has_cath_any"] = df.has_cath_gene3d | df.has_cath_ted
df.to_csv(TAB / "per_protein_mapping.csv", index=False)
c = df[df.clean]
print(len(df), "entries;", len(c), "clean")
print("clean with Gene3D superfamily:", c.has_cath_gene3d.sum(), " TED:", c.has_cath_ted.sum(), " either:", c.has_cath_any.sum())
print(c.org_group.value_counts())
print("no PDB >=50% coverage:", (c.has_pdb_50 == False).sum(), " no overlap:", (c.has_pdb_any_overlap == False).sum())
