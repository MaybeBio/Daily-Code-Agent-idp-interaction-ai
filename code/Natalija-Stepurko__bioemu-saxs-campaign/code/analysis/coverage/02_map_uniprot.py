"""Map each PeptoneDB-SAXS construct to a UniProtKB accession and a UniProt region.

Order of methods: (1) SASBDB molecule uniprot_code; (2) RCSB sequence search (MMseqs2,
>=90% identity) to a PDB polymer entity, then its SIFTS UniProt reference; (3) UniProt
peptide search (exact match of internal 20-mers). The construct is then aligned locally to
the UniProt canonical sequence to define the construct region.
"""
import json, time, re
import requests
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from common import *

tab = pd.read_csv(SAXS_TABLE)
res = pd.read_csv(RESOLV)
clean = set(res[(res.model == "bioemu") & (res.clean == True)].label)
ent = pd.read_csv(ENTRIES).set_index("label")

AL = Align.PairwiseAligner()
AL.mode = "local"; AL.substitution_matrix = substitution_matrices.load("BLOSUM62")
AL.open_gap_score = -11; AL.extend_gap_score = -1


def align(construct, ref):
    aln = AL.align(construct, ref)[0]
    qb, tb = aln.aligned  # blocks
    nid = sum(sum(construct[a] == ref[b] for a, b in zip(range(q0, q1), range(t0, t1)))
              for (q0, q1), (t0, t1) in zip(qb, tb))
    naln = sum(q1 - q0 for q0, q1 in qb)
    return dict(u_first=int(tb[0][0]) + 1, u_last=int(tb[-1][1]), q_first=int(qb[0][0]) + 1,
                q_last=int(qb[-1][1]), aln_identity=nid / max(naln, 1),
                construct_coverage=naln / len(construct))


def uniprot_entry(acc):
    p = DL / "uniprot" / f"{acc}.json"
    if p.exists():
        return json.loads(p.read_text())
    r = get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", ok404=True)
    if r.status_code == 200:
        record_version("UniProtKB", {"release": r.headers.get("X-UniProt-Release"),
                                     "release_date": r.headers.get("X-UniProt-Release-Date"),
                                     "source": "https://rest.uniprot.org/uniprotkb/<ACC>.json"})
    d = r.json() if r.status_code == 200 else {"__status__": r.status_code}
    p.write_text(json.dumps(d))
    return d


def resolve(acc):
    """Follow merged/demerged inactive entries; return (active_acc, entry, note)."""
    d = uniprot_entry(acc)
    note = ""
    if d.get("entryType") == "Inactive":
        inf = d.get("inactiveReason", {})
        to = inf.get("mergeDemergeTo") or []
        note = f"inactive:{inf.get('inactiveReasonType')}->{','.join(to)}"
        if to:
            a2, e2, n2 = resolve(to[0])
            return a2, e2, note + (f"; {n2}" if n2 else "")
        return acc, None, note
    return acc, d, note


def rcsb_search(seq):
    q = {"query": {"type": "terminal", "service": "sequence",
                   "parameters": {"evalue_cutoff": 0.1, "identity_cutoff": 0.9,
                                  "sequence_type": "protein", "value": seq}},
         "return_type": "polymer_entity",
         "request_options": {"scoring_strategy": "sequence", "paginate": {"start": 0, "rows": 25},
                             "results_verbosity": "verbose"}}
    for i in range(4):
        r = S.post("https://search.rcsb.org/rcsbsearch/v2/query", json=q, timeout=120)
        if r.status_code == 204:
            return []
        if r.status_code == 200:
            return r.json().get("result_set", [])
        time.sleep(3)
    r.raise_for_status()


def via_rcsb(lab, seq):
    cache = DL / "rcsb" / f"seqsearch_{lab}.json"
    hits = json.loads(cache.read_text()) if cache.exists() else rcsb_search(seq)
    cache.write_text(json.dumps(hits))
    # rank by query coverage * identity
    best = []
    for h in hits:
        mc = h["services"][0]["nodes"][0]["match_context"][0]
        cov = (mc["query_end"] - mc["query_beg"] + 1) / mc["query_length"]
        ent_id, num = h["identifier"].split("_")
        d = cached_json(DL / "rcsb" / f"entity_{h['identifier']}.json",
                        f"https://data.rcsb.org/rest/v1/core/polymer_entity/{ent_id}/{num}", ok404=True)
        refs = [x for x in (d.get("rcsb_polymer_entity_container_identifiers", {})
                            .get("reference_sequence_identifiers") or [])
                if x.get("database_name") == "UniProt"]
        if len(refs) == 1:  # skip chimeric entities mapping to several UniProt entries
            best.append((cov * mc["sequence_identity"], refs[0]["database_accession"], h["identifier"],
                         refs[0].get("provenance_source")))
    if not best:
        return None
    best.sort(key=lambda x: -x[0])
    return best[0]


def peptide_query(pep):
    for t in range(5):
        try:
            r = S.post("https://peptidesearch.uniprot.org/asyncrest",
                       data={"peps": pep, "lEQi": "off", "spOnly": "off"}, timeout=180, allow_redirects=False)
            loc = r.headers.get("Location")
            if loc:
                break
        except requests.RequestException:
            pass
        time.sleep(20)
    for i in range(90):
        time.sleep(5)
        try:
            g = S.get(loc, timeout=120, allow_redirects=False)
        except requests.RequestException:
            continue
        if g.status_code == 200:
            return [a for a in g.text.strip().split(",") if a]
    raise RuntimeError("peptide search timed out")


def via_peptide(lab, seq):
    """UniProt peptide search: exact matches of three internal 20-mers; keep accessions hit by all."""
    cache = DL / "uniprot" / f"peptide_{lab}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    n = len(seq)
    peps = [seq[int(f * (n - 20)):int(f * (n - 20)) + 20] for f in (0.25, 0.5, 0.75)] if n >= 40 else [seq]
    sets = [peptide_query(p) for p in peps]
    common = [a for a in sets[0] if all(a in s_ for s_ in sets[1:])]
    out = {"peps": peps, "accs": common, "service": "https://peptidesearch.uniprot.org/asyncrest"}
    cache.write_text(json.dumps(out))
    return out


def by_search(lab, seq, synthetic):
    hit = via_rcsb(lab, seq)
    if hit and hit[0] >= 0.5 and resolve(hit[1])[1] is not None:
        return hit[1], "RCSB seq search + SIFTS", f"{hit[2]} score={hit[0]:.2f} ({hit[3]})"
    pep = via_peptide(lab, seq)
    cands = []
    for a in sorted(pep["accs"], key=lambda a: len(a))[:40]:
        e = uniprot_entry(a)
        if "sequence" in e:
            al = align(seq, e["sequence"]["value"])
            cands.append((al["aln_identity"] * al["construct_coverage"], "Swiss-Prot" in e.get("entryType", ""), a))
    cands.sort(key=lambda x: (-x[0], not x[1]))
    if cands and cands[0][0] >= 0.9 and not synthetic:
        return cands[0][2], "UniProt peptide search", f"{len(pep['accs'])} accessions matched all 20-mers"
    if cands and synthetic:
        return None, None, f"synthetic construct; peptide hit {cands[0][2]} not used"
    return None, None, ""


rows = []
for lab, seq, L in zip(tab.label, tab.sequence, tab.length):
    mol = json.loads((DL / "sasbdb" / f"{lab}.json").read_text())["experiment"]["sample"]["molecule"][0]
    row = dict(label=lab, clean=lab in clean, construct_length=len(seq),
               disorder_class=ent.loc[lab, "disorder_class"], disorder_mean=ent.loc[lab, "disorder_mean"],
               sasbdb_uniprot=mol.get("uniprot_code"), sasbdb_range_first=mol.get("uniprot_range_first"),
               sasbdb_range_last=mol.get("uniprot_range_last"), sasbdb_organism=mol.get("organism"),
               molecule_name=mol.get("long_name"))
    synthetic = "synthetic" in (mol.get("organism") or "").lower()
    acc, method, detail = mol.get("uniprot_code"), None, ""
    if acc:
        method = "SASBDB"
        a2, e, note = resolve(acc)
        if e is None:  # deleted UniProtKB entry without a successor: search by sequence
            detail = f"SASBDB accession {acc} {note}; "
            acc, method, d2 = by_search(lab, seq, synthetic)
            detail += d2
        elif "sequence" in e:
            al = align(seq, e["sequence"]["value"])
            if al["aln_identity"] < 0.5 or al["construct_coverage"] < 0.3:
                detail = (f"SASBDB accession {acc} rejected (identity {al['aln_identity']:.2f}, "
                          f"coverage {al['construct_coverage']:.2f}); ")
                acc, method, d2 = by_search(lab, seq, synthetic)
                detail += d2
    else:
        acc, method, detail = by_search(lab, seq, synthetic)
    isoform = None
    if acc:
        acc0 = acc
        acc, e, note = resolve(acc)
        if note:
            detail = (detail + "; " if detail else "") + f"{acc0} {note}"
        if "-" in acc:  # isoform: annotate on the canonical entry (InterPro/SIFTS coordinates)
            isoform = acc
            acc, e, _ = resolve(acc.split("-")[0])
    else:
        e = None
    row.update(uniprot=acc, uniprot_isoform=isoform, map_method=method or "unmapped", map_detail=detail)
    if e and "sequence" in e:
        useq = e["sequence"]["value"]
        al = align(seq, useq)
        row.update(al, uniprot_length=len(useq), reviewed="Swiss-Prot" in e.get("entryType", ""),
                   up_organism=e["organism"]["scientificName"], up_taxid=e["organism"]["taxonId"],
                   up_lineage="|".join(e["organism"].get("lineage", [])))
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(TAB / "mapping_step1.csv", index=False)
print(df.map_method.value_counts())
print(df.groupby("clean").map_method.value_counts())
print((df.aln_identity < 0.9).sum(), "with identity<0.9;", (df.construct_coverage < 0.8).sum(), "with coverage<0.8")
