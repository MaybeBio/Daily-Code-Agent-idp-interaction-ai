"""BioEmu-1's training data: what can be obtained of it, at the size needed (sequences and system identifiers,
no trajectories). Everything is cached under downloads/bioemu_training/ and recorded in versions.json.

1. AFDB pre-training clusters: the cluster definitions published in github.com/microsoft/bioemu
   (dev/afdb_cluster_definitions, pinned to the commit that added them).
2. MD datasets on Zenodo (Lewis et al. 2025, Table S1): CATH1 and CATH2 (record 15629740), octapeptides
   (15641199), MegaSim wild types and mutants (15641184). Dataset cards and record metadata are saved; the zip
   archives are read through HTTP range requests: their member listings, and the topology.pdb of every
   CATH1, CATH2 and MegaSim wild-type system (the sequence that was simulated). Octapeptides and MegaSim
   mutants are counted from the listings only.
3. MEGAscale (Tsuboyama et al. 2023, Zenodo 7992926): the curated Dataset2/Dataset3 table, read out of
   Processed_K50_dG_datasets.zip, for the wild-type and natural-domain sequences.
4. DESRES fast-folding proteins (not redistributable; used by Lewis et al. only for a separate model): the
   twelve named proteins, sequences from the PDB entries of Lindorff-Larsen et al. 2011.
5. The other public MD sets of Table S1 (DDR1 kinases, SETD8, SARS-CoV-2 and related coronavirus proteins,
   MHC class II, barnase-barstar): the named proteins, sequences from UniProtKB (chains of the polyproteins);
   the coronavirus system list from the Folding@home bucket listing.
"""
import csv
import gzip
import io
import json
import re

from concurrent.futures import ThreadPoolExecutor

from common import BT, S, cached_json, get, open_remote_zip, read_remote_member, record_version

BIOEMU_REPO = "microsoft/bioemu"
BIOEMU_COMMIT = "6f059c41017a2c653ef1567ce79d7fa260ad3670"    # "feat: add afdb cluster definitions (#109)"
ZENODO = {15629740: ["ONE_cath1.zip", "MSR_cath2.zip"], 15641199: ["ONE_octapeptides.zip"],
          15641184: ["MSR_megasim_merge.zip", "MSR_megasim_mutants_disp_allatom.zip"]}
TOPOLOGIES = {"ONE_cath1.zip", "MSR_cath2.zip", "MSR_megasim_merge.zip"}
MEGASCALE = 7992926
MEGASCALE_MEMBER = "Tsuboyama2023_Dataset2_Dataset3_20230416.csv"
# DESRES fast folders (Lindorff-Larsen et al. 2011, Science 334:517), PDB entries of the simulated variants
DESRES = {"Chignolin": "5AWL", "Trp-cage": "2JOF", "BBA": "1FME", "Villin": "2F4K", "WW domain": "2F21",
          "NTL9": "2HBA", "BBL": "2WXC", "Protein B": "1PRB", "Homeodomain": "2P6J", "Protein G": "1MI0",
          "alpha3D": "2A3D", "lambda-repressor": "1LMB"}
# Other public MD sets of Lewis et al. Table S1: (dataset, protein, UniProt accession, chain names or None)
OTHER_MD = [("FAH-DDR1", "DDR1", "Q08345", None), ("FAH-SETD8", "SETD8", "Q9NQR1", None),
            ("FUB-barnase-barstar", "barnase", "P00648", None), ("FUB-barnase-barstar", "barstar", "P11540", None),
            ("FUB-MHCII", "HLA-DRA", "P01903", None), ("FUB-MHCII", "HLA-DRB1*01:01", "P04229", None),
            ("FUB-MHCII", "HLA-DMA", "P28067", None), ("FUB-MHCII", "HLA-DMB", "P28068", None),
            ("FAH-sarscov2", "SARS-CoV-2 spike", "P0DTC2", None),
            ("FAH-sarscov2-exascale", "SARS-CoV-2 spike", "P0DTC2", None),
            ("FAH-sarscov2-exascale", "SARS-CoV-2 nucleoprotein", "P0DTC9", None),
            ("FAH-sarscov2-exascale", "SARS-CoV-2 replicase", "P0DTD1",
             ["nsp3", "nsp5", "nsp7", "nsp8", "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16"]),
            ("FAH-sarscov2-exascale", "SARS-CoV spike", "P59594", None),
            ("FAH-sarscov2-exascale", "SARS-CoV nucleoprotein", "P59595", None),
            ("FAH-sarscov2-exascale", "SARS-CoV replicase", "P0C6X7", ["nsp3", "nsp9", "nsp10", "nsp15", "nsp16"]),
            ("FAH-sarscov2-exascale", "HCoV-NL63 spike", "Q6Q1S2", None),
            ("FAH-sarscov2-exascale", "MERS-CoV replicase", "K9N7C7", ["nsp10", "nsp13", "nsp16"]),
            ("FAH-sarscov2-exascale", "ACE2", "Q9BYF1", None), ("FAH-sarscov2-exascale", "IL-6", "P05231", None),
            ("FAH-sarscov2-exascale", "IL-6 receptor", "P08887", None)]

# ------------------------------------------------------------------ 1 AFDB cluster definitions
A = BT / "afdb"
A.mkdir(exist_ok=True)
raw = f"https://raw.githubusercontent.com/{BIOEMU_REPO}/{BIOEMU_COMMIT}/dev/afdb_cluster_definitions"
for f in ("README.md", "afdb_cluster_definitions.json.gz"):
    if not (A / f).exists():
        (A / f).write_bytes(get(f"{raw}/{f}").content)
defs = json.load(gzip.open(A / "afdb_cluster_definitions.json.gz"))
record_version("BioEmu AFDB cluster definitions", {
    "source": f"https://github.com/{BIOEMU_REPO}/tree/{BIOEMU_COMMIT}/dev/afdb_cluster_definitions",
    "version": "file added 2025-04-25", "commit": BIOEMU_COMMIT, "n_clusters": len(defs),
    "n_structure_representatives": sum(len(v) for v in defs.values()), "afdb_release": "v4"})

# ------------------------------------------------------------------ 2 Zenodo MD datasets
Z = BT / "zenodo"
Z.mkdir(exist_ok=True)
for rec, zips in list(ZENODO.items()) + [(MEGASCALE, [])]:
    meta = cached_json(Z / f"record_{rec}.json", f"https://zenodo.org/api/records/{rec}")
    for fi in meta["files"]:
        if fi["key"].endswith((".md", ".txt")) and not (Z / f"{rec}_{fi['key']}").exists():
            (Z / f"{rec}_{fi['key']}").write_bytes(get(fi["links"]["self"]).content)
    lic = (meta["metadata"].get("license") or {}).get("id")
    if not lic and (Z / f"{rec}_LICENSE.txt").exists():           # the BioEmu records carry it as a file
        lic = (Z / f"{rec}_LICENSE.txt").read_text().splitlines()[0].strip()
    record_version(f"Zenodo {rec}", {"title": meta["metadata"]["title"], "doi": meta["doi"],
                                     "version": meta["metadata"].get("version") or f"published {meta['metadata']['publication_date']}",
                                     "updated": meta.get("updated"), "license": lic,
                                     "files": {fi["key"]: fi["size"] for fi in meta["files"]}})
    for zname in zips:
        out = Z / "listings" / f"{zname}.csv.gz"
        sysdir = Z / "systems" / zname.removesuffix(".zip")
        need_top = zname in TOPOLOGIES
        done = sysdir / "DONE"
        if out.exists() and (not need_top or done.exists()):
            continue
        url = f"https://zenodo.org/api/records/{rec}/files/{zname}/content"
        z, h = open_remote_zip(url)
        info = z.infolist()
        out.parent.mkdir(parents=True, exist_ok=True)
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["name", "file_size", "compress_size"])
        w.writerows((i.filename, i.file_size, i.compress_size) for i in info)
        out.write_bytes(gzip.compress(buf.getvalue().encode(), mtime=0))
        print(zname, len(info), "members listed", flush=True)
        if need_top:
            tops = [i for i in info if i.filename.endswith("/topology.pdb") and i.filename.count("/") == 2]
            todo = [i for i in tops if not (sysdir / i.filename.split("/")[1] / "topology.pdb").exists()]

            def one(i, url=url, sysdir=sysdir):
                dst = sysdir / i.filename.split("/")[1] / "topology.pdb"
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_bytes(read_remote_member(url, i))

            with ThreadPoolExecutor(4) as ex:
                for k, _ in enumerate(ex.map(one, todo)):
                    if k % 200 == 0:
                        print(f"  {zname}: {k}/{len(todo)} topologies", flush=True)
            done.write_text(str(len(tops)))

# ------------------------------------------------------------------ 3 MEGAscale curated set
M = BT / "megascale"
M.mkdir(exist_ok=True)
if not (M / f"{MEGASCALE_MEMBER}.gz").exists():
    z, h = open_remote_zip(f"https://zenodo.org/api/records/{MEGASCALE}/files/Processed_K50_dG_datasets.zip/content")
    member = next(i for i in z.infolist() if i.filename.endswith(MEGASCALE_MEMBER))
    (M / f"{MEGASCALE_MEMBER}.gz").write_bytes(gzip.compress(z.read(member), mtime=0))
record_version("MEGAscale (Tsuboyama et al. 2023)", {"source": f"https://doi.org/10.5281/zenodo.{MEGASCALE}",
                                                    "member": f"Processed_K50_dG_datasets.zip/{MEGASCALE_MEMBER}"})

# ------------------------------------------------------------------ 4 DESRES fast folders (PDB sequences)
F = BT / "desres"
F.mkdir(exist_ok=True)
for name, pdb in DESRES.items():
    p = F / f"{pdb}.fasta"
    if not p.exists():
        p.write_text(get(f"https://www.rcsb.org/fasta/entry/{pdb}").text)
record_version("DESRES fast-folder sequences", {"source": "https://www.rcsb.org/fasta/entry/<PDB>",
                                                "pdb_entries": DESRES})

# ------------------------------------------------------------------ 5 other public MD sets (named proteins)
O = BT / "other_md"
O.mkdir(exist_ok=True)
for _, _, acc, _ in OTHER_MD:
    cached_json(O / f"{acc}.json", f"https://rest.uniprot.org/uniprotkb/{acc}.json")
lst = S.get("https://fah-public-data-covid19-cryptic-pockets.s3.amazonaws.com/",
            params={"list-type": "2", "delimiter": "/", "prefix": ""}, timeout=60).text
systems = []
for top in re.findall(r"<Prefix>([^<]+/)</Prefix>", lst):
    sub = S.get("https://fah-public-data-covid19-cryptic-pockets.s3.amazonaws.com/",
                params={"list-type": "2", "delimiter": "/", "prefix": top}, timeout=60).text
    systems += [s for s in re.findall(r"<Prefix>([^<]+/)</Prefix>", sub) if s != top and "final_models" not in s]
(O / "fah_exascale_systems.json").write_text(json.dumps(sorted(systems), indent=1))
(O / "named_proteins.json").write_text(json.dumps(OTHER_MD, indent=1))
record_version("Other public MD sets (named proteins)", {
    "source": "UniProtKB entries named in Lewis et al. 2025 Table S1 and its references; Folding@home bucket "
              "fah-public-data-covid19-cryptic-pockets (system listing)", "n_exascale_systems_listed": len(systems)})
print("exascale systems:", len(systems))
print("done:", sorted(p.name for p in BT.iterdir()))

