---
license: cc-by-4.0
pretty_name: OpenDNA Structures (CATH v4.4.0 S100 domains)
tags:
- biology
- protein
- protein-structure
- cath
size_categories:
- 100K<n<1M
configs:
- config_name: default
  data_files:
  - split: train
    path: data/train.parquet
  - split: validation
    path: data/validation.parquet
  - split: test
    path: data/test.parquet
---

# OpenDNA Structures (CATH v4.4.0 S100 domains)

106910 CATH v4.4.0 S100 domain sequences with CATH class, architecture, topology and superfamily labels and wwPDB entry metadata. Splits are by CATH topology, merged with 30% identity sequence clusters, so no fold or close homolog appears in more than one split. Coordinates are not included.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| CATH v4.4.0 domain classification, S100 non-redundant representatives | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [cath-domain-list-S100-v4_4_0.txt](https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-S100-v4_4_0.txt) | `9a202046a68c80436f76a3631fe0f84fc4f027a89c53657a5dcd4c97231cf562` |
| CATH v4.4.0 domain sequences (from ATOM records), S100 representatives | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [cath-domain-seqs-S100-v4_4_0.fa](https://download.cathdb.info/cath/releases/all-releases/v4_4_0/sequence-data/cath-domain-seqs-S100-v4_4_0.fa) | `d7d1efdf9a169456ad762faca0dc80b64f734ddeaf31481c69eb68ce5a77cd83` |
| wwPDB entries index: deposition date, title, resolution and experimental method per entry | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | [entries.idx](https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx) | `b9d6b90c4b1c4bdfbccb61210276e199aecfdefe73ed754fb6792468e174e55f` |

CATH is CC BY 4.0 and wwPDB data is CC0 1.0; the derived dataset is released under CC BY 4.0.

## Construction

1. Read the CATH v4.4.0 S100 domain list and the matching ATOM-record domain sequences.
2. Drop domains with no sequence or with non-standard residues. Join wwPDB entry metadata (method, deposition date) by PDB id.
3. Cluster domain sequences with `mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 0 -s 7.5`.
4. Build a graph whose nodes are CATH topologies and sequence clusters, with one edge per domain; connected components are the split unit. No topology and no 30% identity cluster spans two splits.
5. Assign components to train/validation/test (80/10/10 of domains, seed 20261114).

## Schema

| Column | Type | Description |
|---|---|---|
| `domain_id` | string | CATH domain id (PDB id, chain, domain number) |
| `pdb_id` | string | PDB entry id (lowercase) |
| `chain_id` | string | PDB chain id |
| `segments` | string | Domain boundaries in PDB residue numbering as given in the CATH FASTA header |
| `sequence` | string | Domain sequence from ATOM records (CATH) |
| `length` | int32 | Sequence length |
| `cath_class` | int32 | CATH class (C) |
| `cath_architecture` | string | CATH architecture code (C.A) |
| `cath_topology` | string | CATH topology code (C.A.T) |
| `cath_superfamily` | string | CATH homologous superfamily code (C.A.T.H) |
| `resolution` | float64 (nullable) | Resolution in Angstrom from CATH; null for NMR and other non-diffraction methods |
| `experimental_method` | string (nullable) | Experimental method from the wwPDB entries index |
| `deposition_date` | string (nullable) | Deposition date (ISO 8601) from the wwPDB entries index |
| `cluster_id` | string | MMseqs2 30% identity cluster (representative domain id) |
| `split_group` | string | Connected component of topologies and clusters; the split unit |

## Splits

| Split | domains | components |
|---|---:|---:|
| train | 85325 | 1006 |
| validation | 10500 | 169 |
| test | 11085 | 245 |
| total | 106910 | |

### Cleaning statistics

| Step | Count |
|---|---:|
| dropped: no CATH sequence | 0 |
| dropped: non-standard residue | 16341 |
| kept without wwPDB metadata (entry obsolete or missing) | 1165 |
| final domains | 106910 |
| CATH topologies | 1465 |
| MMseqs2 clusters | 26748 |
| split groups (components) | 1420 |
| largest split group (domains) | 25433 |

## Known limitations

- Coordinates are not included; this release carries sequences and labels only.
- S100 representatives remove identical sequences only, so near-identical domains remain within a split.
- Splitting by topology means test folds are unseen in training. That is a harder setting than a random split and the realised split fractions can deviate from 80/10/10 when a component is large.
- Domain sequences come from ATOM records, so unobserved residues are absent.
- Domains whose CATH sequence contains `X` (unknown or modified residues) are dropped; this removes about 13% 
             of S100 domains and may under-represent structures with many modified residues.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{sillitoe2021cath,
  title   = {CATH: increased structural coverage of functional space},
  author  = {Sillitoe, Ian and Bordin, Nicola and Dawson, Natalie and others},
  journal = {Nucleic Acids Research},
  volume  = {49},
  number  = {D1},
  pages   = {D266--D273},
  year    = {2021},
  doi     = {10.1093/nar/gkaa1079}
}

@article{wwpdb2019,
  title   = {Protein Data Bank: the single global archive for 3D macromolecular structure data},
  author  = {{wwPDB consortium}},
  journal = {Nucleic Acids Research},
  volume  = {47},
  number  = {D1},
  pages   = {D520--D528},
  year    = {2019},
  doi     = {10.1093/nar/gky949}
}
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
