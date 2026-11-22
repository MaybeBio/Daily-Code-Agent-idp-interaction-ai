---
license: cc-by-4.0
pretty_name: OpenDNA Disorder (DisProt 2025_12 per-residue labels)
tags:
- biology
- protein
- intrinsic-disorder
- token-classification
size_categories:
- 1K<n<10K
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

# OpenDNA Disorder (DisProt 2025_12 per-residue labels)

Per-residue intrinsic disorder labels for 3186 proteins from DisProt release 2025_12, prepared for training a disorder head on frozen protein language model embeddings. The test split holds out the CAID3 Disorder-PDB proteins and every protein in the same 30% identity cluster.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| DisProt release 2025_12, all entries, ambiguous and obsolete annotations excluded | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [disprot_2025_12.json](https://disprot.org/api/search?release=2025_12&show_ambiguous=false&show_obsolete=false&format=json&page_size=100000) | `4dbdde980d29c95c6e14aab9eaad7eaec7cc7bd42c3668c876a1e88b4eae2981` |
| CAID3 Disorder-PDB reference set (BioComputingUP/CAID @ 7bab7e8); only the protein ids and sequences are used, to hold them out | [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/) | [disorder_pdb.fasta](https://raw.githubusercontent.com/BioComputingUP/CAID/7bab7e8880d7f949e480fdb377c5a5c8546ae336/demo-data/references/disorder_pdb.fasta) | `b65848395557c751a818dc38049c89d95534c05fbe6f14651e9f5bfd2bb60469` |

DisProt is distributed under CC BY 4.0, and this derived dataset is released under CC BY 4.0. The CAID3 reference is used only to choose which proteins go to test; none of its labels are redistributed here.

## Construction

1. Download DisProt release 2025_12 through the DisProt API with ambiguous and obsolete annotations excluded.
2. Drop entries whose sequence contains non-standard residues or disagrees with the recorded length.
3. Label each residue 1 if it lies inside any region with term `IDPO:0000002` (Structural state: disorder), using DisProt's 1-based inclusive coordinates, and 0 otherwise. Molten-globule and pre-molten-globule regions are not labelled. Function, transition and binding annotations are not used.
4. Read the CAID3 Disorder-PDB reference (319 proteins) from the CAID repository at a pinned commit.
5. Cluster DisProt sequences together with the CAID3 reference sequences using `mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 0 -s 7.5`.
6. Every cluster containing a CAID3 protein goes to test, so no train or validation protein shares a 30%-identity cluster with a CAID3 target. Non-CAID proteins in those clusters are also in test; filter on `caid3_target` to evaluate on the exact CAID3 subset.
7. The remaining clusters are split 90/10 (by residue count) into train/validation with seed 20261114.
8. Test labels are DisProt's, not CAID's; CAID scores against a different negative definition (residues observed in PDB), so numbers on this test split are not directly comparable to published CAID results.

## Schema

| Column | Type | Description |
|---|---|---|
| `disprot_id` | string | DisProt accession |
| `uniprot_acc` | string | UniProt accession of the DisProt entry |
| `organism` | string (nullable) | Source organism as given by DisProt |
| `ncbi_taxon_id` | int64 (nullable) | NCBI taxonomy id |
| `sequence` | string | Full UniProt sequence used by DisProt |
| `length` | int32 | Sequence length |
| `disorder` | list<int8> | Per-residue label: 1 = inside a DisProt IDPO:0000002 (disorder) region, 0 = not annotated as disordered (not evidence of order) |
| `n_disordered` | int32 | Number of residues labelled 1 |
| `caid3_target` | bool | True if the protein is in the CAID3 Disorder-PDB reference set |
| `cluster_id` | string | MMseqs2 30% identity cluster (representative id); the split unit |

## Splits

| Split | proteins | clusters |
|---|---:|---:|
| train | 2455 | 2031 |
| validation | 327 | 271 |
| test | 404 | 297 |
| total | 3186 | |

### Cleaning statistics

| Step | Count |
|---|---:|
| DisProt 2025_12 entries | 3199 |
| dropped: non-standard residues in sequence | 13 |
| dropped: sequence length disagrees with `length` | 0 |
| IDPO:0000002 regions used | 7100 |
| regions skipped: outside the sequence | 0 |
| regions not labelled: (pre-)molten globule | 152 |
| proteins with no IDPO:0000002 residue (kept, all zeros) | 22 |
| final proteins | 3186 |
| residues / labelled disordered | 1847875 / 319131 |
| CAID3 Disorder-PDB proteins | 319 |
| CAID3 proteins present in this DisProt release | 318 |
| non-CAID proteins moved to test (clustered with a CAID3 protein) | 86 |
| MMseqs2 clusters (including CAID-only) | 2599 |
| residues train / validation / test | 1469314 / 163130 / 215431 |

## Known limitations

- A label of 0 means "not annotated as disordered in DisProt", not "experimentally ordered". DisProt annotates regions with experimental evidence of disorder; the rest of each sequence is mostly unexamined, so negatives are noisy and the positive rate is underestimated.
- Sequences are full UniProt sequences, some longer than 30,000 residues. ESM-2 accepts at most 1,024 residues per pass; training code must crop or tile.
- Coverage is biased toward well-studied proteins and organisms (human, model organisms, viruses).
- Clustering uses 30% identity over 80% bidirectional coverage; a short disordered region shared between otherwise unrelated proteins will not be detected by this criterion.
- Only CAID3 Disorder-PDB is held out. Older CAID rounds are not; their targets may be in train.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{aspromonte2024disprot,
  title   = {DisProt in 2024: improving function annotation of intrinsically disordered proteins},
  author  = {Aspromonte, Maria Cristina and others},
  journal = {Nucleic Acids Research},
  volume  = {52},
  number  = {D1},
  pages   = {D434--D441},
  year    = {2024},
  doi     = {10.1093/nar/gkad928}
}

@misc{caid_repository,
  title        = {CAID: Critical Assessment of protein Intrinsic Disorder, assessment software and reference data},
  author       = {{BioComputingUP}},
  howpublished = {\url{https://github.com/BioComputingUP/CAID}},
  note         = {commit 7bab7e8880d7f949e480fdb377c5a5c8546ae336}
}
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
