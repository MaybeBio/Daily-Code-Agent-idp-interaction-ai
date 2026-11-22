---
license: other
license_name: cc-by-sa-2.1-jp
license_link: https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en
pretty_name: OpenDNA Solubility (eSOL, E. coli cell-free)
tags:
- biology
- protein
- protein-solubility
- escherichia-coli
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

# OpenDNA Solubility (eSOL, E. coli cell-free)

Solubility of 3127 E. coli K-12 proteins expressed in a chaperone-free reconstituted cell-free system (eSOL, Niwa et al. 2009), with sequences from UniProtKB. Prepared for training a solubility regression head on frozen protein language model embeddings. Splits are by 30% sequence identity cluster.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| eSOL (Niwa et al. 2009) as archived by the NBDC LSDB Archive, doi:10.18908/lsdba.nbdc00440-000 | [CC BY-SA 2.1 JP](https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en) | [esol.zip](https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip) | `29afdd2fa5be4947c9cbb1756197f4382e5184f8a2e4ae9162689ad74643d54b` |
| UniProtKB entries of the E. coli K-12 MG1655 reference proteome (UP000000625): accession, review status, ordered locus names, sequence | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [up000000625.tsv](https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000625&fields=accession,reviewed,gene_oln,sequence&format=tsv) | `e4bd795f6189b4816b64defe36cfe708603ef2c5a4f571baa56b515b8bbf5f99` |

eSOL is licensed CC BY-SA 2.1 JP by its maintainers (NBDC archive README). Its share-alike term requires derivatives to carry the same license, so this dataset is released under CC BY-SA 2.1 JP. Required attribution: "eSOL (c) Hideki Taguchi (Tokyo Institute of Technology) licensed under CC Attribution-Share Alike 2.1 Japan". Sequences are from UniProtKB (CC BY 4.0).

### Sources considered and excluded

- **NESG solubility scores (Price et al. 2011)**: The article (doi:10.1186/2042-5783-1-6) is CC BY, but its only supplementary file is a methods document with no per-protein data. The per-target expression and solubility scores used by later work were obtained from the authors on request and have no public release with a stated license.
- **PSI:Biology / TargetTrack**: The PSI Structural Biology Knowledgebase and TargetTrack are retired. No maintained public download with a stated redistribution license was found.

## Construction

1. Read `esol.csv` from the NBDC LSDB Archive copy of eSOL and keep rows with a numeric `Solubility (%)`.
2. Map each b-number to UniProtKB through the ordered-locus names of reference proteome UP000000625. A b-number with several candidate entries is resolved to the single reviewed entry, or dropped if that is not unique.
3. Merge rows that map to the same sequence: average when their values are within 10 percentage points, otherwise drop them all.
4. Cluster sequences with `mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 0 -s 7.5` and assign whole clusters to train/validation/test (80/10/10 of rows, seed 20261114).

## Schema

| Column | Type | Description |
|---|---|---|
| `b_number` | string | Blattner locus tag(s) in E. coli K-12 MG1655; ';'-joined when several genes encode the same sequence |
| `jw_id` | string | ASKA/GenoBase JW id(s) of the expressed clone |
| `gene_name` | string | Gene name(s) as given by eSOL |
| `uniprot_acc` | string | UniProtKB accession the sequence was taken from |
| `sequence` | string | UniProtKB canonical sequence (not the expressed construct; see limitations) |
| `solubility` | float64 | eSOL solubility in percent: supernatant / total protein after cell-free (PURE) synthesis without chaperones. Not clipped; about 7% of upstream values exceed 100 (maximum 147) |
| `n_measurements` | int32 | Number of eSOL rows averaged into solubility |
| `cluster_id` | string | MMseqs2 30% identity cluster (representative sequence hash); the split unit |

## Splits

| Split | proteins | clusters |
|---|---:|---:|
| train | 2501 | 2144 |
| validation | 313 | 284 |
| test | 313 | 255 |
| total | 3127 | |

### Cleaning statistics

| Step | Count |
|---|---:|
| eSOL rows | 4132 |
| dropped: no solubility value (not quantified) | 959 |
| dropped: malformed b-number | 0 |
| dropped: b-number not in UniProt UP000000625 | 40 |
| dropped: b-number maps to several UniProt entries | 0 |
| dropped: non-standard residue | 3 |
| sequences with more than one eSOL row | 2 |
| dropped: conflicting duplicates (spread > 10 points) | 1 groups / 2 rows |
| final rows with solubility > 100% | 226 |
| final rows | 3127 |
| MMseqs2 clusters | 2683 |

## Known limitations

- One organism (E. coli K-12) and one expression system (PURE, no chaperones). Solubility in vivo, with tags, or in other hosts can differ substantially.
- Sequences are UniProt canonical sequences. The expressed constructs came from the ASKA clone library and may differ (tags, start codon choice).
- NESG and PSI:Biology data were requested but are excluded because no public, licensed release exists; the dataset is therefore small.
- About 7% of eSOL values exceed 100% (up to 147%), reflecting measurement noise in the supernatant/total ratio. 
             They are left unclipped; many published solubility models clip to 100%, so do that explicitly if comparing.
- The share-alike license applies to any redistributed derivative of this dataset.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{niwa2009bimodal,
  title   = {Bimodal protein solubility distribution revealed by an aggregation analysis of the entire ensemble of Escherichia coli proteins},
  author  = {Niwa, Tatsuya and Ying, Bei-Wen and Saito, Katsuyo and Jin, WenZhen and Takada, Shoji and Ueda, Takuya and Taguchi, Hideki},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {106},
  number  = {11},
  pages   = {4201--4206},
  year    = {2009},
  doi     = {10.1073/pnas.0811922106}
}

@article{niwa2012chaperone,
  title   = {Global analysis of chaperone effects using a reconstituted cell-free translation system},
  author  = {Niwa, Tatsuya and Kanamori, Takashi and Ueda, Takuya and Taguchi, Hideki},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {109},
  number  = {23},
  pages   = {8937--8942},
  year    = {2012},
  doi     = {10.1073/pnas.1201380109}
}

@misc{esol_nbdc,
  title  = {eSOL, LSDB Archive},
  author = {Taguchi, Hideki},
  doi    = {10.18908/lsdba.nbdc00440-000}
}

@article{uniprot2025,
  title   = {UniProt: the Universal Protein Knowledgebase in 2025},
  author  = {{The UniProt Consortium}},
  journal = {Nucleic Acids Research},
  volume  = {53},
  number  = {D1},
  pages   = {D609--D617},
  year    = {2025},
  doi     = {10.1093/nar/gkae1010}
}
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
