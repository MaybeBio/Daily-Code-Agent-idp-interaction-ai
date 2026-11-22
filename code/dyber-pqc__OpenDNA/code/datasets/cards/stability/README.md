---
license: cc-by-4.0
pretty_name: OpenDNA Stability (Tsuboyama 2023 single substitutions)
tags:
- biology
- protein
- protein-stability
- ddg
- mutation-effect
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

# OpenDNA Stability (Tsuboyama 2023 single substitutions)

Single-substitution folding stability changes (ddG, kcal/mol) for 412 small natural and designed protein domains, measured by cDNA display proteolysis (Tsuboyama et al. 2023). Prepared for training a ddG regression head on frozen protein language model embeddings. Positive ddG means the substitution is stabilising. Splits are by 30% sequence identity cluster, so test domains have no close homolog in train.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| Tsuboyama et al. 2023 mega-scale cDNA display proteolysis, processed dG/ddG tables (Zenodo 10.5281/zenodo.7992926) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [Processed_K50_dG_datasets.zip](https://zenodo.org/api/records/7992926/files/Processed_K50_dG_datasets.zip/content) | `e1c52e14d27b18b846cd87c719eb3f3adae6f736e700962bb7da82fd031488e9` |

The only included source is licensed CC BY 4.0, and this derived dataset is released under CC BY 4.0. Attribution: Tsuboyama et al. 2023 (citation below).

### Sources considered and excluded

- **FireProtDB 2.0**: The FireProtDB terms of use (https://loschmidt.chemi.muni.cz/fireprotdb/terms-of-use/, checked 2026-11-14) grant a non-exclusive, non-transferable license and prohibit sublicensing or transferring rights to the service output and creating derivative works. That does not permit redistribution on Hugging Face. Third-party mirrors label it CC BY 4.0, but that is not the maintainers' statement. FireProtDB 2.0 also already incorporates the Tsuboyama 2023 measurements (its MEGASCALE columns), so including both would double-count them.

## Construction

1. Read `Processed_K50_dG_datasets/Tsuboyama2023_Dataset2_Dataset3_20230416.csv` from the Zenodo archive.
2. Keep rows whose `mut_type` is a single substitution of standard amino acids. Double mutants, insertions and deletions are dropped.
3. Keep rows in Tsuboyama's ML-grade subset (Dataset 3): numeric `ddG_ML` and `dG_ML`. Rows whose dG was outside the assay's dynamic range (`<-1`, `>5`) or was flagged unreliable (`-`) are dropped.
4. Check the sign convention against the source's own `Stabilizing_mut` flag (the build aborts on a contradiction), and convert to kcal/mol with positive = stabilising (a no-op for this source).
5. Verify that each mutant sequence equals the wild-type sequence with the substitution applied; drop rows that do not.
6. Merge measurements of the same substitution in the same wild-type sequence: average when their spread is at most 1 kcal/mol, otherwise drop all of them.
7. Cluster the distinct wild-type sequences with `mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 0 -s 7.5`.
8. Assign whole clusters to train/validation/test (target 80/10/10 of rows) with a seeded deterministic procedure (seed 20261114). All mutations of one domain, and of all its homologs at 30% identity, are in the same split.

## Schema

| Column | Type | Description |
|---|---|---|
| `protein_id` | string | Tsuboyama 2023 wild-type name (WT_name); for merged duplicates, the lexicographically smallest |
| `wt_sequence` | string | Wild-type domain sequence (aa_seq, without the display linkers) |
| `mutation` | string | Substitution as <wt><1-based position><mut>, e.g. D12Q |
| `position` | int32 | 1-based position in wt_sequence |
| `wt_aa` | string | Wild-type residue (one-letter code) |
| `mut_aa` | string | Mutant residue (one-letter code) |
| `ddg` | float64 | dG_unfold(mutant) - dG_unfold(wild type) in kcal/mol; positive = stabilising. Mean over consistent replicates |
| `n_measurements` | int32 | Number of source measurements averaged into ddg |
| `ddg_spread` | float64 | max - min of the averaged measurements in kcal/mol (0 for a single measurement) |
| `source` | string | Upstream source id |
| `cluster_id` | string | MMseqs2 30% identity cluster of the wild-type sequence (representative sequence hash); the split unit |

## Splits

| Split | mutations | clusters |
|---|---:|---:|
| train | 314716 | 123 |
| validation | 37983 | 18 |
| test | 36369 | 24 |
| total | 389068 | |

### Cleaning statistics

| Step | Count |
|---|---:|
| CSV rows read | 776298 |
| wild-type rows | 2348 |
| dropped: not a single substitution (double mutants, insertions, deletions) | 286036 |
| dropped: outside Tsuboyama's ML-grade subset (ddG_ML or dG_ML not numeric) | 98846 |
| dropped: synonymous | 0 |
| dropped: wild type missing or ambiguous | 0 |
| dropped: non-standard residue in wild type | 0 |
| dropped: mutant sequence inconsistent with substitution | 0 |
| measurements kept | 389068 |
| duplicate substitution groups | 0 |
| dropped: conflicting duplicate groups (spread > 1 kcal/mol) | 0 groups / 0 measurements |
| final rows | 389068 |
| distinct wild-type sequences | 412 |
| MMseqs2 clusters | 165 |

## Known limitations

- ddG is inferred from protease resistance through a kinetic model, not measured calorimetrically. The assay's dynamic range is roughly -1 to 5 kcal/mol for dG; mutations that push dG outside it are excluded, which truncates the distribution of strongly destabilising mutations.
- Wild-type domains are 32-74 residues long and many are de novo designs. Performance on large, multi-domain or membrane proteins cannot be inferred from this data.
- Only single substitutions are included.
- FireProtDB was requested as a second source but is excluded for licensing reasons (see above), so all values come from one assay in one laboratory; assay-specific biases are not averaged out.
- Clustering uses 30% identity over 80% bidirectional coverage. Homology below that threshold, or shared structural motifs, can still connect train and test.
- The duplicate-conflict threshold (1 kcal/mol) is a pipeline choice, not a property of the data.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{tsuboyama2023megascale,
  title   = {Mega-scale experimental analysis of protein folding stability in biology and design},
  author  = {Tsuboyama, Kotaro and Dauparas, Justas and Chen, Jonathan and Laine, Elodie and
             Mohseni Behbahani, Yasser and Weinstein, Jonathan J. and Mangan, Niall M. and
             Ovchinnikov, Sergey and Rocklin, Gabriel J.},
  journal = {Nature},
  volume  = {620},
  pages   = {434--444},
  year    = {2023},
  doi     = {10.1038/s41586-023-06328-6}
}

@dataset{tsuboyama2023zenodo,
  title     = {Mega-scale experimental analysis of protein folding stability in biology and design},
  author    = {Tsuboyama, Kotaro and others},
  publisher = {Zenodo},
  year      = {2023},
  doi       = {10.5281/zenodo.7992926}
}
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
