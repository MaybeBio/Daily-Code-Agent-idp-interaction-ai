---
license: cc-by-4.0
pretty_name: OpenDNA Burial (relative solvent accessibility from CATH domains)
tags:
- biology
- protein
- solvent-accessibility
- structure
size_categories:
- n<1K
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

# OpenDNA Burial (relative solvent accessibility from CATH domains)

Relative solvent accessibility for every residue of 939 CATH v4.4.0 domains, computed from deposited coordinates by the Shrake-Rupley method over the theoretical maxima of Tien et al. (2013). 139794 labelled residues. Splits are inherited from the OpenDNA Structures dataset, so no CATH topology and no 30% identity cluster spans two splits.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| CATH v4.4.0 domain classification, S100 non-redundant representatives | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [cath-domain-list-S100-v4_4_0.txt](https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-S100-v4_4_0.txt) | `9a202046a68c80436f76a3631fe0f84fc4f027a89c53657a5dcd4c97231cf562` |
| CATH v4.4.0 domain sequences (from ATOM records), S100 representatives | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [cath-domain-seqs-S100-v4_4_0.fa](https://download.cathdb.info/cath/releases/all-releases/v4_4_0/sequence-data/cath-domain-seqs-S100-v4_4_0.fa) | `d7d1efdf9a169456ad762faca0dc80b64f734ddeaf31481c69eb68ce5a77cd83` |
| wwPDB entries index: deposition date, title, resolution and experimental method per entry | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | [entries.idx](https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx) | `b9d6b90c4b1c4bdfbccb61210276e199aecfdefe73ed754fb6792468e174e55f` |

CATH is CC BY 4.0 and wwPDB coordinate data is CC0 1.0; the derived dataset is released under CC BY 4.0.

## Construction

1. Read the OpenDNA Structures dataset, which indexes CATH v4.4.0 S100 domains with their PDB entry, chain and domain boundaries.
2. Sample 1200 train, 150 validation and 150 test domains, shuffled with seed 20261114 so the sample is not the alphabetically-first and therefore oldest depositions. The sample is drawn within the existing splits, so their no-shared-fold guarantee is inherited.
3. Fetch each domain's coordinate file from RCSB.
4. Compute solvent-accessible area for the domain's whole chain (Shrake-Rupley, 1.4 A probe, 960 points per atom) and divide by the per-residue theoretical maximum.
5. Walk the domain's segments in author numbering, checking each residue name against the expected letter of the CATH sequence. A disagreement drops the domain rather than labelling it from whatever is at those positions.

## Schema

| Column | Type | Description |
|---|---|---|
| `domain_id` | string | CATH domain id (PDB id, chain, domain number) |
| `pdb_id` | string | PDB entry the coordinates came from |
| `chain_id` | string | Chain accessibility was computed on |
| `sequence` | string | Domain sequence, as in the structures dataset |
| `length` | int32 | Sequence length |
| `rsa` | list<float32> | Relative solvent accessibility per residue: Shrake-Rupley area over the theoretical maximum of Tien et al. (2013), computed on the whole chain |
| `buried` | list<int8> | 1 where rsa <= 0.25, the Rost & Sander convention; the rsa column is the real label |

## Splits

| Split | domains |
|---|---:|
| train | 751 |
| validation | 95 |
| test | 93 |
| total | 939 |

Each of these domains is its own split unit, so the no-shared-fold guarantee comes from the grouping described under Construction rather than from grouping rows together here.

### Cleaning statistics

| Step | Count |
|---|---:|
| dropped: domain sequence is not contiguous in the chain | 561 |
| domains annotated | 939 |
| residues labelled | 139794 |
| residues buried (rsa <= 0.25) | 73742 (52.8%) |
| coordinate files read | 1480 |

## Known limitations

- Accessibility is computed on the domain's chain alone. Residues that form an inter-chain interface in an obligate oligomer are therefore labelled exposed, because the partner chain is not present.
- Coordinates are a crystal or NMR model, not the protein in solution. Flexible surface side chains are placed as the depositor modelled them, and unobserved residues are simply absent from the chain.
- This is a sample, not the whole of CATH: 72,522 PDB entries would be hundreds of gigabytes and an unreasonable load on RCSB. Coverage of rare folds is therefore thinner than the structures dataset's.
- Discontinuous domains are excluded. A domain is placed by finding its sequence as a contiguous stretch of its chain, which cannot locate the roughly one in eight CATH domains assembled from separated segments. The kept set is biased toward single-segment domains.
- CATH v4.4.0 dates from 2019 and wwPDB entries are revised, so a domain's sequence sometimes covers residues the current coordinates no longer observe. Where the observed chain lies inside the domain sequence the overlap is labelled and the `sequence` column records what was labelled rather than what CATH listed; where neither contains the other the domain is dropped. This biases the kept set toward entries that have not been substantially revised. Together with the previous point this rejects roughly a third of the sample, counted by reason in the statistics above.
- Domains are aligned by sequence rather than by residue numbering, so insertion codes and renumbered chains do not matter, but a domain whose sequence occurs more than once in its chain is placed at the first occurrence. At these lengths a repeat is rare, and two copies of one domain in a chain would have similar burial.
- Relative accessibility can exceed 1 for a residue whose side chain is unusually extended, because the Tien et al. maxima are theoretical. Values are not clipped.
- The buried column uses a 0.25 threshold by convention. It is a convenience; the rsa column is the measurement.
- Individual coordinate files are pinned by SHA-256 in this dataset's manifest rather than in sources.lock, so a re-run uses the same bytes, but a wwPDB entry revised upstream will be noticed only on a cache miss.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{shrake1973environment,
  title   = {Environment and exposure to solvent of protein atoms. Lysozyme and insulin},
  author  = {Shrake, A. and Rupley, J. A.},
  journal = {Journal of Molecular Biology},
  volume  = {79},
  number  = {2},
  pages   = {351--371},
  year    = {1973},
  doi     = {10.1016/0022-2836(73)90011-9}
}

@article{tien2013maximum,
  title   = {Maximum allowed solvent accessibilites of residues in proteins},
  author  = {Tien, Matthew Z. and Meyer, Austin G. and Sydykova, Dariya K. and
             Spielman, Stephanie J. and Wilke, Claus O.},
  journal = {PLOS ONE},
  volume  = {8},
  number  = {11},
  pages   = {e80635},
  year    = {2013},
  doi     = {10.1371/journal.pone.0080635}
}

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
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
