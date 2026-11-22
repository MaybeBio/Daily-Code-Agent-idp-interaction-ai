---
license: cc-by-4.0
pretty_name: OpenDNA Antibody (human and mouse heavy-chain cohorts, OAS)
tags:
- biology
- protein
- antibody
- immunoglobulin
- repertoire
size_categories:
- 100K<n<1M
configs:
- config_name: default
  data_files:
  - split: human
    path: data/human.parquet
  - split: mouse
    path: data/mouse.parquet
  - split: diagnostic
    path: data/diagnostic.parquet
  - split: out_of_distribution
    path: data/out_of_distribution.parquet
---

# OpenDNA Antibody (human and mouse heavy-chain cohorts, OAS)

Heavy-chain variable regions from human and mouse antibody repertoires in the Observed Antibody Space, labelled by the species of the donor: 176543 rows from 37 donors across 31 studies. Two unlabelled cohorts travel with them -- mice carrying human immunoglobulin loci, where a humanness measure should answer human, and four species with one or two studies each.

These are evaluation cohorts, not a training split. There is no train/validation/test division, because the mouse donors cannot be divided: all 18 of them share identical rearrangements and form one group. Construction says why in full, and `donor_group` carries the grouping so the human side can still be split by anyone who needs to.

## Sources and licenses

| Source | License | File | SHA-256 |
|---|---|---|---|
| Observed Antibody Space, unpaired heavy-chain data units (Olsen, Boyles and Deane 2021); the units read are listed in oas_units.tsv | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [oas_units.tsv](https://opig.stats.ox.ac.uk/webapps/oas/) | (not locked) |

OAS is licensed CC BY 4.0 by its maintainers, stated on <https://opig.stats.ox.ac.uk/webapps/oas/>: "the data contained within OAS is available under a CC-BY 4.0 license". Redistribution is permitted with attribution, so this dataset is CC BY 4.0. Cite Olsen, Boyles and Deane (2021) and the study each row's `study` column names.

### Sources considered and excluded

- **SAbDab / SAbDab-nano**: Not used, and so not relied on. SAbDab holds antibody structures and this is a sequence dataset, so it would add nothing a heavy-chain V region needs. Its licence was therefore not reviewed: third-party pages and its own papers read as CC BY 4.0, but the OPIG site states no terms on a page that renders without JavaScript, and this project does not record a licence it has not seen the maintainers state. If a structural antibody dataset is built later, that review has to happen first.

## Construction

1. Read the 54 data units listed in `oas_units.tsv`, which is committed to the repository and compiled into the builder. OAS has no stable URL for its unit index, so this file is the catalogue; its own SHA-256 is in the statistics above and each unit's is in `manifest.json`. A unit whose digest does not match the catalogue stops the build.
2. Keep rows OAS already marked productive, in frame and free of stop codons, on the heavy chain, whose amino-acid alignment is 90-200 standard residues.
3. Drop any sequence seen under more than one donor, and cap each donor at 5000 sequences chosen by sequence hash.
4. Write four cohorts. `human` and `mouse` carry `is_human`; HIS-Mouse and Kymouse carry human immunoglobulin loci and go to `diagnostic` with no label, where answering "human" is correct; rat, rabbit, rhesus and bactrian camel go to `out_of_distribution`, also unlabelled.
5. Join donors wherever they share an identical rearrangement -- the same V gene, the same J gene and a byte-identical CDR3 -- and record the resulting group in `donor_group`. That is the smallest unit that can be held out without carrying a lineage across the boundary.

**There is no train/validation/test split, because this data does not admit one.** All 18 mouse donors fall into a single donor group, spanning three inbred strains and nine studies, because C57BL/6 and BALB/c animals are genetically near-identical and their repertoires share public rearrangements. One laboratory mouse is not an independent sample of its species the way one person is. The largest group therefore holds 46.2% of the labelled rows, and any grouped split leaves the whole mouse cohort on one side -- which is not a split.

Sequence clustering does not rescue it: clustering the amino-acid V(D)J region at 90% identity and joining donors through shared clusters put all 37 labelled donors into one component. That region is mostly germline-constant, so the threshold groups V-gene families rather than lineages. The shared-rearrangement grouping is the strict alternative, and `donor_group` records it: the human side can be split on it (15 groups), the mouse side cannot (1).

## Schema

| Column | Type | Description |
|---|---|---|
| `sequence` | string | Heavy-chain variable region, amino acids, as OAS aligned it (`sequence_alignment_aa`) |
| `organism` | string | OAS `Species` verbatim, including the strain. Names the animal, not the immunoglobulin: HIS-Mouse and Kymouse are mice with human V genes |
| `study` | string | OAS study the unit belongs to |
| `donor` | string | Donor key. `person:<subject>` when the unit names a subject, which is global because one donor is deposited by several studies; `study:<study>` when it does not |
| `isotype` | string | OAS `Isotype` as rendered to a string |
| `v_call` | string | IMGT V gene call from OAS. Reported so the germline overlap between cohorts can be checked; not a feature |
| `j_call` | string | IMGT J gene call from OAS |
| `cdr3` | string | CDR3 amino acids as OAS delimited them |
| `redundancy` | int32 | Identical reads OAS collapsed into this row |
| `is_human` | bool | True if the sequence came from a human repertoire, false if from a mouse one. The study that deposited it recorded this; it is not a measurement of the molecule and not an immunogenicity label |
| `donor_group` | string | Donors joined where they share an identical rearrangement (same V gene, same J gene, byte-identical CDR3). The smallest unit that can be held out without leaking a lineage. All mouse donors fall in one group, which is why this dataset has no split |

## Cohorts

| Cohort | sequences |
|---|---:|
| human | 94923 |
| mouse | 81620 |
| diagnostic | 30000 |
| out_of_distribution | 40000 |
| total | 246543 |

There is no train/validation/test split. Why not is under Construction, and it is a fact about the data rather than a choice: read it before splitting these cohorts yourself.

### Cleaning statistics

| Step | Count |
|---|---:|
| OAS unit catalogue (oas_units.tsv) sha256 | 083f7ba3a3fe00c47033a72644581ac148e526e0af94e518cb821a27af5bce4b |
| data units read | 54 |
| data units not pinned in the catalogue | 0 |
| rows read from units | 3196883 |
| dropped: not productive | 10828 |
| dropped: stop codon | 0 |
| dropped: V-J out of frame | 0 |
| dropped: not a heavy chain | 0 |
| dropped: empty amino-acid alignment | 0 |
| dropped: non-standard residue or alignment gap | 139505 |
| dropped: outside 90-200 residues | 185566 |
| dropped: duplicate within its unit | 0 |
| dropped: identical sequence under more than one donor | 28027 |
| dropped: over the 5000-per-donor cap | 2586414 |
| labelled rows (human and mouse cohorts) | 176543 |
|   human | 94923 rows, 19 donors, 19 studies |
|   mouse | 81620 rows, 18 donors, 12 studies |
| diagnostic rows (mice with human Ig loci, unlabelled) | 30000 rows, 6 donors |
| out-of-distribution rows (other species, unlabelled) | 40000 rows, 8 donors |
| donor groups (donors joined by a shared rearrangement) | 16 |
|   among human donors | 15 |
|   among mouse donors | 1 |
| largest donor group, share of labelled rows | 46.23% |
| V genes in human rows | 56 |
| V genes in mouse rows | 169 |
| V genes in both | 4 |

## Known limitations

- The label is the species of the donor, recorded by the depositing study. It is not a measurement of the molecule, and it is not an immunogenicity label: a sequence resembling a human repertoire is not thereby safe in humans, and the converse does not hold either.
- Human and mouse V genes differ almost completely -- the card reports how many appear in both cohorts -- so separating the two species is largely a germline lookup. Measured on held-out rows, a lookup table on the V gene call alone reaches 98.5% against a 68.9% majority baseline. An accuracy on this task is therefore not evidence that a model learned anything about antibodies.
- Heavy chains only. Light chains have their own germline repertoire and are a separate task.
- The out-of-distribution cohort is one study per species for rat, rabbit and bactrian camel, and two for rhesus. A score on those species cannot be separated from one laboratory's protocol, which is why those rows carry no label.
- These are cohorts, not folds. There is no train/test split and none can be made: all 18 mouse donors share identical rearrangements and form a single group, so any grouped split puts the whole mouse cohort on one side. The `donor_group` column carries the grouping; the human cohort has 15 groups and can be divided, the mouse cohort has one and cannot.
- Donor keys merge a named subject across studies, because one donor is deposited by several studies -- CAP256 appears in three. Where two studies reuse a placeholder such as Donor-1 for different people, this merges them too, which costs data rather than leaking it.
- Repertoires are not samples of what a B cell can make. They are shaped by the primers used, the cells sorted, the immunisation, and the disease state of the donor, all of which are in the columns and none of which is controlled for here.
- A cap of 5,000 sequences per donor keeps deeply sequenced donors from dominating. It also means row counts per donor say nothing about repertoire size.

## Citation

If you use this dataset, cite the upstream sources:

```bibtex
@article{olsen2022oas,
  title   = {Observed Antibody Space: A diverse database of cleaned, annotated,
             and translated unpaired and paired antibody sequences},
  author  = {Olsen, Tobias H. and Boyles, Fergus and Deane, Charlotte M.},
  journal = {Protein Science},
  volume  = {31},
  number  = {1},
  pages   = {141--146},
  year    = {2022},
  doi     = {10.1002/pro.4205}
}
```

## Reproduction

Built with `opendna-data` from the OpenDNA repository (`datasets/`). The exact upstream files are pinned by SHA-256 in `datasets/sources.lock`; `manifest.json` next to this card records the pipeline parameters and the checksum of every Parquet file.
