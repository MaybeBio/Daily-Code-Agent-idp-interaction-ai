---
license: cc-by-4.0
pretty_name: OpenDNA peptide humanness prevalence table
tags:
- protein
- antibody
- humanness
- immunogenicity
---

# OpenDNA peptide humanness prevalence table

How many human donors hold each 9-residue peptide, over 141 donors from 47 studies of the Observed Antibody Space. Cut an antibody sequence into overlapping 9-mers, look each one up, and every window comes back with the fraction of donors whose repertoire contains it. The per-window list is the useful output: it says *where* a sequence is unlike anything human, which a single score cannot.


## Read this before using the numbers

Each donor contributes one shallow data unit: between 1153 and 7920 sequences survive reading, out of a repertoire of millions. A peptide that donor really has, but that did not appear in the sample, is counted as absent. **Every prevalence here is biased downwards**, and unevenly: the deepest donor contributes 6.9 times as many sequences as the shallowest, so its peptides are the more likely to be seen at all.

Scores from this table are comparable with each other. They are **not** comparable with published OASis values: that work used deep repertoires from 266 subjects, so the same sequence will score lower here. The choice was deliberate -- prevalence is a fraction of donors, so donor breadth decides what thresholds mean at all -- but it is a trade, not a free lunch.

## Sources and licenses

| Source | License | Index | SHA-256 of unit list |
|---|---|---|---|
| Observed Antibody Space, unpaired human heavy-chain data units (Olsen, Boyles and Deane 2021); the units read are listed in `oas_human_donors.tsv` | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [OAS](https://opig.stats.ox.ac.uk/webapps/oas/) | `4225d0a611c07957cf9db5bcf15d82c2b26fc78fe67afc58e41292692c26fb10` |

OAS publishes its data units under CC BY 4.0. This table is a derived work -- counts of peptides, not sequences -- and is released under the same license, with attribution to OAS and to the studies whose repertoires it was built from. Every unit's SHA-256 is recorded in `manifest.json`; the exact list, with donor keys, is `oas_human_donors.tsv` in the pipeline source.

## Construction

1. One OAS unpaired heavy-chain unit per donor, 150 of them, picked round-robin across studies so no single study dominates, from a fixed size band and for breadth of donors rather than depth per donor. 141 of those donors are counted, across 47 studies; the others contributed too little to be evidence, and are named below.
2. Rows kept only where OAS marks them productive, stop-codon-free, in frame and heavy-locus, with an amino-acid V(D)J alignment of 90-200 standard residues. 895955 rows were read and 729494 kept; the drops are itemised below.
3. Sequences grouped by donor key, which is the subject identifier where OAS gives a global one and the study otherwise. Two units from one person would be one repertoire, not two donors.
4. Every overlapping 9-mer of every sequence, counted **once per donor** however many of that donor's sequences contain it. Counting sequences instead would let one deeply sequenced donor make their own peptides look universal.
5. Peptides held by fewer than 2 donors dropped: 6600894 of 7696982 distinct peptides. One donor in 141 is 0.71%, below the loosest threshold anyone scores at, so this changes no threshold answer -- but it is most of the table by count, and the reported prevalence of such a peptide becomes zero rather than 0.71%.

No donor in this table appears in the antibody dataset's catalogue (`oas_units.tsv`), and the build fails if one does. The two donor sets have to be disjoint or an evaluation measures how well the table remembers its own input.

9 donor(s) listed in the catalogue contributed fewer than 500 sequences and are **not** among the 141 counted: person:4B (2), person:4G (0), person:4H (3), person:4K (0), person:Donor-163 (352), person:Patient-1 (102), person:Patient-2 (78), person:Patient-6 (104), person:Patient-8 (109). A donor who contributed almost nothing votes "absent" on almost every peptide, so counting them would divide every prevalence by donors who could not have held it. The floor is not tuned: depth per donor here is bimodal, nine donors fall below 500 and the next has 1,153, so any floor between those two gives the same set.

### Rows read and why they were dropped

| Reason | Rows |
|---|---|
| read | 895955 |
| not productive | 0 |
| stop codon | 0 |
| out of frame | 0 |
| wrong locus | 0 |
| empty sequence | 0 |
| non-standard residue | 47918 |
| outside the length band | 118543 |
| duplicate within the unit | 0 |
| **kept** | 729494 |

## Files

`table.bin` -- 14249180 bytes, SHA-256 `793ac490b8c9754c57089d331eeedafd43be87111f2118c89efd9c9bd69542bd`. Format 1: a 16-byte magic, then format version, peptide length, donor total and entry count as little-endian integers, then one entry per peptide of 9 ASCII residues and a little-endian `u32` donor count, sorted by peptide. Readers refuse a file whose entries are unsorted or whose counts exceed the donor total, because either would produce wrong scores silently rather than failing.

`manifest.json` -- the SHA-256 of every unit read, the catalogue digest, and the counts above.


## What the thresholds can mean

From 141 donors a prevalence takes one of 142 values, so the smallest non-zero prevalence expressible is 0.71%. A threshold finer than that cannot be evaluated at all, however it is reported.

| Stringency | Expressible here | Peptides at or above |
|---|---|---|
| loose (>=1% of donors) | yes | 1096088 |
| medium (>=10% of donors) | yes | 87173 |
| strict (>=50% of donors) | yes | 8512 |

## Limitations

- Prevalence is biased downwards by shallow sampling, and the bias is uneven across donors. Absolute values are not comparable with OASis or with any deep-repertoire table.
- Heavy chains only. Nothing here says anything about a light chain, and a paired antibody scored window-by-window on its heavy chain alone is half measured.
- The donors are whoever OAS has, which is not a sample of humanity: the studies are overwhelmingly from North America and Europe, and several are patient cohorts rather than healthy adults. A peptide common in an unrepresented population reads as rare here.
- A high prevalence means a peptide occurs in human repertoires, not that an antibody carrying it is safe, non-immunogenic or developable. Immunogenicity depends on MHC presentation, dose, route and the recipient, none of which this measures.
- A low prevalence means the peptide was not seen in these donors' sampled sequences. It is not evidence that no human has it.
- The measure is a count over a fixed table, not a model. It will not generalise to a peptide length other than nine or to a different table without being rebuilt, and nothing about it is learned, so there is no accuracy figure to report.

## Citations

The data are OAS; the measure follows OASis, published as
part of BioPhi. This table is built from different, shallower repertoires than
OASis and its values are not interchangeable with those published there.

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

@article{prihoda2022biophi,
  title   = {BioPhi: A platform for antibody design, humanization, and
             humanness evaluation based on natural antibody repertoires and
             deep learning},
  author  = {Prihoda, David and Maamary, Jad and Waight, Andrew and Juan,
             Veronica and Fayadat-Dilman, Laurence and Svozil, Daniel and
             Bitton, Danny A.},
  journal = {mAbs},
  volume  = {14},
  number  = {1},
  pages   = {2020203},
  year    = {2022},
  doi     = {10.1080/19420862.2021.2020203}
}
```
