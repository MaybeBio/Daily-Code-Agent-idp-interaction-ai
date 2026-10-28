# Coverage of known structure space and of BioEmu-1's training data by the SAXS set

This analysis places the 399 clean PeptoneDB-SAXS entries with a BioEmu-1 ensemble in experimentally known
protein structure space. It maps every construct to a UniProtKB region, annotates the region with CATH
superfamilies, organism, disorder and PDB coverage, and compares the set with three references: CATH v4.4.0
at the S35 sequence-cluster level (fold space), the non-redundant PDB (RCSB 30%-identity polymer-entity
clusters, experimental entries only) and UniProtKB/Swiss-Prot (organism, length, disorder). Scripts 11-17
then compare the set with BioEmu-1's own training data and find, for every construct, the nearest training
sequence. Its small outputs are tracked in `results/coverage/`; the `bsc coverage` stage and the project page read only those.

## Steps

| Script | Does |
|---|---|
| `00_download_cath.sh` | CATH latest-release classification files (v4.4.0) and the v4.3.0 S35 list for a sensitivity check |
| `01_fetch_sasbdb.py` | SASBDB entry summaries for the 439 labels |
| `02_map_uniprot.py` | construct → UniProtKB accession: the SASBDB `uniprot_code` where it aligns (identity ≥ 0.5, construct coverage ≥ 0.3), else an RCSB sequence search (≥ 90% identity) with the SIFTS UniProt reference of the hit, else a UniProt peptide search on three internal 20-mers; then a local alignment to the canonical sequence defines the construct region |
| `03_download_reference_sets.sh` | the Swiss-Prot flat file and the RCSB 30%-identity cluster file |
| `04_fetch_annotations.py` | per accession: InterPro CATH-Gene3D and Pfam locations, PDBe SIFTS best structures, TED domains |
| `05_fetch_rcsb_nr.py` | length and source organism of one experimental representative per RCSB 30% cluster |
| `06_parse_swissprot.py` | one row per Swiss-Prot entry: length, organism group, UniProt "Disordered" regions, Gene3D matches, PDB cross-references |
| `07_build_protein_table.py` | the per-entry table: domains whose residues lie at least half inside the construct, SIFTS coverage of the construct by atomic-model chains (`has_pdb_50`: at least 50%), disorder content, organism group |
| `08_compare.py` | all comparisons into `summary.json`: superfamily coverage by rank, class and architecture composition (cluster bootstrap over accessions), over- and under-represented superfamilies (Poisson tails, uncorrected, for ranking), organism, length, domains per protein, disorder, structure availability, redundancy |
| `09_figures.py` | diagnostic figures in the work directory (the page figures are drawn by `bsc coverage`) |
| `10_export.py` | copies the tracked outputs into `results/coverage/` |
| `11_download_bioemu_training.py` | BioEmu-1's training data at the size needed: the AFDB pre-training cluster definitions from github.com/microsoft/bioemu (pinned commit), the Zenodo dataset cards and record metadata, the member listings of the MD archives and the `topology.pdb` of every CATH1, CATH2 and MegaSim wild-type system (read out of the remote zips through HTTP range requests, no trajectories), the curated MEGAscale table, the DESRES fast-folder PDB sequences and the named proteins of the other public MD sets |
| `12_fetch_afdb_sequences.py` | sequence, organism and CATH-Gene3D cross-references of every AFDB cluster member (UniProt ID mapping; entries deleted since AFDB v4 from UniParc) |
| `13_training_sequences.py` | one table of training sequences by component, and a FASTA per component |
| `14_annotate_training.py` | the same three descriptors for every set: CATH superfamilies (Gene3D for the SAXS constructs and AFDB, the domain identifier for the CATH MD systems, the best MMseqs2 hit against CATH S100 for the MEGAscale domains), predicted disorder (metapredict, also for the SAXS constructs and the CATH S35 sequences) and organism group |
| `15_nearest_training.py` | MMseqs2 search of all 439 constructs against each training component: the nearest training sequence covering at least half of the construct, superfamily presence and accession overlap |
| `16_compare_training.py` | composition of the training sets against the SAXS set and CATH (class, superfamilies, length, disorder, organism) and the overlap distribution, into `training_summary.json`; the MEGAscale domains called mostly disordered, by origin (natural, designed) and length |
| `17_export_training.py` | copies the per-construct nearest-neighbour table, the training summary and the versions into `results/coverage/` |
| `18_sasbdb_metadata.py` | one row per entry with the SASBDB record's curve type, size-exclusion coupling, source type, instrument, cell temperature, buffer, pH, additives, concentration range, organisms and molecule names, from the summaries cached by script 01; exported to `results/coverage/sasbdb_metadata.csv` (used by `analysis/arm1/`) |

Units: the domain-level comparisons count unique (UniProt accession, CATH superfamily) pairs inside the
constructs, the protein-family analogue of a CATH S35 cluster. Single proportions carry Wilson 95% intervals.

## BioEmu-1's training data

BioEmu-1 (Lewis et al. 2025) was pre-trained on about 50,000 AFDB sequence clusters, fine-tuned on MD
simulations and then on MEGAscale folding stabilities. Scripts 11-17 obtain what is public of each part and
compare it with the SAXS set:

| Component | Obtained as | Used for |
|---|---|---|
| AFDB pre-training clusters | cluster definitions (50,228 clusters, 179,422 AFDB v4 structure representatives) from the BioEmu repository; sequences from UniProtKB or UniParc | composition (one input sequence per cluster, the member with the highest mean pLDDT) and nearest neighbours (every member) |
| CATH1 and CATH2 MD (Zenodo 15629740) | the simulated topology of each system; CATH domain identifier from the system name | composition and nearest neighbours |
| MegaSim wild types (Zenodo 15641184) | the simulated topology of each system | nearest neighbours; composition together with MEGAscale |
| MegaSim point mutants (Zenodo 15641184) | archive listing only (system count) | counted; each is one substitution or indel away from a wild type |
| Octapeptides (Zenodo 15641199) | archive listing only (system count) | counted; eight residues cannot cover half of any construct |
| MEGAscale (Tsuboyama et al. 2023, Zenodo 7992926) | the curated Dataset2/Dataset3 table; its wild types | composition and nearest neighbours; the 361 wild types used for fine-tuning are a subset the paper does not list |
| Other public MD sets of Table S1 (DDR1, SETD8, coronavirus proteins, MHC class II, barnase-barstar) | the named proteins, full UniProtKB chains or polyprotein chains | nearest neighbours (an upper bound: the simulated constructs can be shorter) |
| DESRES fast folders | the twelve PDB sequences (the trajectories are not redistributable) | reported separately: Lewis et al. used them only for a separate model |

The training data were filtered to under 40% sequence identity to the paper's own test proteins; the SAXS
proteins were not part of that filter. The full MD archives are not downloaded: the scripts read the zip
central directories and single members through HTTP range requests; the training-data downloads take about
0.4 GB, most of it the MEGAscale table. The simulated topologies agree with the CATH domain sequences: of the 576
CATH2 systems whose domain is in the CATH v4.4 S100 set, 575 have the same sequence or one contained in the other.

## Mapping notes

390 of the 399 clean entries map to a UniProtKB entry (353 by the SASBDB accession, 23 by RCSB sequence search,
14 by UniProt peptide search). The SASBDB accession was set aside for 17 entries: 16 pointed to UniProtKB
entries deleted since deposition and one aligned to 6% of the construct; 16 of the 17 were mapped again by
sequence search. Nine entries stay unmapped: four antibody-derived fragments (a Fab, an scFv and two
single-domain antibodies), a synthetic construct, a filamin two-domain fragment, an antimicrobial peptide, a
bacterial harpin and a lanthionine synthetase C-like protein whose UniProtKB entry was deleted. They are left
out of every comparison that needs a mapping. SIFTS residue ranges follow the deposited sequence, so PDB coverage counts residues listed in a
structure record, which can include residues not modelled in the coordinates.

## Database versions

Recorded in `results/coverage/versions.json` when each source was first fetched (October 2026):

| Source | Version |
|---|---|
| CATH | v4.4.0 (601,328 domains, 37,350 S35 clusters) |
| CATH-Gene3D in InterPro | 4.3.0 (InterPro 110.0) |
| Pfam | 38.2 |
| UniProtKB / Swiss-Prot | release 2026_03 (575,748 Swiss-Prot entries) |
| RCSB clusters-by-entity-30 | file of 4 October 2026 (47,069 clusters with an experimental entity) |
| PDBe SIFTS best_structures | live API, weekly update |
| TED | v1 (AFDB v4 models with CATH labels) |
| SASBDB | REST API entry summaries |
| BioEmu-1 AFDB cluster definitions | github.com/microsoft/bioemu, file added 2025-04-25 (commit recorded in `versions.json`) |
| BioEmu-1 MD datasets | Zenodo 15629740 (CATH1, CATH2), 15641199 (octapeptides), 15641184 (MegaSim), CDLA-Permissive-2.0 |
| MEGAscale | Zenodo 7992926, version v2_230420, CC-BY 4.0 |
| UniParc | for the 129,202 AFDB cluster members deleted from UniProtKB since AFDB v4 |
| metapredict | 3.1.1 |
| MMseqs2 | release 18-8cc5c |

The live APIs (SIFTS, InterPro, TED, RCSB search) move between releases, so a rerun after deleting the download
cache can differ slightly from the tracked outputs; with the cache in place the rerun reproduces them byte for
byte.

## Rerun

```bash
uv venv /path/to/venv-cov && uv pip install --python /path/to/venv-cov -r analysis/coverage/requirements.txt \
  --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match
curl -sSfL -o mmseqs.tar.gz https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-avx2.tar.gz
COVERAGE_PY=/path/to/venv-cov/bin/python BSC_ROOT=<where bsc fetch put data/> MMSEQS=/path/to/mmseqs/bin/mmseqs \
  COVERAGE_WORK=/path/with/2GB/free analysis/coverage/run_all.sh
uv run bsc coverage      # structure split, closeness to training data, categories and the page figures
```

Scripts 14 and 15 need the MMseqs2 binary (`MMSEQS`; the static `mmseqs-linux-avx2` release build was used,
kept under `$COVERAGE_WORK/tools/`). Needs network access, about 1 GB of disk for the downloads (the Swiss-Prot flat file is 0.7 GB) and the
extra requirements in `requirements.txt`, which the main environment does not carry. The SAXS table comes from
`bsc fetch` (`BSC_ROOT/data/PeptoneDB-SAXS/`) and the clean-entry flags from `results/resolvability.csv`.
Downloads are cached in `$COVERAGE_WORK/downloads` (`COVERAGE_DOWNLOADS` moves them); delete that directory
to refetch everything, which also starts a new `versions.json`.
