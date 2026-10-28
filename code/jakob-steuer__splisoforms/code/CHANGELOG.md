# Changelog

All notable changes to the SPLISOFORMS project will be documented in this file.

## [1.4.0] - 2026-10-27

### Methods-review rebuild

A full methods review of the pipeline against the preprint led to a database rebuild. The
public resource now serves the rebuilt data; the preprint is linked from the site and cited
in `CITATION.cff` (bioRxiv, doi:10.64898/2026.08.06.743218).

#### 🧬 Ingestion & canonical assignment
- One isoform per protein sequence: when identical proteins collide, the survivor is chosen
  by priority (MANE Select > Ensembl canonical > non-readthrough > named gene > has structure)
  instead of file order; every MANE protein is now flagged and is its gene's canonical.
- PacBio isoform IDs are resolved per sample; the same ID naming different proteins in
  different samples is split (`<id>@<sample>`), and per-sample expression is joined in the
  right sample file.
- Gene `canonical_sequence` is the designated canonical isoform's own sequence, the baseline
  for every canonical-relative layer.

#### 🔬 Analysis layers
- **Pfam**: 1-based coordinates straight from pyhmmer; same-clan overlapping hits collapsed
  to the best E-value (pfam_scan convention) before the canonical diff.
- **TED**: full human consensus table (42,599 domains); Chainsaw segments shifted to 1-based.
- **NMD**: the 50-nt rule on real CDS coordinates; GENCODE `cds_end_NF`/`mRNA_end_NF` →
  unknown, `NMD_exception` → escape; a CDS that runs past the transcript end is unknown.
- **Disorder**: documented as AIUPred + AIUPred-binding (the method actually used).
- **PTMs**: sites whose residue or position disagrees with the canonical sequence are dropped.
- **Structure**: PAE rigid-block scoring rewritten (summed-area table); parallel metrics.
- **Interactions**: homodimer layer repaired (12,006 genes, was 2,893).
- **Neoantigens**: every row carries MHCflurry, NetMHCpan 4.2 and BigMHC scores (the earlier
  release had NetMHCpan missing on ~40 % of rows); P1–P4 is a true three-tool vote; new
  `proteome_uniqueness` flag (unique / same-gene / other-gene) with atlas filters and chips;
  per-row tool availability.
- **Search**: exact gene symbols and transcript names rank first; protein names are searched.

#### 🛠 Pipeline & operations
- `backfill_neoantigen_uniqueness.py` (new), `hpc/build_isoform_id_map.py` (new), resumable
  allele-parallel NetMHCpan runner, sparse-SLURM-array guard, cross-shard dedup in the loader.
- Backfill scripts exit non-zero when they load nothing; retired scripts removed.
- Atlas partial indexes declared on the model (migration `r1b6e3d8c5a2`) with a covering
  stats index and a short-lived stats cache.
- Serving image pins `SQLAlchemy<2.1`.

## [1.3.0] - 2026-07-15

### Major Release: Drug Layer, Neoantigen Rebuild & Data Integrity

#### 💊 Drug Interactions (new)
- **Co-folded drug pockets**: 20,091 AlphaFold 3 drug–protein complexes across 3,725 compounds and 2,889 target genes. Each drug's 4.5 Å contact pocket is carried onto every splice isoform and scored as preserved / disrupted / lost, with a pLDDT≥70 no-call gate.
- **Drug pages & lens**: `/explore/drugs` ranked lens, per-gene drug-pocket pages with paired structure viewers, and a drug summary panel (2D skeleton, bound conformer, PubChem chemistry). Drug search resolves target–drug pairs.

#### 🧬 Data Integrity
- **Gene re-attribution**: 20,420 novel long-read isoforms re-homed to their true gene using each dataset's own ground truth (GTF / SQANTI) instead of a UniProt best-hit fallback. Canonical is now MANE-first, one per gene.
- **Neoantigen ingest**: 5.17M predictions loaded and corrected (removed peptides that occur in the current canonical; Gold Medal tier applied). Discovery index rebuilt.
- **Readable transcript names**: every isoform gets a display name (e.g. `MET-201`, `MET-206-NOV-a6d00a`) alongside its stable source ID.

#### 🔗 Protein Interactions & Docs
- Homodimer **and** heterodimer interface-disruption analysis surfaced and documented; STRING seven-channel evidence.
- Landing page reworked into a 3×3 feature tour; documentation completed (Drug Interactions written, tooltips clarified).

## [1.1.0] - 2026-04-15

### Major Release: Structural & Genomic Analytics

#### 🧬 Genomic Visualization
- **Multi-Exon Architecture**: Implemented an inheritance-based genomic coordinate system. Novel transcripts now inherit detailed exon boundaries from the human transcriptome cache, enabling high-fidelity visualization of FSM/ISM variants.
- **Transcriptome Cache**: Ingested the entire Ensembl human transcriptome (~200k records) to ensure 100% genomic coverage for novel clinical isoforms.

#### 💎 Structural Analytics
- **AlphaFold 3 (AF3) Support**: Enabled custom AF3 structural predictions for novel clinical isoforms. The pipeline now supports sub-angstrom accuracy for non-canonical sequences.
- **SASA-based PTM Accessibility**: Refactored the PTM analysis engine to use **Gemmi** and **FreeSASA** for precise solvent accessibility calculations, allowing for nuanced before/after comparison between canonical and novel structures.

#### 📊 Functional Impact Dashboard
- **Granular Pfam Monitoring**: Activated a live functional assessment panel that tracks Domain Loss, Gain, and Truncation relative to reference proteins.
- **NMD & Neoantigen Fidelity**: Enhanced prediction reliability for Nonsense-Mediated Decay (NMD) and peptide binding categories.

#### 🎨 User Interface & Experience
- **Nomenclature Standardization**: Updated scientific terminology across the app (e.g., "Loop" replaced "Coil" in structural legends).
- **Branding & Data Provenance**: Redesigned the footer and integration panels for better integration with the Human Protein Atlas (HPA) and AlphaFoldDB.

---

## [1.0.2] - 2026-04-10
### Added
- Initial support for Mol* interactive structure viewing.
- Basic Pfam domain scanning via pyHMMER.

## [1.0.0] - 2026-04-01
### Added
- Initial release of the SPLISOFORMS platform with canonical AlphaFoldDB coverage.
