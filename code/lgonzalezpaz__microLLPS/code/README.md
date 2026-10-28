# microLLPS — neuronal microexons as phase-separation switches

Code and processed data accompanying the article
**"Neuronal microexons are tuned switches of protein phase separation across the human proteome"**
([AUTHORS], 2026; DOI: [ARTICLE DOI]). Archived at Zenodo: `10.5281/zenodo.[RECORD]`.

The pipeline reconstructs paired protein isoforms that differ only by inclusion (ME+) or skipping
(ME−) of each catalogued human microexon, locates the microexon in its structural context
(disorder, AlphaFold confidence) and quantifies how its inclusion changes the phase-separation
propensity of the host intrinsically disordered region (IDR), relative to a composition- and
position-matched null model. All steps run on a laptop CPU (≈5–6 h on 2 cores, 8 GB RAM, no GPU).

## Quick start

```bash
conda env create -f environment.yml && conda activate microllps
bash run_all.sh          # downloads public inputs, runs scripts 00–07, writes results/ and figures/
```

## Repository layout

| Path | Content |
|---|---|
| `scripts/download_data.sh` | Downloads all public inputs (Ensembl 115, microexon catalogue, PhaSePro, CD-CODE, SFARI Gene) and third-party code |
| `scripts/00_fetch_me_sequences.py` | Extracts microexon nucleotide sequences from GRCh38 |
| `scripts/01_build_isoforms.py` | Maps microexons to Ensembl CDS and builds ME+/ME− protein pairs |
| `scripts/02_llps_features.py` | Disorder, predicted pLDDT, IDR context, ΔG/c_sat (CALVADOS-trained predictor), FINCHES ε, null insertions |
| `scripts/03_alphafold_plddt.py` | AlphaFold DB v6 per-residue pLDDT at the microexon |
| `scripts/04_benchmark_predictors.py` | Experimental c_sat benchmark and PhaSePro-vs-proteome AUROC |
| `scripts/05_neural_psi.py` | Neural vs non-neural ΔPSI from 454 RNA-seq conditions |
| `scripts/06_analysis.py` | Statistics, sensitivity analyses, tables and figures |
| `scripts/08_orthogonal_llphyscore.py` | Orthogonal validation with LLPhyScore (locally executable ML predictor) |
| `scripts/09_switch_permutation.py` | 1,000-null permutation and scrambled-peptide tests for called switches |
| `results/` | `stats_summary.json`, Tables 1, S1–S3, benchmark outputs |
| `figures/` | Figures 1–7 (SVG vector, PDF, 600-dpi PNG) |
| `data/processed/isoform_pairs.tsv` | All reconstructed ME+/ME− isoforms |

## Inputs (all public)

* Microexon catalogue and neural PSI matrix: Parada *et al.*, *Nature* (2026), github.com/geparada/microexon-regulatory-code (MIT).
* Ensembl release 115 (GRCh38): GTF, CDS FASTA, UniProt cross-references, primary assembly.
* AlphaFold Protein Structure Database v6 (CC-BY 4.0), retrieved through its public API.
* PhaSePro (phasepro.elte.hu), CD-CODE (cd-code.org), SFARI Gene (gene.sfari.org).
* Phase-separation predictor: von Bülow *et al.*, *PNAS* 122, e2417920122 (2025) (GPL-3.0; downloaded, not redistributed).
* LLPhyScore standalone 2022-06-27 (github.com/julie-forman-kay-lab/LLPhyScore).
* FINCHES 1.0 (Ginell *et al.*, *Science* 2025) and metapredict 3.1 (Emenecker *et al.*, *Biophys. J.* 2021).

## Notes on reproducibility

* Random seeds are fixed (null insertions: per-event streams seeded with (2026, CRC32(event ID)), independent of sharding or resumption; benchmark: 7; figures: 42).
* `scikit-learn==1.6.1` is required to load the published ΔG/c_sat models; `afrc==0.3.7` is required by FINCHES 1.0.
* Orthogonal validation uses LLPhyScore (Cai, Vernon & Forman-Kay, *Biomolecules* 2022), run locally on every tested microexon.
* Ensembl REST was unavailable during development; sequences are extracted locally from the genome FASTA instead and agree 100% with the catalogue's 344 published microexon sequences.

## Licence

MIT (see `LICENSE`). Third-party data and code keep their own licences.

## Citation

See `CITATION.cff`.
