# mdSAM: Discovery of self-assembly motifs regulating membrane-protein organization and function

Analysis notebooks and figure panels for the manuscript *Discovery of self-assembly motifs regulating membrane-protein organization and function* (Elazar, Trivedi, et al., Babu lab, St. Jude Children's Research Hospital).

## Abstract

While cell-surface clustering of membrane proteins is central to numerous processes, its molecular determinants remain elusive. We postulated that peptide motifs within disordered regions of transmembrane proteins may drive self-association and clustering and developed an approach to identify them. Analysis of the human proteome revealed that one-fifth of membrane proteins (945 proteins; >1,400 motifs), including receptors and transporters, contain motifs that we term mdSAMs (membrane-disordered Self-Assembly Motifs). Biophysical, biochemical, structural, and cell-based analyses of mdSAMs in model and native membrane proteins demonstrate that mdSAMs can drive self-assembly and clustering. In the serotonin receptor HTR7, a single mdSAM-disrupting substitution alters nanoscale organization, attenuates basal and agonist-induced signaling, and shifts the transcriptional state of the cell. Notably, mdSAM sequences are unique, which might enable discrimination of self from non‑self during self‑assembly, yet conserved among orthologs, underscoring their relevance. These findings reveal a widespread mechanism linking membrane-protein self-assembly to cell-surface organization and function.

## What is an mdSAM?

An **mdSAM** (membrane-disordered Self-Assembly Motif) is a short peptide motif, usually 7-9 residues (median 8), located in an intrinsically disordered, solvent-exposed region of a transmembrane protein, that is predicted to self-assemble through beta-strand-mediated, amyloid-like interactions. The working hypothesis is that mdSAMs provide weak, sequence-specific self-association that can cluster membrane proteins at the cell surface.

**How mdSAMs are identified**

1. A random-forest classifier is trained on WALTZ-DB 2.0 (1,402 hexapeptides: 510 amyloid-forming, 892 non-forming) using the one-hot-encoded hexapeptide sequence plus the summed beta-strand free energy (F1 = 0.77, stratified 5-fold cross-validation).
2. Every hexapeptide window in human membrane proteins is scored; positive windows are Self-Assembly Motifs (SAMs).
3. A SAM is called an mdSAM if it lies in a disordered region: AlphaFold pLDDT < 70, non-helical and solvent exposed, excluding transmembrane segments and signal peptides. A ROC-derived threshold on consecutive positive windows sets the high-confidence set (Fig. S1I).

**What was found.** 1,425 high-confidence mdSAMs in 945 membrane proteins (about 20% of the 4,886 human membrane proteins analysed), 91% of them in N- or C-terminal tails, and 72% of mdSAM-containing proteins carry a single motif. Of 44 synthesised mdSAM peptides, 36 (81%) gave Thioflavin-T signal above buffer; in the eight peptides tested, proline substitutions reduced the signal, abolishing it in seven; cryo-EM showed ordered fibrils. In cells, mdSAMs increase surface puncta density of a synthetic scaffold and of native full-length proteins. In the serotonin receptor HTR7, the V451P substitution in the C-terminal mdSAM (VSVWVLQS) changes nanoscale cluster organisation, lowers basal and agonist-induced cAMP signalling (Emax reduced by 25%) and shifts the transcriptional state of the cell. mdSAM sequences are diverse and unique, and conserved among orthologs.

> Naming: older notebooks call these motifs **dSAM** (disordered SAM) and refer to the structured-region counterpart as **mSAM/sSAM/msSAM**. dSAM in a file name or variable means mdSAM.

## Repository contents

```
notebooks/   the 23 analysis notebooks listed below, numbered in pipeline order (outputs cleared)
files/       shared input tables read by the notebooks (see files/README.md)
FIGURE_MAP.csv   machine-readable map: paper panel -> notebook -> cell
```

## What each notebook does

Cell numbers are 0-based indices into the notebook (as in the `cells` list of the `.ipynb`).

### Data preparation and classifier

| Notebook | Task | Paper panels |
|---|---|---|
| `01_train_waltzdb_classifier.ipynb` | Trains and evaluates the random-forest classifier on WALTZ-DB 2.0 hexapeptides (one-hot sequence + summed beta-strand energy); reports test-set counts, stratified K-fold ROC and the beta-strand energy distributions. Saves the trained model. | Fig 1D; S1E-G |
| `02_alphafold_residue_table_and_early_plots.ipynb` | Builds the shared AlphaFold-derived residue table (`AF2_with_ENSG.pkl`, cell 205) used by most other notebooks; remaining cells hold earlier versions of exploratory plots. | (input table) |
| `03_mdsam_sequence_logos.ipynb` | Sequence logos and composition of mdSAMs and SAMs in membrane proteins (logomaker); also builds `hTMs_uniprot_v1.pkl` and `AF_IDR_v3.pkl`. Contains additional gnomAD analyses that are not part of the paper figures. | S1H |

### Discovery and diversity (Figs 1-2, S1-S2)

| Notebook | Task | Paper panels |
|---|---|---|
| `04_pasta_aggregation_comparison.ipynb` | Compares mdSAM predictions with PASTA energies and an aggregation threshold on windows around WALTZ-DB hexapeptides; ROC for the number of consecutive positive windows. | S1I-J |
| `05_uniprot_keyword_enrichment.ipynb` | UniProt keyword enrichment (observed minus expected, FDR) for mdSAM-containing membrane proteins: biological process, molecular function, disease. | Fig 1G; S1K |
| `06_mdsam_sequence_space_embedding.ipynb` | Embeds the mdSAM sequence space from physicochemical features, colours it by composition, clusters it (GMM) with keyword enrichment, draws the charge-archetype ternary plot, the dot matrix of mdSAM-containing proteins and the pairwise Levenshtein-distance analysis. | Fig 1F, 2F-I, 3A; S2E-G |
| `07_hyperlopit_organelle_enrichment.ipynb` | Re-analyses the hyperLOPIT subcellular fractionation (U2OS): fraction profiles, t-SNE and Gaussian-mixture organelle assignment, and mdSAM enrichment per organelle. | Fig 2A; S2A-C |
| `08_mdsam_topology_and_membrane_distance.ipynb` | Topology of mdSAM-containing proteins: number of TM segments (bitopic vs polytopic), N-/C-terminus vs loop, distance to the nearest membrane leaflet, number and length of mdSAMs. | Fig 2B-F |
| `09_mdsams_per_protein.ipynb` | Histogram of mdSAMs per protein (1, 2, 3 or more). | Fig 2E |

### In vitro self-assembly (Fig 3, S3)

| Notebook | Task | Paper panels |
|---|---|---|
| `10_tht_peptide_screen_96well.ipynb` | Thioflavin-T screen of mdSAM peptides in 96-well plates: per-peptide intensities, spectra, wild-type vs proline-mutant plates, overlay on the charge-archetype ternary plot. | Fig 3A-B; S3A, S3C |
| `11_tht_gpcr_peptides_proline_mutants.ipynb` | ThT of GPCR mdSAM peptides (HTR7, CCR4, DRD1) and their proline mutants, including HTR7 V451P. | Fig 4C, 5B; S3C |
| `12_amytracker_peptides.ipynb` | Amytracker 520 fluorescence (100 uM) of wild-type and proline-mutant peptides. | Fig 4C; S3D |

### Cell-based clustering (Figs 4-5, S5-S6)

| Notebook | Task | Paper panels |
|---|---|---|
| `13_tirf_puncta_poisson_glm.ipynb` | Poisson generalised linear model on CellProfiler cell and puncta tables: puncta per cell area and rate ratio (WT/mutant) per protein and cell line (HEK, U2OS), SNAPtm motif swaps, full-length APP and HTR7. | Fig 4C-D, 5E; S5C |
| `14_flow_cytometry_total_vs_surface.ipynb` | Flow cytometry (FCS files): total and surface SNAP signal for wild-type vs mutant constructs by DNA amount, with percentile gating against mock cells. | S6A |
| `15_biogrid_interactions_by_mdsam.ipynb` | BioGRID analysis of membrane-protein interactions by mdSAM status: self-interactions and plasma-membrane interaction networks with hypergeometric and binomial tests. Contains further interaction analyses that are not part of the paper figures. | S6B-C |
| `16_dstorm_htr7_cluster_analysis.ipynb` | dSTORM cluster analysis: ThunderSTORM localisation tables, Voronoi thresholding against randomised points, HDBSCAN clustering, cluster size and density stratified by local density, bootstrap of WT-vs-V451P skewness (tail-heaviness) for HTR7. Also holds earlier SNAPtm/APP analyses that are not part of the paper figures. | Fig 5F-G |

### HTR7 signalling and transcription (Fig 6)

| Notebook | Task | Paper panels |
|---|---|---|
| `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` | HTR7 cAMP Emax and EC50 comparison between WT and V451P; per-protein total- and surface-expression plots with a mixed model adjusted for total expression. | Fig 6D; S6A |
| `18_htr7_rnaseq_pca_volcano.ipynb` | RNA-seq PCA and volcano plots from DESeq2 tables (HTR7 WT, V451P, empty vector). | Fig 6E |
| `19_htr7_rnaseq_gsea.ipynb` | Re-runs differential expression (pyDESeq2) and gene-set enrichment (gseapy) for the HTR7 RNA-seq data; heatmaps and volcano plots. | supports Fig 6E |

### Evolution (Fig 7, S7)

| Notebook | Task | Paper panels |
|---|---|---|
| `20_conservation_alignment_based.ipynb` | Alignment-based conservation (Kullback-Leibler, Shannon, Lockless; pycanal) of mdSAMs vs disordered and structured residues; retention of mdSAMs in one-to-one orthologs vs alignment depth. | Fig 7B; S7A-B |
| `21_conservation_esm_embedding.ipynb` | Alignment-free conservation from protein language-model embeddings for mdSAMs, random IDR segments, ELM motifs and structured residues. | Fig 7B |
| `22_ortholog_peptide_tht.ipynb` | ThT of NGFR and TNFRSF13C ortholog mdSAM peptides (wild type and proline mutant) across species. | Fig 7D; S7C |
| `23_overlap_with_larks_wirs_finches.ipynb` | Overlap of mdSAM-containing membrane proteins with LARKS, WIRS and FINCHES-predicted fuzzy interactions (UpSet plot). | S7D |

## Where the paper figures come from

Each row points to the notebook and cell that draws the panel. Cell numbers are 0-based indices into the notebook. Notebooks are distributed with outputs cleared; run a notebook to regenerate the plot. Final figures are in the paper.

### Figure 1: Discovery of mdSAMs

| Panel | Content | Notebook (cell) |
|---|---|---|
| 1D | RF classifier test-set performance (amyloid vs non-amyloid) | `01_train_waltzdb_classifier.ipynb` (cell 26) |
| 1F | Dot matrix: membrane proteins with predicted mdSAMs | `06_mdsam_sequence_space_embedding.ipynb` (cell 14) |
| 1G | UniProt keyword categories (donut) | `05_uniprot_keyword_enrichment.ipynb` (cell 32) |

### Figure 2: Diversity of mdSAMs

| Panel | Content | Notebook (cell) |
|---|---|---|
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 79) |
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 87) |
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 89) |
| 2B | Bitopic vs polytopic enrichment | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 70) |
| 2B | Bitopic vs polytopic enrichment | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 73) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 146) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 148) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 168) |
| 2D | Distance to nearest membrane leaflet | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 45) |
| 2D | Distance to nearest membrane leaflet | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 55) |
| 2E | mdSAMs per protein | `09_mdsams_per_protein.ipynb` (cell 2) |
| 2F | mdSAM length distribution | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 76) |
| 2G | Pairwise Levenshtein distance distribution / matrix | `06_mdsam_sequence_space_embedding.ipynb` (cell 207) |
| 2G | Pairwise Levenshtein distance distribution / matrix | `06_mdsam_sequence_space_embedding.ipynb` (cell 211) |
| 2H-I | Physicochemical embedding coloured by composition | `06_mdsam_sequence_space_embedding.ipynb` (cell 28) |
| 2H-I | Physicochemical embedding coloured by composition | `06_mdsam_sequence_space_embedding.ipynb` (cell 37) |

### Figure 3: mdSAM peptides self-assemble

| Panel | Content | Notebook (cell) |
|---|---|---|
| 3A | Charge archetypes (ternary plot) | `10_tht_peptide_screen_96well.ipynb` (cell 25) |
| 3A | Charge archetypes (ternary plot) | `10_tht_peptide_screen_96well.ipynb` (cell 26) |
| 3B | ThT fluorescence per peptide with buffer/MAPT lines | `10_tht_peptide_screen_96well.ipynb` (cell 54) |
| 3B | ThT fluorescence per peptide with buffer/MAPT lines | `10_tht_peptide_screen_96well.ipynb` (cell 73) |

### Figure 4: mdSAMs promote clustering in cells

| Panel | Content | Notebook (cell) |
|---|---|---|
| 4C | Amytracker WT vs proline peptides | `12_amytracker_peptides.ipynb` (cell 3) |
| 4C | SNAPtm motif swaps: puncta density, HTR7 and KCNE4 mdSAM | `13_tirf_puncta_poisson_glm.ipynb` (cell 42) |
| 4C | SNAPtm motif swaps: puncta density, HTR7 and KCNE4 mdSAM | `13_tirf_puncta_poisson_glm.ipynb` (cell 43) |
| 4C / 5B | ThT of HTR7 WT vs V451P peptide | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 16) |
| 4C / 5B | ThT of HTR7 WT vs V451P peptide | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 17) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 12) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 15) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 19) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 22) |

### Figure 5: HTR7 mdSAM and nanoscale organisation

| Panel | Content | Notebook (cell) |
|---|---|---|
| 5B | HTR7 and other GPCR peptide ThT spectra | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 32) |
| 5E | HTR7 WT vs V451P TIRF puncta | `13_tirf_puncta_poisson_glm.ipynb` (cell 25) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 409) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 410) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 440) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 451) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 452) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 453) |

### Figure 6: HTR7 mdSAM, signalling and transcription

| Panel | Content | Notebook (cell) |
|---|---|---|
| 6D | HTR7 cAMP Emax and EC50 summary statistics | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 18) |
| 6E | RNA-seq PCA and V451P vs WT volcano | `18_htr7_rnaseq_pca_volcano.ipynb` (cell 3) |
| 6E | RNA-seq PCA and V451P vs WT volcano | `18_htr7_rnaseq_pca_volcano.ipynb` (cell 8) |

### Figure 7: Evolutionary conservation

| Panel | Content | Notebook (cell) |
|---|---|---|
| 7B | ESM2 embedding conservation by segment class | `21_conservation_esm_embedding.ipynb` (cell 15) |
| 7B | Conservation by segment class (related analyses) | `20_conservation_alignment_based.ipynb` (cell 22) |
| 7B | Conservation by segment class (related analyses) | `20_conservation_alignment_based.ipynb` (cell 28) |
| 7D | Ortholog ThT: WT vs proline mutant per species | `22_ortholog_peptide_tht.ipynb` (cell 6) |
| 7D | Ortholog ThT: WT vs proline mutant per species | `22_ortholog_peptide_tht.ipynb` (cell 16) |

### Supplementary Figure S1: Validation of mdSAM prediction

| Panel | Content | Notebook (cell) |
|---|---|---|
| S1F | ROC, stratified 5-fold CV | `01_train_waltzdb_classifier.ipynb` (cell 15) |
| S1G | Beta-strand energetics | `01_train_waltzdb_classifier.ipynb` (cell 12) |
| S1H | Sequence logos, mdSAMs vs SAMs | `03_mdsam_sequence_logos.ipynb` (cell 35) |
| S1H | Sequence logos, mdSAMs vs SAMs | `03_mdsam_sequence_logos.ipynb` (cell 37) |
| S1I | ROC for number of consecutive windows classified as amyloid | `04_pasta_aggregation_comparison.ipynb` (cell 17) |
| S1I | ROC for number of consecutive windows classified as amyloid | `04_pasta_aggregation_comparison.ipynb` (cell 21) |
| S1J | PASTA energy of mdSAMs | `04_pasta_aggregation_comparison.ipynb` (cell 14) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 47) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 73) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 74) |

### Supplementary Figure S2: mdSAM diversity and localisation

| Panel | Content | Notebook (cell) |
|---|---|---|
| S2A | hyperLOPIT profiles per compartment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 26) |
| S2A | hyperLOPIT profiles per compartment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 29) |
| S2B-C | hyperLOPIT t-SNE and GMM organelle assignment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 48) |
| S2B-C | hyperLOPIT t-SNE and GMM organelle assignment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 68) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 101) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 122) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 124) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 49) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 55) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 65) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 67) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 69) |

### Supplementary Figure S3: Peptide assembly

| Panel | Content | Notebook (cell) |
|---|---|---|
| S3A | ThT spectra, all peptides | `10_tht_peptide_screen_96well.ipynb` (cell 86) |
| S3C | ThT spectra, WT vs proline mutants | `10_tht_peptide_screen_96well.ipynb` (cell 80) |
| S3D | Amytracker WT vs proline | `12_amytracker_peptides.ipynb` (cell 3) |

### Supplementary Figure S6: Expression controls and BioGRID analysis

| Panel | Content | Notebook (cell) |
|---|---|---|
| S6A | Total expression bar plots, per protein | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 13) |
| S6A | Total expression bar plots, per protein | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 17) |
| S6A | Total expression, WT vs proline, by DNA amount | `14_flow_cytometry_total_vs_surface.ipynb` (cell 45) |
| S6A | Total expression, WT vs proline, by DNA amount | `14_flow_cytometry_total_vs_surface.ipynb` (cell 98) |
| S6B | BioGRID self-interactions, mdSAM+ (n=658) and mdSAM- (n=2,642) | `15_biogrid_interactions_by_mdsam.ipynb` (cell 171) |
| S6B | BioGRID self-interactions, mdSAM+ (n=658) and mdSAM- (n=2,642) | `15_biogrid_interactions_by_mdsam.ipynb` (cell 172) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 78) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 81) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 84) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 86) |

### Supplementary Figure S7: Conservation

| Panel | Content | Notebook (cell) |
|---|---|---|
| S7A | Alignment-based conservation metrics by residue class | `20_conservation_alignment_based.ipynb` (cell 24) |
| S7A | Alignment-based conservation metrics by residue class | `20_conservation_alignment_based.ipynb` (cell 25) |
| S7B | mdSAM retention in orthologs vs alignment depth | `20_conservation_alignment_based.ipynb` (cell 53) |
| S7B | mdSAM retention in orthologs vs alignment depth | `20_conservation_alignment_based.ipynb` (cell 71) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 20) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 21) |
| S7D | UpSet: overlap with LARKs, WIRS, FINCHES | `23_overlap_with_larks_wirs_finches.ipynb` (cell 16) |

Panels that are schematics, microscopy or EM images, NMR or mass-photometry data, or were produced outside these notebooks are not listed.

## Running the notebooks

The notebooks are distributed with their outputs cleared.

The notebooks were written for an interactive session and read three shared tables from `files/` (`dSAM_MPs_v_paper_v2_new_locs.pkl`, `hTMs_uniprot_v1.pkl`, `AF2_with_ENSG.pkl`); see `files/README.md`. Run each notebook with `notebooks/` as the working directory. Some notebooks still read other intermediate files that are not in this repository. Public inputs are WALTZ-DB 2.0, the AlphaFold human-proteome table of N. Davey (github.com/normandavey/ProcessedAlphafold), UniProt, the hyperLOPIT data (Thul et al., 2017), BioGRID 4.4.x, PANTHER 17 and GTEx v8. Supplementary data files accompany the paper.

Python packages imported by these notebooks (pip names can differ from import names): `Bio`, `Levenshtein`, `PIL`, `PyPDF2`, `adjustText`, `bokeh`, `cmasher`, `fcsparser`, `geopandas`, `gseapy`, `holoviews`, `intersect`, `iupred2a_lib`, `joblib`, `localcider`, `logomaker`, `matplotlib`, `matplotlib_venn`, `mycolorpy`, `networkx`, `numpy`, `pandas`, `plotly`, `pyaaisc`, `pycanal`, `pydeseq2`, `radialtree`, `rdfpy`, `ripleyk`, `sanbomics`, `scanpy`, `scipy`, `seaborn`, `shapely`, `sklearn`, `statsmodels`, `ternary`, `tqdm`, `umap`, `uncertainties`, `upsetplot`. `iupred2a_lib`, `pyaaisc` and `intersect` are local modules rather than PyPI packages.

## Citation and contact

Citation: to be added after publication. Correspondence: M. Madan Babu (madan.babu@stjude.org) and Assaf Elazar (assaf.elazar@stjude.org).
