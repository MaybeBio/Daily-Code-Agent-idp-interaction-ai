# files/

Shared intermediate tables used by the notebooks in `notebooks/`. They are pickled pandas DataFrames, read with `pd.read_pickle`. The notebooks reference them as `../files/<name>`, so run each notebook with `notebooks/` as the working directory.

| File | Contents | Created by | Read by notebooks |
|---|---|---|---|
| `dSAM_MPs_v_paper_v2_new_locs.pkl` | The mdSAM table loaded by the analysis notebooks; the motif sequence is in the `AAseq` column. | The mdSAM identification pipeline (STAR Methods: *Self-assembly prediction*, *Intrinsically disordered region prediction*) | 05, 06, 07, 08, 09, 15, 18, 20, 21, 23 |
| `hTMs_uniprot_v1.pkl` | Table of human transmembrane proteins with UniProt topology annotation. | `03_mdsam_sequence_logos.ipynb` (`TMs.to_pickle`) | 03, 05, 06, 07, 08, 09, 15, 18, 20, 21, 23 |
| `AF2_with_ENSG.pkl` | Table of AlphaFold-derived per-protein annotations with Ensembl gene identifiers (ENSG) added. | `02_alphafold_residue_table_and_early_plots.ipynb` (`AF2.to_pickle`) | 02, 03, 05, 06, 07, 08, 11, 15, 18, 19, 20, 21, 23 |

Notebook numbers refer to the file-name prefixes in `notebooks/` (for example `05` is `05_uniprot_keyword_enrichment.ipynb`).

Pickle files can execute code when loaded and depend on the pandas version used to write them. Load them only from this repository, and use a pandas version compatible with the one used to create them.
