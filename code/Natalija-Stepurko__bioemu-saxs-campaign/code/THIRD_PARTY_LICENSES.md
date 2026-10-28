# Third-party licences

The software, models, data and services this project uses, with their licences or terms, how each is
used, and what of each this repository holds. Each curated row names the page or file where I read the
licence and the date I read it; the Python tables are generated from the installed package metadata.

## Summary

**This repository.** The code, documents and derived results written for this project are under the MIT
licence in [`LICENSE`](https://github.com/Natalija-Stepurko/bioemu-saxs-campaign/blob/HEAD/LICENSE). It does not cover the third-party tools, models and data below, nor the
third-party material copied into some tables (listed under *Copied into tracked files*); those keep their own
licences.

**Not redistributed.** The repository holds no third-party program, binary, package, model weight or raw data
file: no SAXS profiles, computed curves, ensembles, trajectories, structures, database dumps or flat files.
This was checked against the tracked files and against every file ever added in the history.

**Downloaded at run time** and kept in git-ignored folders: the PeptoneBench archive (`bsc fetch`, `data/`);
the SASBDB entry summaries; CATH, UniProt, RCSB, InterPro, SIFTS and TED records and BioEmu's training-data
pieces (coverage analysis, `COVERAGE_WORK`); Pepsi-SAXS, PULCHRA, the CALVADOS 3 parameters and the CALVADOS
training lists (simulations, `SIM_WORK`); the BioEmu checkpoints and the AlphaFold2 parameters, and the MSAs from
the ColabFold server (Colab notebook). MMseqs2 is downloaded by hand.

**Copied into tracked files**, as table columns, with attribution given in the data table below:

- construct sequences and SASBDB entry titles from PeptoneDB-SAXS (CC-BY-4.0) in `results/entries.csv`,
  `results/simulations/inputs.csv` and `notebooks/resample_list.csv`;
- SASBDB entry metadata, including project titles and contributor names, in `results/coverage/sasbdb_metadata.csv`
  (SASBDB data carry no copyright restrictions; attribution to the depositors is asked);
- CATH superfamily names and S35 counts (CC-BY-4.0), UniProt accessions and lineages (CC-BY-4.0), InterPro and Pfam
  identifiers (CC0-1.0) and TED labels (CC-BY-4.0) in `results/coverage/`;
- names and sequences of the CALVADOS training and validation proteins in `results/simulations/calvados_sets.csv`
  (CALVADOS 1 lists CC-BY-4.0; CALVADOS 3 lists MIT, Copyright (c) 2024 Linderstrøm-Lang Centre for Protein
  Science, University of Copenhagen; CALVADOS 2 lists: the GitHub copy carries no licence file, the Zenodo
  deposit of the same folder is CC-BY-4.0).

**Restrictions that matter in practice.**

- **Pepsi-SAXS is academic use only** ("All rights reserved. The academic version is free."; the program prints
  "For academic use only!"). It computes the curves of the simulations (`analysis/simulations/06_saxs.py`), and
  the PeptoneBench authors used it for the archive's curves. Rerunning that step requires being an academic user;
  the terms say nothing about the computed curves.
- **Copyleft in the toolchain, none of it shipped.** CALVADOS (GPL-3.0) is imported by
  `analysis/simulations/05_run.py`; DENSS (GPL-3.0) is run as a separate program; mdtraj (LGPL-2.1-or-later) is
  imported by the main package; MDAnalysis and GridDataFormats (LGPL-3.0-or-later), OpenMM's CUDA and OpenCL
  platforms (LGPL) and certifi and tqdm (MPL-2.0) are installed dependencies. These terms bind whoever
  distributes those packages; the repository distributes none of them.
- **No non-commercial model terms were found.** The BioEmu code and all three checkpoints are MIT; the AlphaFold2
  parameters BioEmu downloads are CC-BY-4.0.
- **The ColabFold MSA server has no published terms**; its README asks for serial queries from one IP. Sampling
  with BioEmu sends each sequence to that server.
- **Attribution** is required for the CC-BY-4.0 sources (PeptoneBench, MEGAscale, CATH, UniProt, AlphaFold DB,
  TED, the AlphaFold2 parameters, the CALVADOS 1 and 2 lists) and asked for by SASBDB, EMBL-EBI and RCSB; the
  attribution column of each table states what each asks.

**Not confirmed:** the terms of the ColabFold MSA server, of Kaggle, of the SASBDB and TED APIs as services, a
licence for SIFTS data and for CATH-Gene3D as served by InterPro, and that the CALVADOS 2 Zenodo deposit holds
the two list files copied from GitHub. Each row says what was found.

**Regenerate** after changing an environment or `docs/licenses/curated.json`: `python3 scripts/licenses.py`
(`make licenses`). It rewrites the generated blocks below, `docs/licenses/python-packages.md` (every installed
package of every environment) and `docs/licenses/python-packages.csv`; `--venv ID=PATH` points at an environment
installed elsewhere.

## Tools, models and services

<!-- generated:tools:start -->
| Tool | How the project uses it | Obtained | In this repo | Licence or terms | Restrictions in practice | Attribution | Status | Checked at |
|---|---|---|---|---|---|---|---|---|
| **Pepsi-SAXS** 3.0 (May 2020), Linux binary | Computes the SAXS curve of each simulated conformer (`analysis/simulations/06_saxs.py`). The PeptoneBench authors used it to compute the archive's curves. | Downloaded at run time by `analysis/simulations/02_fetch_sources.py` from the authors' server into `SIM_WORK` (git-ignored); sha256 in `results/simulations/sources.json` | No | Proprietary. Web page: "All rights reserved. The academic version is free." Program banner: "For academic use only!" | Academic use only. The page publishes no terms for other use. | Cite Grudinin, Garkavenko & Kazennov, Acta Cryst. D73 (2017) 449-464 | confirmed | [Pepsi-SAXS page, section License](https://team.inria.fr/nano-d/software/pepsi-saxs/) (2026-10-16); [banner printed by the downloaded binary (--help)](https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip) (2026-10-16) |
| **PULCHRA** 3.06 | Rebuilds all-atom models from the CALVADOS Cα frames before Pepsi-SAXS (`06_saxs.py`) | C source downloaded at run time from KULL-Centre/papers at a pinned commit and compiled with gcc into `SIM_WORK` (git-ignored) | No | MIT (LICENSE file in pulchra_306.tgz, Copyright (c) 2000-2009 Piotr Rotkiewicz) | None beyond keeping the notice in copies | Keep the MIT notice in copies; cite Rotkiewicz & Skolnick, J Comput Chem 2008 | confirmed | [pulchra_306.tgz, file pulchra_306/LICENSE](https://raw.githubusercontent.com/KULL-Centre/papers/85979b1b4123b6b5391b617d16551969eda9f56e/2021/CG-IDPs-Tesei-et-al/pulchra_306.tgz) (2026-10-16) |
| **MMseqs2** release 18-8cc5c, static linux-avx2 binary | Sequence searches of the constructs against the training sets and CATH S100 (`analysis/coverage/14_annotate_training.py`, `15_nearest_training.py`) | Downloaded by hand from the GitHub release into `COVERAGE_WORK/tools` (command in `analysis/coverage/README.md`) | No | MIT (LICENSE.md at tag 18-8cc5c) | None | Cite Steinegger & Söding, Nat Biotechnol 2017 | confirmed | [MMseqs2 LICENSE.md at 18-8cc5c](https://github.com/soedinglab/MMseqs2/blob/18-8cc5c/LICENSE.md) (2026-10-16) |
| **DENSS** 1.8.7 | Ab initio envelopes of the three worked examples; `bsc envelopes` runs the `denss-all` program as a separate process | Python package in its own environment (`/scratch/.venv-denss`); maps written to `results/denss/` (git-ignored) | No; only fit statistics (`results/denss_stats.json`) and the figure drawn from the maps are tracked | GPL-3.0 (package metadata 'GPLv3'; repository LICENSE is GPL version 3) | Copyleft applies to anyone who distributes DENSS or a work based on it; the repo calls it as a program and ships none of it | Cite Grant, Nat Methods 2018 | confirmed | [tdgrant1/denss LICENSE](https://github.com/tdgrant1/denss/blob/master/LICENSE) (2026-10-16) |
| **CALVADOS (code)** commit 6981542 (package 0.9.0) | Coarse-grained simulations; imported by `analysis/simulations/05_run.py` | Installed into the simulation environment from GitHub at the pinned commit (`analysis/simulations/requirements.txt`) | No | GPL-3.0 (LICENSE at the pinned commit; package metadata GPL-3.0-only) | Copyleft applies to anyone who distributes CALVADOS or a work based on it. The repo's scripts import it but ship none of its code | Cite Tesei et al. PNAS 2021; Tesei & Lindorff-Larsen, Open Res Europe 2022; Cao et al. Protein Sci 2024 | confirmed | [KULL-Centre/CALVADOS LICENSE at 6981542](https://github.com/KULL-Centre/CALVADOS/blob/69815423632cb9b9b08f4b0aea34d8ba43a057c4/LICENSE) (2026-10-16) |
| **CALVADOS 3 residue parameters** residues_CALVADOS3.csv at commit 6981542 | Force-field parameters of the simulations; `07_analyse.py` reads the residue masses from it | Downloaded at run time by `02_fetch_sources.py` into `SIM_WORK` (git-ignored); sha256 in `results/simulations/sources.json` | No | GPL-3.0 (the file is part of the CALVADOS repository, whose LICENSE covers it; no separate licence for the parameters was found) | As for the CALVADOS code | Cite Cao et al. Protein Sci 2024 | confirmed | [KULL-Centre/CALVADOS repository at 6981542 (LICENSE at the root, no other licence file)](https://github.com/KULL-Centre/CALVADOS/tree/69815423632cb9b9b08f4b0aea34d8ba43a057c4) (2026-10-16) |
| **OpenMM** 8.4.0.post2, CPU platform | Simulation engine under CALVADOS | Python package in the simulation environment | No | MIT for the API, Reference and CPU platforms; LGPL for the CUDA and OpenCL platforms (the wheel ships the OpenCL plugins) | None for use; LGPL terms bind redistribution of the CUDA and OpenCL platforms | Cite Eastman et al. J Phys Chem B 2024 (OpenMM 8) | confirmed | [OpenMM Licenses.txt](https://github.com/openmm/openmm/blob/master/docs-source/licenses/Licenses.txt) (2026-10-16) |
| **metapredict** 3.1.1 | Predicted disorder of the constructs, training sequences and CATH S35 sequences (`14_annotate_training.py`) | Python package in the coverage environment; its network weights ship inside the package | No | MIT (package metadata and repository LICENSE) | None | Cite Emenecker et al. Biophys J 2021 | confirmed | [idptools/metapredict LICENSE](https://github.com/idptools/metapredict/blob/master/LICENSE) (2026-10-16) |
| **BioEmu (code)** archive ensembles: commit ac7455d (by the PeptoneBench authors); Colab re-sample: unpinned `bioemu[cuda]`; local: 1.4.1 | The model under study. The archive's ensembles were generated by the PeptoneBench authors; the notebook re-samples with bioemu-v1.2 on Colab; the local environment converts the sampled batches | Python package (Colab and `/scratch/.venv-bioemu`) | No | MIT. The package bundles code derived from ColabFold (MIT), a patched subset of AlphaFold2 (Apache-2.0) and OpenFold code (Apache-2.0), as listed in its README | None | Cite Lewis et al. Science 2025 | confirmed | [microsoft/bioemu LICENSE](https://github.com/microsoft/bioemu/blob/main/LICENSE) (2026-10-16); [microsoft/bioemu README, third-party code](https://github.com/microsoft/bioemu#readme) (2026-10-16) |
| **BioEmu checkpoints** bioemu-v1.0, v1.1, v1.2 | Model weights; the notebook samples with bioemu-v1.2. The archive does not record which checkpoint it used | Downloaded at run time from Hugging Face by `bioemu.sample` | No | MIT. The Hugging Face model repository, which holds all three checkpoints, declares `license: mit` and its card states 'License: MIT License'; no separate terms for any checkpoint were found | None | Cite Lewis et al. Science 2025 | confirmed | [Hugging Face model card microsoft/bioemu](https://huggingface.co/microsoft/bioemu) (2026-10-16); [Hugging Face API record (cardData.license, file list)](https://huggingface.co/api/models/microsoft/bioemu) (2026-10-16) |
| **AlphaFold2 parameters** alphafold_params_2021-07-14.tar | BioEmu computes its sequence embeddings with AlphaFold2 and downloads these weights on first use | Downloaded at run time by BioEmu (Colab) from storage.googleapis.com | No | CC-BY-4.0 for the AlphaFold parameters; Apache-2.0 for the AlphaFold code. The README's statement covers 'the AlphaFold parameters' without naming release files | None beyond attribution | Cite Jumper et al. Nature 2021 | confirmed | [google-deepmind/alphafold README, Model Parameters License](https://github.com/google-deepmind/alphafold#model-parameters-license) (2026-10-16) |
| **ColabFold MMseqs2 MSA server** api.colabfold.com | BioEmu sends each query sequence to this server for its MSA when sampling (Colab notebook) | External service called at run time | No; the MSAs stay in the sampling output | No licence or terms of service for the server were found. ColabFold's README asks that queries from a local computer be serial from a single IP and not spread over several computers. ColabFold code is MIT | Usage limits only, as stated in the README | Cite Mirdita et al. Nat Methods 2022 | not confirmed | [ColabFold README, FAQ on the MSA server](https://github.com/sokrypton/ColabFold#faq) (2026-10-16); [ColabFold LICENSE](https://github.com/sokrypton/ColabFold/blob/main/LICENSE) (2026-10-16) |
| **Google Colab** | Runs `notebooks/bioemu_resample.ipynb` on a free GPU; output written to the user's Google Drive | External service; the notebook only | Not applicable | Google Terms of Service plus the Colab Additional Terms of Service (last modified 14 May 2024) | No sharing of accounts and no resale of access; resource limits as set by Google |  | confirmed | [Colab Additional Terms of Service](https://research.google.com/colaboratory/tos_v5.html) (2026-10-16) |
| **Kaggle notebooks** | The notebook also detects Kaggle and writes its output there; an alternative runtime only | External service | Not applicable | Kaggle Terms of Use; the page renders with JavaScript and its text could not be read |  |  | not confirmed | [Kaggle Terms of Use](https://www.kaggle.com/terms) (2026-10-16) |
| **GitHub, GitHub Actions and GitHub Pages** | Repository hosting; CI (`.github/workflows/ci.yml`: actions/checkout, astral-sh/setup-uv, ruff, pytest); the project page served from `docs/` | External service | Not applicable | GitHub Terms of Service (effective 27 April 2026) and the additional terms for Actions and Pages. The two actions are MIT | Actions only for developing and testing the project; Pages not for commercial hosting; Pages sites at most 1 GB, soft limit 100 GB per month bandwidth. A public repository grants other users a licence to view and fork it through GitHub |  | confirmed | [GitHub Terms of Service](https://docs.github.com/en/site-policy/github-terms/github-terms-of-service) (2026-10-16); [Additional product terms (Actions, Pages)](https://docs.github.com/en/site-policy/github-terms/github-terms-for-additional-products-and-features) (2026-10-16); [GitHub Pages limits](https://docs.github.com/en/pages/getting-started-with-github-pages/github-pages-limits) (2026-10-16); [actions/checkout LICENSE](https://github.com/actions/checkout/blob/main/LICENSE) (2026-10-16); [astral-sh/setup-uv LICENSE](https://github.com/astral-sh/setup-uv/blob/main/LICENSE) (2026-10-16) |
| **uv** 0.11 | Creates the environments and runs the pipeline and CI | Installed by hand; in CI by astral-sh/setup-uv | No | Apache-2.0 OR MIT (the repository carries LICENSE-APACHE and LICENSE-MIT) | None |  | confirmed | [astral-sh/uv LICENSE-APACHE](https://github.com/astral-sh/uv/blob/main/LICENSE-APACHE) (2026-10-16) |
| **CPython** 3.11 (main, coverage, DENSS, BioEmu), 3.12 (simulation, screenshots) | Interpreter of every environment | Installed by uv | No | PSF-2.0 | None |  | confirmed | [Python licence page](https://docs.python.org/3/license.html) (2026-10-16) |
| **GCC** system compiler | Compiles PULCHRA (`02_fetch_sources.py`) | System package | No; the compiled PULCHRA binary stays in `SIM_WORK` | GPL-3.0 with the GCC Runtime Library Exception 3.1 | None for compiling; the exception leaves the compiled program under its own licence (MIT for PULCHRA) |  | confirmed | [GCC runtime library licence](https://gcc.gnu.org/onlinedocs/libstdc++/manual/license.html) (2026-10-16) |
| **Playwright and Chromium** playwright 1.63.0, Chromium build 1243 | Screenshots of the project page while checking it; not part of the pipeline | Python package in `/scratch/scfm-work/pw`; Playwright downloads its browser builds to a local folder | No; no screenshots are tracked | Playwright: Apache-2.0. Chromium: BSD-3-Clause, with the third-party notices shipped in the build's LICENSE.headless_shell | None |  | confirmed | [microsoft/playwright-python LICENSE](https://github.com/microsoft/playwright-python/blob/main/LICENSE) (2026-10-16); [Chromium LICENSE](https://chromium.googlesource.com/chromium/src/+/main/LICENSE) (2026-10-16) |
| **PeptoneBench code** | Not installed or imported. The score conventions in `src/bsc/saxs.py` and the Pepsi-SAXS flag in `analysis/simulations/common.py` follow it | Not obtained | No | Apache-2.0 | None | Cited as Invernizzi et al. 2025 | confirmed | [PeptoneLtd/peptonebench LICENSE](https://github.com/PeptoneLtd/peptonebench/blob/main/LICENSE) (2026-10-16) |
<!-- generated:tools:end -->

## Data sources

<!-- generated:data:start -->
| Source | How the project uses it | Licence or terms | Attribution | What the repo tracks from it | Raw data redistributed | Status | Checked at |
|---|---|---|---|---|---|---|---|
| **PeptoneBench archive, Zenodo 17306061** (published 2025-10-09; PeptoneDBs.tar.gz, Predictions.tar.gz) | The 439 SAXS profiles with sequences, pH and disorder scores (PeptoneDB-SAXS, curated from SASBDB), and the ensembles and Pepsi-SAXS curves of BioEmu-1 and six other models (AlphaFold2, ESMFold, Boltz-2, IDPFold, idpSAM, PepTron) | CC-BY-4.0 | Invernizzi et al. 2025 (bioRxiv 10.1101/2025.10.18.680935) and the Zenodo record; the profiles also carry SASBDB's attribution to their depositors | Derived results: scores, reweighting paths, radii of gyration, classes, figures. Also copied as columns: the construct sequences (`results/entries.csv`, `results/simulations/inputs.csv`, `notebooks/resample_list.csv`) and the SASBDB entry titles | No profiles, curves or ensembles. Construct sequences are copied as table columns, which CC-BY-4.0 permits with attribution | confirmed | [Zenodo record 17306061 (API, rights field)](https://zenodo.org/records/17306061) (2026-10-16) |
| **SASBDB data and REST API** (entry summaries fetched 2026-10-13) | Entry summaries of the 439 profiles: measurement conditions, molecular weight, Dmax, organisms, UniProt codes (`bsc fetch`, `analysis/coverage/01_fetch_sasbdb.py`, `18_sasbdb_metadata.py`) | "The data and models deposited in SASBDB are free of all copyright restrictions and made fully and freely available for both non-commercial and commercial use." No terms specific to the REST API were found | "Users of the data should attribute the original authors." Each entry is identified by its SASBDB code and the page links every code to its record; cite Kikhney et al. Protein Sci 2020 | `results/coverage/sasbdb_metadata.csv` (curve type, instrument, buffer, pH, concentration, organisms, molecule names, project titles and contributor names per entry) and entry titles in `results/candidates.csv` and `results/entries.csv` | No scattering profiles. Entry metadata are copied as tables, which the SASBDB statement permits with attribution | confirmed (data); not confirmed (API-specific terms) | [SASBDB About page](https://www.sasbdb.org/about/) (2026-10-16); [SASBDB policies (retraction policy only)](https://www.sasbdb.org/policies/) (2026-10-16) |
| **BioEmu MD training sets, Zenodo 15629740, 15641199, 15641184** (CATH1 and CATH2; octapeptides; MegaSim (published June 2025)) | Zip listings and the `topology.pdb` of each CATH and MegaSim wild-type system, read through HTTP range requests; no trajectories | CDLA-Permissive-2.0 | None required for Results (section 3.1: no restriction or obligation on the use, modification or sharing of Results); whoever shares the Data must include the agreement text. Cite Lewis et al. Science 2025 | Derived results only: nearest-training identifiers and identities (`results/coverage/training_nearest.csv`), counts and compositions (`training_summary.json`) | No | confirmed | [Zenodo record 15629740 (rights field)](https://zenodo.org/records/15629740) (2026-10-16); [Zenodo record 15641199 (rights field)](https://zenodo.org/records/15641199) (2026-10-16); [Zenodo record 15641184 (rights field)](https://zenodo.org/records/15641184) (2026-10-16); [CDLA-Permissive-2.0 text](https://cdla.dev/permissive-2-0/) (2026-10-16) |
| **BioEmu AFDB cluster definitions** (microsoft/bioemu dev/afdb_cluster_definitions at commit 6f059c4) | The 50,228 pre-training clusters and their 179,422 AFDB v4 members (`11_download_bioemu_training.py`) | MIT (the repository LICENSE at that commit; the folder's README states no other licence). The members are AlphaFold DB entries (CC-BY-4.0) | Cite Lewis et al. Science 2025 | Derived results only (cluster counts, compositions, nearest-member identifiers) | No | confirmed | [microsoft/bioemu LICENSE at 6f059c4](https://github.com/microsoft/bioemu/blob/6f059c41017a2c653ef1567ce79d7fa260ad3670/LICENSE) (2026-10-16); [dev/afdb_cluster_definitions/README.md at 6f059c4](https://github.com/microsoft/bioemu/tree/6f059c41017a2c653ef1567ce79d7fa260ad3670/dev/afdb_cluster_definitions) (2026-10-16) |
| **AlphaFold DB** (v4 entries) | Not downloaded. The BioEmu clusters are AFDB v4 entries identified by UniProt accession, and TED's domains are cut from AFDB models | CC-BY-4.0 ("Data is available for academic and commercial use, under a CC-BY-4.0 license"); use of the site is subject to the EMBL-EBI Terms of Use | Cite Jumper et al. Nature 2021 (the AlphaFold methods paper, as the licence page asks) | AFDB-derived identifiers and counts only | No | confirmed | [AlphaFold DB License and Disclaimer (PDF)](https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf) (2026-10-16) |
| **MEGAscale, Zenodo 7992926** (v2_230420; Tsuboyama2023_Dataset2_Dataset3_20230416.csv) | The curated stability table and its wild-type domains (one zip member read through HTTP range requests) | CC-BY-4.0 | Cite Tsuboyama et al. Nature 2023 and the Zenodo record | Derived results only | No | confirmed | [Zenodo record 7992926 (rights field)](https://zenodo.org/records/7992926) (2026-10-16) |
| **CATH** (v4.4.0 classification files; v4.3.0 S35 list (sensitivity check)) | Superfamily classification, S35 clusters and S100 sequences (`00_download_cath.sh`, `08_compare.py`, `14_annotate_training.py`) | CC-BY-4.0 for the data. The CATH wiki text itself is CC BY-NC-SA 4.0 | Attribute to I. Sillitoe, N. Dawson, T. Lewis, D. Lee, J. Lees, C. Orengo, based on work at cath.biochem.ucl.ac.uk | `results/coverage/cath_superfamily_rank.csv.gz` (superfamily identifiers, S35 cluster counts, rank and CATH superfamily names) and derived comparisons | No classification files. Superfamily names and cluster counts are copied as table columns, which CC-BY-4.0 permits with attribution | confirmed | [CATH licence page](https://www.cathdb.info/wiki?id=data:licence) (2026-10-16) |
| **InterPro, Pfam and CATH-Gene3D (InterPro API)** (InterPro 110.0, Pfam 38.2, CATH-Gene3D 4.3.0) | Domain locations per UniProt accession (`04_fetch_annotations.py`) | InterPro and Pfam data: CC0-1.0 ("All of the InterPro, Pfam, PRINTS and SFLD downloadable data ... CC0 1.0"). CATH-Gene3D is not named in that statement; CATH's own data are CC-BY-4.0 | Cite InterPro and the member databases (requested by InterPro); EMBL-EBI expects attribution | Gene3D superfamily identifiers and Pfam accessions per construct (`results/coverage/per_protein_mapping.csv`, `results/arm1/relatedness.csv`) | No | confirmed (InterPro, Pfam); not confirmed (CATH-Gene3D as served by InterPro) | [InterPro licence page](https://interpro-documentation.readthedocs.io/en/latest/license.html) (2026-10-16); [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/) (2026-10-16) |
| **UniProtKB, Swiss-Prot, UniParc, ID mapping, peptide search, taxonomy** (release 2026_03) | Mapping constructs to accessions, sequences of the AFDB cluster members and named MD proteins, the Swiss-Prot reference set, organism lineages | CC-BY-4.0 for all copyrightable parts of the databases; no separate terms for the REST services were found | Cite UniProt as on its citation page | Accessions, isoforms, ranges, organism and lineage per construct (`results/coverage/per_protein_mapping.csv`) and derived comparisons | No flat files or sequences from UniProt | confirmed | [UniProt License & disclaimer (help page)](https://rest.uniprot.org/help/license) (2026-10-16) |
| **RCSB PDB: search API, data API, clusters-by-entity-30, FASTA** (cluster file of 4 October 2026) | Sequence search for unmapped constructs, polymer-entity data, the 30% clusters with their representatives, and the FASTA of the twelve DESRES fast-folder PDB entries | CC0-1.0 for PDB archive files and data from the RCSB PDB APIs | Encouraged, not required: the original authors of the structures | PDB identifiers and derived counts | No | confirmed | [RCSB PDB usage policy](https://www.rcsb.org/pages/usage-policy) (2026-10-16) |
| **PDBe SIFTS (best_structures API)** (live API, weekly update) | Structure coverage of each construct (`04_fetch_annotations.py`, `07_build_protein_table.py`) | No SIFTS-specific licence was found. The EMBL-EBI Terms of Use place no restrictions beyond those of the original data owners and expect attribution | The SIFTS page asks users to cite Dana et al. NAR 2019 (or Velankar et al. NAR 2013) | Derived coverage per construct (`has_pdb_50` and coverage fractions in `results/coverage/per_protein_mapping.csv`) | No | not confirmed | [SIFTS page](https://www.ebi.ac.uk/pdbe/docs/sifts/) (2026-10-16); [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/) (2026-10-16) |
| **TED (The Encyclopedia of Domains)** (v1, ted.cathdb.info API) | CATH labels of domains in the AFDB model of each accession (`04_fetch_annotations.py`) | The TED data release on Zenodo (record 10848710) is CC-BY-4.0. The API site renders with JavaScript and no terms could be read from it | Cite Lau et al. Science 2024 and the Zenodo record | TED superfamily labels per construct (`ted_superfamilies` in `results/coverage/per_protein_mapping.csv`) | No | confirmed (Zenodo release); not confirmed (API terms) | [Zenodo record 10848710 (rights field)](https://zenodo.org/records/10848710) (2026-10-16); [TED site](https://ted.cathdb.info/) (2026-10-16) |
| **CALVADOS training and validation protein lists** (KULL-Centre/papers at 85979b1; KULL-Centre/_2024_Cao_CALVADOSCOM at 51564e0) | The proteins the CALVADOS 1, 2 and 3 parameters were trained and validated on, for the circularity check (`02_fetch_sources.py`, `03_circularity.py`) | CALVADOS 1 lists (papers/2021/CG-IDPs-Tesei-et-al): CC-BY-4.0, from the LICENSE file in that folder. CALVADOS 2 lists (papers/2022/CG-cutoffs-Tesei-et-al): no licence file in the folder or at the repository root; the folder's README links Zenodo DOI 10.5281/zenodo.6815067, whose current version (record 7437501, an archive of the folder) is CC-BY-4.0; I did not check that the archive holds the same two files. CALVADOS 3 lists (_2024_Cao_CALVADOSCOM): MIT, Copyright (c) 2024 Linderstrøm-Lang Centre for Protein Science, University of Copenhagen | CC-BY-4.0: cite Tesei et al. PNAS 2021 and Tesei & Lindorff-Larsen 2022. MIT: keep the copyright and permission notice with copies of the files. Cite Cao et al. Protein Sci 2024 | `results/simulations/calvados_sets.csv`: set, role, file, protein name and sequence copied from the eight list files; `circularity.csv` derived | Yes, names and sequences. Permitted for the CC-BY-4.0 and MIT files with attribution (given here); for the CALVADOS 2 files it rests on the Zenodo deposit | confirmed (CALVADOS 1, 3); not confirmed (CALVADOS 2 GitHub copy) | [papers/2021/CG-IDPs-Tesei-et-al/LICENSE at 85979b1](https://github.com/KULL-Centre/papers/blob/85979b1b4123b6b5391b617d16551969eda9f56e/2021/CG-IDPs-Tesei-et-al/LICENSE) (2026-10-16); [papers/2022/CG-cutoffs-Tesei-et-al at 85979b1 (README, no LICENSE)](https://github.com/KULL-Centre/papers/tree/85979b1b4123b6b5391b617d16551969eda9f56e/2022/CG-cutoffs-Tesei-et-al) (2026-10-16); [Zenodo record 7437501 (rights field)](https://zenodo.org/records/7437501) (2026-10-16); [_2024_Cao_CALVADOSCOM LICENSE at 51564e0](https://github.com/KULL-Centre/_2024_Cao_CALVADOSCOM/blob/51564e0edd126c44d6a7299fbe4414c579394ef4/LICENSE) (2026-10-16) |
| **Folding@home COVID-19 cryptic-pockets bucket** (s3://fah-public-data-covid19-cryptic-pockets) | Only the object listing, to count the simulated systems of one of BioEmu's MD sets (`11_download_bioemu_training.py`) | CC0 (AWS Open Data Registry entry for the Folding@home COVID-19 datasets) |  | One count (`n_exascale_systems_listed` in `results/coverage/versions.json`) | No | confirmed | [awslabs/open-data-registry foldingathome-covid19.yaml](https://github.com/awslabs/open-data-registry/blob/main/datasets/foldingathome-covid19.yaml) (2026-10-16) |
| **DESRES fast-folder simulations** | Cited only. No DESRES data are downloaded or used; the twelve proteins' sequences come from their PDB entries at RCSB (row above) | Not reviewed, since nothing is obtained from DESRES |  | Nothing from DESRES | No | not applicable | [the code that fetches the twelve PDB sequences](https://github.com/Natalija-Stepurko/bioemu-saxs-campaign/blob/HEAD/analysis/coverage/11_download_bioemu_training.py) (2026-10-16) |
<!-- generated:data:end -->

## Python packages

### Direct dependencies

<!-- generated:python-direct:start -->
| Package | Environment | Version | Licence | Review | Homepage |
|---|---|---|---|---|---|
| joblib | main | 1.6.0 | BSD-3-Clause |  | https://joblib.readthedocs.io |
| matplotlib | main | 3.11.2 | LicenseRef-Matplotlib (PSF-based, permissive) |  | https://matplotlib.org |
| mdtraj | main | 1.11.1.post2 | LGPL-2.1-or-later | copyleft | http://mdtraj.org |
| mrcfile | main | 1.5.4 | BSD-3-Clause |  | https://github.com/ccpem/mrcfile |
| numpy | main | 2.4.6 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |  | https://numpy.org |
| pandas | main | 3.0.6 | BSD-3-Clause |  | https://pandas.pydata.org |
| pyarrow | main | 25.0.1 | Apache-2.0 |  | https://arrow.apache.org/ |
| pytest | main | 9.1.1 | MIT |  | https://docs.pytest.org/en/latest/ |
| requests | main | 2.34.2 | Apache-2.0 |  | https://github.com/psf/requests |
| ruff | main | 0.16.10 | MIT |  | https://github.com/astral-sh/ruff |
| scikit-learn | main | 1.9.1 | BSD-3-Clause |  | https://scikit-learn.org |
| scipy | main | 1.17.1 | BSD-3-Clause |  | https://scipy.org/ |
| biopython | coverage | 1.88 | LicenseRef-Biopython-License-Agreement |  | https://biopython.org/ |
| matplotlib | coverage | 3.11.2 | LicenseRef-Matplotlib (PSF-based, permissive) |  | https://matplotlib.org |
| metapredict | coverage | 3.1.1 | MIT |  | https://metapredict.net |
| numpy | coverage | 2.4.6 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |  | https://numpy.org |
| pandas | coverage | 3.0.6 | BSD-3-Clause |  | https://pandas.pydata.org |
| requests | coverage | 2.34.2 | Apache-2.0 |  | https://github.com/psf/requests |
| scipy | coverage | 1.17.1 | BSD-3-Clause |  | https://scipy.org/ |
| torch | coverage | 2.14.1+cpu | Apache-2.0 AND Apache-2.0 WITH LLVM-exception AND BSD-2-Clause AND BSD-3-Clause AND BSL-1.0 AND MIT |  | https://pytorch.org |
| biopython | simulation | 1.88 | LicenseRef-Biopython-License-Agreement |  | https://biopython.org/ |
| calvados | simulation | 0.9.0 | GPL-3.0-only | copyleft | https://github.com/KULL-Centre/CALVADOS |
| Jinja2 | simulation | 3.1.6 | BSD-3-Clause |  | https://github.com/pallets/jinja/ |
| matplotlib | simulation | 3.11.2 | LicenseRef-Matplotlib (PSF-based, permissive) |  | https://matplotlib.org |
| MDAnalysis | simulation | 2.10.0 | LGPL-3.0-or-later | copyleft | https://www.mdanalysis.org |
| mdtraj | simulation | 1.11.2 | LGPL-2.1-or-later | copyleft | https://mdtraj.org |
| numba | simulation | 0.68.0 | BSD-2-Clause |  | https://numba.pydata.org |
| numpy | simulation | 2.5.4 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |  | https://numpy.org |
| OpenMM | simulation | 8.4.0.post2 | MIT (API, Reference and CPU platforms) AND LGPL (CUDA and OpenCL platforms) | copyleft | https://openmm.org |
| pandas | simulation | 3.0.6 | BSD-3-Clause |  | https://pandas.pydata.org |
| pydantic | simulation | 2.14.1 | MIT |  | https://github.com/pydantic/pydantic |
| PyYAML | simulation | 6.0.3 | MIT |  | https://github.com/yaml/pyyaml |
| requests | simulation | 2.34.2 | Apache-2.0 |  | https://github.com/psf/requests |
| scipy | simulation | 1.18.1 | BSD-3-Clause |  | https://scipy.org/ |
| statsmodels | simulation | 0.15.0 | BSD-3-Clause |  | https://www.statsmodels.org |
| tqdm | simulation | 4.70.1 | MPL-2.0 AND MIT | copyleft | https://tqdm.github.io |
| denss | denss | 1.8.7 | GPL-3.0 | copyleft |  |
| numpy | denss | 1.26.4 | BSD-3-Clause |  | https://numpy.org |
| bioemu | bioemu | 1.4.1 | MIT |  | https://github.com/microsoft/bioemu |
| pandas | bioemu | 3.0.6 | BSD-3-Clause |  | https://pandas.pydata.org |
| playwright | screenshots | 1.63.0 | Apache-2.0 |  | https://github.com/Microsoft/playwright-python |
<!-- generated:python-direct:end -->

### Packages for review in any environment

<!-- generated:python-review:start -->
| Package | Licence | Review | Environments (version) |
|---|---|---|---|
| calvados | GPL-3.0-only | copyleft | simulation 0.9.0 (direct) |
| certifi | MPL-2.0 | copyleft | main 2026.7.22, coverage 2026.7.22, simulation 2026.7.22, bioemu 2026.7.22 |
| denss | GPL-3.0 | copyleft | denss 1.8.7 (direct) |
| GridDataFormats | LGPL-3.0-or-later | copyleft | simulation 1.2.0 |
| MDAnalysis | LGPL-3.0-or-later | copyleft | simulation 2.10.0 (direct) |
| mdtraj | LGPL-2.1-or-later | copyleft | main 1.11.1.post2 (direct), simulation 1.11.2 (direct), bioemu 1.11.1.post2 |
| OpenMM | MIT (API, Reference and CPU platforms) AND LGPL (CUDA and OpenCL platforms) | copyleft | simulation 8.4.0.post2 (direct) |
| tqdm | MPL-2.0 AND MIT | copyleft | coverage 4.70.1, simulation 4.70.1 (direct), bioemu 4.70.1 |
<!-- generated:python-review:end -->
