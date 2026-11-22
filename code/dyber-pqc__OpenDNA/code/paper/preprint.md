---
title: "OpenDNA: A Free, Open-Source Protein Engineering Platform for Consumer Hardware"
authors:
  - name: Zachary Kleckner
    affiliation: Dyber, Inc.
    email: zkleckner@dyber.org
    orcid: ""
date: 2026-04-06
keywords:
  - protein engineering
  - structure prediction
  - inverse folding
  - drug discovery
  - open source
  - democratization
  - ESMFold
  - ESM-IF1
license: CC-BY-4.0
---

# OpenDNA: A Free, Open-Source Protein Engineering Platform for Consumer Hardware

## Abstract

We present **OpenDNA**, a free, open-source protein engineering platform that combines structure prediction, inverse folding sequence design, comprehensive molecular property analysis, mutation effect prediction, molecular dynamics, and 3D visualization in a single integrated desktop application. Unlike commercial protein engineering suites that cost ~$200,000 per seat per year and require workstation-class hardware, OpenDNA runs on a gaming laptop and is accessible to anyone with a web browser. The platform integrates state-of-the-art open models — ESMFold for structure prediction, ESM-IF1 for inverse folding, Boltz for multimer prediction, and OpenMM for molecular dynamics — alongside ~50 sequence- and structure-based analyses that replicate the functionality of commercial QikProp, MolProbity, PROPKA, TANGO, and SignalP. We also introduce a **language-model agent framework** with tool calling that lets users issue natural-language goals (e.g. *"Analyze ubiquitin and predict the effect of K48R"*) and have the system plan, execute, and report multi-step protein engineering workflows. OpenDNA is licensed Apache 2.0 with Commons Clause and is available at https://github.com/dyber-pqc/OpenDNA.

## 1. Introduction

Computational protein engineering has been transformed in the last five years by deep learning. AlphaFold2 [^jumper2021], ESMFold [^lin2023], and RoseTTAFold [^baek2021] now produce experimentally-accurate structure predictions in minutes; ProteinMPNN [^dauparas2022] and ESM-IF1 [^hsu2022] generate novel sequences that fold into desired backbones; and RFdiffusion [^watson2023] enables the de novo design of binders, enzymes, and symmetric assemblies.

These advances should be in everyone's hands, but they are not. The dominant commercial protein engineering platforms — Schrödinger Maestro/Prime/BioLuminate, Discovery Studio, MOE — cost ~$50,000–$200,000 per seat per year, lock features behind paywalls, run only on workstations or HPC clusters, and require days of training to use productively. Academic alternatives like Rosetta exist, but they are command-line driven, written in C++, and steep enough to learn that most undergraduates never get past the first tutorial.

The result is a stark imbalance: a high school student in Lagos curious about protein design has fundamentally different tools available than a postdoc at Genentech. Even at academic institutions, lab affiliation determines whether a graduate student can use the best methods. For citizen scientists, open educators, and underfunded labs, the gap is effectively unbridgeable.

We built **OpenDNA** to close that gap.

## 2. Design Principles

OpenDNA was designed around five principles:

1. **Consumer hardware first.** Every default workflow must run on a $1,000 gaming laptop. GPU is optional, not required.
2. **Plain English where possible.** Non-experts should be able to issue natural-language commands and get actionable results.
3. **Transparent and explainable.** Every prediction shows confidence intervals, citations, and the reasoning behind it.
4. **Local-first and private.** No data leaves the user's machine unless they explicitly opt in.
5. **Free and open.** Apache 2.0 licensed core; no feature gating; no telemetry; no surveillance.

## 3. Architecture

OpenDNA is a multi-process system:

- A **Python backend** (FastAPI + Uvicorn) hosts ~20 analysis engines and orchestrates ML model inference.
- A **React + TypeScript frontend** (Vite + Molstar) provides a desktop-class web UI with 3D visualization.
- A **Rust core** (PyO3 bindings) handles data models, PDB parsing, version control, and SQLite metadata storage.
- An optional **Tauri shell** packages the entire system as a native installer for Windows, macOS, and Linux.
- A **Python SDK** (`opendna.sdk.Client`) provides a typed client for programmatic use.

The platform exposes ~38 REST API endpoints and a YAML workflow engine for declarative reproducible pipelines.

### 3.1 Compute engines

| Engine | Method | Reference |
|---|---|---|
| Folding | ESMFold v1 | Lin et al. 2023 [^lin2023] |
| Inverse folding | ESM-IF1 | Hsu et al. 2022 [^hsu2022] |
| Multimer | Boltz-1 (optional) | Wohlwend et al. 2024 [^boltz] |
| MD | OpenMM (AMBER14 + TIP3P) | Eastman et al. 2017 [^openmm] |
| Docking | Heuristic + DiffDock optional | Corso et al. 2023 [^diffdock] |
| Quantum chem | xTB optional | Bannwarth et al. 2019 [^xtb] |

### 3.2 Sequence-based analyses

OpenDNA implements heuristic equivalents of commercial property suites:

- **Sequence properties**: molecular weight, isoelectric point (bisection on Henderson-Hasselbalch), GRAVY hydropathy, instability index (Guruprasad et al. 1990), aliphatic index (Ikai 1980), N-terminal half-life via the N-end rule (Bachmair et al. 1986), aromaticity, extinction coefficient (Pace et al.).
- **Drug-likeness**: Lipinski's Rule of Five (Lipinski et al. 1997).
- **Hydropathy**: Kyte-Doolittle sliding window profile (Kyte & Doolittle 1982).
- **Disorder**: IUPred-style propensity scoring (Linding et al. 2003).
- **Transmembrane prediction**: TMHMM-style hydrophobicity scanning.
- **Signal peptide**: SignalP-style n/h/c-region scoring.
- **Aggregation**: TANGO/Aggrescan-style propensity (Pawar et al.).
- **PTM sites**: kinase consensus motif matching (PKA/PKC/CK2/CDK/GSK3/MAPK), N-/O-glycosylation sequon scanning.
- **Mutation effects**: heuristic ΔΔG estimation.
- **Conservation**: per-residue ESM-2 perplexity (mask each position, query the model).
- **Pairwise alignment**: Needleman-Wunsch with BLOSUM62.

### 3.3 Structure-based analyses

When a 3D structure is available, OpenDNA computes:

- Secondary structure (DSSP-like phi/psi assignment)
- Ramachandran plot (full phi/psi calculation)
- Radius of gyration
- Solvent-accessible surface area (Shrake-Rupley simplified)
- Binding pocket detection (cavity scoring)
- Bond network: hydrogen bonds, salt bridges, disulfides
- pKa shifts (PROPKA-style, with desolvation and Coulombic contributions)
- Structure validation (MolProbity-style: Ramachandran outliers, steric clashes, bond geometry)
- Pharmacophore feature extraction (donors, acceptors, charges, aromatics, hydrophobics)
- Antibody numbering and CDR detection (Kabat/Chothia heuristic)
- RMSD via Kabsch superposition
- MMGBSA-style binding energy estimation

### 3.4 Iterative design loop

A core differentiating feature is the **iterative design loop**, which automates the cycle that protein engineers normally do manually:

1. Fold the input sequence with ESMFold
2. Generate N candidate alternatives via ESM-IF1
3. Score each candidate with the composite quality scorer
4. Fold the highest-scoring candidate to confirm
5. If confirmed better, keep as new best
6. Repeat for N rounds

This automates what is typically a multi-day workflow in commercial tools and produces a clear optimization trajectory with provenance.

### 3.5 LLM agent framework

OpenDNA includes a unified LLM interface that auto-detects available providers:

1. **Ollama** (local, free, private) — preferred when running tool-capable models like llama3.2:3b
2. **Anthropic API** (if `ANTHROPIC_API_KEY` is set)
3. **OpenAI API** (if `OPENAI_API_KEY` is set)
4. **Heuristic fallback** (always available, regex-based intent parser)

The agent framework defines 13 tools that map high-level capabilities to OpenDNA functions. Given a natural-language goal, the agent plans a sequence of tool calls, executes them, and returns a final answer with full trace. For models that do not support tool calling (e.g., Mistral 7B base), the system falls back to plain LLM chat for general questions.

This is the first protein engineering platform we are aware of that lets users issue commands like *"Import ubiquitin, predict the effect of K48R, and tell me if the mutation is stabilizing"* and have the system actually execute the multi-step workflow.

## 4. Validation

OpenDNA ships with a **self-benchmarking framework** that compares its ESMFold predictions against reference AlphaFold DB structures (which use the more-accurate AlphaFold2 with multiple sequence alignments). The reference set includes ubiquitin, insulin, lysozyme, myoglobin, and GFP.

Metrics computed:

- **RMSD** via Kabsch superposition over Cα atoms
- **TM-score proxy** with the standard length-dependent d0 normalization
- **Mean pLDDT** confidence
- **Inference time**

Preliminary results on consumer hardware (RTX 3070, 32 GB RAM) show RMSDs of 1.5–4.5 Å against AlphaFold2 references for proteins under 250 residues, with inference times of 5–60 seconds. This is consistent with published ESMFold benchmarks (~5–10% accuracy gap behind AlphaFold2 with MSA).

All benchmarks are reproducible: the runner is part of the platform and can be invoked via `opendna benchmark run` or the `/v1/benchmark` API endpoint.

## 5. Differentiation

| Feature | Schrödinger Suite | Rosetta | **OpenDNA** |
|---|---|---|---|
| Annual cost | ~$200,000 | Free for academia | **Free, all users** |
| Structure prediction | Prime | RoseTTAFold | ESMFold (built-in) |
| Sequence design | Prime | ProteinMPNN | ESM-IF1 (built-in) |
| Property suite | QikProp | Manual | ~50 analyses (built-in) |
| MD | Desmond | RosettaMD | OpenMM (built-in) |
| Docking | Glide | RosettaLigand | DiffDock-ready (heuristic) |
| Drug-likeness | ✅ | Manual | ✅ |
| LLM agent | ❌ | ❌ | **✅ (first in class)** |
| Plain English UI | ❌ | ❌ | ✅ |
| Runs on a laptop | ❌ | ❌ | ✅ |
| 3D viewer | Maestro | PyMOL | Molstar (built-in) |
| Open source | ❌ | Partial | **✅ Apache 2.0** |
| Workflow YAML | ✅ | Manual | ✅ |
| Python SDK | ❌ | Limited | ✅ |

## 6. Limitations and Future Work

OpenDNA is in active development. Current limitations:

- **Single-sequence ESMFold** is ~5–10% less accurate than full AlphaFold2 with MSA for novel sequences. For known proteins, OpenDNA fetches the AlphaFold DB structure directly, sidestepping this gap.
- **DiffDock integration is currently heuristic.** Real DiffDock requires installing the diffdock pip package and ~5 GB of model weights.
- **Multimer prediction via Boltz** is supported but optional; the default fallback folds chains independently.
- **No multi-user authentication yet.** The API listens on localhost; networked deployments require user-supplied auth.
- **No cloud burst.** All compute is local in the current version.
- **Tested mostly on monomeric proteins under 500 residues.** Larger systems and complex multimers need further validation.

Planned in v0.5+:
- Real DiffDock as a default install option
- AlphaFold-Multimer / Boltz as the default multimer predictor
- Multi-user authentication and SSO
- Reproducibility benchmarks against published CASP15/CAMEO results
- LTS release line for production users
- Tauri desktop installer with bundled Python sidecar

## 7. Availability

OpenDNA is freely available at **https://github.com/dyber-pqc/OpenDNA** under Apache 2.0 + Commons Clause. Documentation, tutorials, and recipes are at https://dyber-pqc.github.io/OpenDNA/. The Python SDK is `pip install opendna`.

Tested platforms: Windows 10/11, macOS 11+, Ubuntu 20.04+. Requires Python 3.10+ and a modern web browser.

## 8. Author Contributions

Z.K. designed the platform, wrote all the source code, and authored the manuscript.

## 9. Acknowledgments

OpenDNA stands on the shoulders of the open-source ML community. We thank the developers of ESMFold and ESM-IF1 (Meta AI), Boltz (MIT/Wohlwend lab), OpenMM, Molstar, Biotite, FastAPI, Tauri, and many others.

Special thanks to the maintainers of AlphaFold DB for making 200 million high-quality predicted structures freely available — these are the gold standard against which OpenDNA is measured.

## 10. Competing Interests

The author declares no competing financial interests.

---

## References

[^jumper2021]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Jumper%2C J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.

[^lin2023]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Lin%2C Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure. *Science* 379, 1123–1130.

[^baek2021]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Baek%2C M. et al. (2021). Accurate prediction of protein structures and interactions using a three-track neural network. *Science* 373, 871–876.

[^dauparas2022]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Dauparas%2C J. et al. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science* 378, 49–56.

[^hsu2022]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Hsu%2C C. et al. (2022). Learning inverse folding from millions of predicted structures. *ICML 2022*.

[^watson2023]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Watson%2C J.L. et al. (2023). De novo design of protein structure and function with RFdiffusion. *Nature* 620, 1089–1100.

[^boltz]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Wohlwend%2C J. et al. (2024). Boltz-1: Democratizing biomolecular interaction modeling. *bioRxiv* 2024.11.19.624167.

[^openmm]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Eastman%2C P. et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Comp. Bio.* 13, e1005659.

[^diffdock]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Corso%2C G. et al. (2023). DiffDock: Diffusion steps, twists, and turns for molecular docking. *ICLR 2023*.

[^xtb]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/Bannwarth%2C C., Ehlert, S., & Grimme, S. (2019). GFN2-xTB. *J. Chem. Theory Comput.* 15, 1652–1671.
