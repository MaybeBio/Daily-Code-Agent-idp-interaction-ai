---
title: "OpenDNA v0.5.0: An Open, Reproducible, Post-Quantum-Secured Workbench for Computational Biology"
author:
  - name: OpenDNA Contributors
    affiliation: dyber-pqc
date: 2026-04-08
bibliography: references.bib
---

## Abstract

Computational structural biology has made extraordinary progress since
AlphaFold2 [@jumper2021alphafold], but the infrastructure that lets a working
scientist turn a protein sequence into a structure, into a design, into an
ordered DNA plasmid still lives across a dozen incompatible command-line
tools, notebooks, and web servers. We present **OpenDNA v0.5.0**, an
open-source desktop and server workbench that composes best-in-class
folding, design, docking, and simulation engines behind a single
provenance-tracking, post-quantum-secured API. v0.5.0 ships 19 coordinated
subsystems ("phases") covering component discovery, post-quantum
authentication, encrypted workspaces, a visual workflow graph, full run
provenance with blame and bisect, hardware selection, external-vendor
integration, a scientist-facing notebook, an academy, compliance reporting,
and an in-process self-benchmarking suite validated against CASP15 targets.
We describe the architecture, the methods that let us compose heterogeneous
engines while preserving reproducibility, and a benchmark showing that our
default ESMFold path reaches within 0.05 TM-score of AlphaFold2 on easy
CASP15 targets at an order of magnitude less wall-clock. OpenDNA is
released under a permissive licence, with a published software bill of
materials, signed builds, and a post-quantum identity layer aligned with
recent NIST standardisation [@nist2022pqc].

## 1 Introduction

Structural biology tooling has been transformed by deep learning. Predicting
a protein structure, once a multi-month campaign, is now a single call to
AlphaFold [@jumper2021alphafold], ESMFold [@lin2023esmfold], RoseTTAFold
[@baek2021rosettafold], or a hosted service such as ColabFold
[@mirdita2022colabfold]. De novo backbone generation with RFdiffusion
[@watson2023rfdiffusion], sequence design with ProteinMPNN
[@dauparas2022proteinmpnn] or ESM-IF1 [@hsu2022esmif1], and blind docking
with DiffDock [@corso2023diffdock] put capabilities that used to require a
dedicated methods lab on the desktop of any graduate student with a decent
GPU. Open complex-prediction models such as Boltz [@wohlwend2024boltz] and
AlphaFold 3 [@abramson2024alphafold3] continue that trajectory.

And yet the **experience** of doing science with these tools remains
disjointed. The typical project touches several of: a Colab notebook for
folding, a Rosetta [@leaver2011rosetta3] script for relaxation, a PyMOL
[@pymol2015] session to look at the result, a spreadsheet to track
mutations, a manual email to a DNA synthesis vendor, and a backup drive
whose provenance nobody trusts. The field has tools but no *workbench*.

OpenDNA is an attempt at a workbench. We treat folding, design, docking,
and simulation as composable **engines** behind a uniform API; attach a
cryptographically-signed **provenance graph** to every run; expose a
**visual workflow editor** for multi-step pipelines; gate everything
behind a **post-quantum identity layer**; and ship a **self-benchmarking**
suite that validates the platform against community-standard targets each
release. The v0.5.0 release collects this work into a single stable version
consisting of 19 tightly integrated phases.

This paper describes the architecture (Section 3), the methods that make
composition across engines tractable (Section 4), the validation we carry
out against CASP15 targets (Section 5), and the things we explicitly did
not solve (Section 6).

## 2 Related Work

The space OpenDNA sits in is crowded at the component level and almost
empty at the integration level.

**Folding and prediction.** ColabFold [@mirdita2022colabfold] wraps
AlphaFold2 with an MMseqs2 [@steinegger2017mmseqs2] server and a friendly
Colab interface, but it is a notebook rather than a workbench. ESMFold
[@lin2023esmfold] removes the MSA dependency and runs at a latency
compatible with an interactive GUI. RoseTTAFold [@baek2021rosettafold] and
AlphaFold 3 [@abramson2024alphafold3] extend the model family. Boltz
[@wohlwend2024boltz] is our default choice for open complex prediction.

**Design.** RFdiffusion [@watson2023rfdiffusion] generates backbones and is
typically paired with ProteinMPNN [@dauparas2022proteinmpnn] or ESM-IF1
[@hsu2022esmif1] for sequence design. These tools are distributed as
command-line scripts with conflicting dependency pins, which is one of the
problems OpenDNA addresses.

**Simulation and chemistry.** OpenMM [@eastman2017openmm] provides
GPU-accelerated classical MD. xTB [@bannwarth2019gfn2xtb] offers
semi-empirical QM, and ANI-2x [@devereux2020ani2x] provides a neural
potential for organic molecules.

**Viewers and UI.** PyMOL [@pymol2015] and Mol\* [@sehnal2021molstar] remain
the structural-biology viewers of record. OpenDNA embeds Mol\* in its
Tauri-based desktop frontend. The workflow editor uses React Flow
[@reactflow2019], and collaborative editing is built on Yjs [@yjs2020].

**Security.** Post-quantum cryptographic algorithms are provided via
liboqs [@liboqs2023], tracking the NIST standardisation timeline
[@nist2022pqc]. Our SBOM follows the CycloneDX specification
[@cyclonedx2023].

No existing open-source project we are aware of combines all of these
pieces behind a single provenance-tracked, reproducible API. That gap is
what OpenDNA tries to fill.

## 3 Architecture

OpenDNA v0.5.0 is organised as 19 phases that evolved during the release.
They are not runtime layers; they are subsystems that share a single
process boundary and a single component registry. The phases are:

1. **Component registry.** A typed, introspectable registry of every model
   or tool the workbench can call, including licence, size, and citation.
2. **Post-quantum auth.** Identity, token issuance, and validation using
   ML-KEM / ML-DSA primitives via liboqs [@liboqs2023].
3. **Encrypted workspaces.** Per-user, password-derived keys protect
   project files on disk.
4. **Priority queue.** Async job runner that serialises access to scarce
   GPU resources.
5. **Reliability.** Retries with jitter, crash reporting, and secret
   scrubbing.
6. **Provenance.** A content-addressed DAG of every step, with diff, blame,
   and bisect operators.
7. **Workflow graph runner.** Execution engine for the React Flow editor.
8. **External vendors.** Abstract DNA-synthesis quoting and webhook
   delivery.
9. **Notebook.** Markdown entries plus SVG/figure/GIF exporters.
10. **Data sources.** Unified access to PDB [@berman2000pdb], UniProt, and
    AlphaFold DB [@varadi2022alphafolddb; @tunyasuvunakool2021alphafold_human].
11. **Academy.** Structured learning levels and badges.
12. **Compliance.** CycloneDX SBOM [@cyclonedx2023], privacy, HIPAA and
    GDPR self-checks.
13. **Benchmarking.** In-process comparison against CASP15
    [@kryshtafovych2023casp15].
14. **Hardware advisor.** Chooses CPU/GPU targets per job.
15. **LLM gateway.** Local Ollama integration for natural-language help.
16. **Collaboration.** Yjs-backed shared editing [@yjs2020].
17. **Runtime.** Tauri desktop + FastAPI sidecar with health checks.
18. **CLI / SDK.** Thin clients over the same API.
19. **Docs and science.** Auto-generated model cards and science notes.

Each phase exposes a small Python module under `python/opendna/`. The
server process (`opendna.api.server`) is the only surface the desktop
frontend and the CLI call; every feature in the paper is accessible
through its REST endpoints.

## 4 Methods

### 4.1 Engine composition

Every engine (folding, design, docking, MD) implements a narrow Python
interface: a callable that takes a typed input, a progress callback, and a
hardware preference, and returns a typed output plus a provenance record.
Because each engine declares its inputs and outputs in the component
registry, the workflow graph runner can statically check that a
backbone-generation node (RFdiffusion) feeds a sequence-design node
(ESM-IF1 or ProteinMPNN) only along compatible edges. This turns the
"run a script, rename files, run another script" pattern into a typed DAG.

### 4.2 Provenance

Every engine call emits a provenance node with a content hash of its
inputs, outputs, and parameters. Nodes form a DAG per project. On top of
that we expose three operators:

- **diff**, which enumerates residue mutations between two structures and
  surfaces parameter changes;
- **blame**, which walks the DAG from a residue in the current structure
  back to the node(s) that introduced it;
- **bisect**, which finds the first node whose score dropped below a
  threshold, analogous to `git bisect`.

These operators mean a scientist can ask "why does residue 42 look wrong?"
and get a chain of steps, not just a file.

### 4.3 Post-quantum authentication

Credentials are stored in an encrypted SQLite database keyed by a
password-derived value. Session tokens are signed with ML-DSA (Dilithium);
the server validates every request against the public key stored in the
user's identity record. We chose post-quantum primitives now so that the
workbench's long-lived project files and signed workflow exports remain
verifiable against future adversaries [@nist2022pqc]. Key exchange for
encrypted workspace sync uses ML-KEM (Kyber). Implementation is via liboqs
[@liboqs2023] with a thin Python wrapper.

### 4.4 Benchmark methodology

Section 5 describes the results; the method is simple. For each target the
suite folds the sequence with ESMFold, downloads a reference PDB, aligns
Cα atoms by Kabsch superposition [@kabsch1976rmsd], and reports
RMSD and a length-normalised TM-score [@zhang2004tmscore]. Reference
TM-scores are curated from CASP15 public-server results
[@kryshtafovych2023casp15] and stored alongside the test case so that any
reader can update them.

## 5 Validation

We run the benchmark suite against a curated subset of 10 CASP15 monomer
targets, split evenly between easy single-domain cases (~80–120 aa) and
hard free-modelling targets (~180 aa). Table 1 summarises the aggregate
result and compares it to the AlphaFold2 / Rosetta CASP15 medians. Detailed
per-target numbers are generated on every release by
`python -m opendna.benchmarks.casp` and published at
`docs/BENCHMARKS.md`.

| Metric                     | OpenDNA (ESMFold) | AF2 / Rosetta median |
|----------------------------|-------------------|----------------------|
| Mean TM-score (easy)       | ~0.90             | ~0.92                |
| Mean TM-score (hard)       | ~0.55             | ~0.63                |
| Mean Cα RMSD               | 1.5 – 5 Å         | 1.0 – 4 Å            |
| Time-to-solution (GPU)     | seconds           | minutes – hours      |

*Table 1: Aggregate CASP15 results. Reference numbers are conservative
medians of CASP15 public-server submissions and are carried in the
benchmark source so they can be refreshed per release.*

The pattern is the one the community reported for ESMFold at launch
[@lin2023esmfold]: https://github.com/dyber-pqc/OpenDNA/blob/HEAD/paper/single-sequence folding gives up a small amount of
accuracy on easy targets and a larger amount on hard free-modelling cases
in exchange for two orders of magnitude less wall-clock. The benefit of
putting this inside a workbench is that the scientist can decide, per
target, whether the latency-quality trade is worth running ColabFold
[@mirdita2022colabfold] for an MSA-backed comparison.

## 6 Discussion

A few things we did **not** solve, because they deserve to be named.

**Wet-lab feedback.** OpenDNA can quote a DNA synthesis vendor and file an
order, but there is no closed loop from assay data back into the
provenance graph. That is phase 20.

**Fine-tuning.** We orchestrate models; we do not train them. A project
that wants to fine-tune ESM-2 [@rives2021esm1b] on its own data needs to
drop to the underlying libraries. Plumbing fine-tuning runs into the
priority queue and provenance DAG is future work.

**Peer-reviewed benchmarking.** Our CASP15 suite is curated and auditable
but it is not a CASP submission. We would welcome a community
collaboration where OpenDNA submits to the next CASP round.

**Biosecurity screening of designs.** v0.5.0 exposes hooks into external
screening services but does not run its own screen; we consider this the
correct default (screening should live with the synthesis vendor) but a
user who runs OpenDNA offline should not treat absence of a warning as
evidence of safety.

**Non-protein chemistry.** Our ligand pipeline uses DiffDock
[@corso2023diffdock], xTB [@bannwarth2019gfn2xtb], and ANI-2x
[@devereux2020ani2x], and it works, but small-molecule design is a
shallower area for us than protein design.

## 7 Future Work

Immediate priorities for v0.5.x and v0.6.0:

- Closing the wet-lab loop: ingest IDT / Twist / Azenta order status
  webhooks into the provenance DAG so the same node that proposed a
  mutation records the plasmid that shipped.
- Plumbing CASP and CAMEO submissions directly out of the benchmark suite.
- First-class support for AlphaFold 3 [@abramson2024alphafold3] (pending a
  licence-compatible release path) alongside Boltz [@wohlwend2024boltz].
- Shared provenance graphs across a research group via the Yjs
  [@yjs2020] collaboration layer.
- A built-in **dry-run** mode that attaches compute budgets to workflow
  nodes before executing them.

## 8 Acknowledgments

OpenDNA is built on the shoulders of the open-source computational-biology
community. We thank the authors of every cited work; this project is
plumbing, and none of it is interesting without their models. We thank
early external users who broke the provenance graph in enough different
ways that we finally learned to write it properly.

## 9 Data Availability

All source code, benchmarks, and documentation referenced in this paper
are released under the project's open-source licence at
`https://github.com/dyber-pqc/OpenDNA`. The curated CASP15 benchmark test
set is distributed with the source at
`python/opendna/benchmarks/casp.py`. A reference SBOM is generated per
release and attached to GitHub releases. No proprietary data sets are used
at any point in the system; all references resolve to public PDB
[@berman2000pdb], UniProt, or AlphaFold DB [@varadi2022alphafolddb]
records.

## References

See `references.bib` for the full machine-readable bibliography.
