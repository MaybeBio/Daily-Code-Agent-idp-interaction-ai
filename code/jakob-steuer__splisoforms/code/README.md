# SPLISOFORMS

![Version](https://img.shields.io/badge/version-1.4.0-blue.svg)
![Code license: MIT](https://img.shields.io/badge/code-MIT-blue.svg)
![Data license: CC BY-NC 4.0](https://img.shields.io/badge/data-CC%20BY--NC%204.0-lightgrey.svg)
[![Website](https://img.shields.io/badge/web-splisoforms.org-2b6cb0.svg)](https://splisoforms.org)

**SPLISOFORMS** is a structure-resolved knowledge base of alternative-splicing isoforms. It
integrates long-read cancer transcriptomes with **AlphaFold 3** predictions to map the
structural and functional consequences of alternative splicing — comparing every novel
isoform to its canonical counterpart, residue by residue. It closes the gap between
sequence-level transcriptomic discovery and protein function, where the largest structure
databases model only one (canonical) protein per gene.

Live resource: **https://splisoforms.org** · public REST API under `/api/v1` · no login required.

> Accompanies: Steuer J & Kahraman A, *SPLISOFORMS: a Structure-Resolved Knowledge Base of
> Alternative Splicing Isoforms.* bioRxiv 2026, [doi:10.64898/2026.08.06.743218](https://doi.org/10.64898/2026.08.06.743218).
> See [`CITATION.cff`](https://github.com/jakob-steuer/splisoforms/blob/HEAD/CITATION.cff).

---

## What's inside

Current release (October 2026): **22,075 genes · 131,554 isoforms · 124,734 AlphaFold 3 structures.**

| Layer | Coverage |
|---|---|
| Structures | 124,734 (94.8% of isoforms) |
| Pfam domains | 216,364 clan-de-duplicated hits on 104,135 isoforms (plus canonical-domain loss records) |
| TED (structure-based) domains | 42,599 consensus domains on 21,154 isoforms; 193,290 Chainsaw segments |
| Intrinsic disorder (AIUPred + AIUPred-binding) | 130,501 (99.2%) |
| NMD predictions | 131,554 (all isoforms; 9,262 predicted targets) |
| PTM sites (re-scored on structure) | 2,087,743 sites, 107,603 isoforms |
| Neoantigens | 7,930,448 peptide–HLA (1,941,536 peptides, 66,747 isoforms), all three tools on every row |
| STRING partners | 176,930 gene–partner pairs |
| Homo- / heterodimer interfaces (AFDB) | 12,006 / 79,628 |
| Drug complexes | 20,091 (3,725 drugs) |
| Expression | 86,331 cell-type (HPA) · 39,768 per-sample (long-read) |

---

## Datasets

| Source | Platform | Isoforms | Structures |
|---|---|---:|---:|
| GENCODE v47 (reference) | — | 75,895 | 69,075 |
| ccRCC (Karakulak 2025) | PacBio MAS-seq | 32,511 | 32,511 |
| Breast cancer (Head 2026) | Nanopore ESPRESSO | 23,148 | 23,148 |

Counts are unique protein sequences after cross-source de-duplication on `protein_md5`. All
isoforms are anchored to GRCh38 and reconciled to GENCODE v47 / Ensembl 114, with MANE Select
transcripts flagged as the per-gene ground truth.

---

## Analysis layers

Every isoform is compared against its gene's canonical protein:

- **Structure & structural impact** — AF3 model, per-residue pLDDT/PAE, DSSP secondary
  structure, and a residue-level canonical-vs-isoform diff; a composite discovery index
  ranks the most consequential variants.
- **Domains** — Pfam (HMMER, clan-de-duplicated) + TED/Chainsaw, classified as present /
  truncated / lost / gained.
- **Intrinsic disorder** — AIUPred + AIUPred-binding (disorder + binding regions).
- **Nonsense-mediated decay** — the 50-nt rule on annotated stop codons, GENCODE NMD tags
  honoured; transcripts without an annotated stop are reported as unknown.
- **PTMs** — curated modification sites re-scored on structure (intact / exposed / buried /
  lost, from the change in solvent accessibility); sites that do not match the canonical
  sequence are dropped.
- **Neoantigens** — junction-spanning 8–11-mers scored by an MHCflurry + NetMHCpan + BigMHC
  ensemble across 35 HLA class-I alleles, ranked into P1–P4 priority tiers, and flagged by
  proteome uniqueness (absent from every reference protein / present in another isoform of
  the same gene / present in another gene).
- **Interactions** — AFDB homo-/heterodimer interface disruption + STRING partners.
- **Drug pockets** — AF3 co-folded drug–target complexes, each pocket mapped onto sibling
  isoforms as preserved / disrupted / lost.
- **Expression** — Human Protein Atlas single-cell + long-read per-sample.

---

## Tech stack

- **Backend** — Python, FastAPI, SQLAlchemy 2 + Alembic, Pydantic, PostgreSQL.
- **Frontend** — Next.js 15.5, React 19, TypeScript, Mol\* (3D viewer), Plotly.
- **Structure prediction** — AlphaFold 3, run on the EuroHPC LEONARDO cluster.

---

## Getting started (development)

Requires Docker + Docker Compose.

```bash
docker compose up -d                              # PostgreSQL + backend + frontend
docker compose exec backend alembic upgrade head  # apply the schema
```

Backend → http://localhost:8000 · frontend → http://localhost:3000.

Predicted structure files (`.cif`) are large and **not tracked in git** (see `.gitignore`);
provision them into `backend/static/data/` to serve structures locally.

---

## REST API

Open and read-only — **no key required** for standard use (anonymous requests are
rate-limited per IP; an API key raises the quota for bulk access). Base URL
`https://splisoforms.org/api/v1`; interactive docs at https://splisoforms.org/api-docs.

```bash
curl https://splisoforms.org/api/v1/genes/TP53
curl https://splisoforms.org/api/v1/isoforms/TP53-208            # transcript name or id
curl -L https://splisoforms.org/api/v1/isoforms/TP53-208/structure -o model.cif
```

---

## Repository layout

```
backend/    FastAPI app (api/, analysis/, models/, services/), Alembic migrations, pipeline scripts
frontend/   Next.js app (gene / isoform / explore / documentation / api-docs)
scripts/    database construction & backfill pipeline
database/   base schema
docs/       supplementary documentation
docker-compose.yml
```

---

## Citation

If you use SPLISOFORMS, please cite the preprint (see [`CITATION.cff`](https://github.com/jakob-steuer/splisoforms/blob/HEAD/CITATION.cff)):

> Steuer J, Kahraman A. *SPLISOFORMS: a Structure-Resolved Knowledge Base of Alternative
> Splicing Isoforms.* bioRxiv 2026. https://doi.org/10.64898/2026.08.06.743218

The journal reference will be added on publication.

---

## License

- **Code** — the source in this repository is released under the **MIT License** (see [LICENSE](https://github.com/jakob-steuer/splisoforms/blob/HEAD/LICENSE)).
- **Database content & annotations** — the SPLISOFORMS data and derived annotations are provided under **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)**.
- **Predicted structures** — all AlphaFold 3-predicted structures are additionally subject to the **AlphaFold 3 Model Parameters Terms of Use**; any ongoing use or adaptation of the structural models remains governed by those terms.

---

## Funding & acknowledgements

Developed by the **Data Science in Life Sciences Group** (Kahraman Lab), Institute for
Chemistry and Bioanalytics, **FHNW School of Life Sciences**, and the **Swiss Institute of
Bioinformatics**. Supported by the **Swiss National Science Foundation** (SNSF grant 222878),
with structure prediction on the EuroHPC **LEONARDO** supercomputer (CINECA).

The database-construction and analysis pipeline was developed with the assistance of Claude
Code (Anthropic); all code and analyses were reviewed and verified by the authors.
