<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/banner-dark.png">
  <img alt="Abubakar Saeed — interpretable physicochemical graph encoding for protein and RNA" src="https://github.com/Abubakar-Saeed/Abubakar-Saeed/raw/HEAD/assets/banner-light.png">
</picture>

<p align="center">
  <a href="https://doi.org/10.1093/bioadv/vbag181"><img alt="Bioinformatics Advances" src="https://img.shields.io/badge/Bioinformatics_Advances-2026-0f766e?style=flat-square"></a>
  <a href="https://doi.org/10.1186/s12859-026-06619-5"><img alt="BMC Bioinformatics" src="https://img.shields.io/badge/BMC_Bioinformatics-2026-0f766e?style=flat-square"></a>
  <img alt="Under review" src="https://img.shields.io/badge/under_review-2_manuscripts-64748b?style=flat-square">
  <a href="https://orcid.org/0009-0008-9936-4707"><img alt="ORCID" src="https://img.shields.io/badge/ORCID-0009--0008--9936--4707-a6ce39?style=flat-square&logo=orcid&logoColor=white"></a>
  <a href="https://scholar.google.com/citations?hl=en&user=Kp1VU0QAAAAJ"><img alt="Google Scholar" src="https://img.shields.io/badge/Google_Scholar-profile-4285f4?style=flat-square&logo=googlescholar&logoColor=white"></a>
  <a href="https://linkedin.com/in/abubakarsaeed915"><img alt="LinkedIn" src="https://img.shields.io/badge/LinkedIn-abubakarsaeed915-0a66c2?style=flat-square&logo=linkedin&logoColor=white"></a>
  <a href="mailto:abubakarsaeedgcuf@gmail.com"><img alt="Email" src="https://img.shields.io/badge/Email-contact-334155?style=flat-square&logo=gmail&logoColor=white"></a>
</p>

Hi — I'm Abubakar. I work on **interpretable machine learning for molecular
biology**, from Chiniot, Pakistan.

I came to this from the biology side first: a pre-medical background before
computer science, which is probably why my work tends to start from a
biophysical mechanism and end at a model, rather than the other way round.
Since 2025 I have been doing that independently — four papers, two of them
published — all circling one question: *how much of a molecule's behaviour can
you predict from its sequence alone, if you encode the physics carefully
enough?*

The answer so far is: more than I expected, and the encoding stays readable
while doing it. A residue or nucleotide is described by the interactions it
*could* form given its sequence neighbourhood — a deterministic graph, no
predicted structure anywhere in the pipeline, every feature traceable back to
the biophysical rule that produced it.

I am currently applying for graduate study in computational biology.

<sub>BS Computer Science (2026) · GC University Faisalabad · CGPA 3.81/4.0 · pre-medical background · independent researcher since 2025</sub>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/pipeline-dark.gif">
    <img alt="Sequence to constraint graph to per-residue feature profile" src="https://github.com/Abubakar-Saeed/Abubakar-Saeed/raw/HEAD/assets/pipeline-light.gif" width="100%">
  </picture>
</p>

<p align="center"><sub>Sequence in, profile out. Nodes are residues coloured by chemistry; edges appear where two residues are chemically compatible and close enough along the chain.</sub></p>

<p align="center">
  <a href="https://abubakar-saeed.github.io"><b>Try it in your browser →</b></a><br>
  <sub>Paste any protein sequence and build the graph live. Hover a residue to see its features and which rules fired.</sub>
</p>

---

## Published

<table>
<tr>
<td width="50%" valign="top">

**BioGraphX**
*Bridging the Sequence–Structure Gap via Physicochemical Graph Encoding for Interpretable Subcellular Localization Prediction*

Saeed A., Abbas W.
**Bioinformatics Advances** (Oxford University Press), 2026

State of the art on DeepLoc 2.0 from sequence alone.

[`10.1093/bioadv/vbag181`](https://doi.org/10.1093/bioadv/vbag181)

</td>
<td width="50%" valign="top">

**BioGraphX-RNA**
*A Universal Physicochemical Graph Encoding for Interpretable RNA Subcellular Localization Prediction*

Saeed A., Abbas W.
**BMC Bioinformatics**, 2026

Graph-only miRNA AUROC **0.9396** · 2.05M parameters.

[`10.1186/s12859-026-06619-5`](https://doi.org/10.1186/s12859-026-06619-5)

</td>
</tr>
</table>

## Under review

<table>
<tr>
<td width="50%" valign="top">

**BioGraphX-IDP** — intrinsic disorder
*Bioinformatics* (Oxford)

61 biophysical features gated per residue against frozen ESM-2.
Blind CAID-3: AUC **0.899** (PDB), **0.886** (NOX).

</td>
<td width="50%" valign="top">

**BioGraphX-Sol** — protein solubility
*BMC Bioinformatics*

Physics-only R² **0.499** exceeds FGNNSol (AlphaFold3-based, 0.487).
Ensemble matches ProtSATT at 0.545 — 1.44M parameters against 6M.

</td>
</tr>
</table>

---

## One encoding, four domains

| Domain | What it encodes | Headline |
|:--|:--|:--|
| **Protein localisation** | Residue interaction graph, 158 features | State of the art on DeepLoc 2.0 |
| **RNA localisation** | Base pairing, stacking, stem-loops, pseudoknots | Graph-only miRNA AUROC 0.9396 |
| **Intrinsic disorder** | Per-residue gating against ESM-2 | Blind CAID-3: 0.899 PDB / 0.886 NOX |
| **Solubility** | 200-dimensional physicochemical encoding | Beats a structure-based model, without structures |

<details>
<summary><b>The finding that ties them together</b></summary>

<br>

Each paper measures how much of the signal is biophysics and how much is
evolutionary information, by sweeping the blend weight between the two streams.
The answer turns out to depend on **how the property is defined**, not on the
molecule.

Solubility — physically measured — settles at 59% physics. Intrinsic disorder
settles at 49% on the structurally observed CAID-3 PDB track but 13% on the
functionally annotated NOX track: a 3.8-fold difference between two tracks of
the *same* task. Properties that are measured physically want physics;
properties that are annotated by function want evolution.

</details>

<details>
<summary><b>What the encoding is, and is not</b></summary>

<br>

Interaction capacity, not geometry. Edges are placed by sequence separation and
residue chemistry, so the representation is **complementary to a contact map
rather than an approximation of one**. That is why the method needs no
predicted structure, no multiple sequence alignment, and no solved structure
for the target — and why every feature stays traceable to the rule that
produced it.

</details>

---

## Interests

`interpretable AI for biomedicine` · `protein & RNA structure–function` ·
`graph representations of sequence` · `parameter-efficient & Green AI` ·
`generative models for RNA and protein design` · `bioinformatics tooling`

## Tools

![Python](https://img.shields.io/badge/Python-1e293b?style=flat-square&logo=python&logoColor=white)
![PyTorch](https://img.shields.io/badge/PyTorch-1e293b?style=flat-square&logo=pytorch&logoColor=white)
![XGBoost](https://img.shields.io/badge/XGBoost-1e293b?style=flat-square)
![SHAP](https://img.shields.io/badge/SHAP-1e293b?style=flat-square)
![ESM&#8209;2](https://img.shields.io/badge/ESM--2-1e293b?style=flat-square)
![RiNALMo](https://img.shields.io/badge/RiNALMo-1e293b?style=flat-square)
![BioPython](https://img.shields.io/badge/BioPython-1e293b?style=flat-square)
![LaTeX](https://img.shields.io/badge/LaTeX-1e293b?style=flat-square&logo=latex&logoColor=white)

---

> [!NOTE]
> **Open to collaboration.** I am always glad to discuss joint work, reviews,
> or applying this encoding to a new molecule or property — enzyme engineering,
> aptamer design, protein property prediction, or something I have not thought
> of. Reach me at
> **[abubakarsaeedgcuf@gmail.com](mailto:abubakarsaeedgcuf@gmail.com)**.

<sub>Chiniot, Punjab, Pakistan · [ORCID](https://orcid.org/0009-0008-9936-4707) · [Google Scholar](https://scholar.google.com/citations?hl=en&user=Kp1VU0QAAAAJ) · [LinkedIn](https://linkedin.com/in/abubakarsaeed915)</sub>
