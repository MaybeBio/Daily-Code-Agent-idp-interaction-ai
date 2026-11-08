

Disobind accepts protein pair inputs in two formats—**CSV** and **FASTA**—each designed for a different use case. The CSV format is the preferred and more flexible option, supporting batch prediction of multiple protein pairs by referencing UniProt accessions. The FASTA format enables prediction for custom sequences that lack a UniProt entry, but is limited to a single pair per file. Both formats also support the optional **Disobind+AF2** mode, where you supply AlphaFold2 structural predictions alongside the sequence information for refined results.

## Input Requirements

Before preparing your input files, keep these constraints in mind:

| Requirement | Detail |
|:---|:---|
| **Binary complexes only** | Disobind operates on pairs (AB). For non-binary complexes (ABC), decompose into binary pairs (AB, BC, AC) and run each separately. |
| **Protein 1 must be an IDR** | The first protein in the pair is assumed to be an intrinsically disordered region. Protein 2 may or may not be disordered. |
| **Interacting pair assumed** | Disobind predicts *where* the proteins interact (contact maps and interface residues), not *whether* they interact. |
| **Residue range validation** | The specified `start`–`end` range must not exceed the length of the corresponding protein sequence. |

Sources: [run_disobind.py](/run_disobind.py#L4-L14), [run_disobind.py](/run_disobind.py#L435-L453), [README.md](/README.md#L43-L46)

## CSV Format (Preferred)

The CSV format is the recommended way to provide input. Each row represents one protein fragment pair for which Disobind predictions are required. There are **two variants** depending on whether you want Disobind-only or Disobind+AF2 predictions.

### Disobind-Only (6 fields)

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2
```

| Field | Description | Example |
|:---|:---|:---|
| `UniProt_ID1` | UniProt accession for the IDR (Protein 1) | `P04273` |
| `start1` | Start residue position of Protein 1 fragment (UniProt numbering) | `95` |
| `end1` | End residue position of Protein 1 fragment (UniProt numbering) | `193` |
| `UniProt_ID2` | UniProt accession for the binding partner (Protein 2) | `P04273` |
| `start2` | Start residue position of Protein 2 fragment | `95` |
| `end2` | End residue position of Protein 2 fragment | `192` |

When using this format, Disobind automatically downloads the full UniProt sequences for both accessions via the UniProt REST API, then extracts the specified residue fragments for prediction.

### Disobind+AF2 (12 fields)

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2,AF2_struct_file,AF2_pae_file,chain1,chain2,offset1,offset2
```

| Field | Description | Example |
|:---|:---|:---|
| `AF2_struct_file` | Path to the AF2 predicted structure file (`.pdb` or `.cif`) | `./example/unrelaxed_model_4_multimer_v3_pred_4.pdb` |
| `AF2_pae_file` | Path to the AF2 PAE data file (`.json`) | `./example/pae_model_4_multimer_v3_pred_4.json` |
| `chain1` | Chain ID in the AF2 structure corresponding to Protein 1 | `B` |
| `chain2` | Chain ID in the AF2 structure corresponding to Protein 2 | `C` |
| `offset1` | Residue position offset for Protein 1 (0 if AF2 structure matches UniProt range) | `0` |
| `offset2` | Residue position offset for Protein 2 (0 if AF2 structure matches UniProt range) | `0` |

> [!TIP]
> The **offset** fields account for differences between AF2 residue numbering and UniProt numbering. Set both to `0` if the AF2 structure corresponds to the full UniProt sequence or exactly to the specified fragment range.

### CSV Example

The repository includes a ready-to-run example at `example/test.csv`:

```csv
P04273,95,193,P04273,95,192
P04273,95,193,P04273,95,193,./example/unrelaxed_model_4_multimer_v3_pred_4.pdb,./example/pae_model_4_multimer_v3_pred_4.json,B,C,0,0
```

- **Row 1**: A Disobind-only prediction for the prion protein (P04273) fragment 95–193 interacting with itself at fragment 95–192.
- **Row 2**: A Disobind+AF2 prediction for the same prion protein fragment 95–193, augmented with an AF2 multimer model (chains B and C, no residue offset).

Sources: [example/test.csv](/example/test.csv#L1-L2), [run_disobind.py](/run_disobind.py#L312-L359), [README.md](/README.md#L52-L65)

## FASTA Format

The FASTA format is for proteins **without a UniProt accession**. You provide the amino acid sequences directly. Each file must contain **exactly two entries**—one for each protein in the pair.

### Disobind-Only (3 header fields)

```
>Protein_ID1,start1,end1
SEQUENCE1
>Protein_ID2,start2,end2
SEQUENCE2
```

### Disobind+AF2 (7 header fields)

```
>Protein_ID1,start1,end1,AF2_struct_file,AF2_pae_file,chain1,offset1
SEQUENCE1
>Protein_ID2,start2,end2,AF2_struct_file,AF2_pae_file,chain2,offset2
SEQUENCE2
```

| Header Field | Disobind-Only | Disobind+AF2 | Description |
|:---|:---:|:---:|:---|
| `Protein_ID` | ✓ | ✓ | An identifier for the protein (used in output labeling) |
| `start` | ✓ | ✓ | Start residue position of the fragment |
| `end` | ✓ | ✓ | End residue position of the fragment |
| `AF2_struct_file` | — | ✓ | Path to the AF2 structure file (`.pdb`/`.cif`) |
| `AF2_pae_file` | — | ✓ | Path to the AF2 PAE file (`.json`) |
| `chain` | — | ✓ | Chain ID in the AF2 structure for this protein |
| `offset` | — | ✓ | Residue position offset for this protein |

> [!TIP]
> In FASTA mode, the provided sequence is assumed to be the **full protein sequence**. Disobind uses this sequence directly instead of fetching from UniProt. The `start`/`end` values define the fragment to extract from this full sequence for prediction.

### FASTA Example

The repository includes an example at `example/test.fasta`:

```
>prion,95,193,./example/unrelaxed_model_4_multimer_v3_pred_4.pdb,./example/pae_model_4_multimer_v3_pred_4.json,B,0
MANLSYWLLALFVAMWTDVGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGTWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHNQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPMMHFGNDWEDRYYRENMNRYPNQVYYRPVDQYNNQNNFVHDCVNITIKQHTVTTTTKGENFTETDIKIMERVVEQMCTTQYQKESQAYYDGRRSSAVLFSSPPVILLISFLIFLMVG

>prion,95,193,./example/unrelaxed_model_4_multimer_v3_pred_4.pdb,./example/pae_model_4_multimer_v3_pred_4.json,C,0
MANLSYWLLALFVAMWTDVGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGTWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHNQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPMMHFGNDWEDRYYRENMNRYPNQVYYRPVDQYNNQNNFVHDCVNITIKQHTVTTTTKGENFTETDIKIMERVVEQMCTTQYQKESQAYYDGRRSSAVLFSSPPVILLISFLIFLMVG
```

This defines a Disobind+AF2 prediction for the prion protein homodimer, where both chains map to the same AF2 multimer structure (chains B and C respectively).

Sources: [example/test.fasta](/example/test.fasta#L1-L5), [run_disobind.py](/run_disobind.py#L240-L308), [README.md](/README.md#L67-L72)

## AF2 Supplementary Files

When using Disobind+AF2 mode, you must provide two additional files from your AlphaFold2 multimer prediction:

| File | Format | Content | Key Field |
|:---|:---|:---|:---|
| **Structure file** | `.pdb` or `.cif` | 3D coordinates with per-residue pLDDT scores in the B-factor column | Atom coordinates, B-factor (pLDDT) |
| **PAE file** | `.json` | Predicted Aligned Error matrix | `predicted_aligned_error` |

The PAE JSON file must be a JSON array containing an object with the key `"predicted_aligned_error"`, which stores the asymmetric PAE matrix as a 2D array. Disobind symmetrizes this internally by computing `(PAE + PAE^T) / 2`. The structure file is parsed using Biopython's `PDBParser` (for `.pdb`) or `MMCIFParser` (for `.cif`), and per-residue pLDDT scores are read from the Cα B-factor column.

### Confidence Filtering

Disobind applies three internal confidence thresholds when combining AF2 predictions:

| Parameter | Default | Meaning |
|:---|:---:|:---|
| `dist_threshold` | 8 Å | Maximum Cα–Cα distance to define a residue contact |
| `plddt_threshold` | 70 | Minimum pLDDT score to consider a residue confidently predicted |
| `pae_threshold` | 5 | Maximum PAE value to consider a residue pair confidently aligned |

Sources: [run_disobind.py](/run_disobind.py#L966-L1090), [run_disobind.py](/run_disobind.py#L70-L76)

## Choosing Between CSV and FASTA

| Feature | CSV | FASTA |
|:---|:---:|:---:|
| Multiple pairs per file | ✓ | ✗ |
| UniProt accession required | ✓ | ✗ |
| Custom sequences supported | ✗ | ✓ |
| Disobind-only mode | ✓ | ✓ |
| Disobind+AF2 mode | ✓ | ✓ |
| Auto sequence download | ✓ | ✗ |

Use **CSV** when your proteins have UniProt entries—you can batch dozens of pairs in a single file and let Disobind fetch the sequences automatically. Use **FASTA** when working with custom or synthetic sequences that are not deposited in UniProt.

Sources: [run_disobind.py](/run_disobind.py#L212-L237), [README.md](/README.md#L48-L51)

## Running the Example

Execute Disobind with the provided example input using:

```bash
python run_disobind.py -i csv -f ./example/test.csv
```

For the FASTA example:

```bash
python run_disobind.py -i fasta -f ./example/test.fasta
```

Both commands use default settings: interface residue prediction at coarse-grained resolution 1 on CPU. See [Quick Start](2-quick-start) for the full list of command-line flags and [Output Interpretation](10-output-interpretation) for understanding the results.

Sources: [run_disobind.py](/run_disobind.py#L1312-L1357), [README.md](/README.md#L76-L80)

## Format Validation

Disobind performs strict validation on your input. The following errors will be raised if the format is incorrect:

| Condition | Error Type | Message |
|:---|:---|:---|
| CSV row has neither 6 nor 12 comma-separated fields | `ValueError` | `Incorrect input format...` |
| FASTA file contains ≠ 2 entries | `ValueError` | `Expected exactly 2 FASTA entries, found N...` |
| FASTA header has neither 3 nor 7 comma-separated fields | `ValueError` | `Incorrect input header format...` |
| Specified residue range exceeds sequence length | `ValueError` | `Sequence length (L residues) is shorter than the specified residue range S-E...` |
| AF2 structure file is not `.pdb` or `.cif` | `ValueError` | `Incorrect file format.. Supported .pdb/.cif only.` |
| PAE file is not `.json` | `ValueError` | `Incorrect file format.. Supported .json only.` |
| PAE JSON missing `predicted_aligned_error` key | `ValueError` | `PAE matrix not found...` |
