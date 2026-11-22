# OpenDNA Benchmark Methodology

OpenDNA's self-benchmarking framework compares its predictions against reference structures from AlphaFold DB. This document describes the methodology so anyone can reproduce the results.

## What we measure

For each protein in the reference set:

1. **Sequence is fetched from UniProt** using the standard accession.
2. **Reference structure is fetched from AlphaFold DB** (the AlphaFold2 prediction with multiple sequence alignment).
3. **OpenDNA folds the sequence** using ESMFold (single-sequence, no MSA).
4. **RMSD is computed** between the OpenDNA prediction and the AlphaFold reference using Kabsch superposition over Cα atoms.
5. **TM-score proxy** is computed as `1 / (1 + (RMSD/d0)^2)` where `d0 = 1.24 * (n - 15)^(1/3) - 1.8` is the standard length-dependent normalization.
6. **Mean pLDDT confidence** is reported.
7. **Inference time** is measured.

## Reference set (default)

The default benchmark set contains five well-characterized monomeric proteins:

| Protein | UniProt ID | Length | Class |
|---|---|---|---|
| Ubiquitin | P0CG48 | 76 (single chain) | Beta-grasp |
| Insulin | P01308 | 110 | Mostly helix |
| Lysozyme C | P00698 | 147 | All alpha enzyme |
| Myoglobin | P02185 | 154 | All alpha |
| GFP | P42212 | 238 | Beta barrel |

These cover the main fold classes and the small-to-medium size range that ESMFold handles well on consumer hardware.

## How to run

### Via the CLI

```bash
opendna benchmark run --output benchmark_results.json
```

### Via the API

```bash
curl -X POST http://localhost:8765/v1/benchmark
# Returns a job_id; poll /v1/jobs/{job_id} for results
```

### Via the Python SDK

```python
from opendna.benchmarks import run_benchmark_suite

suite = run_benchmark_suite(output_path="results.json")
print(f"Mean RMSD: {suite.mean_rmsd:.2f} Å")
print(f"Mean TM-score proxy: {suite.mean_tm_score:.3f}")
```

### Customizing the reference set

```python
from opendna.benchmarks import run_benchmark_suite

custom_proteins = [
    {"uniprot_id": "P00533", "name": "EGFR", "expected_length": 1210},
    {"uniprot_id": "P04637", "name": "p53", "expected_length": 393},
]

suite = run_benchmark_suite(proteins=custom_proteins)
```

## Interpretation

- **RMSD < 2 Å**: Effectively identical to the AlphaFold reference (very accurate)
- **RMSD 2–4 Å**: Same overall fold, minor local differences (good)
- **RMSD 4–6 Å**: Similar topology but with significant variation (acceptable for many uses)
- **RMSD > 6 Å**: Different fold (failure case for novel sequences)

- **TM-score proxy > 0.8**: Same fold (high confidence)
- **TM-score proxy 0.5–0.8**: Similar fold
- **TM-score proxy < 0.5**: Different fold

ESMFold's published performance is ~5–10% behind AlphaFold2 with MSA on average across CASP14, so we expect RMSDs in the 1.5–4.5 Å range for proteins similar to the training set.

## Limitations

- **AlphaFold DB is the reference, not the ground truth.** Real ground truth would be experimental crystal/cryo-EM structures from RCSB PDB.
- **Single-sequence ESMFold has known weaknesses on novel folds and orphans.** It cannot use evolutionary information.
- **Inference time varies by hardware.** GPU is ~50× faster than CPU. Results from `benchmark run` always include hardware info in the output JSON.
- **Larger proteins (>500 residues) are excluded by default** because they exceed CPU memory budgets on a typical laptop.

## Why we benchmark against AlphaFold DB instead of CASP

1. **Reproducibility**: AlphaFold DB has 200 million pre-computed predictions. Anyone can re-run our benchmarks immediately.
2. **No expensive infrastructure**: CASP comparisons require either submitting to the assessors or running AlphaFold2 yourself, which takes hours per protein and needs the full forward pass with MSA.
3. **Continuous monitoring**: We can run our benchmarks on every commit in CI.

That said, **comparing OpenDNA against published CASP15/CAMEO results is on our v0.5 roadmap.** It requires more infrastructure and we want to do it correctly.

## Citing the benchmark

If you use OpenDNA's benchmark framework in published work, please cite:

```
Kleckner, Z. (2026). OpenDNA: A Free, Open-Source Protein Engineering Platform
for Consumer Hardware. https://github.com/dyber-pqc/OpenDNA
```

And the upstream methods:

- ESMFold: Lin, Z. et al. (2023). *Science* 379, 1123–1130.
- AlphaFold DB: Varadi, M. et al. (2022). *NAR* 50, D439–D444.
