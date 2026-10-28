# Checks on the Arm 1 replication list

Arm 1 of the campaign lists the twelve eligible systems with the largest observed raw χ². Two questions about
that list are answered here, as code with tracked outputs; the numbers are not on the project page.

1. **Precision.** Does a large raw χ² reflect the precision of the deposited data? Reduced χ² weights each
   point by its reported error, so profiles with small errors give large χ² for the same deviation.
2. **Relatedness.** Are the twelve systems related to one another (organism, sequence, superfamily)?

| Script | Does |
|---|---|
| `arm1_checks.py` | the experimental relative error of every entry (median σ/I over the lower half of the measured q range); over the 399 clean BioEmu-1 profiles, its Spearman correlation with log raw χ² (and with NRMSD), the median raw χ² by SASBDB curve type and by source type (Kruskal–Wallis; Mann–Whitney for in-house against synchrotron); for the twelve listed systems against the eligible pool, the medians of relative error, NRMSD and raw χ² with each system's percentile in the pool and a Mann–Whitney test against the rest of the pool; and the organisms, UniProt accessions, CATH-Gene3D superfamilies, near-duplicate clusters, sequence groups and largest pairwise 5-mer containment of the twelve |

Outputs, tracked in `results/arm1/`:

| File | Contents |
|---|---|
| `relative_error.csv` | per entry (all 439): relative error and the number of points it is taken over |
| `precision.json` | the correlation, the medians by curve and source type, and the twelve against the pool |
| `relatedness.csv` | one row per listed system: organism, molecule, UniProt accession, organism group, superfamilies, near-duplicate cluster, sequence group |
| `relatedness.json` | counts of species, accessions and sequence groups, the largest pairwise containment, superfamilies shared by two or more systems |

The SASBDB metadata (curve type, source type, instrument, sample conditions) come from
`results/coverage/sasbdb_metadata.csv`, written by `analysis/coverage/18_sasbdb_metadata.py` from the
entry summaries that `analysis/coverage/01_fetch_sasbdb.py` caches. The eligible pool is
`bsc.report.arm1_eligible`, the rule the campaign uses; the script stops if the twelve it derives differ from
the listed ones in `results/campaign_arms.json`.

## Rerun

```bash
COVERAGE_PY=/path/to/venv-cov/bin/python   # as for analysis/coverage
(cd analysis/coverage && BSC_ROOT=<where bsc fetch put data/> COVERAGE_WORK=/path/with/2GB/free $COVERAGE_PY 18_sasbdb_metadata.py)
BSC_ROOT=<where bsc fetch put data/> uv run python analysis/arm1/arm1_checks.py
```

`arm1_checks.py` runs in the main environment; it reads the measured curves from `BSC_ROOT/data` and writes
to the repository's `results/arm1/`.
