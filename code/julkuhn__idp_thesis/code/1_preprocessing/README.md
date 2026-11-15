# Target Preprocessing

1. Download PDB and PAE data for each target from uniprot
2. Create domain annotation based on PAE matrix using https://github.com/tristanic/pae_to_domains.
3. Segment target PDBs into domain PDBs using domain annotation.
4. Compute per-domain structural statistics (`04_segment_and_compute_domain_stats.py`) → `domain_stats.csv`.
5. Filter domains and build the structured-domain sweep index (`05_filter_domains_and_create_sweep_index.py`,
   structured mode — the default, no `--idp` flag) ->  `filtered_domain_stats.csv` + `sweep_index*.tsv`.
6. Generate IDP crops (`04b_generate_idp_crops.py`) from the residues NOT covered by an
   accepted structured domain (step 5's output) and not membrane-embedded — dynamic-length
   (30-300 aa by default) crops with hydrophilic-capped termini, allowed to overlap each
   other but never an accepted domain → `idp_crops_stats.csv`.
7. Filter IDP crops and build the IDP sweep index (`05_filter_domains_and_create_sweep_index.py --idp`)
   → `filtered_idp_crops_stats.csv` + `idp_sweep_index*.tsv`.

```bash
python 05_filter_domains_and_create_sweep_index.py

python 04b_generate_idp_crops.py \
    --filtered-domain-stats-csv Data/Filtered/filtered_domain_stats.csv \
    --design-root /path/to/output_idp \
    --bindcraft-advanced-settings /path/to/idp_advanced_settings.json \
    --max-crops-per-segment 3   # bounds total crop/job volume; tune against a sample first (see below)

python 05_filter_domains_and_create_sweep_index.py --idp \
    --min-domain-length 30 --max-domain-length 300 --idp-max-plddt 0
```