# UniPPI-OWL v1.8 — Snakemake High-Throughput Structure-to-Affinity Workflow
# Usage: snakemake --cores 1

INPUT_DIR      = config.get("input_dir", "demo_structures")
OUTPUT_DIR     = config.get("output_dir", "demo_results")
TEMP_K         = config.get("temp", 298.15)
LAMBDA_PAE     = config.get("lambda_pae", 8.0)
IONIC_STRENGTH = config.get("ionic_strength", 150.0)
PH_VAL         = config.get("ph", 7.2)
PTM_STATE      = config.get("ptm_state", "native")
MUTATIONS      = config.get("mutate", "A:LEU5PRO,A:PRO5LEU,A:ASP2ALA,A:TYR5ALA")

rule all:
    input:
        f"{OUTPUT_DIR}/master_interaction_ranking.csv",
        f"{OUTPUT_DIR}/unippi_full_report.json",
        f"{OUTPUT_DIR}/REGRESSION_VERIFICATION_REPORT.tsv"

rule unippi_batch_screen:
    input:
        folder = INPUT_DIR
    output:
        csv  = f"{OUTPUT_DIR}/master_interaction_ranking.csv",
        json = f"{OUTPUT_DIR}/unippi_full_report.json"
    shell:
        """
        python3 unippi_cli.py \
            -i {input.folder} \
            -o {OUTPUT_DIR} \
            --temp {TEMP_K} \
            --lambda-pae {LAMBDA_PAE} \
            --ionic-strength {IONIC_STRENGTH} \
            --ph {PH_VAL} \
            --ptm-state {PTM_STATE} \
            --mutate "{MUTATIONS}"
        """

rule unippi_regression_verify:
    input:
        csv = f"{OUTPUT_DIR}/master_interaction_ranking.csv"
    output:
        tsv = f"{OUTPUT_DIR}/REGRESSION_VERIFICATION_REPORT.tsv"
    shell:
        "python3 test_unippi_regression.py"
