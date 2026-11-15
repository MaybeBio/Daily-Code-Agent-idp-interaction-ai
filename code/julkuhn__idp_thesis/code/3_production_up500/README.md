# idp_bindcraft — the `actif` production configuration

The configuration used for the **500-crop production run** (`up500`), and how to run it.

BindCraft is driven by **three JSON files**; the code is identical across every arm of this
thesis, only these files change. That is what an "ablation version" is here.

```
python -u bindcraft.py --settings <target>.json --filters <filters>.json --advanced <advanced>.json
```

## Files

```
advanced_actif_uncapped.json     the actif production config
advanced_full_uncapped.json      the IDP baseline it was compared against
filters_idp.json                 acceptance filters, identical for both arms
EXAMPLE_target_settings.json     one per-target file, as a template
crops_500.txt, crops_100_full.txt, jobs_actif.txt, jobs_full.txt
run_up500_actif.slurm, run_up500_actif_requeue.slurm, run_up500_full.slurm
manifest.csv                     one row per run
binder_structural_features.csv   per-design structural measurements
output/                          per-crop design and trajectory statistics
```

`jobs_actif.txt` is **tab-separated**: `<settings>\t<advanced>\t<filters>`, one line per crop,
so every run carries its own arm config. The SLURM script reads the line matching its task index.

## What makes it "actif"

`advanced_actif_uncapped.json` differs from the baseline in exactly these keys — the design
objective, nothing else:

| Key | actif | baseline (`full`) |
|---|---|---|
| `use_i_ptm_loss` | `false` | `true` |
| `weights_iptm` | `0.0` | `0.15` |
| `use_actifptm_loss` | `true` | — |
| `use_actif_i_ptm_energy_loss` | `true` | — |
| `weights_actif_iptmenergy_loss` | `1.0` | — |
| `weights_actifptm` | `0.0` | — |
| `iptmenergy_lambda` | `1.0` | — |
| `actifptm_contact_dist` | `8.0` | — |

The standard i_pTM interface loss is switched **off** and replaced by the contact-weighted
actif-pTM energy loss. Everything else — multimer design + validation, the IDP interface
losses, 4-stage design, recycles, filters — is unchanged. `max_trajectories` is **false**
(uncapped): runs end on walltime via the resubmit chain, never on a trajectory cap.

## Running it

Single crop (interactive, one GPU):

```bash
source /miniforge3/etc/profile.d/conda.sh
conda activate BindCraft                       # inside the jax_env container
cd <idp_bindcraft>                       # ../bindcraft_code/idp_bindcraft
python -u bindcraft.py \
  --settings  EXAMPLE_target_settings.json \
  --advanced  advanced_actif_uncapped.json \
  --filters   filters_idp.json
```

To target your own crop, copy `EXAMPLE_target_settings.json` and set `starting_pdb`
(the cropped target PDB), `design_path` (output dir), `binder_name`, and
`lengths` (`[130,130]` was used throughout).

At scale (`run_up500_actif.slurm`): 1 node = **16 crops on 4 GPUs** (4 per GPU,
`XLA_PYTHON_CLIENT_MEM_FRACTION=0.22`), array `1-32` to cover 500 crops. The script traps
`USR1` 120 s before walltime and resubmits itself with `--dependency=afterany`, up to
`MAX_LEGS=4` legs, so a campaign survives the 12 h limit. Edit `ROOT`, `BINDCRAFT_DIR`,
`--account` and the log paths before reuse.

## What it produced

500 crops, uncapped: **38,354 accepted designs from 117,610 trajectories**. The paired
baseline (100 crops run both ways) gave **2,901 accepted from 20,053 trajectories** — the
actif objective is the single change responsible for the difference.

Per-crop numbers: `output/<crop>__actif/final_design_stats.csv`.

## Outputs per run

Inside each `design_path`: `final_design_stats.csv` (one row per **accepted** design, all
confidence + PyRosetta metrics), `trajectory_stats.csv` (one row per trajectory — the
denominator for acceptance rate), and `Accepted/` / `Rejected/` PDB complexes (binder = chain B).
