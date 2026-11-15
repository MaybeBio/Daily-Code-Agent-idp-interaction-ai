# BindCraft code — the ablation versions

Both codebases used in this thesis, source only.

```
idp_bindcraft/        BindCraft1 + the IDP modifications  
bindcraft2_v1.0.4/          BindCraft2 baseline                
ablation_arm_configs/       what actually defines each arm
```
```
python bindcraft.py --settings <target>.json --advanced <arm>.json --filters <filters>.json
```

`ablation_arm_configs/`:

| Folder | Contents |
|---|---|
| `panelA_advanced/` | 13 arm configs — losses (`no_losses`, `no_pae_loss`, `no_plddt_loss`), template (`with_template`), objectives (`actif_w1_lam0`, `becraft_w1_lam0`, `energy_*`), filter arms (`ipsae`, `ipsae_replace`, `vanilla_filters`), and the `full` baseline |
| `panelB_advanced/` | 5 arm configs — crop/binder length and validation model (`multimer`, `monomer_validation_only`, `actif_w1_lam1`) |
| `filters/` | 5 acceptance-filter JSONs: `idp_final_filters_iface_ipae` (default), `default_filters` (vanilla BindCraft), `filters_ipsae`, `filters_ipsae_replace`, `filters_idp` (production) |
| `cmp50_advanced/` | the budget-matched BC1-vs-BC2 comparison config |

Each arm changes exactly one factor against its baseline — verified by diffing the JSONs.
**One documented exception:** the Panel B length/crop/validation arms use the standard
i_pTM objective (`multimer.json`) while the Panel B reference `actif_w1_lam1` uses the actif
objective, so comparisons against that reference carry an objective change as well as the
intended one. See the ablation notes in the reproducibility bundle.

## Which engine ran what

- **Panel A, Panel B, the 500-crop production, the conformation runs** → `idp_bindcraft`
- **BC2 arms** (cmp50 V1/V2/V3, the V7/V8 IDP-head ports, the FASTA-IDR mode) → `bindcraft2_v1.0.4`

## Environment

BC1 runs inside the `jax_env` container (`conda activate BindCraft`); BC2 has its own
`.venv` (excluded here — recreate with `pyproject.toml` / `install.sh`).
