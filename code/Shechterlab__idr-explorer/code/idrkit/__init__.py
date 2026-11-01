"""idrkit: the sequence definitions shared by the Shechter lab's IDR work (glu_proteomics, MicroDomains, PTM_Pipeline
and the IDR Explorer), so each is defined once."""
from .charge import (CHARGE_TYPES, FCR_CUT, HIS, NET_CUT, basic_patches, charge_type, charge_vector, delta, delta_z, fcr,
                     ncpr, runs, scd, scd_z)
from .ptm import MODIFIED_CHARGE, charge_shift, modified_charge
from .regions import idr_intervals, merge_intervals
from .tracts import (acidic_runs, acidic_stretches, ck2_motif_count, ck2_motif_sites, e_tract_count, longest_acidic_run,
                     longest_run, residue_runs, rg_clusters, rg_motif_count)

__version__ = "0.1.0"
__all__ = ["CHARGE_TYPES", "FCR_CUT", "HIS", "NET_CUT", "charge_type", "MODIFIED_CHARGE", "charge_shift", "modified_charge", "acidic_runs", "acidic_stretches", "basic_patches", "charge_vector", "ck2_motif_count",
           "ck2_motif_sites", "delta", "delta_z", "e_tract_count", "fcr", "idr_intervals", "longest_acidic_run", "longest_run",
           "merge_intervals", "ncpr", "residue_runs", "rg_clusters", "rg_motif_count", "runs", "scd", "scd_z"]
