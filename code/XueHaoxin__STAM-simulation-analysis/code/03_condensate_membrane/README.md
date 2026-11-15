# Condensate–membrane interactions

This folder contains analysis scripts about the interaction between the condensate and membrane, corresponding to Fig. 3g, Fig. 3k, and Extended Data Fig. 6a,b in the manuscript.

- **Fig. 3g:**
1. Run `condensate_membrane_contacts.py` to calculate bead–bead contacts between the condensate and membrane.
2. Run `residue_membrane_contacts.py` to group contacts by residue ID and average over time.
3. Use `resname-beads_number.xlsx` to divide each residue's contact count by its number of beads.
4. Divide by 9 to average over the nine protein chains.

Keep `get_chainID_resid.py` alongside `condensate_membrane_contacts.py`.

- **Extended Data Fig. 6a,b:** `membrane_z-height_curvature_after_interaction.py` calculates membrane height and mean curvature maps and exports the maximum mean curvature for each leaflet.
- **Fig. 3k:** `max_curvature_summary.py` summarizes curvature results across systems and replicates.

Update input file paths before running. The residue_membrane_contacts.py reads `contact_details.csv` from each immediate subfolder. The curvature summary script expects `<replicate_folder>/mean_curvature/<system_name>_mean_curvature_max.csv`.

Requirements we used: MDAnalysis(version 2.10.0), MembraneCurvature(version 1.1.2)

Example data are available in the [Releases](https://github.com/XueHaoxin/STAM-simulation-analysis/releases) section.
