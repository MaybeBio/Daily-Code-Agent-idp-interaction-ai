#!/usr/bin/env python
# coding: utf-8

# In[1]:


from datetime import datetime
start_time = datetime.now()
print(f"start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")


# In[2]:


import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda

import warnings
warnings.filterwarnings("ignore")

from lipyphilic.leaflets.assign_leaflets import AssignLeaflets, AssignCurvedLeaflets
from lipyphilic.analysis import AreaPerLipid, MembThickness
import glob
import os


# In[3]:


def run_membrane_analysis(dir):
    origin_cwd = os.getcwd()
    try:
        os.chdir(dir)
        print(f"================ Analyzing {dir} ================")
        
        u = mda.Universe("md.tpr", "md_pbcwhole.xtc")
        print(f"Total frames{len(u.trajectory)} dt={u.trajectory.dt} ps ")
        lipid_num=0
        for seg in u.segments:
            print(f"seg{seg.segid}:\t{set(seg.residues.resnames)}-{seg.residues.n_residues}")
            if ('POPC' in seg.segid) or ('POPE' in seg.segid) or ('POPA' in seg.segid) or ('CHOL' in seg.segid):
                lipid_num=lipid_num+seg.residues.n_residues
        print(f'Total number of lipids:{lipid_num}')

        # 1. Assign lipids to the upper and lower leaflets
        leaflets = AssignLeaflets(
            universe=u,
            lipid_sel="name GL1 GL2 ROH"  # atom selection covering all lipids in the bilayer
        )
        leaflets.run(start=1000,stop=None,step=2,verbose=True)

        # Save leaflet assignments with shape (n_lipids, n_frames); leaflet labels are -1 and +1
        np.save(f"../{results_directory}/{dir}_leaflets.npy", leaflets.leaflets)

        # 2. Calculate membrane thickness
        memb_thickness = MembThickness(
          universe=u,
          leaflets=leaflets.filter_leaflets("resname POPC POPE POPA"),  # exclude cholesterol from thickness calculation
          lipid_sel="resname POPC POPE POPA and name PO4",
          #n_bins=10,
          return_surface=True
        )
        memb_thickness.run(start=1000,stop=None,step=2,verbose=True)
        print(f"Thickness_final:{memb_thickness.memb_thickness.shape}")
        np.save(f"../{results_directory}/{dir}_membrane_thickness.npy", memb_thickness.memb_thickness)

        # 3. Calculate area per lipid
        areas = AreaPerLipid(
          universe=u,
          lipid_sel="name GL1 GL2 ROH",  # assuming we're using the MARTINI forcefield
          leaflets=leaflets.leaflets
        )
        areas.run(start=1000, stop=None, step=2,verbose=True)
        print(areas.areas.shape)
        # Save per-lipid areas with shape (n_lipids, n_frames)
        np.save(f"../{results_directory}/{dir}_APL.npy", areas.areas)

        # Calculate the mean APL across lipids for each frame
        frame_avg_area = np.mean(areas.areas, axis=0)
        print(f"APL_final:{frame_avg_area.shape}")
        np.save(f"../{results_directory}/{dir}_APL_avg_per_frame.npy", frame_avg_area)


        mean_thickness = float(np.mean(memb_thickness.memb_thickness))
        std_thickness = float(np.std(memb_thickness.memb_thickness))
        mean_APL = float(np.mean(frame_avg_area))
        std_APL = float(np.std(frame_avg_area))
        
        print(f"{dir} analysis complete | Mean thickness: {mean_thickness:.4f} | Mean APL: {mean_APL:.4f}")
        return mean_thickness, std_thickness, mean_APL, std_APL
        
    finally:
        os.chdir(origin_cwd)


# In[4]:


#pwd


# In[5]:


results_directory = pathlib.Path(rf"results")
results_directory.resolve().mkdir(exist_ok=True, parents=True)# Create the directory if it doesn't already exist


# In[6]:

target_dirs = [
    "noPA_ca0", "noPA_ca25", "noPA_ca50",
    "PA_ca0", "PA_ca25", "PA_ca50"
]

root_path = "."
sub_folders = glob.glob(os.path.join(root_path, "*"))
sub_folders = [d for d in sub_folders if os.path.isdir(d)]

summary_list = []

for folder in sub_folders:
    dir_name = os.path.basename(folder)
    if dir_name not in target_dirs:
        continue
    
    print(f"\n===== Processing {dir_name} =====")
    mean_t, std_t, mean_a, std_a = run_membrane_analysis(dir_name)
    summary_list.append({
        "dir": dir_name,
        "mean_thickness": mean_t,
        "std_thickness": std_t,
        "mean_APL": mean_a,
        "std_APL": std_a
    })

summary_df = pd.DataFrame(summary_list)
summary_df = summary_df.set_index("dir").reindex(target_dirs).reset_index()
summary_df.to_csv(results_directory / "membrane_summary.csv", index=False)

print("\nAll analyses complete. Summary table saved.")
print(summary_df)

