#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def analysis_contact_data(work_dir):
    """
    Analyze contact data in a single subdirectory.
    :param work_dir: Full path to the subdirectory
    :return: Tuple of total-contact and lipid-group statistics dictionaries; (None, None) if the input file is missing.
    """
    f_detail = os.path.join(work_dir, "contact_details.csv")
    
    if not os.path.exists(f_detail):
        print(f"[Skipping {os.path.basename(work_dir)}] Missing file: contact_details.csv")
        return None, None

    folder_name = os.path.basename(work_dir)

    df = pd.read_csv(f_detail)
    df_filter = df[df['time'] >= 800].copy()

    if df_filter.empty:
        print(f"Warning: {folder_name} contains no data with time >= 800 ns.")

    ALL_RESIDS = np.arange(1, 541)          # Protein residues 1–540
    all_frames = df_filter['frame'].unique()

    groups = {
        "all": None,
        "PO4": ["PO4"],
        "pos_head": ["NC3", "NH3"],
        "polar": ["GL1", "GL2", "ROH"],
        "nopolar": [
            "C1", "C1A", "C1B", "C2", "C2B", "C3A", "C3B", "C4A", "C4B",
            "D2A", "R1", "R2", "R3", "R4", "R5", "R6"
        ]
    }

    def compute_avg_contact(sub_df, atom_list):
        """Return the mean contact count per residue over the frames in all_frames."""
        if atom_list is not None:
            sub_df = sub_df[sub_df['lipid_atom_name'].isin(atom_list)]
        # Count contacts by frame and residue
        contact_series = sub_df.groupby(['frame', 'stam_resid']).size()
        # Build the complete frame–residue index and fill missing entries with zero
        full_idx = pd.MultiIndex.from_product(
            [all_frames, ALL_RESIDS], names=['frame', 'stam_resid']
        )
        contact_full = contact_series.reindex(full_idx, fill_value=0).reset_index(name='count')
        avg = contact_full.groupby('stam_resid')['count'].mean().reset_index(name='avg')
        return avg

    avg_dict = {}
    for gname, atoms in groups.items():
        avg_df = compute_avg_contact(df_filter, atoms)
        avg_df = avg_df.rename(columns={'avg': f'avg_{gname}'})
        avg_dict[gname] = avg_df

    merged = avg_dict['all']
    for g in ['PO4', 'pos_head', 'polar', 'nopolar']:
        merged = merged.merge(avg_dict[g], on='stam_resid', how='outer')
    merged = merged[['stam_resid', 'avg_all', 'avg_PO4', 'avg_pos_head', 'avg_polar', 'avg_nopolar']]

    csv_path = os.path.join(work_dir, "lipid_avg_contacts.csv")
    merged.to_csv(csv_path, index=False, float_format='%.4f')

    plt.figure(figsize=(14, 6))
    plt.plot(merged['stam_resid'], merged['avg_all'],      label='All atoms',       color='black')
    plt.plot(merged['stam_resid'], merged['avg_PO4'],      label='PO4 head',        color='red')
    plt.plot(merged['stam_resid'], merged['avg_pos_head'], label='Charged head (N)', color='blue')
    plt.plot(merged['stam_resid'], merged['avg_polar'],    label='Polar (GL/ROH)',  color='green')
    plt.plot(merged['stam_resid'], merged['avg_nopolar'],  label='Nonpolar tails',  color='orange')
    plt.xlabel('Protein residue ID')
    plt.ylabel('Average contact count per frame')
    plt.title(f'Lipid-protein average contacts ({folder_name}, time ≥ 800)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    fig_path = os.path.join(work_dir, "lipid_contact_profile.png")
    plt.savefig(fig_path, dpi=150)
    plt.close()

    contact_stats = dict(zip(merged['stam_resid'], merged['avg_all']))
    
    lipid_stats = {}
    for _, row in merged.iterrows():
        resid = int(row['stam_resid'])
        lipid_stats[resid] = {
            'PO4': row['avg_PO4'],
            'pos_head': row['avg_pos_head'],
            'polar': row['avg_polar'],
            'nopolar': row['avg_nopolar']
        }

    return contact_stats, lipid_stats


# In[2]:


if __name__ == "__main__":

    main_folder = "."

    all_subdirs = []
    for entry in os.listdir(main_folder):
        full_path = os.path.join(main_folder, entry)
        if os.path.isdir(full_path):
            all_subdirs.append(full_path)

    print(f"Found {len(all_subdirs)} subdirectories. Starting batch analysis...")
    total_summary_list = [] 
    lipid_summary_list = []  

    for sub_dir in all_subdirs:
        try:
            total_dict, lipid_dict = analysis_contact_data(sub_dir)
            if total_dict is not None and lipid_dict is not None:
                total_summary_list.append(total_dict)
                lipid_summary_list.append(lipid_dict)
        except Exception as e:
            print(f"[Error] Failed to process directory {sub_dir}: {str(e)}")
            continue


# In[ ]:




