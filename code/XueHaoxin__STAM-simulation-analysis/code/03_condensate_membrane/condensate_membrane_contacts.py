#!/usr/bin/env python
# coding: utf-8

# In[1]:


from datetime import datetime

start_time = datetime.now()
print(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")


# In[2]:


import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import contacts
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import csv
from get_chainID_resid import *


# In[3]:


#pwd


# In[4]:


#cd PA_Ca50  #stam_cluster_withPA_Ca50


# In[5]:


TPR=r"md.tpr"  
XTC=r"md_pbcwhole.xtc"


# In[6]:


u=mda.Universe(TPR,XTC)

dt=u.trajectory.dt/1000 #ns
print(f"{dt}ns, {len(u.trajectory)}frame")


# In[7]:


membrane=u.select_atoms("resname POPC POPE POPA CHOL")
stam=u.select_atoms("protein")

print(membrane.segments.segids)
print(stam.segments.segids)


# In[8]:


stam.resindices #The first STAM residue index is 5406 in this system

# In[9]:


with open("contact_lipid_number.csv", "w", newline="") as f_count, \
     open("contact_details.csv", "w", newline="") as f_detail, \
     open("contact_timeseries.csv", "w", newline="") as f_time:

    count_writer = csv.writer(f_count)
    detail_writer = csv.writer(f_detail)
    time_writer = csv.writer(f_time)

    count_writer.writerow(["frame", "time", "n_total", "n_POPC", "n_POPE", "n_POPA", "n_CHOL"])
    detail_writer.writerow(["frame", "time",
                            "stam_atom_name", "stam_resname", 
                            "stam_resindex", "stam_chainID", "stam_resid",
                            "lipid_atom_name", "lipid_resname",
                            "contact_distance"])
    time_writer.writerow(["frame", "time", "n_contacts"])

    #Start contact analysis
    for ts in tqdm(u.trajectory, desc="	Calculating condensate–membrane contacts", unit="frame"):
    #for ts in u.trajectory[::2000]:
        frame=ts.frame
        time=dt*frame  #ns

        contact_membrane=membrane.select_atoms("around 6 global group topsel", topsel=stam, updating=True)
        n_total=contact_membrane.n_residues
        n_popc=contact_membrane.select_atoms("resname POPC").n_residues
        n_pope=contact_membrane.select_atoms("resname POPE").n_residues
        n_popa=contact_membrane.select_atoms("resname POPA").n_residues
        n_chol=contact_membrane.select_atoms("resname CHOL").n_residues
        #print(f"Total number of contacting lipids {n_total}\t:POPC {n_popc}\tPOPE {n_pope}\tPOPA {n_popa}\tCHOL {n_chol}")
        count_writer.writerow([frame, time, n_total, n_popc, n_pope, n_popa, n_chol])

        #Calculate pairwise contact distances
        dist_arr=distances.distance_array(stam.positions,contact_membrane,box=ts.dimensions)
    
        contact_matrix = contacts.contact_matrix(dist_arr, radius=6.0)
        n_contact=contact_matrix.sum()
        time_writer.writerow([frame, time, n_contact])
    
        for stam_idx, lipid_idx in np.argwhere(contact_matrix):

            contact_distance = dist_arr[stam_idx, lipid_idx]
        
            stam_atom  = stam.atoms[stam_idx]
            stam_atom_name  = stam_atom.name
            stam_res  = stam_atom.residue
            stam_resindex = stam_res.resindex
            stam_resname  = stam_res.resname

            stam_chainID=get_chainID(stam_resindex-5406,9,540,stam) #The first STAM residue index is 5406 in this system
            stam_resid=get_local_resid(stam_resindex-5406,9,540)
        
            lipid_atom = contact_membrane.atoms[lipid_idx]    
            lipid_atom_name = lipid_atom.name      
            lipid_res = lipid_atom.residue            
            lipid_resname  = lipid_res.resname
        
            detail_writer.writerow([
                frame, time,
                stam_atom_name, stam_resname, 
                stam_resindex,stam_chainID, stam_resid,
                lipid_atom_name, lipid_resname,
                round(contact_distance, 3)
                ])


# In[10]:


end_time = datetime.now()
print(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
