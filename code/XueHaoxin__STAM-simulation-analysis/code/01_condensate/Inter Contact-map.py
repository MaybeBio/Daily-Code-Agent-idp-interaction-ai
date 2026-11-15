#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS
import matplotlib.pyplot as plt
from tqdm import tqdm


# In[2]:


from vdw_def import *


# In[3]:


# Open the folder containing the input data files.
# cd 

# In[4]:


TPR=r'protein.tpr'
XTC=r'protein_only.xtc'


# In[5]:


u = mda.Universe(TPR,XTC)


# In[6]:


u.residues.resindices


# In[7]:


# ================= 1. Analysis parameters =================
cutoff=6  # Contact cutoff 6 angstroms (0.6 nm)
nchain=10  
start=15000  #analyze final 5001 frame
end=len(u.trajectory)
step=1


# In[8]:


protein = u.select_atoms("protein")


# In[9]:


nres=protein.n_residues
resname=protein.residues.resnames
nres=len(resname) 
atom_name=protein.names 
resindex = protein.resindices 


# In[10]:


# ================= 2. Initialize contact matrices =================
intra = np.zeros((nres, nres))
inter = np.zeros((nres, nres))
vdwall = np.zeros((nres, nres), dtype=int)


# In[11]:

frame_foldfold = []
frame_foldidr = []
frame_idridr = []
len_chain=540


# In[12]:


for tdx, ts in enumerate(tqdm(u.trajectory[start:end:step], desc="Analyzing trajectory")):
    status_intra = np.zeros((nres, nres), dtype=int)
    status_inter = np.zeros((nres, nres), dtype=int)
    status_vdw = np.zeros((nres, nres), dtype=int)
    
    nsr = FastNS(cutoff=cutoff, coords=protein.positions, box=ts.dimensions, pbc=True)
    result = nsr.self_search()
    nb = result.get_pairs()

    for pair in nb:
        res1 = resindex[pair[0]]
        res2 = resindex[pair[1]]
        name1 = atom_name[pair[0]]
        name2 = atom_name[pair[1]]
        chain1 = get_chainID(res1, nchain, u)
        chain2 = get_chainID(res2, nchain, u)

        status_vdw[res1, res2] = 1
        status_vdw[res2, res1] = 1

        # Classify intrachain and interchain contacts
        if chain1 == chain2:
            status_intra[res1, res2] = 1           
            status_intra[res2, res1] = 1
        else:
            status_inter[res1, res2] = 1
            status_inter[res2, res1] = 1

    intra += status_intra
    inter += status_inter
    vdwall += status_vdw

    # Calculate regional contact sums for the current frame
    M_frame = aggregate_by_local_resid(status_inter, nchain, len_chain, u)
    if M_frame.shape == (540,540):
        sum_A = np.sum(M_frame[0:390, 0:390])
        sum_B1 = np.sum(M_frame[0:390, 390:540])
        sum_B2 = np.sum(M_frame[390:540, 0:390])
        sum_C = np.sum(M_frame[390:540, 390:540])

        ff = sum_A / 2
        fid = (sum_B1 + sum_B2) / 2
        ii = sum_C / 2

        frame_foldfold.append(ff)
        frame_foldidr.append(fid)
        frame_idridr.append(ii)
    else:
        print(f"Unexpected contact matrix shape at frame {ts.frame}")
    # =====================================================
    
print("\nAnalysis complete. Saving results...")
print(intra.shape)
print(inter.shape)
print(vdwall.shape)


# In[13]:


import os 
os.makedirs("Inter_contacts-frame", exist_ok=True)
output_dir='Inter_contacts-frame'


# In[14]:

np.save(os.path.join(output_dir, 'vdwall.npy'), vdwall)
np.save(os.path.join(output_dir, 'inter.npy'), inter)
np.save(os.path.join(output_dir, 'intra.npy'), intra)


# In[15]:


# =================3. Post-processing and output=================
calc_frame = (end - start - 1) // step + 1  # Number of analyzed frames


# In[16]:



inter_avg = inter / calc_frame

inter_resid = aggregate_by_local_resid(inter_avg, nchain, len_chain, u)

pd.DataFrame(inter_resid).to_csv(os.path.join(output_dir, "inter_resid.csv"), header=False, index=False)


# In[17]:

n_foldfold = inter_resid[0:390,0:390].size / 2
n_foldidr  = inter_resid[0:390,390:540].size
n_idridr   = inter_resid[390:540,390:540].size / 2

perframe_avg_ff = np.array(frame_foldfold) / n_foldfold
perframe_avg_fid = np.array(frame_foldidr) / n_foldidr
perframe_avg_ii = np.array(frame_idridr) / n_idridr

avg_foldfold = np.mean(perframe_avg_ff)
std_foldfold = np.std(perframe_avg_ff, ddof=1) # Sample SD across analyzed frames
avg_foldidr  = np.mean(perframe_avg_fid)
std_foldidr  = np.std(perframe_avg_fid, ddof=1)
avg_idridr   = np.mean(perframe_avg_ii)
std_idridr   = np.std(perframe_avg_ii, ddof=1)

print(f"avg fold-fold: {avg_foldfold:.4f} ± {std_foldfold:.4f}")
print(f"avg fold-idr:  {avg_foldidr:.4f} ± {std_foldidr:.4f}")
print(f"avg idr-idr:   {avg_idridr:.4f} ± {std_idridr:.4f}")

df_frames = pd.DataFrame({
    "fold_fold": frame_foldfold,
    "fold_idr": frame_foldidr,
    "idr_idr": frame_idridr,
    "avg_foldfold": perframe_avg_ff,
    "avg_foldidr": perframe_avg_fid,
    "avg_idridr": perframe_avg_ii
})
df_frames.to_csv(os.path.join(output_dir, "per_frame_contact.csv"), index=False)
# =====================================================


# In[18]:


# ====================== 4. Plot the interchain contact map ======================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

plt.figure(figsize=(12,10))
matrix = inter_resid

from matplotlib.colors import LinearSegmentedColormap
color_list = [
    "#ffffff",
    "#87ceeb",
    "#2060c0",
    "#401060"
]
custom_cmap = LinearSegmentedColormap.from_list("white_blue_purple", color_list, N=200)

im = plt.imshow(matrix, cmap=custom_cmap, interpolation='nearest',
                aspect='equal', origin='lower'
                #,vmin=0.0,vmax=0.6
               ) 

# residue index = 390 
plt.axvline(x=390, linestyle='--', color='black', lw=1.2)
plt.axhline(y=390, linestyle='--', color='black', lw=1.2)

# colorbar 
cbar = plt.colorbar(im, shrink=0.82, pad=0.03)
cbar.set_label('Contact Frequency', fontsize=16)
cbar.ax.tick_params(labelsize=14)

plt.title("Inter-chain Contact Map", fontsize=20, pad=15)
plt.xlabel('Residue Index', fontsize=16)
plt.ylabel('Residue Index', fontsize=16)

plt.tick_params(axis='both', labelsize=16)

n = matrix.shape[0]

save_path = os.path.join(output_dir, 'Fig1_Inter_Contact_Map.pdf')
plt.savefig(save_path, bbox_inches='tight', dpi=600)
plt.show()


