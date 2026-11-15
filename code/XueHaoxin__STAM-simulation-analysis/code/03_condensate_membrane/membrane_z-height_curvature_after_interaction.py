#!/usr/bin/env python
# coding: utf-8

# In[1]:


from datetime import datetime

start_time = datetime.now()
print(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")


# In[2]:


#pwd


# In[3]:


#cd PA_Ca50  enter data dir


# In[4]:


dir='PA_Ca50'  #stam_cluster_withPA_Ca50


# In[5]:


TPR=r"non-Solvent_pbcwhole_conect.pdb"
XTC=r"non-Solvent_pbcwhole.xtc"


# In[6]:


import MDAnalysis as mda
import numpy as np
import pandas as pd

from membrane_curvature.base import MembraneCurvature
import matplotlib.pyplot as plt
import nglview as nv
from scipy import ndimage
from tqdm.notebook import tqdm


# In[7]:


u=mda.Universe(TPR,XTC)


# In[8]:


u.dimensions


# In[9]:


protein = u.select_atoms('protein')
not_protein = u.select_atoms('not protein')
for ts in tqdm(u.trajectory, desc="Centering the condensate in the xy plane"):
    protein.unwrap(compound='fragments')
    
    protein_center = protein.center_of_geometry(wrap=True)
    #print(protein_center)
    dim = ts.triclinic_dimensions
    box_center = np.sum(dim, axis=0) / 2
    #print(box_center)
    translation = box_center - protein_center
    translation[2] = 0.0
    #print(translation)
    u.atoms.translate(translation)
    not_protein.wrap(compound='residues')


# ## 2.identify the upper and lower leaflets

# In[10]:


protein  = u.select_atoms("protein")
print("\nProtein has {} residues." \
      "".format(protein.n_residues))

lipids  = u.select_atoms("resname POPC POPA POPE CHOL")
print("\nLipids molecules include {} residues and {} atoms." \
      "".format(lipids.n_residues, lipids.n_atoms))


# In[11]:


P_headgroups = lipids.select_atoms('name PO4')


# In[12]:


from MDAnalysis.analysis.leaflet import LeafletFinder

L = LeafletFinder(u, 'name PO4', cutoff=20)

upper_leaflet = L.groups(0) # upper leaflet
lower_leaflet = L.groups(1) # lower leafet


# In[13]:


type(upper_leaflet)


# In[14]:


leaflets = ['Lower', 'Upper']

upper_leaflet_P = upper_leaflet.select_atoms("name PO4")
lower_leaflet_P = lower_leaflet.select_atoms("name PO4")


for name, new_lf in zip(leaflets,[lower_leaflet_P, upper_leaflet_P]):
    print("{} leaflet includes {} elements.".format(name, len(new_lf)))


# In[15]:


sel_upper = " ".join([str(r) for r in upper_leaflet.residues.resids])
sel_lower = " ".join([str(r) for r in lower_leaflet.residues.resids])

upper_string = "resid {} and name PO4".format(sel_upper)
lower_string = "resid {} and name PO4".format(sel_lower)


# Calculate membrane curvature

# In[16]:


curvature_upper_leaflet_stop = MembraneCurvature(u,
                                                    select=upper_string,
                                                    surface_method='binning',
                                                    n_x_bins=19,
                                                    n_y_bins=19,
                                                    wrap=True
                                                    ).run(start=1600,verbose=True)


# In[17]:


curvature_lower_leaflet_stop = MembraneCurvature(u,
                                                    select=lower_string,
                                                    surface_method='binning',
                                                    n_x_bins=19,
                                                    n_y_bins=19,
                                                    wrap=True
                                                    ).run(start=1600,verbose=True)


# In[18]:


surface_upper_leaflet_stop=curvature_upper_leaflet_stop.results.average_z_surface


# In[19]:


normalize_surface_upper=surface_upper_leaflet_stop-surface_upper_leaflet_stop.mean()


# In[20]:


surface_upper_leaflet_stop.mean()


# In[21]:


surface_lower_leaflet_stop=curvature_lower_leaflet_stop.results.average_z_surface


# In[22]:


normalize_surface_lower=surface_lower_leaflet_stop-surface_lower_leaflet_stop.mean()


# In[23]:


surface_lower_leaflet_stop.mean()


# In[24]:


surface_upper_leaflet_stop_meancurvature=curvature_upper_leaflet_stop.results.average_mean


# In[25]:


surface_lower_leaflet_stop_meancurvature=curvature_lower_leaflet_stop.results.average_mean

# In[26]:


def plot_contours(title, results, label, levels_, cm, 
                  vmin=None, vmax=None, save_path=None,):
    """
    Plot contours of MembraneCurvature results.

    Additional controls for colorbar range (vmin, vmax) and saving the figure.

    Parameters
    ----------
    results : list
        List with results by leaflets as elements [lower_leaflet, upper_leaflet]
    label : str
        Label for the colorbar.
    levels_ : int
        Number of contour levels.
    cm : str
        Colormap.
    vmin, vmax : float or None
        Fixed color scale limits. If None, automatic symmetric or data range is used.
    save_path : str or None
        If provided, the figure is saved to this path.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage 

    leaflets = ['Upper','Lower']

    fig, [ax1, ax2] = plt.subplots(ncols=2, figsize=(4, 3.5), dpi=200)

    finite = [rs[np.isfinite(rs)] for rs in results]
    max_ = max((np.max(np.abs(v)) for v in finite if v.size), default=0.0)

    plt.suptitle(title, fontsize=8, y=0.8)

    for ax, rs, lf in zip((ax1, ax2), results, leaflets):
        rs = np.asarray(rs).T
        valid_mask = np.isfinite(rs)

        if vmin is None and vmax is None:
            if np.min(rs[valid_mask]) < 0 < np.max(rs[valid_mask]):
                use_vmin, use_vmax = -max_, max_
                sym = True
            else:
                use_vmin, use_vmax = np.min(rs[valid_mask]), np.max(rs[valid_mask])
                sym = False
        else:

            if vmin is None or vmax is None:
                if np.min(rs[valid_mask]) < 0 < np.max(rs[valid_mask]):
                    use_vmin, use_vmax = -max_, max_
                    sym = True
                else:
                    use_vmin, use_vmax = np.min(rs[valid_mask]), np.max(rs[valid_mask])
                    sym = False
            else:
                use_vmin, use_vmax = vmin, vmax
                sym = use_vmin < 0 < use_vmax

        levs = np.linspace(use_vmin, use_vmax, levels_)

        im = ax.contourf(
            rs, cmap=cm, origin='lower', levels=levs, alpha=0.95,
            vmin=use_vmin, vmax=use_vmax
        )

        tcs = [use_vmin, 0, use_vmax] if sym else [use_vmin, use_vmax]

        ax.set_aspect('equal')
        ax.set_title(f'{lf} Leaflet', fontsize=6)
        ax.axis('off')
        
        cbar = plt.colorbar(
            im, ticks=tcs, orientation='horizontal', ax=ax,
            shrink=0.7, aspect=10, pad=0.05
        )

        cbar.ax.tick_params(labelsize=4, width=0.5)
        cbar.set_label(label, fontsize=6, labelpad=2)

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    
    plt.show()


# In[27]:


plot_contours(
    title=f"{dir}",
    results=[normalize_surface_upper, normalize_surface_lower],
    label="Å",
    levels_=20,
    cm='YlGnBu_r',
    vmin=-50,  # min
    vmax=30,   # max
    save_path="z_surface_membrane.png"  
)


# In[28]:


plot_contours(
    title=f"{dir}",
    results=[surface_upper_leaflet_stop_meancurvature, surface_lower_leaflet_stop_meancurvature],
    label=r"$H$ ($\AA^{-1}$)",
    levels_=50,
    cm='bwr',
    vmin=-2E-2,  # min
    vmax=2E-2,   # max
    save_path="mean_curvature.png"
)


# In[30]:


mean_max_up = surface_upper_leaflet_stop_meancurvature.max()
mean_max_low = surface_lower_leaflet_stop_meancurvature.max()


# In[31]:


import os
import csv

with open("mean_curvature_max.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Leaflet", "Mean_Curvature_Max", "Unit"])
    writer.writerow(["Upper Leaflet", mean_max_up, "Å^-1"])
    writer.writerow(["Lower Leaflet", mean_max_low, "Å^-1"])

print("Results saved to: mean_curvature_max.csv")


# In[ ]:




