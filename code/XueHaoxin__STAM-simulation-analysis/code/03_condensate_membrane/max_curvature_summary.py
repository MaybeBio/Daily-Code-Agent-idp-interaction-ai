#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Open the folder containing the input data files.
# Run from the directory containing: <replicate_folder>/mean_curvature/<system_name>_mean_curvature_max.csv

# In[8]:


import os
import glob
import pandas as pd

pattern = "*/mean_curvature/*_mean_curvature_max.csv"  
out_file = "all_mean_curvature_summary.csv"
rows = []

for csv_path in sorted(glob.glob(pattern, recursive=True)):
    parts = csv_path.split(os.sep)
    if len(parts) < 3:
        continue
    first_sub = parts[0]       
    fname = os.path.basename(csv_path)
    sys_name = fname.replace("_mean_curvature_max.csv", "")

    df = pd.read_csv(csv_path)
    value_map = dict(zip(df["Leaflet"].astype(str).str.strip(), df["Mean_Curvature_Max"].astype(float)))

    upper_ang = value_map.get("Upper Leaflet")
    lower_ang = value_map.get("Lower Leaflet")

    rows.append({
        "first_subdir": first_sub,     
        "dir": sys_name,
        "upper-Mean_Curvature_Max(Å⁻¹)": upper_ang,
        "lower-Mean_Curvature_Max(Å⁻¹)": lower_ang,
        "upper-Mean_Curvature_Max(nm⁻¹)": upper_ang * 10 if upper_ang is not None else None,
        "lower-Mean_Curvature_Max(nm⁻¹)": lower_ang * 10 if lower_ang is not None else None,
    })

summary = pd.DataFrame(rows)
summary.to_csv(out_file, index=False, float_format="%.15g")
print(summary.to_string(index=False))


# In[10]:


import pandas as pd

unit = "Å⁻¹" 

df = pd.read_csv("all_mean_curvature_summary.csv")

def pivot_metric(col):
    p = df.pivot_table(index="dir", columns="first_subdir",
                       values=col, aggfunc="mean")   # Equivalent to direct lookup when no duplicate entries exist
    p.columns = [c.replace("Analysis-", "") for c in p.columns]
    p = p.reindex(columns=sorted(p.columns)) 
    return p

upper = pivot_metric(f"upper-Mean_Curvature_Max({unit})").add_prefix("upper_")
lower = pivot_metric(f"lower-Mean_Curvature_Max({unit})").add_prefix("lower_")

out = upper.join(lower).reset_index() 

out.to_csv("pivoted_summary_upper_lower.csv", index=False, float_format="%.15g")
print(out.to_string(index=False))


# In[ ]:




