

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage


#Make clustered distance tree


def gen_clustered_tree2(df, color_dict):
    """
    Generates a dendrogram from a given all-by-all distance matrix using UPGMA (Unweighted Pair Group Method with Arithmetic Mean).

    Parameters:
    df (pandas.DataFrame): A pandas DataFrame representing the all-by-all distance matrix.
    color_dict (dict): A dictionary mapping species names to colors for the x-tick labels.

    The dendrogram visualizes the hierarchical clustering of species based on the provided distance matrix.
    """
    # Convert the distance matrix to a NumPy array
    distance_matrix = np.asarray(df)

    # Species names corresponding to the distance matrix
    species_names = list(df.index)

    # Convert the distance matrix to a condensed form
    condensed_matrix = distance_matrix[np.triu_indices(len(distance_matrix), 1)]

    # Perform hierarchical clustering using UPGMA
    Z = linkage(condensed_matrix, method='average')

    # Plot the dendrogram
    fig, ax = plt.subplots()
    #fig.set_size_inches(15,10)#(12, 8)
    fig.set_size_inches(10,7)#(12, 8)
    dendrogram(Z, labels=species_names, leaf_font_size=10, ax=ax)

    # Set x-tick labels and apply colors
    plt.xticks(rotation=90)
    
    for tick_label in ax.get_xticklabels():
        tick_value = tick_label.get_text()
        tick_label.set_color(color_dict.get(tick_value, 'black'))  # Default to black if not in dictionary
    
    plt.title('UPGMA Nearest Neighbour Distances Between CDK Human IDRs')
    plt.xlabel('Species')
    plt.ylabel('Distance')
    plt.show()

# Example usage
# df = dataframe_here
# color_dict = {'species1': 'red', 'species2': 'blue', ...}
# gen_clustered_tree(df, color_dict)


#Set theme
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
sns.set(style='whitegrid', context='talk', font_scale=1.2)
plt.rcParams['pdf.fonttype'] = 42  # For editable text in Illustrator
plt.rcParams['ps.fonttype'] = 42
font_manager.fontManager.addfont("/home/moseslab/denise/Thesis/misc/arial.ttf") 
plt.rcParams['font.family'] = 'Arial'



import pandas as pd
import matplotlib.pyplot as plt

def make_pairwise_plot(df, x_col_a, x_col_b, 
                       y_col_a, y_col_b,
                       log_scale=False,
                       a_name= 'CDK IDRs',
                       b_name = 'CDK PFAMs',
                       select_idx=None, #or list
                       xlabel = "Human-Species Paralogs",
                        ylabel = "Human-Species Orthologs",
                        title= "Comparison of all CDK pairs vs Unrelated IDRs and domains between species Distances (n=100)"):
    
    if log_scale:
        x_a =  np.log10(df[x_col_a])
        y_a = np.log10(df[y_col_a])
        x_b = np.log10(df[x_col_b])
        y_b = np.log10(df[y_col_b])

    else:
        x_a = df[x_col_a]
        y_a = df[y_col_a]
        x_b = df[x_col_b]
        y_b = df[y_col_b]
    
    min_val = min(x_a.min(),y_a.min(),x_b.min(), y_b.min())
    max_val = max(x_a.max(),y_a.max(),x_b.max(), y_b.max())

    plt.figure(figsize=(6, 7))
    plt.scatter(x=x_a, y= y_a, alpha=0.5, label = a_name, color='#1e90ff',edgecolors='black', linewidths=0.5)
    plt.scatter(x=x_b, y= y_b, alpha=0.5,color='orange', label=b_name,edgecolors='black', linewidths=0.5)
    
    if select_idx: #spp name
        for idx in select_idx:
            x_select = x_a[select_idx]
            y_select = y_a[select_idx]
            plt.scatter(
                x=x_select, 
                y=y_select, 
                s=80, 
                edgecolors='red', 
                facecolors='none', 
                linewidth=2, 
                label=None,  # Avoid legend clutter,
                zorder=3
            )
            
            x_select = x_b[select_idx]
            y_select = y_b[select_idx]
            plt.scatter(
                x=x_select, 
                y=y_select, 
                s=80, 
                edgecolors='red', 
                facecolors='none', 
                linewidth=2, 
                label=None,  # Avoid legend clutter,
                zorder=3
            )
        
    
    #Annotate specific labels
    
    #Where differences between orthologs and paralogs are equal
    plt.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', linewidth=3)

    # Log scale change ticks to real numbers
    if log_scale:
        ax = plt.gca()
        ticks = ax.get_yticks()
        ax.set_yticklabels([f"{10**t:.0f}" for t in ticks])
        
        ticks = ax.get_xticks()
        ax.set_xticklabels([f"{10**t:.0f}" for t in ticks])

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.legend()
    #plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.grid(True)
    #plt.tight_layout()
    plt.show()