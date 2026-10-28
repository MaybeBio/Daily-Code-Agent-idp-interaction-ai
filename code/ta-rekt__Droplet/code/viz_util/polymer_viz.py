"""
Visualize polymers as horizontal colored bars.

Each bead is represented as a segment in a horizontal bar, colored according
to its monomer type. This allows quick visual inspection of polymer composition
and patterning.

Usage:
    from polymer_viz import polymer_bar
    import json
    
    with open('path/to/data.json') as f:
        data = json.load(f)
    
    x, y, z, lab = data["final_coordinates_poly"]
    
    mono_colors = {
        1: "yellow",
        2: "red",
        3: "blue",
    }
    
    mono_labels = {
        1: "A",
        2: "B",
        3: "C",
    }
    
    fig = polymer_bar(lab, mono_colors, mono_labels=mono_labels, title="My Polymer")
    fig[0].show()
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def polymer_bar(monomer_types, mono_colors, mono_labels=None, 
                title="Polymer", figsize=(14, 2), bar_height=0.5, ax=None):
    """
    Visualize a polymer as a horizontal colored bar.
    
    Parameters
    ----------
    monomer_types : list or array
        List of monomer type labels for each bead (integers or strings).
    mono_colors : dict
        Mapping from monomer type to color name or hex code.
        e.g., {1: "yellow", 2: "red", 3: "blue"}
    mono_labels : dict, optional
        Mapping from monomer type to label strings to display below each region.
        e.g., {1: "A", 2: "B", 3: "C"}
        One label is displayed per contiguous region of each monomer type.
        If None, no labels are displayed.
    title : str, optional
        Title for the plot.
    figsize : tuple, optional
        Figure size (width, height). Default (14, 2).
    bar_height : float, optional
        Height of the bar (in data units). Default 0.5.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, creates a new figure.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    ax : matplotlib.axes.Axes
        The axes object.
    """
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    n_beads = len(monomer_types)
    
    # Create rectangles for each bead
    y_center = 0.5
    for i, mono_type in enumerate(monomer_types):
        color = mono_colors.get(mono_type, "gray")
        rect = mpatches.Rectangle(
            (i, y_center - bar_height/2),
            1, bar_height,
            linewidth=0,
            edgecolor='none',
            facecolor=color
        )
        ax.add_patch(rect)
    
    # Add labels for contiguous regions (one label per monomer type region)
    if mono_labels is not None:
        i = 0
        while i < n_beads:
            current_type = monomer_types[i]
            if current_type in mono_labels:
                # Find the end of this contiguous region
                j = i + 1
                while j < n_beads and monomer_types[j] == current_type:
                    j += 1
                # Place label at the center of the region
                center_pos = (i + j) / 2.0
                label = mono_labels[current_type]
                ax.text(center_pos, y_center - bar_height/2 - 0.05, label,
                       ha='center', va='top', fontsize=14)
                i = j
            else:
                i += 1
    
    # Set axis limits and remove axes
    ax.set_xlim(0, n_beads)
    ax.set_ylim(0, bar_height)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Add title
    ax.set_title(f'n = {n_beads} beads', ha='center', va='top',
            fontsize=13, style='italic', color='black')
    
    return fig, ax


def polymer_bar_with_legend(monomer_types, mono_colors, mono_labels=None, 
                            title="Polymer", figsize=(14, 2.5), bar_height=0.5):
    """
    Visualize a polymer as a horizontal colored bar with a legend.
    
    Parameters
    ----------
    monomer_types : list or array
        List of monomer type labels for each bead.
    mono_colors : dict
        Mapping from monomer type to color.
    mono_labels : dict, optional
        Mapping from monomer type to label strings.
    title : str, optional
        Title for the plot.
    figsize : tuple, optional
        Figure size (width, height).
    bar_height : float, optional
        Height of the bar.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    (ax_bar, ax_legend) : tuple of matplotlib.axes.Axes
        Bar axes and legend axes.
    """
    
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.3)
    
    ax_bar = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])
    
    # Plot polymer bar
    polymer_bar(monomer_types, mono_colors, mono_labels=mono_labels, 
                title=title, bar_height=bar_height, ax=ax_bar)
    
    # Create legend
    unique_types = sorted(set(monomer_types))
    patches = [mpatches.Patch(facecolor=mono_colors.get(t, "gray"), 
                              label=f'Type {t}')
               for t in unique_types]
    ax_legend.legend(handles=patches, loc='center left', fontsize=10)
    ax_legend.axis('off')
    
    return fig, (ax_bar, ax_legend)


def polymer_comparison(polymer_list, mono_colors, mono_labels=None, 
                       titles=None, figsize=None, bar_height=0.4):
    """
    Visualize multiple polymers for comparison.
    
    Parameters
    ----------
    polymer_list : list of arrays
        List of monomer type arrays.
    mono_colors : dict
        Mapping from monomer type to color.
    mono_labels : dict, optional
        Mapping from monomer type to label strings.
    titles : list of str, optional
        Titles for each polymer.
    figsize : tuple, optional
        Figure size. Default auto-scales based on max polymer length.
    bar_height : float, optional
        Height of each bar.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axes : list of matplotlib.axes.Axes
        One axes per polymer.
    """
    
    n_polymers = len(polymer_list)
    max_length = max(len(p) for p in polymer_list)
    
    if figsize is None:
        figsize = (max(12, max_length / 3), max(2, n_polymers * 1.5))
    
    if titles is None:
        titles = [f'Polymer {i+1}' for i in range(n_polymers)]
    
    fig, axes = plt.subplots(n_polymers, 1, figsize=figsize, 
                              gridspec_kw={'hspace': 0.5})
    if n_polymers == 1:
        axes = [axes]
    
    for ax, poly, title in zip(axes, polymer_list, titles):
        polymer_bar(poly, mono_colors, mono_labels=mono_labels, 
                   title=title, bar_height=bar_height, ax=ax)
        ax.set_xlim(0, max_length)
    
    return fig, axes


if __name__ == "__main__":
    # Example usage
    import numpy as np
    
    # Create example polymer with 3 monomer types
    monomer_types = np.array([1, 1, 1, 2, 2, 1, 3, 3, 1, 1, 2, 2, 2, 1, 3])
    
    mono_colors = {
        1: "yellow",
        2: "red",
        3: "blue",
    }
    
    mono_labels = {
        1: "A",
        2: "B",
        3: "C",
    }
    
    # Single polymer visualization
    fig, ax = polymer_bar(monomer_types, mono_colors, mono_labels=mono_labels,
                         title="Example Polymer - Single View")
    plt.tight_layout()
    plt.savefig("polymer_bar_single.png", dpi=150, bbox_inches='tight')
    plt.show()
    
    # Polymer with legend
    fig, (ax_bar, ax_leg) = polymer_bar_with_legend(
        monomer_types, mono_colors, mono_labels=mono_labels,
        title="Example Polymer - With Legend")
    plt.tight_layout()
    plt.savefig("polymer_bar_legend.png", dpi=150, bbox_inches='tight')
    plt.show()
    
    # Multiple polymer comparison
    poly1 = np.array([1, 1, 1, 2, 2, 1, 3, 3, 1, 1])
    poly2 = np.array([2, 2, 1, 1, 3, 3, 3, 2, 1, 1, 1, 1])
    poly3 = np.array([1, 2, 3, 1, 2, 3])
    
    fig, axes = polymer_comparison(
        [poly1, poly2, poly3],
        mono_colors,
        mono_labels=mono_labels,
        titles=["Polymer A", "Polymer B", "Polymer C"]
    )
    plt.tight_layout()
    plt.savefig("polymer_bar_comparison.png", dpi=150, bbox_inches='tight')
    plt.show()
