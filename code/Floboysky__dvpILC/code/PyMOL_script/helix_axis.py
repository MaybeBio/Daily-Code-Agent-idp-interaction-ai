#!/usr/bin/env python3

import math
import numpy as np
from pymol import cmd
from pymol.cgo import BEGIN, END, LINEWIDTH, COLOR, VERTEX, LINES


"""
A program that places pseudo-atoms at the centroid of each helix turn (protein/DNA), then fits a 3D line passing as close as possible (using PCA). Can be used directly in PyMOL.
1) In PyMOL, navigate to the directory of the program and run "run helix_axis.py".
2) Select the residues of the protein or DNA (helix).
3) "helix_centroids" (for the protein), or "helix_centroids_dna" (for the DNA), to create the pseudo-atoms.
4) "fit_line" to fit the 3D line.
"""


def helix_centroids(selection="sele", step=3.6, name="helix_point"):
    """
    Creates pseudo-atoms at the centroid of each helix turn (~step residues per turn) for protein.
    Can be used directly in PyMOL.
    
    Parameters:
        selection (str): PyMOL selection of the protein
        step (float): Number of base pairs per helical turn (~3.6 Å)
    Returns:
        name (str): Prefix of the created pseudo-atoms (helix_point_1, helix_point_2...)
    """
    
    # Set the type of step
    try:
        step = float(step)
        if step <= 0:
            raise ValueError
    except Exception:
        print("Invalid 'step' argument (must be a positive number).")
        return

    # Retrieve alpha carbons (CA) in the selection
    model = cmd.get_model(f"{selection} and name CA")
    if not model.atom:
        print("No CA atoms found for the provided selection.")
        return

    # Get the sorted list of residues (converted to integers for more reliable sorting)
    residues_int = sorted({int(atom.resi) for atom in model.atom})

    # Construct the list of average coordinates (one point per residue)
    coords = []
    for r in residues_int:
        pts = [atom.coord for atom in model.atom if int(atom.resi) == r]
        if pts:
            coords.append(np.mean(pts, axis=0))
    coords = np.array(coords)
    n = len(coords)
    if n < 2:
        print("Not enough residues/CA atoms to calculate the centroid.")
        return

    # Starting indices for each block
    starts = np.arange(0.0, float(n), step)
    pseudo_coords = []
    for i, s in enumerate(starts):
        start = int(np.floor(s))
        end = int(min(np.floor(s + step), n))
        if end <= start:
            # If the block is too small (rounded), enlarge it by at least +1 residue
            end = min(start + 1, n)
        centroid = coords[start:end].mean(axis=0)
        pseudo_name = f"{name}_{i+1}"
        # Delete if already exists
        try:
            cmd.delete(pseudo_name)
        except Exception:
            pass
        cmd.pseudoatom(pseudo_name, pos=list(centroid))
        pseudo_coords.append((pseudo_name, centroid))

    print(f"{len(pseudo_coords)} pseudo-atoms created.")
    

def helix_centroids_dna(selection="sele", step=3.6, name="helix_point"):
    """
    Creates pseudo-atoms at the centroid of each helix turn (DNA).
    Can be used directly in PyMOL.

    Parameters:
        selection (str): PyMOL selection of the DNA
        step (float): Number of base pairs per helical turn (~3.6 Å)
    Returns:
        name (str): Prefix of the created pseudo-atoms (helix_point_1, helix_point_2...)
    """

    # Retrieve the sorted list of residues (by segi/chain/resi)
    model = cmd.get_model(selection)
    residues = []
    current = None
    atoms_buffer = []
    
    # Group atoms by residue (using segi, chain, resi as key)
    for atom in model.atom:
        key = (atom.segi, atom.chain, atom.resi)
        if key != current:
            if atoms_buffer:
                residues.append(atoms_buffer)
            atoms_buffer = [atom]
            current = key
        else:
            atoms_buffer.append(atom)

    # Add the last buffered residue
    if atoms_buffer:
        residues.append(atoms_buffer)

    # Check if any residues were found
    if len(residues) == 0:
        print("No results found in the selection.")
        return

    residues_per_turn = int(round(step))
    n_turns = int(math.ceil(len(residues) / residues_per_turn))
    count = 1
    
    # For each turn, calculate the centroid of the corresponding residues and create a pseudo-atom
    for i in range(n_turns):
        chunk = residues[i*residues_per_turn:(i+1)*residues_per_turn]
        if not chunk:
            continue

        coords = []
        for res in chunk:
            for atom in res:
                coords.append(atom.coord)

        if len(coords) == 0:
            continue

        # Centroid
        x = sum(c[0] for c in coords) / len(coords)
        y = sum(c[1] for c in coords) / len(coords)
        z = sum(c[2] for c in coords) / len(coords)

        pseudo_name = f"{name}_{count}"
        cmd.pseudoatom(object=pseudo_name, pos=[x, y, z])
        count += 1

    print(f"{count-1} pseudo-atoms created.")


def fit_line(selection="helix_point*", name="helix_axis", extend=5.0, couleur="red"):
    """
    Fit a 3D line (linear regression via PCA) that passes as close as possible to the pseudo-atoms previously created by 'helix_centroids' or 'helix_centroids_dna'.
    Can be used directly in PyMOL.
    
    Parameters:
        selection (str): PyMOL selection of the pseudo-atoms to fit
        name (str): Name of the PyMOL object created
        extend (float): Margin added on the line (in Å)
        couleur (str): Line color. Options: red, green, blue, yellow... Default: red (1.0, 0.2, 0.2)
    """
    
    if couleur == "red":
        color = (1.0, 0.2, 0.2)
    elif couleur == "green":
        color = (0.2, 1.0, 0.2)
    elif couleur == "blue":
        color = (0.2, 0.2, 1.0)
    elif couleur == "yellow":
        color = (1.0, 1.0, 0.0)
    elif couleur == "cyan":
        color = (0.0, 1.0, 1.0)
    elif couleur == "magenta":
        color = (1.0, 0.0, 1.0)
    elif couleur == "orange":   
        color = (1.0, 0.5, 0.0)
    elif couleur == "purple":
        color = (0.55, 0.15, 0.55)
    elif couleur == "aquamarine":
        color = (0.498, 1.000, 0.831)
    elif couleur == "sand":
        color = (0.961, 0.871, 0.702)
    elif couleur == "gray":
        color = (0.5, 0.5, 0.5)
    else:
        print(f"Color '{couleur}' not recognized, using red as default.")
        print("List of available colors: red, green, blue, yellow, cyan, magenta, orange, purple, aquamarine, sand and gray.")
        color = (1.0, 0.2, 0.2)
        
    model = cmd.get_model(selection)
    if not model.atom:
        print(f"No atoms found for '{selection}'.")
        return

    coords = np.array([atom.coord for atom in model.atom])
    centroid = coords.mean(axis=0)

    # PCA
    coords_centered = coords - centroid
    u, s, vh = np.linalg.svd(coords_centered)
    direction = vh[0]

    # Normalize the direction vector
    direction /= np.linalg.norm(direction)

    # Project all the points onto this line
    projections = coords_centered.dot(direction)
    min_proj, max_proj = projections.min(), projections.max()

    # Extreme points of the line in space
    p1 = centroid + (min_proj - extend) * direction
    p2 = centroid + (max_proj + extend) * direction

    # CGO construction
    obj = [
        BEGIN, LINES,
        LINEWIDTH, 10.0,
        COLOR, *color,
        VERTEX, *p1,
        VERTEX, *p2,
        END
    ]

    cmd.delete(name)
    cmd.load_cgo(obj, name)
    cmd.zoom(selection)
    print(f"Line '{name}' created between {p1.round(2)} and {p2.round(2)}.")


# Main
cmd.extend("helix_centroids", helix_centroids)
cmd.extend("helix_centroids_dna", helix_centroids_dna)
cmd.extend("fit_line", fit_line)
