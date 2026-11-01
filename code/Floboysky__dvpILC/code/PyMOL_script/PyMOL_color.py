#!/usr/bin/env python3

import os
import json
import pandas as pd
from pymol import cmd


def color_by_residues(file, key="Conservation", sele="sele", grad="rainbow_rev"):
    """
    Color the residues based on the values in the jalview_output.json or data_all.csv files.
    Used directly in PyMOL.
    
    Parameters:
        file: Path to the file containing the values to be used for coloring.
        key: Key in the file containing the values to be used for coloring (default: "Conservation").
        sele: PyMOL selection to color (default: “sele”).
        grad: List of colors defining the gradient (e.g., “blue_white_red”).
    """
    
    # Check the file extension
    ext = os.path.splitext(file)[1].lower()
    
    if ext == ".json":
        # Load JSON file
        with open(file, "r") as f:
            data = json.load(f)
            
        nums = data["num"]
        values = data[key]
    
        # Ungroup the lists
        if isinstance(nums, list) and isinstance(nums[0], list):
            nums = nums[0]
        if isinstance(values, list) and isinstance(values[0], list):
            values = values[0]
            
    elif ext == ".csv":
        # Load CSV file
        data = pd.read_csv(file, sep=";", encoding="utf-8")
        
        nums = data["Position"].tolist()
        values = data[f"{key}"].tolist()
    
    else:
        raise ValueError("Unsupported file format. Please use the 'jalview_output.json' or 'data_all.csv' file.")
    
    # Find the minimum and maximum values to define the gradient scale
    min_val = min(values)
    max_val = max(values)
    if max_val == min_val:  # Avoid division by zero if all values are identical
        raise ValueError(f"All values for '{key}' are identical, impossible to generate a gradient.")

    # Apply the colors to the residues based on the key value
    for resi, val in zip(nums, values):
        cmd.alter(f"resi {resi} and {sele}", f"b = {val}")

    # Apply a predefined gradient to the selection
    cmd.spectrum("b", grad, sele)


# Main
cmd.extend("color_by_residues", color_by_residues)
