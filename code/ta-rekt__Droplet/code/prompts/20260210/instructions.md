I would like to be able to initialize a simulation from a lattice state, output by a previous simulation. I'd like the API for the init function to be as follows, with two additional optional keyword arguments:
coords_keys = ["final_coordinates_poly", "final_coordinates_tf", "final_coordinates_rnap"]
file_path = "path/to/output_cpu1.json"
init!(l, [poly, tf, rnap], file_path, coords_keys)

the function would then set the lattice and the full state of each species (pos, nMonomerTypes, keys_vec, customPattern, patternTails, and every other relevant state variable)