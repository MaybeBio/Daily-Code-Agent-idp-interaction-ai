import chimerax
from chimerax.core.commands import run

from typing import List
import os, json
import numpy as np


def ensemble_in_background( pdb_file: str, chain1: str, chain2: str,
							predicted_interface1: List, predicted_interface2: List,
							out_path: str, name: str ):
	"""
	Given a PDB containing an NMR structure, show one representative
		model with the other models having high transparency in the background.
	"""
	# Base ribbon representation.
	run( session, f"open {pdb_file} format mmcif" )

	# Set ChimeraX preset 1
	run( session, "preset 1" )

	# Hide the complete structure - ribbon and atom representations.
	run( session, "hide #1" )
	run( session, "hide #1 cartoon" )

	run( session, f"show #1/{chain1},{chain2} cartoon" )

	# Base colour for the chains.
	run( session, f"color #1/{chain1} #87CEFA" )
	run( session, f"color #1/{chain2} #F08080" )
	if pdb_id in ["7lna", "8cmk"]:
		model = 1
		run( session, "turn x 50" )
	else:
		model = 1.1
		# Reduce opacity of all chains.
		run( session, f"transparency #1/{chain1} 95 target c" )
		run( session, f"transparency #1/{chain2} 95 target c" )

		# Reduce opacity of all chains.
		run( session, f"transparency #{model}/{chain1} 0 target c" )
		run( session, f"transparency #{model}/{chain2} 0 target c" )

	if pdb_id == "2kqs":
		run( session, "turn y 50" )
		run( session, "turn z -100" )
		run( session, "turn x -20" )

	# Colour interface residues and increase opacity.
	p1_res = ",".join( predicted_interface1 )
	run( session, f"color #{model}/{chain1}:{p1_res} blue" )
	p2_res = ",".join( predicted_interface2 )
	run( session, f"color #{model}/{chain2}:{p2_res} red" )
	# Save image.
	file_name = os.path.join( out_path, f"{pdb_id}_{name}_bg_ensemble.png" )
	run( session, f"save ./{file_name} format png width 3600 height 2400" )
	run( session, "close")



######################################################

with open( "..datbase/Misc/Summary.json", "r" ) as f:
	entry_dict = json.load( f )
with open( "..datbase/Misc/Misc_dict_19.json", "r" ) as f:
	preds_dict = json.load( f )


def get_interface( predicted_interface: List,
					pdb_positions1: np.array,
					pdb_positions2: np.array ):
	p1_interface = predicted_interface[:200]
	p2_interface = predicted_interface[200:]
	
	contact_idx1 = np.where( np.array( p1_interface ) >= 0.5 )
	contact_idx2 = np.where( np.array( p2_interface ) >= 0.5 )

	predicted_interface1 = list( map( str, pdb_positions1[contact_idx1] ) )
	predicted_interface2 = list( map( str, pdb_positions2[contact_idx2] ) )
	return predicted_interface1, predicted_interface2


def get_pdb_positions( pdb_res1,pdb_res2 ):
	start, end = list( map( int, pdb_res1 ) )
	pdb_positions1 = np.arange( start, end+1, 1 )
	start, end = list( map( int, pdb_res2 ) )
	pdb_positions2 = np.arange( start, end+1, 1 )
	return pdb_positions1, pdb_positions2


af2_png_dir = "./png_files/AF2+Diso/"
target_dir = "./png_files/target/"
os.makedirs( "png_files/", exist_ok = True )
os.makedirs( af2_png_dir, exist_ok = True )
os.makedirs( target_dir, exist_ok = True )

for entry_id in entry_dict:
	pdb_id = entry_dict[entry_id]["pdb_id"]

	if pdb_id != "2kqs":
		continue
	_, chain1, chain2 = list( entry_dict[entry_id].keys() )
	pdb_res1 = entry_dict[entry_id][chain1]["pdb_res"]
	pdb_res2 = entry_dict[entry_id][chain2]["pdb_res"]
	pdb_positions1, pdb_positions2 = get_pdb_positions( pdb_res1, pdb_res2 )

	for model_name in ["AF2+Diso"]:
		( predicted_interface1,
			predicted_interface2 ) = get_interface( preds_dict[entry_id][model_name],
														pdb_positions1, pdb_positions2 )
		out_path = af2_png_dir

		ensemble_in_background( f"./pdb_files/{pdb_id}.cif",
								chain1, chain2,
								predicted_interface1,
								predicted_interface2,
								out_path,
								model_name )
		# break

	( predicted_interface1,
		predicted_interface2 ) = get_interface( preds_dict[entry_id]["Target"],
													pdb_positions1, pdb_positions2 )
	# Get the target interface image.
	ensemble_in_background( f"./pdb_files/{pdb_id}.cif",
							chain1, chain2,
							predicted_interface1,
							predicted_interface2,
							target_dir,
							"Target" )
