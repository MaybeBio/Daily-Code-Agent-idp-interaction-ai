"""
Assuming predictions from AF2 and AF3 already exist, obtain contact maps
	and interface residues across all CG (1, 5, 10) for evaluation with Disobind.
Select confidently predicted interactions based on pLDDT and PAE.

Used for OOD test and Misc dataset.
"""
import os, json, pickle as pkl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
from torch import nn
import Bio
from Bio.PDB import PDBParser, MMCIFParser
import h5py

from src.metrics import torch_metrics
from dataset.utility import get_contact_map

import warnings
warnings.filterwarnings("ignore")


MAX_LEN_DICT = {"ood": {19: 100}, "misc": {19: 200}}

class AF2MPredictions():
	def __init__( self ):
		# Dataset version.
		self.version = 19
		# path for the dataset dir.
		self.base_path = f"../database/v_{self.version}/"

		# Distance threshold for creating the binary contact map.
		self.dist_threshold = 8
		# pLDDt cutoff for confident predictions.
		self.plddt_threshold = 70
		# PAE cutoff for confident predictions.
		self.pae_threshold = 5
		# ipTM cutoff for confident predictions.
		self.iptm_cutoff = 0.0
		# Datset type.
		self.mode = "ood" # "ood" or "misc"
		# Max prot1/2 lengths.
		self.max_len = MAX_LEN_DICT[self.mode][self.version]
		self.pad = True
		# If True, get contact map from struct, else from contact_probs.
		self.use_af3_struct = True
		# Resolution of coarse graining.
		self.coarse_grain_sizes = [1, 5, 10]
		self.device = "cuda"
		# AF2 or AF3.
		self.af_model = "AF2"

		if self.mode == "ood":
			# Dir containing AF2/3 predictions.
			self.af_test_pred_dir = f"{self.af_model}_preds"
			# OOD set target contact maps.
			self.cmap = h5py.File( f"{self.base_path}Target_bcmap_test_v_{self.version}.h5", "r" )		
			# File to store AF2/3 derived contact maps for OOD set.
			self.af_pred_dir = f"./AF_preds_v{self.version}/"

		elif self.mode == "misc":
			# Dir containing AF2/3 predictions.
			self.af_test_pred_dir = f"{self.af_model}_misc_preds"
			# OOD set target contact maps.
			self.cmap = h5py.File( f"../database/Misc/misc_test_target.h5", "r" )
			# File to store AF2/3 derived contact maps for OOD set.
			self.af_pred_dir = f"./AF_misc_preds_v{self.version}/"

		else:
			raise ValueError( "Incorrect mode specified (ood/peds supported)..." )

		# Dict to store prediction results.
		self.predictions = {}


	def forward( self ):
		if not os.path.exists( self.af_pred_dir ):
			os.makedirs( self.af_pred_dir )

		if self.af_model == "AF2":
			self.get_af2m_prediction()
			np.save( f"{self.af_pred_dir}Predictions_af2m_results_{self.iptm_cutoff}.npy", self.predictions )
		
		elif self.af_model == "AF3":
			self.get_af3_prediction()
			np.save( f"{self.af_pred_dir}Predictions_af3_results_{self.iptm_cutoff}.npy", self.predictions )

		else:
			raise Exception( f"Incorrect AF model - {self.af_model} specified. (AF2/Af3 supported)..." )



	def get_best_model( self, model_path, header ):
		"""
		For using the relaxed models, AF2 provides the iptm+ptm values 
				for each and hence their ranking.
		This can be obtained from ranking_debug.json file.
		If using ranked model, just take the ranked_0.pdb.

		Input:
		----------
		model_path --> path to the AF2/3 predicted structure.
		header --> (str) OOD entry identifier.

		Returns:
		----------
		best_model --> (str) file for the best model.
		score --> list of AF2/3 confidence metrics - [ipTM, pTM, ranking confidence/score]
		"""
		if self.af_model == "AF2":
			with open( f"{model_path}ranking_debug.json", "r" ) as f:
				data = json.load( f )
			best_model = data["order"][0]

			score = [data["iptm+ptm"][best_model]]
			with open( f"{model_path}/result_{best_model}.pkl", "rb" ) as f:
				result = pkl.load( f )
			score = [float( result["iptm"] ), float( result["ptm"] ), result["ranking_confidence"]]

		elif self.af_model == "AF3":
			with open( f"{model_path}fold_{header}_summary_confidences_0.json", "r" ) as f:
				data = json.load( f )
			best_model = "0"
			score = [data["iptm"], data["ptm"], data["ranking_score"]]

		else:
			raise Exception( f"Incorrect AF model - {self.af_model} specified. (AF2/Af3 supported)..." )

		return best_model, score


	def get_PAE_matrix( self, data_file_path ):
		"""
		PAE provides inter-domain confidence in the predicted structure.

		Input:
		----------
		path --> path for the AF2/3 PAE file.

		Returns:
		----------		
		pae --> (np.array) PAE matrix.
		"""
		if self.af_model == "AF2":
			with open( data_file_path, "rb" ) as f:
				data = pkl.load( f )
			pae = data["predicted_aligned_error"]
		
		elif self.af_model == "AF3":
			with open( data_file_path, "rb" ) as f:
				data = json.load( f )
			pae = np.array( data["pae"] )

		else:
			raise Exception( f"Incorrect AF model - {self.af_model} specified. (AF2/Af3 supported)..." )

		return pae


	def get_coordinates( self, model_file ):
		"""
		Get the Calpha coordinates for each residue for both chains from the AF2M predicted structure.

		Input:
		----------
		model_file --> AF2/3 predicted structure file.

		Returns:
		----------
		coords_dict --> dict containing coordinates for all chains.
		plddt_dict --> dict containing Ca-plddt for all residues in all chains.
		"""
		coords_dict = {}
		plddt_dict = {}
		if self.af_model == "AF2":
			models = PDBParser().get_structure( "pdb", model_file )
		else:
			models = MMCIFParser().get_structure( "cif", model_file )

		for model in models:
			for chain in model:
				coords_dict[chain.id[0]] = []
				plddt_dict[chain.id[0]] = []
				for residue in chain:
					# Take only the ATOM entries and ignore the 
					# 		HETATM entries (which contain "w" instead of " ").
					if residue.id[0] == " ":
						coords_dict[chain.id[0]].append( residue['CA'].coord )
						plddt_dict[chain.id[0]].append( residue["CA"].get_bfactor() )
				coords_dict[chain.id[0]] = np.array( coords_dict[chain.id[0]] )

		return coords_dict, plddt_dict


	def create_contact_map( self, coords_dict ):
		"""
		Obtain the inter-residue distance map from the AF2/3 predicted structure.
		Create a binary contact map if residues are within threshold distance.

		Input:
		----------
		coords_dict --> dict containing coords for chains A and B.

		Returns:
		----------
		contact_map --> binary contact map.
		"""
		chainA, chainB = coords_dict.keys()
		contact_map = get_contact_map( coords_dict[chainA], coords_dict[chainB], self.dist_threshold )

		return contact_map


	def get_af3_contact_map_from_probs( self, data_file_path, uni_pos1, uni_pos2 ):
		"""
		Obtain contact probabilities from the AF3 JSON file.
		"""
		start1, end1 = list( map( int, uni_pos1 ) )
		start2, end2 = list( map( int, uni_pos2 ) )

		p1_len = end1 - start1 + 1
		p2_len = end2 - start2 + 1

		with open( data_file_path, "rb" ) as f:
			data = json.load( f )

		contact_probs = np.array( data["contact_probs"] )

		ur = contact_probs[:p1_len, p1_len:]
		ll = contact_probs[p1_len:, :p1_len]
		avg_probs = ( ur + ll.T )/2

		contact_map = np.where( avg_probs >= 0.5, 1, 0 )
		contact_map = contact_map.reshape( p1_len, p2_len )
		return contact_map



	def get_interface_from_cmap( self, pred ):
		"""
		Given a contact map, obtain the interface for prot1 and prot2.
			i.e. Obtain interacting residues for both proteins.

		Input:
		----------
		pred --> AF2/3 predicted contact map.

		Returns:
		----------
		interface --> (np.array) containing interfaces for prot1/2.
		"""
		idx1, idx2 = np.where( pred == 1.0 )
		idx1 = np.unique( idx1 )
		idx2 = np.unique( idx2 )
		p1 = np.zeros( ( pred.shape[0], 1 ) )
		p1[idx1,:] = 1
		p2 = np.zeros( ( pred.shape[1], 1 ) )
		p2[idx2,:] = 1
		interface = np.concatenate( ( p1, p2 ), axis = 0 )

		return interface


	def coarse_grain( self, pred, key, plddt_mat, pae_mat ):
		"""
		Given the AF2/3 predicted contact map, obtain interaction/interface
			predictions for all CG values specified.

		Input:
		----------
		pred --> AF2/3 predicted contact map.
		key --> OOD set entry identifier.
		plddt_mat --> binary mask indicating residues with high pLDDT.
		pae_mat --> binary mask indicating residues with high PAE.

		Returns:
		----------
		None
		"""
		if self.pad:
			padded_pred = np.zeros( ( self.max_len, self.max_len ) )
			h, w = pred.shape
			padded_pred[:h,:w] = pred
			pred = padded_pred
		
		plddt_pae_pred = pred*plddt_mat*pae_mat

		plddt_pae_pred = torch.from_numpy( plddt_pae_pred ).to( self.device ).float()
		plddt_pae_pred = plddt_pae_pred.unsqueeze( 0 )
		plddt_pae_pred = plddt_pae_pred.unsqueeze( 1 )

		for cg in self.coarse_grain_sizes:
			if cg != 1:
				m = nn.MaxPool2d( kernel_size = cg, stride = cg )
				plddt_pae_pred_ = m( plddt_pae_pred )

			else:
				plddt_pae_pred_ = plddt_pae_pred

			plddt_pae_pred_ = torch.squeeze( plddt_pae_pred_, ( 0, 1 ) ).detach().cpu().numpy()
			plddt_pae_interaction = plddt_pae_pred_
			plddt_pae_interface = self.get_interface_from_cmap( plddt_pae_pred_ )

			self.predictions[key][f"interaction_{cg}"] = {
														f"{self.af_model}_pLDDT_PAE": plddt_pae_interaction
														}
			self.predictions[key][f"interface_{cg}"] = {
														f"{self.af_model}_pLDDT_PAE": plddt_pae_interface
														}


	def get_plddt_pae_mat( self, plddt_dict, pae ):
		"""
		Create a binary matrix based on:
			plDDT values of prot1/2 for AF2 structures.
				plDDT >= 70 --> 1 else 0
			PAE matrix
				plDDT <= 5 --> 1 else 0

		Input:
		----------
		plddt_dict --> dict containing Ca-plddt for all residues in all chains.
		pae --> AF2/3 predicted PAE matrix.

		Returns:
		----------
		plddt_mat --> binary mask indicating residues with high pLDDT.
		pae_mat --> binary mask indicating residues with high PAE.
		"""
		chain1, chain2 = list( plddt_dict.keys() )
		plddt1 = np.array( plddt_dict[chain1] ).reshape( -1, 1 )
		plddt2 = np.array( plddt_dict[chain2] ).reshape( 1, -1 )
		plddt1 = np.where( plddt1 >= self.plddt_threshold, 1, 0 )
		plddt2 = np.where( plddt2 >= self.plddt_threshold, 1, 0 )
		plddt_mat = plddt1*plddt2

		m, n = plddt_mat.shape
		padded_mat = np.zeros( ( self.max_len, self.max_len ) )

		padded_mat[:m, :n] = plddt_mat
		plddt_mat = padded_mat

		# Taking the upper-right and lower-left quadrants.
		pae_ur = pae[:m, m:]
		pae_ll = pae[m:, :m]
		pae = ( pae_ur + pae_ll.T )/2

		if m != pae.shape[0] or n != pae.shape[1]:
			raise Exception( "Incorrect PAE matrix dimensions..." )

		m, n = pae.shape
		pae_mat = np.zeros( ( self.max_len, self.max_len ) )
		pae_mat[:m, :n] = pae

		# Create a binary PAE matrix to highlight confident predictions.
		pae_mat = np.where( pae_mat <= self.pae_threshold, 1, 0 )

		return plddt_mat, pae_mat


	def get_af2m_prediction( self ):
		"""
		Obtain the contact map from predicted AF2M structures.
		Output a dict containing -
			Interactions and Interface predictions for all CG values.
			Best model.
			pLDDT for both chains.
			PAE matrix.
			Scores for best model.
		
		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		if not os.path.exists( self.af_pred_dir ):
			os.makedirs( self.af_pred_dir )
		
		for i, key in enumerate( self.cmap.keys() ):
			if key == "P0DTD1--P0DTD1_1":
				continue
			self.predictions[key] = {}

			head1, head2 = key.split( "--" )
			head2, num = head2.split( "_" )
			uni_id1, _, _ = head1.split( ":" )
			uni_id2, _, _ = head2.split( ":" )
			header = f"{uni_id1}--{uni_id2}_{num}"

			# Ignoring this entry, as AF2-multimer crashed for this.
			if header == "P0DTD1--P0DTD1_1":
				continue

			# AF2 fasta file names consist of: "UniID1--UniID2_num"
			path = f"{self.base_path}{self.af_test_pred_dir}/{header}/"

			best_model, score = self.get_best_model( path, header )
			print( f"Entry {i} Best model --> {key} \t\t {best_model} \t\t {score}" )
			
			pae = self.get_PAE_matrix( f"{path}/result_{best_model}.pkl" )
			
			
			model_file = f"{path}unrelaxed_{best_model}.pdb"
			coords_dict, plddt_dict = self.get_coordinates( model_file )
			contact_map = self.create_contact_map( coords_dict )

			if score[0] <= self.iptm_cutoff:
				contact_map = np.zeros( ( contact_map.shape ) )
			print( score[0], "  ", np.count_nonzero( contact_map ), "\n" )

			self.predictions[key]["best_model"] = best_model
			self.predictions[key]["plddt"] = plddt_dict
			self.predictions[key]["pae"] = pae
			self.predictions[key]["scores"] = score

			plddt_mat, pae_mat = self.get_plddt_pae_mat( plddt_dict, pae )

			self.coarse_grain( contact_map, key, plddt_mat, pae_mat )



	def get_af3_prediction( self ):
		"""
		Obtain the contact map from predicted AF3 structures.
		Output a dict containing -
			Interactions and Interface predictions for all CG values.
			Best model.
			pLDDT for both chains.
			PAE matrix.
			Scores for best model.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		if not os.path.exists( self.af_pred_dir ):
			os.makedirs( self.af_pred_dir )
		
		for i, key in enumerate( self.cmap.keys() ):
			self.predictions[key] = {}

			head1, head2 = key.split( "--" )
			head2, num = head2.split( "_" )
			uni_id1, start1, end1 = head1.split( ":" )
			uni_id2, start2, end2 = head2.split( ":" )
			header = f"{uni_id1}_{uni_id2}_{num}".lower()

			# AF3 fasta file names consist of: "UniID1--UniID2_num"
			if os.path.exists( f"{self.base_path}{self.af_test_pred_dir}/fold_{header}/" ):
				path = f"{self.base_path}{self.af_test_pred_dir}/fold_{header}/"
			elif os.path.exists( f"{self.base_path}{self.af_test_pred_dir}/{header}/" ):
				path = f"{self.base_path}{self.af_test_pred_dir}/{header}/"
			else:
				raise Exception( f"AF3 output directory not found for {header}..." )

			best_model, score = self.get_best_model( path, header )
			print( f"Entry {i} Best model --> {key} \t\t {best_model} \t\t {score}" )
			
			pae = self.get_PAE_matrix( f"{path}/fold_{header}_full_data_0.json" )
			
			model_file = f"{path}fold_{header}_model_0.cif"
			coords_dict, plddt_dict = self.get_coordinates( model_file )
			
			if self.use_af3_struct:
				contact_map = self.create_contact_map( coords_dict )
			else:
				contact_map = self.get_af3_contact_map_from_probs(
															f"{path}/fold_{header}_full_data_0.json",
															( start1, end1 ), ( start2, end2 )
															)

			if score[0] <= self.iptm_cutoff:
				contact_map = np.zeros( ( contact_map.shape ) )

			self.predictions[key]["best_model"] = best_model
			self.predictions[key]["plddt"] = plddt_dict
			self.predictions[key]["pae"] = pae
			self.predictions[key]["scores"] = score

			plddt_mat, pae_mat = self.get_plddt_pae_mat( plddt_dict, pae )

			self.coarse_grain( contact_map, key, plddt_mat, pae_mat )


if __name__ == "__main__":
	AF2MPredictions().forward()
	print( "May the Force be with you..." )
