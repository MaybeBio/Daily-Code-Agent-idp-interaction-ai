"""
Create Disobind input files for the IDPPI dataset.
DOI: https://doi.org/10.1038/s41598-018-28815-x
"""
from typing import List, Tuple, Dict
import os, json
import numpy as np
import pandas as pd
from multiprocessing import Pool
import tqdm

from dataset.from_APIs_with_love import get_uniprot_seq


class IdppiInput():
	"""
	Evaluate Disobind on the IDPPI dataset using interface_1 model.
	"""
	def __init__( self ):
		self.base_dir = "../database/"
		self.idppi_output_dir = os.path.join( self.base_dir, "./idppi/" )

		self.idppi_file = os.path.join( self.base_dir, "input_files/41598_2018_28815_MOESM2_ESM.xlsx" )
		self.diso_uni_seq_file = os.path.join( self.base_dir, "v_19/Uniprot_seq.json" )

		self.cores = 100
		self.max_seq_len = 10000
		# If True, ignores protein pairs which are present in Disobind dataset.
		self.remove_diso_seq = True
		
		self.uniprot_seq_dict = {}
		self.logger = {}

		# .json file containing Uniprot sequences for IDPPI dataset.
		self.uniprot_seq_file = os.path.join( self.idppi_output_dir, "Uniprot_seq_idppi.json" )
		# .csv file containing IDPPI entry_ids for Disobind.
		self.diso_input_file = os.path.join( self.idppi_output_dir, "IDPPI_input_diso.csv" )
		self.idppi_target_file = os.path.join( self.idppi_output_dir, "IDPPI_target.json" )
		self.logs_file = os.path.join( self.idppi_output_dir, "Logs.txt" )


	def forward( self ):
		"""
		Parse the IDPPI .xslx file.
		Download the unique UniProt seq.
		Remove protein pairs for which UniProt seq was not obtained.
		"""
		self.create_dir()
		idppi_pairs, unique_uni_ids = self.parse_idppi_file()
		self.logger["idppi_pairs"] = len( idppi_pairs )
		self.logger["unique_uni_ids"] = len( unique_uni_ids )

		print( f"Total protein pairs obtained = {len( idppi_pairs )}" )
		print( f"Unique UniProt IDs obtained = {len( unique_uni_ids )}" )

		self.get_uniprot_seq_dict( unique_uni_ids )
		idppi_test_dict = self.filter_idppi_pairs( idppi_pairs )

		self.logger["selected_idppi_pairs"] = len( idppi_test_dict )
		print( f"Selected IDPPI protein pairs = {len( idppi_test_dict )}" )

		self.create_input_for_disobind( idppi_test_dict )
		self.write_logs()


	def create_dir( self ):
		"""
		Create the required directories.
		"""
		os.makedirs( self.idppi_output_dir, exist_ok = True )


	def parse_idppi_file( self ) -> Tuple[Dict, List]:
		"""
		Obtain all protein pairs and the corresponding labels from the 5 IDPPI Test sets.
		Obtain all entries in the Test sets (1-5)
			in sheets 17-21.
		"""
		with open( self.diso_uni_seq_file, "r" ) as f:
			diso_uni_seq_dict = json.load( f )

		idppi_pairs = {}
		unique_uni_ids = set()
		for sheet_name in ["Table S17-TestSet1", "Table S18-TestSet2 ",
							"Table S19-TestSet3 ", "Table S20-TestSet4 ",
							"Table S21-TestSet5"]:
			df = pd.read_excel( self.idppi_file, sheet_name = sheet_name, header = None )
			print( sheet_name, "  ", df.shape[0] )
			for row in df[0].str.split( " " ):
				pair_id = f"{row[1]}--{row[2]}"
				if self.remove_diso_seq:
					# Remove protein  pairs if they are present in Disobind dataset.
					if row[1] in diso_uni_seq_dict or row[2] in diso_uni_seq_dict:
						continue
				idppi_pairs[pair_id] = row[0]
				unique_uni_ids.update( row[1:] )
		
		return idppi_pairs, sorted( unique_uni_ids )


	def download_uniprot_seq( self, uni_id: str ) -> Tuple[str, str, bool]:
		"""
		Get the UniProt seq given the UniProt ID.
		"""
		uni_seq = get_uniprot_seq( uni_id, max_trials = 5, wait_time = 5, return_id = False )
		if len( uni_seq ) == 0:
			success = False
		else:
			success = True
		return uni_id, uni_seq, success


	def get_uniprot_seq_dict( self, unique_uni_ids: List ):
		"""
		Download all UniProt IDs.
		"""
		if os.path.exists( self.uniprot_seq_file ):
			print( "Loading pre-downloaded UniProt seq..." )
			with open( self.uniprot_seq_file, "r" ) as f:
				self.uniprot_seq_dict = json.load( f )
		else:
			with Pool( self.cores ) as p:
				for result in tqdm.tqdm( p.imap_unordered( self.download_uniprot_seq, unique_uni_ids ),
															total = len( unique_uni_ids ) ):

					uni_id, uni_seq, success = result
					if success:
						self.uniprot_seq_dict[uni_id] = uni_seq

			with open( self.uniprot_seq_file, "w" ) as w:
				json.dump( self.uniprot_seq_dict, w )

		self.logger["uni_seq_dwnld"] = len( self.uniprot_seq_dict )
		print( f"UniProt sequences obtained = {len( self.uniprot_seq_dict )}" )


	def filter_idppi_pairs( self, idppi_pairs: Dict ) -> Dict:
		"""
		Ignore protein pairs for which:
			UniProt seq could not be downloaded.
			UniProt seq len > max_seq_len.
		"""
		idppi_test_dict = {}
		for pair_id in idppi_pairs:
			uni_id1, uni_id2 = pair_id.split( "--" )
			if uni_id1 not in self.uniprot_seq_dict or uni_id2 not in self.uniprot_seq_dict:
				continue
			elif ( 
					( len( self.uniprot_seq_dict[uni_id1] ) > self.max_seq_len ) or 
					( len( self.uniprot_seq_dict[uni_id1] ) > self.max_seq_len )
					):
				continue
			else:
				idppi_test_dict[pair_id] = idppi_pairs[pair_id]
		return idppi_test_dict


	def create_input_for_disobind( self, idppi_test_dict: Dict ):
		"""
		For all selected protein pairs from IDPPI test set,
			create entry_id for input to Disobind.
		"""
		disobind_input_pairs = []
		for pair_id in idppi_test_dict:
			uni_id1, uni_id2 = pair_id.split( "--" )

			seq_len1 = len( self.uniprot_seq_dict[uni_id1] )
			seq_len2 = len( self.uniprot_seq_dict[uni_id2] )

			copy_num = 0
			entry_id = f"{uni_id1}:{1}:{seq_len1}--{uni_id2}:{1}:{seq_len2}_{copy_num}"

			disobind_input_pairs.append( entry_id )

		print( len( disobind_input_pairs ) )
		self.logger["total_idpi_entry_ids"] = len( disobind_input_pairs )
		with open( self.diso_input_file, "w" ) as w:
			w.writelines( ",".join( disobind_input_pairs ) )
		with open( self.idppi_target_file, "w" ) as w:
			json.dump( idppi_test_dict, w )

		target_labels = list( idppi_test_dict.values() )
		total = len( target_labels )
		pos_pairs = np.count_nonzero( np.array( target_labels ).astype( int ) )
		self.logger["interacting_pairs_count"] = pos_pairs
		self.logger["non_interacting_pairs_count"] = total-pos_pairs


	def write_logs( self ):
		"""
		Write logs to a .txt file.
		"""
		with open( self.logs_file, "w" ) as w:
			w.writelines( "------------------- IDPPI Logs -------------------\n" )
			w.writelines( "Configs -----\n" )
			w.writelines( f"Cores: {self.cores}\n" )
			w.writelines( f"Max sequence len: {self.max_seq_len}\n" )
			w.writelines( f"Redundancy reduce with v_21 UniProt seq: {self.remove_diso_seq}\n" )
			w.writelines( "\nStats -----\n" )
			w.writelines( f"Total IDPPI protein pairs = {self.logger['idppi_pairs']}\n" )
			w.writelines( f"Unique IDPI UniProt IDs obtained = {self.logger['unique_uni_ids']}\n" )
			w.writelines( f"Total UniProt seq obtained = {self.logger['uni_seq_dwnld']}\n" )
			w.writelines( f"Selected IDPPI pairs = {self.logger['selected_idppi_pairs']}\n" )
			w.writelines( f"Total IDPPI entry_ids obtained = {self.logger['total_idpi_entry_ids']}\n" )
			w.writelines( f"\tTotal Interacting pairs = {self.logger['interacting_pairs_count']}\n" )
			w.writelines( f"\tTotal non-Interacting pairs = {self.logger['non_interacting_pairs_count']}\n" )


if __name__ == "__main__":
	IdppiInput().forward()

