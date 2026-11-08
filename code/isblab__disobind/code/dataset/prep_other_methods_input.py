"""
Create input files for running AIUPred, DeepDISOBind, and MORFchibi.
"""
############ ------>"May the Force serve u well..." <------##############
#########################################################################
import math
import json, os, subprocess, json
from typing import List, Tuple, Dict
import numpy as np
import os, json


class CreateInput():
    """
    Get input files for using AIUPred, DeepDISOBind, and MORFchibi on the OOD set.
    """
    def __init__(self):
        self.base_dir = "../database/"
        self.data_version = 19
        self.meth_dir = os.path.join( self.base_dir, "other_methods/" )

        self.ood_csv_file = os.path.join( self.base_dir, f"v_{self.data_version}/prot_1-2_test_v_{self.data_version}.csv" )
        self.uni_seq_file = os.path.join( self.base_dir, f"v_{self.data_version}/Uniprot_seq.json" )
        self.deepdiso_fasta_file = os.path.join( self.meth_dir, "deepdisobind_fasta" )
        self.aiupred_input_file = os.path.join( self.meth_dir, "aiupred_input.json" )
        self.morfchibi_fasta_file = os.path.join( self.meth_dir, "morfchibi_fasta.fasta" )
        # self.entry_ids_map_file = os.path.join( self.meth_dir, "entry_id_map.json" )

        # Dict to store all protein sequences in OOD set.
        self.ood_seq_dict = {}
        # # Dict to mpa old entry_ids to new entry_ids.
        # self.old_to_new = {}


    def forward( self ):
        with open( self.uni_seq_file, "r" ) as f:
            self.uniprot_seq = json.load( f )
        self.create_dir()
        entry_ids = self.get_entry_ids()
        self.get_seq_for_ood_entry( entry_ids )
        self.write_deepdiso_fasta_file()
        self.write_aiupred_input_file()
        self.write_morfchibi_fasta_file()

        # with open( self.entry_ids_map_file, "w" ) as w:
        #     json.dump( self.old_to_new, w )


    def create_dir( self ):
        """
        Create the required directories if they don't already exist.
        """
        os.makedirs( self.meth_dir, exist_ok = True )


    def split_entry_id( self, entry_id: str ):
        """
        entry_id --> "{uni_id1}:{start1}:{end1}--{uni_id2}:{start2}:{end2}_{copy_num}"
        Split and return uni_id1, start1, end1, uni_id2, start2, end2, copy_num
        """
        uni_id1, uni_id2 = entry_id.split( "--" )
        uni_id2, copy_num = uni_id2.split( "_" )
        uni_id1, start1, end1 = uni_id1.split( ":" )
        uni_id2, start2, end2 = uni_id2.split( ":" )

        return uni_id1, start1, end1, uni_id2, start2, end2, copy_num


    def get_entry_ids( self ):
        """
        Get all OOD set entry_id's.
        """
        with open( self.ood_csv_file, "r" ) as f:
            entry_ids = f.readlines()[0].split( "," )
            if entry_ids[-1] == "":
                entry_ids = entry_ids[:-1]
        return entry_ids


    def select_uniprot_seq( self, uni_id:str, start1: str, end: str ):
        """
        Select the required Uniprot seq from the full sequence.
        """
        seq = self.uniprot_seq[uni_id][int( start1 )-1:int( end )]
        return seq


    def get_seq_for_ood_entry( self, entry_ids: List ):
        """
        Both AIUPred and DeepDisoBind provide partner-independent,
            interface residue predictions.
        For all OOD entries, both proteins will be considered as
            separate inputs for the above two methods.
        """
        for entry_id in entry_ids:
            uni_id1, start1, end1, uni_id2, start2, end2, copy_num = self.split_entry_id( entry_id )

            seq1 = self.select_uniprot_seq( uni_id1, start1, end1 )
            seq2 = self.select_uniprot_seq( uni_id2, start2, end2 )
            new_p1_id = f"{uni_id1}_{start1}_{end1}"
            new_p2_id = f"{uni_id2}_{start2}_{end2}"

            self.ood_seq_dict[new_p1_id] = {"ood_entry_id": f"{uni_id1}--{uni_id2}_{copy_num}",
                                            "seq": seq1}
            self.ood_seq_dict[new_p2_id] = {"ood_entry_id": f"{uni_id1}--{uni_id2}_{copy_num}",
                                            "seq": seq2}


    def write_deepdiso_fasta_file( self ):
        """
        Create a FASTA file for input to DeepDisoBind.
        DepDisoBind server accepts only 20 seq per job.
            So create batches of 20 seq each.
        """
        all_ids = list( self.ood_seq_dict.keys() )
        for s in np.arange( 0, len( all_ids ), 20 ):
            e = s+20 if s+20 < len( all_ids ) else len( all_ids )
            with open( f"{self.deepdiso_fasta_file}_{s}-{e}.fasta", "w" ) as w:
                for id_ in all_ids[s:e]:
                    seq = self.ood_seq_dict[id_]["seq"]
                    w.writelines( f">{id_}\n" )
                    w.writelines( f"{seq}\n\n" )


    def write_aiupred_input_file( self ):
        """
        AIUPred will be used from python, so just saving the ood_seq_dict as is.
        """
        with open( self.aiupred_input_file, "w" ) as w:
            json.dump( self.ood_seq_dict, w )


    def write_morfchibi_fasta_file( self ):
        """
        Create a FASTA file for input to MORFchibi web.
        """
        # all_ids = list( self.ood_seq_dict.keys() )
        # for s in np.arange( 0, len( all_ids ), 20 ):
        #     e = s+20 if s+20 < len( all_ids ) else len( all_ids )
        with open( self.morfchibi_fasta_file, "w" ) as w:
            for id_ in self.ood_seq_dict:
                seq = self.ood_seq_dict[id_]["seq"]
                w.writelines( f">{id_}\n" )
                w.writelines( f"{seq}\n\n" )



################################################
if __name__ == "__main__":
    CreateInput().forward()

