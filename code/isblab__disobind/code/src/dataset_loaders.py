##### Preparing dataset and creating batchs for model training #####
########## ------>"May the Force serve u well..." <------###########
####################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import os

import torch
from torch import nn
from torch.utils.data.dataset import random_split


def create_residue_pairs( prot1_emb, prot2_emb ):
	"""
	Create all vs all residue pairs.
	Input must be residue wise embeddings for prot1 and prot2.
	"""
	prot1 = np.array( prot1_emb )
	prot2 = np.array( prot2_emb )
	
	num_res1, num_res2 = prot1.shape[0], prot2.shape[0]

	prot1 = np.repeat( prot1, num_res2, axis = 0 )
	prot2 = np.tile( prot2, ( num_res1, 1 ) )

	return prot1, prot2, num_res1, num_res2


def create_dataset( prot1_input, prot2_input, output = None, dtype = 16 ):
	prot1_input = torch.from_numpy( np.array( prot1_input, dtype = np.float32 ) ).squeeze( 1 )
	prot2_input = torch.from_numpy( np.array( prot2_input, dtype = np.float32 ) ).squeeze( 1 )

	dataset = [prot1_input, prot2_input]

	if output is not None:
		dataset.append( torch.from_numpy( np.array( output, dtype = np.float32 ) ) )

	dataset = torch.cat( dataset, axis = 1 )
	return dataset


def split_dataset( dataset, partitions ):
	"""
	Randomly split the dataset into a Train:Dev:Test set.
	"""
	print( "Spliting the dataset..." )
	train_set, dev_set, test_set = random_split( dataset, partitions )

	return train_set, dev_set, test_set



def dataloader( dataset, batch_size, shuffle = False ):
	print("Loading dataset")

	dataset = torch.utils.data.DataLoader( dataset,
										batch_size = batch_size,
										shuffle = shuffle
										)

	return dataset



#########################################################################
#########################################################################
class DatasetLoader ( nn.Module ):
	def __init__ ( self, config, seed_worker, seed ):
		super(DatasetLoader, self).__init__()
		self.config = config


	def load_dataset( self ):
		dev_set = np.load( os.path.join( self.config.input_files, self.config.dev_file ) )
		print( "Loaded Dev set..." )

		test_set = np.load( os.path.join( self.config.input_files, self.config.test_file ) )
		print( "Loaded Test set..." )
		
		train_set = np.load( os.path.join( self.config.input_files, self.config.train_file ) )
		print( "Loaded Train set..." )

		print( "Train\t\t Dev\t\t Test \n", train_set.shape, "   ", dev_set.shape, "  ", test_set.shape )
		dev_set = torch.from_numpy( dev_set )
		test_set = torch.from_numpy( test_set )
		train_set = torch.from_numpy( train_set )

		return train_set, dev_set, test_set


	def prepare_dataset( self ):
		dataset = self.load_dataset()

		partitions = [self.config.train_set_size,
					self.config.dev_set_size,
					self.config.test_set_size]

		train_set, dev_set, test_set = split_dataset( dataset, partitions )

		return train_set, dev_set, test_set



	def create_dataloaders( self, train_set, dev_set, test_set ):
		"""
		Create dataloader objects for the train, dev, test sets.
		"""
		train_set = dataloader( train_set, self.config.batch_size, self.config.batch_shuffle[0] )
		dev_set = dataloader( dev_set, self.config.batch_size, self.config.batch_shuffle[1] )
		test_set = dataloader( test_set, self.config.batch_size, self.config.batch_shuffle[2] )

		return train_set, dev_set, test_set

