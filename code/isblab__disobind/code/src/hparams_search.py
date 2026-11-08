####################### A Quest for Hyperparameters #######################
############# ------>"May the Force serve u well..." <------###############
###########################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import argparse
from omegaconf import OmegaConf
import random
import os
import time
import joblib

# import optuna
# from optuna.trial import TrialState

import torch
from torch import nn

from models.get_model import get_model
from dataset_loaders import DatasetLoader
from build_model import Trainer
from utils import ( dump_metrics, create_plots )

import warnings
warnings.filterwarnings( "ignore" )

"""
Following script executes the model training using configs from the config_file.
Implements hyper-parameter tuning process in multiple ways:
	manual --> specify hparams to be used in a grid search like fashion. -- Preferred
	optuna --> use optuna for tuning hparams (deprecated).
"""


##Load Training details##############################
#####################################################
parser = argparse.ArgumentParser(description="Get the Version details")
parser.add_argument('--Version_file', '-f', dest="f", 
					help="*.yml file containing the parameters and hyperparameters for the model", 
					required = True )
parser.add_argument('--mode', '-m', dest="m", 
					help="manual/optuna -- optuna has been deprecated.", 
					required = True )
version = parser.parse_args().f

with open( version, 'r' ) as f:
    config_file = OmegaConf.load( f )


class HparamSearch( nn.Module ):
	def __init__( self ):
		super( HparamSearch, self ).__init__()
		self.global_seed = config_file.Global_seed
		# self.g = torch.Generator().manual_seed( self.global_seed )
		self.train_set, self.dev_set, self.train_set = 0, 0, 0
		self.num_metrics = config_file.conf.train_params.num_metrics[0]
		
		# For hyperparameter optimization using Optuna.
		# self.activation_func1 = config_file.conf.model_params.activation1
		# self.projection_layer = config_file.conf.model_params.projection_layer
		# self.num_hid_layers = config_file.conf.model_params.num_hid_layers
		# self.log_weight = config_file.conf.train_params.log_weight
		# self.weight_decay = config_file.conf.train_params.weight_decay
		# self.lr = config_file.conf.train_params.learning_rate
		# self.contact_threshold = config_file.conf.train_params.contact_threshold


	def seed_worker( self, seed ):
		# Set the seeds for PRNG.
		torch.manual_seed( seed )
		# torch.cuda.manual_seed( worker_seed )
		torch.cuda.manual_seed_all( seed )
		np.random.seed( seed )
		random.seed( seed )


	def call_of_duty ( self ):
		self.seed_worker( self.global_seed )
		dataset_config = config_file.conf.dataset

		_dataset = DatasetLoader( config_file.conf.dataset, self.seed_worker, self.global_seed )

		train_set, dev_set, test_set = _dataset.load_dataset()
		self.train_set, self.dev_set, self.test_set = _dataset.create_dataloaders( train_set, dev_set, test_set )
		
		print( "No. of batches created -- Train:Dev:Test = {}: {}: {}\n".format( len( self.train_set ),
																				len( self.dev_set ),
																				len( self.test_set ) ) )


	# ## Using Optuna#######################################
	# #####################################################
	# def objective( self, trial ):
	# 	# Generate the parameters to be optimized.
	# 	# optimizer_name = trial.suggest_categorical("optimizer", ["AdamW", "RMSprop", "SGD"])
	# 	model_config = config_file.conf.model_params
	# 	train_config = config_file.conf.train_params

	# 	model_config.activation1 = trial.suggest_categorical( "activation1", self.activation_func1 )
	# 	model_config.projection_layer = trial.suggest_categorical( "projection_layer", self.projection_layer )
	# 	model_config.op_od = trial.suggest_categorical( "op_od", self.op_od )
	# 	model_config.aggregate = trial.suggest_categorical( "aggregate", self.aggregate )
	# 	model_config.drop1_first = trial.suggest_categorical("drop1_first", self.drop1_first )
	# 	model_config.dropout1_prob = trial.suggest_categorical("dropout1_prob", self.dropout1_prob )
	# 	model_config.dropout2_prob = trial.suggest_categorical("dropout2_prob", self.dropout2_prob )
	# 	train_config.contact_threshold = trial.suggest_categorical( "threshold", self.contact_threshold )
	# 	train_config.weight_decay = trial.suggest_float( "weight_decay", self.weight_decay[0], self.weight_decay[1], log = True )
	# 	train_config.learning_rate = trial.suggest_float( "lr", self.lr[0], self.lr[1], log = True )

	#     # Initialize the model.
	# 	device = model_config.device
	# 	model = get_model( model_config )
	# 	# Using multiple GPUs
	# 	# model = nn.DataParallel( model, device_ids = [0,1] )
	# 	model.to( device )

	# 	# Initialize the trainer.
	# 	trainer = Trainer(  train_config, device )

	# 	# Training the model.
	# 	for epoch in range( train_config.max_epochs ):
	# 		print( "\nEpoch ---------->", epoch )
	# 		# Get F1 score for dev set.
	# 		dev_f1 = trainer.optuna_mode( model, self.train_set, self.dev_set )

	# 		trial.report( dev_f1, epoch )

	# 		# Handle pruning based on the intermediate value.
	# 		if trial.should_prune():
	# 			raise optuna.exceptions.TrialPruned()

	# 	return dev_f1


	# def optuna_search ( self, config_file ):
	# 	# Using optuna for hparam search.
	# 	trials = config_file.conf.train_params.optuna_trials

	# 	study = optuna.create_study(
	# 	    direction = "maximize",
	# 	    study_name = "Optimize",
	# 	    pruner=optuna.pruners.MedianPruner(),
	# 	)
	# 	study.optimize( self.objective, n_trials = trials )  #, timeout=1000
	# 	pruned_trials = study.get_trials( deepcopy = False, states = [TrialState.PRUNED] )
	# 	complete_trials = study.get_trials( deepcopy = False, states = [TrialState.COMPLETE] )

	# 	print( "\n -----------> Best params: ", study.best_params )
	# 	print( "\n -----------> Best trial: ", study.best_trial )
		
	# 	with open( "Optimization.txt", "w" ) as w:
	# 		w.writelines( "<<-------------Summary------------->>\n" )
	# 		w.writelines( "Best params: " + str(study.best_params) + "\n" )
	# 		w.writelines( "No. of finished trials = " + str( len( study.trials ) ) + "\n" )
	# 		w.writelines( "No. of pruned trials = " + str( len( pruned_trials ) ) + "\n" )
	# 		w.writelines( "No. of complete trials = " + str( len( complete_trials ) ) + "\n" )

	# 		trial = study.best_trial
	# 		w.writelines("\n\nBest trial:" )
	# 		w.writelines("\n  Value: " + str( trial.value ) )
	# 		w.writelines("\n  Params: ")
	# 		for key, value in trial.params.items():
	# 			w.writelines("    {}: {}".format( str( key ), str( value ) ) )

	# 	print( "Alas, the journey comes to an end..." )



	## Manual hparam search###############################
	######################################################
	def executor( self, model_config, train_config, file_name ):
		# Execute the training process
		model = get_model( model_config )
		device = model_config.device
		# Using multiple GPUs
		# model = nn.DataParallel( model, device_ids = [0,1] )
		model.to( device )
		print( "Using learning_rate = ", train_config.learning_rate )
		trainer = Trainer(  train_config, device )
		model, cal_model, train_logs, dev_logs, test_logs = trainer.forward( model, 
																self.train_set, 
																self.dev_set, 
																self.test_set, 
																file_name )
		return model, cal_model, train_logs, dev_logs, test_logs


	def manual_search( self, config_file ):
		# Manual hparam search.
		# Initialize dict to store logs for all hparam searched.
		train_logs_dict, dev_logs_dict,test_logs_dict = {}, {}, {}
		total_time = {}

		version = config_file.Version
		device = config_file.conf.model_params.device

		model_config = config_file.conf.model_params
		train_config = config_file.conf.train_params
		
		# Create combinations of tunable hyperparameters.
		# Currently tuning - projection_dim, num_hid_layers, activation1, learning rate, 
		# 					 contact threshold, log-weight, dropouts, weight_decay, gamma.
		hparam_comb = []
		for player in model_config.projection_layer:
			for nlayers in model_config.num_hid_layers:
				for activ1 in model_config.activation1:
					for lr in train_config.learning_rate:
						for thresh in train_config.contact_threshold:
							for w1 in train_config.log_weight:
								for drop in model_config.dropouts:
									for wd in train_config.weight_decay:
										for g in train_config.scheduler.gamma:
											hparam_comb.append( ( player, nlayers, activ1, lr, thresh, w1, drop, wd, g ) )


		num_runs = train_config.Nruns
		trial_seeds = [self.global_seed + i*1111  for i in range( num_runs )]
		
		for comb in hparam_comb:
			# key represents the hparam combination used.
			key = ""
			# key = "_".join( comb )
			for param_idx, param in enumerate( comb):
				if param_idx == len( comb ) - 1:
					key += f"{param}"
				else:
					key += f"{param}_"
			# key = f"{comb[0]}_{comb[1]}_{comb[2]}_{comb[3]}_{comb[4]}_{comb[5]}_{comb[6]}_{comb[7]}_{comb[8]}"
			
			# Train set also has loss, epoch and time apart from metrics.
			train_logs_dict[key] = np.zeros( ( train_config.max_epochs, self.num_metrics + 3 ) )
			# Dev set has loss and epoch apart from metrics.
			dev_logs_dict[key] = np.zeros( ( train_config.max_epochs, self.num_metrics + 2 ) )
			# test set only has the metrics.
			test_logs_dict[key] = np.zeros( ( train_config.max_epochs, self.num_metrics ) )

			model_config.projection_layer = comb[0]
			model_config.num_hid_layers = comb[1]
			model_config.activation1 = comb[2]
			train_config.learning_rate = comb[3]
			train_config.contact_threshold = comb[4]
			train_config.log_weight = comb[5]
			model_config.dropouts = comb[6]
			train_config.weight_decay = comb[7]
			train_config.scheduler.gamma = comb[8]

			tic = time.time()
			print( "----------- Starting for: " )
			print( "projection_layer = ", comb[0], "\tnum_hid_layers = ", comb[1], 
					"\tactivation1 = ", comb[2], 
					"\tlr = ", comb[3], "\tcontact threshold = ", comb[4],
					"\tlog_weight = ", comb[5],
					"\tdropouts = ", comb[6],
					"\tweight_decay = ", comb[7], "\tgamma = ", comb[8] )

			for run, seed in enumerate( trial_seeds ):
				print( f"------->> Run {run}" )
				self.seed_worker( seed )

				model, cal_model, train_logs, dev_logs, test_logs = self.executor( 
																	model_config,
																	train_config,
																	key )

				train_logs_dict[key] += train_logs
				dev_logs_dict[key] += dev_logs
				test_logs_dict[key] += test_logs

				if train_config.save_model:
					torch.save( model.state_dict(), f"./{train_config.model_path}-{key}__{run}.pth" )

					if cal_model != None:
						joblib.dump( cal_model, f"./{train_config.model_path}-{key}__{run}.pkl" )
				
				del model
				print( f"Aah...Some much needed rest..." )
				time.sleep( 5 )

			toc = time.time()
			time_ = ( toc-tic )/3600
			total_time[key] = time_
			# Averaging predictions across all runs.
			train_logs_dict[key] = train_logs_dict[key]/train_config.Nruns
			dev_logs_dict[key] = dev_logs_dict[key]/train_config.Nruns
			test_logs_dict[key] = test_logs_dict[key]/train_config.Nruns
			
			# Plot train and dev set metric and loss.
			create_plots( train_logs_dict, dev_logs_dict, key )
			print( "----------- Completed for: " )
			print( "projection_layer = ", comb[0], "\tnum_hid_layers = ", comb[1], 
					"\tactivation1 = ", comb[2], 
					"\tlr = ", comb[3], "\tcontact threshold = ", comb[4],
					"\tlog_weight = ", comb[5], "\tdropouts = ", comb[6],
					"\tweight_decay = ", comb[7], "\tgamma = ", comb[8] )


		OmegaConf.save( config = total_time, f = "Total_time.yml" )
		np.save( "logs_train_{}_{}.npy".format( model_config.Model, version ), train_logs_dict )
		np.save( "logs_dev_{}_{}.npy".format( model_config.Model, version ), dev_logs_dict )

		# Write out the summary report.
		dump_metrics( train_logs_dict, 
						dev_logs_dict, 
						test_logs_dict, 
						total_time, "_{}_{}".format( model_config.Model, version ), self.num_metrics )

		print("May the Force b with u...")


	def forward( self, args ):
		##Begin the search for hyperparameters#################
		#######################################################
		# Load the datasets.
		# Split data and create dataloaders or use presaved dataloader.
		t1 = time.time()
		self.call_of_duty()
		t2 = time.time()
		print( f"Time taken to load the dataset = {( t2 - t1 )/60} minutes" )

		# Output path directory.
		version = config_file.Version
		PATH = config_file.conf.dataset.output_path + f"/Version_{version}/"

		# Make the directory if it does not exist.
		if not os.path.isdir( PATH ):
			os.makedirs( PATH )
		else:
			print( "/n Directory already exists..." )
			reply = input( "Wanna overwrite? (Y/n) = " )
			if reply == "Y":
				print( "Overwriting\n" )
			else:
				print( "Terminating the process..." )
				exit()

		# Move the versions file to the respective directory.
		os.rename( f"./version_{config_file.Version}.yml", PATH + f"version_{version}.yml")

		# Change to output directory.
		os.chdir( PATH )

		# Check the mode of optimization.
		if args.m == "optuna":
			self.optuna_search( config_file )

		elif args.m == "manual":
			self.manual_search( config_file )



# The Showdown ####################################
###################################################
args = parser.parse_args()
search = HparamSearch().forward( args )
print( "May the Force be with you..." )
###################################################
###################################################