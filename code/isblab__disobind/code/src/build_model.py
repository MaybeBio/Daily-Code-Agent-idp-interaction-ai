############ Creating employment for another downstream script ###########
## U can think of it as running the model with specified configurations ##
############# ------>"May the Force serve u well..." <------##############
##########################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import time
import random
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
from betacal import BetaCalibration

import torch
from torch import nn

from metrics import torch_metrics
from loss import get_loss_function
from utils import ( plot_output, prepare_input, plot_reliabity_diagram )

import warnings
warnings.filterwarnings("ignore")

class Trainer (nn.Module):
	def __init__ (self, config, device):
		super(Trainer, self).__init__()
		self.objective = config.objective
		self.emb_type = config.emb
		self.optim = config.optimizer
		self.amsgrad = config.amsgrad
		self.wd = config.weight_decay
		self.scheduler_config = config.scheduler
		self.lr = config.learning_rate
		self.max_norm = config.max_norm
		self.mask = config.mask
		self.loss = config.loss
		self.loss_func = get_loss_function( self.loss )
		self.weight = config.log_weight
		self.num_metrics = config.num_metrics[0]
		self.multidim_avg = config.num_metrics[1]
		self.method = config.calibration
		self.max_epochs = config.max_epochs
		self.threshold = config.contact_threshold
		self.device = device
		self.max_len = config.max_seq_len
		self.prec = 4    # Significant figures.


	def optimizer( self ):
		if self.optim == "Adam":
			return torch.optim.Adam( self.model1.parameters(), lr = self.lr, weight_decay = self.wd, amsgrad = self.amsgrad )
		
		elif self.optim == "SGD":
			return torch.optim.SGD( self.model1.parameters(), lr = self.lr, weight_decay = self.wd, momentum = 0.9 )
		
		elif self.optim == "AdamW":
			return torch.optim.AdamW( self.model1.parameters(), lr = self.lr, weight_decay = self.wd, amsgrad = self.amsgrad )


	def scheduler( self, model, optimizer ):
		if self.scheduler_config.apply:
			name = self.scheduler_config.name
			if name == "swa":
				print( "Using Stochastic Weight Averaging" )
				self.swa_model = torch.optim.swa_utils.AveragedModel( model )
				self.swa_start = scheduler_config.swa_start
				self.swa_scheduler = torch.optim.swa_utils.SWALR( optimizer, swa_lr = self.swa_lr )
			elif name == "cycliclr":
				print( "Using cyclical learning rate scheduler..." )
				base_lr = self.scheduler.base_lr
				step_size_up = self.scheduler.step_size_up
				step_size_down = self.scheduler.step_size_down
				scheduler = torch.optim.lr_scheduler.CyclicLR( optimizer, base_lr = base_lr, max_lr = self.lr, 
																step_size_up = step_size_up, step_size_down = step_size_down )
			elif name == "multistep":
				print( "Using MultiStepLR scheduler..." )
				milestone = self.scheduler_config.milestone
				gamma = self.scheduler_config.gamma
				scheduler = torch.optim.lr_scheduler.MultiStepLR( optimizer, milestones = milestone, gamma = gamma )

			elif name == "exp":
				print( "Using ExponentialLR scheduler..." )
				gamma = self.scheduler_config.gamma
				scheduler = torch.optim.lr_scheduler.ExponentialLR( optimizer, gamma = gamma )

			elif name == "linear":
				print( "Using LinearLR scheduler..." )
				start = self.scheduler_config.start_factor
				end = self.scheduler_config.end_factor
				total_iters = self.scheduler_config.total_iters
				scheduler = torch.optim.lr_scheduler.LinearLR( optimizer, start_factor = start, 
																end_factor = end, 
																total_iters = total_iters )

		else:
			print( "Not using any scheduler.." )
			scheduler = None

		return scheduler


	def get_inputs( self, batch ):
		"""
		Create input and target for different kinds of training objective.
		Interaction prediction
			Given prot1 and prot2 embeddings, predict the interactions for all vs all residues.
				e.g. ( L1, C ) and ( L2, C ) --> ( L1, L2 )
		Interface prediction
			Given prot1 and prot2 embeddings, predict the residues in volved in interaction
				for prot1 and prot2.
				e.g. ( L1, C ) and ( L2, C ) --> ( L1 + L2 )
		Binned interaction prediction
			Same as interaction prediction, but only for binned (coarse grained) residues.
				Slide a kernel over the target to get binned output.
				No need to bin the input.
		Binned interface prediction
			Same as interface prediction, but only for binned (coarse grained) residues.
				Slide a kernel over the target to get binned output.
				No need to bin the input.
		target --> tensor of shape (N, L1, L2)
		p1_emb --> tensor of shape (N, L1, C )
		p2_emb --> tensor of shape (N, L2, C ) 
		"""
		target = batch[:,:,-2*self.max_len:]
		# Separate the target and mask.
		target, target_mask = target[:,:,:self.max_len], target[:,:,self.max_len:]

		if self.emb_type == "ProSE":
			p1_emb, p2_emb = batch[:,:,:6165], batch[:,:,6165:-2*self.max_len]
		elif self.emb_type in ["T5", "ProstT5", "BERT"]:
			p1_emb, p2_emb = batch[:,:,:1024], batch[:,:,1024:-2*self.max_len]
		elif self.emb_type == "ESM2-650M":
			p1_emb, p2_emb = batch[:,:,:1280], batch[:,:,1280:-2*self.max_len]
		else:
			raise Exception( "Incorrect embedding used..." )

		( p1_emb, p2_emb, target,
			target_mask, interaction_mask ) = prepare_input( prot1 = p1_emb, prot2 = p2_emb, target = target, 
															target_mask = [self.mask[1], target_mask], 
															objective = self.objective[0], 
															bin_size = self.objective[1], 
															bin_input = self.objective[4], 
															single_output = self.objective[5] )

		p1_emb = p1_emb.to( self.device ).float()
		p2_emb = p2_emb.to( self.device ).float()
		target = target.to( self.device ).float()
		target_mask = target_mask.to( self.device ).float()
		interaction_mask = interaction_mask.to( self.device ).float()

		return p1_emb, p2_emb, target, target_mask, interaction_mask



	def predict( self, batch, train = False ):
		"""
		Perform forward pass.

		Input:
		----------
		batch --> minibatch for the forward pass.
		train --> (bool) flag to highlight the training setp.

		Return:
		----------
		preds --> model predictions.
		target --> target for the minibatch.
		target_mask --> binary mask for ignoring contribution from padding.
		"""
		p1_emb, p2_emb, target, target_mask, interaction_mask = self.get_inputs( batch )

		preds = self.model1( p1_emb, p2_emb, interaction_mask )

		del p1_emb, p2_emb
		return preds, target, target_mask


	def calculate_loss_n_metrics( self, preds, target, target_mask, train = False ):
		"""
		Calculate loss and evaluation metrocs.

		Input:
		----------
		preds --> model predictions.
		target --> target for the minibatch.
		target_mask --> binary mask for ignoring contribution from padding.
		train --> (bool) flag to highlight the training setp.

		Return:
		----------
		loss --> loss calculated for the minibatch.
		metrics --> all evaluation metrics for the minibatch.
		preds --> model predictions.
		target --> target for the minibatch.
		target_mask --> binary mask for ignoring contribution from padding.
		"""
		loss = self.loss_func.forward( preds, target, 
											self.device, self.weight, [self.mask[1], target_mask] )

		if self.loss in ["bce_with_logits", "representation_loss"]:
			preds = torch.sigmoid( preds )

		elif self.loss == "interface":
			preds = torch.cat( preds, axis = 1 )
			target_mask = torch.cat( target_mask, axis = 1 )
			target = torch.cat( target, axis = 1 )

		if self.mask[1]:
			preds = preds*target_mask
		metrics = torch_metrics( preds, target, self.threshold, self.multidim_avg, self.device )

		preds = preds.detach().cpu().numpy()
		target = target.detach().cpu().numpy()
		return loss, metrics, preds, target, target_mask



	def training_step( self, train_set, epoch ):
		"""
		Perform forward pass for all mini-batches in training set.

		Input:
		----------
		train_set --> dataloader object for train set.
		epoch --> (int) training epoch.

		Return:
		----------
		batch_dict --> dict containing metrics for all mini-batches in train set.
		"""
		##############------------------##############
		print( "\nTraining step..." )
		##############------------------##############
		# Dictionary to store batch loss and metric values
		batch_dict = np.array( [] )
		batch_size = 0

		uncal_preds, target, mask = [], [], []

		for batch in train_set:
			batch_size += 1
			
			self.optimizer1.zero_grad()

			train_pred, train_target, train_mask = self.predict( batch, train = True )

			train_loss, train_metrics, train_pred, train_target, train_mask = self.calculate_loss_n_metrics( train_pred, train_target, train_mask, train = True )

			batch_dict = np.append( batch_dict, train_loss.item() )
			batch_dict = np.append( batch_dict, [train_metrics[0].item(), 
												train_metrics[1].item(), 
												train_metrics[2].item(),
												train_metrics[3].item(),
												train_metrics[4].item(),
												train_metrics[5].item(),
												train_metrics[6].item()] )

			# Calculate gradients per batch.
			train_loss.backward()

			if epoch == self.max_epochs - 1:
				if len( uncal_preds ) == 0:
					uncal_preds = train_pred
					target = train_target
					mask = train_mask.detach().cpu().numpy()
				else:
					uncal_preds = np.concatenate( ( uncal_preds, train_pred ), axis = 0 )
					target = np.concatenate( ( target, train_target ), axis = 0 )
					mask = np.concatenate( ( mask, train_mask.detach().cpu().numpy() ), axis = 0 )

			if self.max_norm != None:
				torch.nn.utils.clip_grad_norm_( self.model1.parameters(), max_norm = 1 )

			# Update params for every batch.
			self.optimizer1.step()
			
			if self.scheduler_config.name == "swa" and epoch > self.swa_start:
				self.swa_model.update_parameters( self.model1 )
				self.swa_scheduler.step()
		
		if self.scheduler1 != None:
			self.scheduler1.step()

		batch_dict = batch_dict.reshape( batch_size, self.num_metrics + 1 )

		if epoch == self.max_epochs - 1:
			uncal_preds = uncal_preds.flatten()
			target = target.flatten()
			mask = mask.flatten()
			self.calibrate_model( uncal_preds*mask, target*mask )
				

		return batch_dict



	def calibrate_model( self, preds, target ):
		"""
		Train a calibration model using the specified method.
		"""
		if self.method == "platt":
			print( "Using Platt scaling..." )
			self.cal_model = LogisticRegression()
			self.cal_model.fit( preds.reshape( -1, 1 ), target.reshape( -1, 1 ) )

		elif self.method == "iso":
			print( "Using Isotonic regression..." )
			self.cal_model = IsotonicRegression().fit( preds, target )

		elif "beta" in self.method:
			print( "Using beta calibration..." )
			_, params = self.method.split( "-" )
			self.cal_model = BetaCalibration( parameters = params )
			self.cal_model.fit( preds, target )

		elif self.method in ["temp", "None"]:
			self.cal_model = None



	def get_calibrated_preds( self, preds ):
		"""
		Get calibrated predictions.
		"""
		if self.method == "platt":
			cal_preds = self.cal_model.predict_proba( preds.reshape( -1, 1 ) )[:,1]

		elif self.method == "iso":
			cal_preds = self.cal_model.predict( preds )

		elif "beta" in self.method:
			cal_preds = self.cal_model.predict( preds )

		elif self.method in ["temp", "None"]:
			cal_preds = preds

		return cal_preds



	def validation_step ( self, dev_set, epoch ):
		"""
		Perform forward pass for all mini-batches in dev set.

		Input:
		----------
		dev_set --> dataloader object for dev set.
		epoch --> (int) training epoch.

		Return:
		----------
		batch_dict --> dict containing metrics for all mini-batches in dev set.
		"""
		##############------------------##############
		print("Validation step")
		##############------------------##############
		batch_dict = np.array( [] )

		uncal_preds, target, mask = [], [], []

		batch_size = 0
		with torch.no_grad():
			for batch in dev_set:
				batch_size += 1

				dev_pred, dev_target, dev_mask = self.predict( batch )

				dev_loss, dev_metrics, dev_pred, dev_target, dev_mask = self.calculate_loss_n_metrics( dev_pred, dev_target, dev_mask )

				if epoch == self.max_epochs - 1:
					if len( uncal_preds ) == 0:
						uncal_preds = dev_pred
						target = dev_target
						mask = dev_mask.detach().cpu().numpy()
					else:
						uncal_preds = np.concatenate( ( uncal_preds, dev_pred ), axis = 0 )
						target = np.concatenate( ( target, dev_target ), axis = 0 )
						mask = np.concatenate( ( mask, dev_mask.detach().cpu().numpy() ), axis = 0 )
				
				# Keep logs for all the batches
				batch_dict = np.append( batch_dict, dev_loss.item() )
				batch_dict = np.append( batch_dict, [dev_metrics[0].item(), 
													dev_metrics[1].item(), 
													dev_metrics[2].item(),
													dev_metrics[3].item(),
													dev_metrics[4].item(),
													dev_metrics[5].item(),
													dev_metrics[6].item()] )

		batch_dict = batch_dict.reshape( batch_size, self.num_metrics + 1 )

		if epoch == self.max_epochs - 1:
			uncal_preds = uncal_preds.flatten()
			target = target.flatten()
			mask = mask.flatten()
			self.calibrate_model( uncal_preds*mask, target*mask )


		return batch_dict



	def test_step( self, test_set, file_name ):
		"""
		Perform forward pass for all mini-batches in test set.

		Input:
		----------
		test_set --> dataloader object for test set.
		epoch --> (int) training epoch.

		Return:
		----------
		test_met
		test_target
		test_target
		batch_dict --> dict containing metrics for all mini-batches in test set.
		"""
		##############------------------##############
		print("Predicting for Test set")
		##############------------------##############
		batch_dict = np.array( [] )

		uncal_preds, cal_preds = [], []

		batch_size = 0
		with torch.no_grad():
			for batch in test_set:
				batch_size += 1

				test_pred, test_target, test_mask = self.predict( batch )

				test_pred = test_pred.detach().cpu().numpy()
				test_target = test_target.detach().cpu().numpy()
				test_mask = test_mask.detach().cpu().numpy()
				shape = test_target.shape

				if len( uncal_preds ) == 0:
					uncal_preds = test_pred*test_mask
					target = test_target
				
				else:
					uncal_preds = np.concatenate( ( uncal_preds, test_pred*test_mask ), axis = 0 )
					target = np.concatenate( ( target, test_target ), axis = 0 )
				
				test_pred = self.get_calibrated_preds( test_pred.flatten() )
				test_pred = test_pred.reshape( shape )

				if len( cal_preds ) == 0:
					cal_preds = test_pred*test_mask
				
				else:
					cal_preds = np.concatenate( ( cal_preds, test_pred*test_mask ), axis = 0 )
                
				test_pred = torch.from_numpy( test_pred ).to( self.device )
				test_target = torch.from_numpy( test_target ).to( self.device )
				test_mask = torch.from_numpy( test_mask ).to( self.device )

				_, test_metrics, test_pred, test_target, test_mask = self.calculate_loss_n_metrics( test_pred, test_target, test_mask )

				# Keep logs for all the batches
				batch_dict = np.append( batch_dict, [test_metrics[0].item(), 
													test_metrics[1].item(), 
													test_metrics[2].item(),
													test_metrics[3].item(),
													test_metrics[4].item(),
													test_metrics[5].item(),
													test_metrics[6].item()] )

		batch_dict = batch_dict.reshape( batch_size, self.num_metrics )
	
		# Stored values - Recall, Precision, F1score, AvgPreciison, MCC, AUROC, Accuracy.
		test_met = np.round( np.mean( batch_dict, axis = 0 ), self.prec )

		# For self.method = temp and None, uncal_preds = cal_preds.
		if self.method == "temp":
			plot_reliabity_diagram( [], cal_preds, target, file_name )

		elif self.method == "None":
			pass
			# plot_reliabity_diagram( uncal_preds, [], target, file_name )
		
		else:
			plot_reliabity_diagram( uncal_preds, cal_preds, target, file_name )

		# Get calibrated predictions for test set.
		print( f"Recall = {test_met[0]} \t", f"Precision = {test_met[1]} \t",
				f"F1score = {test_met[2]}\t", "\n" )

		return test_met, test_pred, test_target


	# def optuna_mode( self, model, train_set, dev_set ):
	# 	# Using optuna for optimizing hyperparameters.
	# 	optimizer1 = self.optimizer( model )

	# 	tic = time.time()
	# 	model.train()
	# 	model, _ = self.training_step( model, train_set, optimizer, None )
	# 	model.eval()
	# 	dev_batch_dict, _, _ = self.validation_step( model, dev_set )
	# 	dev_met = np.round( np.mean( dev_batch_dict, axis = 0 ), self.prec )
	# 	print( dev_met )
	# 	f1_score = dev_met[3]
	# 	toc = time.time()
	# 	print( f"F1: {f1_score} \t Time: {toc-tic}")

	# 	return f1_score


	def forward ( self, model, train_set, dev_set, test_set, file_name ):
		"""
		Perform model training.

		Input:
		----------
		model --> instance of the model class to be trained.
		train_set --> dataloader object for train set.
		dev_set --> dataloader object for dev set.
		test_set --> dataloader object for test set.
		file_name --> (str) identifier for output files.

		Returns:
		----------
		model1 --> trained base model.
		cal_model --> trained calibration model.
		train_mat --> loss, metrics, epochs, time for training set.
		dev_mat --> loss, metrics, epochs for dev set.
		test_mat --> loss, metrics for test set.
		"""
		##############------------------##############
		# Creating the logs dictionary
		train_mat = np.array( [] )
		dev_mat = np.array( [] )

		print( "Using device: ", torch.cuda.get_device_name( 0 ) )
		print( f"Input mask: {self.mask[0]} \t Output mask: {self.mask[1]}" )

		##############------------------##############
		self.model1 = model
		self.optimizer1 = self.optimizer()
		self.scheduler1 = self.scheduler( self.model1, self.optimizer1 )

		for epoch in range( self.max_epochs ):
			tic = time.time()

			self.model1.train()
			train_batch_dict = self.training_step( train_set, epoch )

			self.model1.eval()
			dev_batch_dict = self.validation_step( dev_set, epoch )

			toc = time.time()
			time_ = toc - tic
						
			# Stored values - Loss, Recall, Precision, F1score, AvgPrecision, MatthewsCorrCoef, AUROC, Accuracy, Epoch, Time.
			epoch_met = np.round( np.mean( train_batch_dict, axis = 0 ), self.prec )
			# print( epoch_met )
			# exit()
			train_mat = np.append( train_mat, epoch_met )
			train_mat = np.append( train_mat, [epoch, ( time_ )] )

			# Stored values - Loss, Recall, Precision, F1score, AUROC, MatthewsCorrCoef, Epoch.
			dev_met = np.round( np.mean( dev_batch_dict, axis = 0 ), self.prec )
			dev_mat = np.append( dev_mat, dev_met )
			dev_mat = np.append( dev_mat, epoch )

			percent_dif = ( epoch_met[0] - dev_met[0] )/epoch_met[0]

			print("Loss:      {} ---- {} :: {} %".format( epoch_met[0]   , dev_met[0], percent_dif*100 ) )
			print("Recall:    {} ---- {}".format( epoch_met[1]   , dev_met[1] ) )
			print("Precision: {} ---- {}".format( epoch_met[2]   , dev_met[2] ) )
			print("F1_Score:  {} ---- {}".format( epoch_met[3]   , dev_met[3] ) )
			print("Time taken per Epoch: {} seconds".format( time_ ) )
			print("<<--Completed for epoch {} -->>".format( epoch ), "\n" )

		# Predicting for the dev set with swa_model.
		if self.scheduler_config.apply and self.scheduler_config.name == "swa":
			model = self.swa_model
			print("Predicting for Dev set")
			_, _, _ = self.test_step( dev_set, file_name )

		# Predicting for the test set.
		test_mat,  test_pred, test_target = self.test_step( test_set, file_name )

		if "interaction" in self.objective[0]:
			if "bin" in self.objective[0]:
				length = self.max_len//self.objective[1]
			else:
				length = self.max_len
			test_pred = test_pred.reshape( -1, length, length )
			test_target = test_target.reshape( -1, length, length )

			plot_output(
				pred = test_pred, 
				target = test_target,
				threshold = self.threshold,
				file_name = file_name, interface = False
				)
		elif "interface" in self.objective[0]:
			plot_output(
				pred = test_pred, 
				target = test_target,
				threshold = self.threshold,
				file_name = file_name, interface = True
				)

		else:
			raise Exception( "Incorrect objective specified..." )

		del test_pred, test_target
		train_mat = train_mat.reshape( self.max_epochs, self.num_metrics + 3 )
		dev_mat = dev_mat.reshape( self.max_epochs, self.num_metrics + 2 )

		del train_set, dev_set, test_set
		torch.cuda.empty_cache()
		return self.model1, self.cal_model, train_mat, dev_mat, test_mat

