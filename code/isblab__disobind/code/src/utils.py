####################### Accessory functions #######################
########## ------>"May the Force serve u well..." <------##########
###################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean
from imblearn.over_sampling import SMOTE, ADASYN
from sklearn.calibration import calibration_curve

import torch
from torch import nn



## Plot Reliability diagram########################
###################################################
def plot_reliabity_diagram( uncal_preds, cal_preds, target, file_name ):
	target = target.flatten()
	fig, axis = plt.subplots( 1, 1, figsize = ( 10, 10 ) )

	df1 = pd.DataFrame()
	df2 = pd.DataFrame()

	if len( cal_preds ) != 0:
		cal_preds = cal_preds.flatten()
		target_prob1, cal_pred_prob = calibration_curve( target, cal_preds, n_bins = 10 )
		df1["Cal_preds"] = cal_preds
		# Mean predicted probability
		df2["Cal_pred_prob"] = cal_pred_prob
		# Fraction of positives
		df2["Target_prob1"] = target_prob1

		axis.plot( cal_pred_prob, target_prob1, marker = "o", color = "blue" )

	if len( uncal_preds ) != 0:
		uncal_preds = uncal_preds.flatten()
		target_prob2, uncal_pred_prob = calibration_curve( target, uncal_preds, n_bins = 10 )
		df1["UnCal_preds"] = uncal_preds
		# Mean predicted probability
		df2["UnCal_pred_prob"] = uncal_pred_prob
		# Fraction of positives
		df2["Target_prob2"] = target_prob2

		axis.plot( uncal_pred_prob, target_prob2, marker = "v", color = "red" )

	axis.plot( [0, 1], [0, 1], color = "gray" )
	axis.xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
	axis.yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
	axis.set_xlabel( "Mean predicted probability", fontsize = 16 )
	axis.set_ylabel( "Fraction of Positives", fontsize = 16 )

	plt.savefig( f"Reliability_diagram_{file_name}.png", dpi = 300 )
	plt.close()
	
	df1["Target"] = target

	df1.to_csv( f"Test_preds_{file_name}.csv" )
	df2.to_csv( f"Calibration_probs_{file_name}.csv" )



## Oversampling using MSOTE and ADASYN#############
###################################################
def oversample( self, sampler, entry, output_labels ):
	"""
	Create synthetic data for oversampling the minority class.
		Currently using SMOTE and ADASYN.
	"""
	if sampler == None:
		return entry

	else:		
		if sampler == "smote":
			sampler = SMOTE( random_state = 42, k_neighbors = 6 )
		
		elif sampler == "adasyn":
			sampler = ADASYN( random_state = 42, n_neighbors = 5 )
		
		entry, output_labels = sampler.fit_resample( entry, output_labels )
		output_labels = output_labels.reshape( len( output_labels ), 1 )
		entry = np.hstack( ( entry, output_labels ) )
		# np.random.shuffle( entry )
		return entry


## Prepare input tensors for specific objective####
###################################################
def prepare_input( prot1, prot2, target, target_mask, objective, bin_size, bin_input, single_output ):
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

		L1 = L2 if input is padded.
	prot1 --> tensor of shape (N, L1, C )
	prot2 --> tensor of shape (N, L2, C ) 
	target --> tensor of shape (N, L1, L2)
	target_mask --> tensor of shape (N, L1, L2)
	objective --> one of [interaction, interaction_bin, interface, interface_bin]
	bin_size --> for coarse grained prediction, this specifies level of coarse graining.
	bin_input --> if True average the input embeddings.
	single_output --> convert multi-output prediction ro single-output prediction.
	"""

	# print( prot1.shape, "\t", prot2.shape )
	apply_mask, target_mask = target_mask
	if apply_mask:
		max_len = [100, 100]
	else:
		max_len = [prot1.shape[1], prot2.shape[1]]
	
	if bin_size > 1:
		# MaxPool2d requires input of shape ( N, C, L1, L2 ).
		# 	target: ( N, L1, L2 ) --> ( N, 1, L1, L2 )
		eff_len = [max_len[0]//bin_size, max_len[1]//bin_size]
		
		with torch.no_grad():
			if bin_input:
				a = nn.AvgPool1d( kernel_size = bin_size, stride = bin_size )
				prot1 = torch.permute( prot1, ( 0, 2, 1 ) )
				prot2 = torch.permute( prot2, ( 0, 2, 1 ) )
				prot1 = a( prot1 )
				prot2 = a( prot2 )
				prot1 = torch.permute( prot1, ( 0, 2, 1 ) )
				prot2 = torch.permute( prot2, ( 0, 2, 1 ) )

			m = nn.MaxPool2d( kernel_size = bin_size, stride = bin_size )
			if target != None:
				target = target.unsqueeze( 1 )
				target = m( target )
				# Remove the extra C dimension.
				target = target.squeeze( 1 )

			if apply_mask:
				target_mask = target_mask
				target_mask = m( target_mask )
				target_mask = target_mask.squeeze( 1 )
	else:
		eff_len = max_len

	# For masking pads in the interaction_tensor in model.
	interaction_mask = target_mask.clone().unsqueeze( -1 )

	if "interaction" in objective:
		# target: ( N, L1, L2 ) --> ( N, L1xL2 )
		if target != None:
			target = target.reshape( target.shape[0], -1 )

		if apply_mask:
			target_mask = target_mask.reshape( target_mask.shape[0], -1 )

	elif "interface" in objective:
		# Find indexes of conatct elements.
		if target != None:
			idx = torch.where( target == 1 )
			# Get interacting residues for prot1.
			p1_target = torch.zeros( ( prot1.shape[0], target.shape[1] ) )
			p1_target[idx[0],idx[1]] = 1

			# Get interacting residues for prot2.
			p2_target = torch.zeros( ( prot2.shape[0], target.shape[2] ) )
			p2_target[idx[0], idx[2]] = 1

			target = torch.cat( ( p1_target, p2_target ), axis = 1 )

		if apply_mask:
			# Project the mask for interface prediction
			idx_mask = torch.where( target_mask == 1 )
			p1_target_mask = torch.zeros( ( prot1.shape[0], target_mask.shape[1] ) )
			p1_target_mask[idx_mask[0],idx_mask[1]] = 1
			
			p2_target_mask = torch.zeros( ( prot2.shape[0], target_mask.shape[2] ) )
			p2_target_mask[idx_mask[0],idx_mask[2]] = 1
			
			# target: ( N, L1, L2 ) --> ( N, L1+L2 )
			target_mask = torch.cat( ( p1_target_mask, p2_target_mask ), axis = 1 )

	# prot1 = torch.cat( ( prot1, prot2 ), axis = 0 )

	# print( prot1.size(), "\t", target.size(), "\t", target_mask.size() )
	# exit()
	if single_output:
		prot1 = prot1.numpy()
		prot1 = np.repeat( prot1, eff_len[1], axis = 1 )
		prot1 = torch.from_numpy( prot1 )
		prot2 = torch.tile( prot2, ( eff_len[0], 1 ) )

	return prot1, prot2, target, target_mask, interaction_mask



def plot_output( pred, target, threshold, file_name, interface = False ):
	# path = "./Output_Plots/"
	# os.makedirs( path )
	# Plot contact maps of IDR and R
	pred = np.where( pred > threshold, 1, 0 )
	target = np.where( target > threshold, 1, 0 )
	
	# if write_counts:
	# 	counts = {}
	# 	counts["Actual_counts"] = [np.count_nonzero(y_dev[i]) for i in range(y_dev.shape[0])] 
	# 	counts["Predicted_counts"] = [np.count_nonzero(cm_dev[i]) for i in range(cm_dev.shape[0])]
	# 	OmegaConf.save( config = counts, f = "contact_counts_{}.yml".format( file_name ) )
	if interface:
		x = 1
		fig, axis = plt.subplots( x, 2, figsize = ( 20, 15 ) )
		# cm = plt.cm.gray.reversed()
		axis[0].imshow( pred[:10,:], aspect = "auto" )
		axis[1].imshow( target[:10,:], aspect = "auto"  )
		axis[0].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis[1].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis[0].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis[1].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis[0].set_xlabel( "Predicted_contact_map", fontsize = 16 )
		axis[1].set_xlabel( "Target_contact_map", fontsize = 16 )

	else:
		if pred.shape[0] == 1:
			x = 1
			fig, axis = plt.subplots( x, 2, figsize = ( 20, 22 ) )
			# cm = plt.cm.gray.reversed()
			axis[0].imshow( pred[0,:], aspect = "auto" )
			axis[1].imshow( target[0,:], aspect = "auto"  )
			axis[0].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			axis[1].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			axis[0].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			axis[1].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			axis[0].set_xlabel( "Predicted_contact_map", fontsize = 16 )
			axis[1].set_xlabel( "Target_contact_map", fontsize = 16 )
		else:
			x = cm_dev.shape[0] if pred.shape[0] < 8 else 8
			fig, axis = plt.subplots( x, 2, figsize = ( 20, 22 ) )
			for i in range( x ):
				# cm = plt.cm.gray.reversed()
				axis[i,0].imshow( pred[i,:,:], aspect = "auto" )
				axis[i,1].imshow( target[i,:,:], aspect = "auto"  )
				axis[i,0].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
				axis[i,1].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
				axis[i,0].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
				axis[i,1].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
				axis[i,0].set_xlabel( "Predicted_contact_map", fontsize = 16 )
				axis[i,1].set_xlabel( "Target_contact_map", fontsize = 16 )
	
	plt.savefig( "Predictions_cm_dev_{}.png".format( file_name ), dpi = 300 )
	# plt.savefig( "{}Predictions_cm_{}_{}_{}.png".format( path, set_, file_name, i ), dpi = 300 )
	plt.close()


# Plot Metrics for train and dev set with 3 hparams
###################################################
def create_plots( train, dev, key_ = None ):
	# Plot metric values for the train and dev set
	for key in train.keys():
		fig, axis = plt.subplots( 4, 2, figsize = ( 14, 16 ) )
		x, y = 0, 0
		for label, idx in zip( ["Loss", "Recall", "Precision", "F1score", "AvgPrecision", "MatthewsCorrCoef", "AUROC", "Accuracy"], range( 8 ) ):
			if idx % 2 == 0 and idx != 0:
				x += 1
			
			axis[x,y].plot( train[key][:,8], train[key][:,idx], c = "blue", label = "Train" )
			axis[x,y].plot( dev[key][:,8], dev[key][:,idx], c = "red", label = "Dev" )
			axis[x,y].set_xlim( [0, train[key][:,8][-1]] )
			axis[x,y].set_xlabel( "Epochs", fontsize = 18 )
			axis[x,y].xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			if label != "Loss":
				axis[x,y].set_ylim( [-0.05, 1.1] )
			axis[x,y].set_ylabel( f"{label} per Epoch", fontsize = 18 )
			axis[x,y].yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
			axis[x,y].legend()

			y = 0 if y == 1 else 1

		plt.subplots_adjust( left = 0.1,
		                    bottom = 0.1, 
		                    right = 0.9, 
		                    top = 0.9, 
		                    wspace = 0.3, 
		                    hspace = 0.3 )

		plt.savefig( "Metrics_perEpoch_{}.png".format( key ), dpi = 300 )
		plt.close()


def dump_metrics( train, dev, test, time_dict, version, num_metrics, epoch = -1 ):
	# Write the metric, loss, time to file.
	w = open ( f"summary_report{version}.txt", "w" )
	w.writelines("-------------- Summary Report -------------- \n\n")

	for key in train.keys():
		projection_layer, hid_layers, activ1, lr, thresh, w1, drop, wd, gamma = key.split( "_" )
		w.writelines( f"\n Projection layer = {str( projection_layer )}" )
		w.writelines( f"\n No. of hidden layers = {str( hid_layers )}" )
		w.writelines( f"\n Activation1 = {str( activ1 )}" )
		w.writelines( f"\n Learning rate = {str( lr )}" )
		w.writelines( f"\n Contact threshold = {str( thresh )}" )
		w.writelines( f"\n Log weight = {str( w1 )}" )
		# w.writelines( f"\n Activation1 param = {str( alpha )}" )
		w.writelines( f"\n Dropout probability = {str( drop )}" )
		# w.writelines( f"\n Dropout2 probability = {str( d2 )}" )
		w.writelines( f"\n Weight decay = {str( wd )}" )
		w.writelines( f"\n Gamma (scheduler param) = {str( gamma )}" )
		# w.writelines( f"\n Max norm = {str( max_norm )}" )
				
		w.writelines( "\n--------- Train set metrics ---------\n" )
		w.writelines( "Recall = " + str( train[key][-1,1] ) + "\t" + \
			"Precision = " + str( train[key][-1,2] ) + "\t" + \
			"F1 = " + str( train[key][-1,3] ) + "\t" + \
			"AvgPrecision = " + str( train[key][-1,4] ) + "\t" + \
			"MatthewsCorrCoef = " + str( train[key][-1,5] ) + "\t" + \
			"AUROC = " + str( train[key][-1,6] ) + "\t" + \
			"Accuracy = " + str( train[key][-1,7] ) + "\t" + \
			"Loss = " + str( train[key][-1,0] ) )
		w.writelines("\nAverage time = " + str( round( mean( train[key][:, 9] ), 2 ) ) )
		w.writelines("\nTotal time = " + str( round( sum( train[key][:, 9] ), 2 ) ) )

		w.writelines("\n --------- Dev set metrics ---------\n")
		w.writelines("Recall = " + str( dev[key][-1,1] ) + "\t" + \
			"Precision = " + str( dev[key][-1,2] ) + "\t" + \
			"F1 = " + str( dev[key][-1,3] ) + "\t" + \
			"AvgPrecision = " + str( dev[key][-1,4] ) + "\t" + \
			"MatthewsCorrCoef = " + str( dev[key][-1,5] ) + "\t" + \
			"AUROC = " + str( dev[key][-1,6] ) + "\t" + \
			"Accuracy = " + str( dev[key][-1,7] ) + "\t" + \
			"Loss = " + str( dev[key][-1,0] ) + "\n\n" )

		w.writelines("\n --------- Test set metrics ---------\n")
		w.writelines("Recall = " + str( test[key][-1,0] ) + "\t" + \
					"Precision = " + str( test[key][-1,1] ) + "\t" + \
					"F1 = " + str( test[key][-1,2] ) + "\t" + \
					"AvgPrecision = " + str( test[key][-1,3] ) + "\t" + \
					"MatthewsCorrCoef = " + str( test[key][-1,4] ) + "\t" + \
					"AUROC = " + str( test[key][-1,5] ) + "\t" + \
					"Accuracy = " + str( test[key][-1,6] ) + "\n\n" )

	w.close()

	df = pd.DataFrame()
	params = list( train.keys() )
	df["HParams"] = params

	tmp = []
	epochs = None
	for key in train.keys():
		if epochs == None:
			epochs = train[key].shape[0]
			num = int( 0.2*epochs ) if epochs > 5 else epochs
			print( num )
		dif = ( train[key][-num:,0] - dev[key][-num:,0] )/train[key][-num:,0]
		tmp.append( np.mean( dif )*100 )

	df["Over/Under-fitting"] = tmp

	for set_, name in zip( [test, dev, train], ["Test", "Dev", "Train"] ):
		idx = 0 if name == "Test" else 1
		for heads in ["Recall", "Precision", "F1score", "AvgPrecision", "MatthewsCorrCoef", "AUROC", "Accuracy"]:
			tmp = []
			for key in train.keys():
				# print( name, "\t", heads, "\t", idx, "\t", key )
				tmp.append( set_[key][-1,idx] )
			df[f"{name}_{heads}"] = tmp
			idx += 1

	tmp = [time_dict[key] for key in time_dict.keys()]
	df["Total_time_(hrs)"] = tmp
	
	df.to_csv( f"summary_metrics{version}.csv" )
