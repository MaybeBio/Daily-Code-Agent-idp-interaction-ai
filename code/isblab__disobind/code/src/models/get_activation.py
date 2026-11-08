import torch
from torch import nn


#Logistic Activation##########################
##############################################
class LogisticActivation(nn.Module):
	# A eneralized Sigmoid function applied elementwise.
	def __init__( self, x0 = 0, eta = 1, train = False ):
		super( LogisticActivation, self ).__init__()
		self.x0 = x0
		self.eta = nn.Parameter( torch.FloatTensor( [float( eta )] ) )
		self.eta.requiresGrad = train

	def forward( self, x ):
		o = torch.clamp(
		    1 / ( 1 + torch.exp( -self.eta * ( x - self.x0 ) ) ), min = 0, max = 1 )
		return o

	def clip( self ):
		# Restricts the sigmoid slope (k) to be >= 0 if k is trained.
		self.eta.data.clamp_( min = 0 )


#get Activation function##########################
##################################################
def get_activation( activation_, device, param = None ):

	if activation_ == "sigmoid":
		return nn.Sigmoid()

	elif activation_ == "tanh":
		return nn.Tanh()

	elif activation_ == "softmax":
		return nn.Softmax( dim = 1 )

	elif activation_ == "relu":
		return nn.ReLU()

	elif activation_ == "leakyrelu":
		if param == None:
			return nn.LeakyReLU()
		else:
			return nn.LeakyReLU( negative_slope = param )

	elif activation_ == "elu":
		if param == None:
			return nn.ELU()
		else:
			return nn.ELU( alpha = param )

	elif activation_ == "gelu":
		return nn.GELU()

	elif activation_ == "celu":
		return nn.CELU()

	elif activation_ == "prelu":
		if param == None:
			return nn.PReLU( device = device )
		else:
			return nn.PReLU( num_parameters = param, device = device )

	elif activation_ == "silu":
		return nn.SiLU()

	elif activation_ == "mish":
		return nn.Mish()

	elif activation_ == "logistic_activation":
		x0, eta, train = param
		return LogisticActivation( x0 = x0, eta = eta, train = train )
	
	else:
		return None

