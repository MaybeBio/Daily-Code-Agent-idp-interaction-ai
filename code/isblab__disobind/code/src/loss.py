import numpy as np

import torch
from torch import nn
import torch.utils.data as Data
import torch.nn.functional as F
from torch.autograd import Variable



# MSE Loss##########################################
class MSE():
    def __init__( self ):
        super().__init__()
        print( "Using MSE as loss." )
    
    def forward( self, pred, target, device, weight1 = None, apply_pos_weight = None, smooth = None ):
        criterion = nn.MSELoss()
        loss = criterion( pred, target )
        return loss


# MSE Loss##########################################
class RMSE():
    def __init__( self ):
        super().__init__()
        print( "Using RMSE as loss." )
    
    def forward( self, pred, target, device, weight1 = None, apply_pos_weight = None, smooth = None ):
        criterion = nn.MSELoss()
        loss = criterion( pred, target )
        return torch.sqrt( loss )


# Distribution based Loss functions#########################
############################################################
# Binary Cross Entropy Loss#########################
####################################################
class CE():
    def __init__( self ):
        super().__init__()
        print( "Using Cross Entropy loss." )
    
    def forward( self, pred, target, device, weight1 = None ):
        # Calculates the unweighted BCE loss.
        target = target.type( torch.LongTensor ).to( device )
        criterion = nn.CrossEntropyLoss()

        loss = criterion( pred, target )
        return loss


# Binary Cross Entropy Loss#########################
####################################################
class BCE():
    def __init__( self ):
        super().__init__()
        print( "Using BCE loss." )
    
    def forward( self, pred, target, device, weight, target_mask ):
        # Calculates the unweighted BCE loss.
        apply_mask, target_mask = target_mask

        criterion = nn.BCELoss( reduction = "none" )

        if apply_mask:
            loss = criterion( pred, target )*target_mask
        else:
            loss = criterion( pred, target )
        loss = torch.mean( loss )
        # exit()

        return loss


# BCE + MSE Loss####################################
####################################################
class BCE_MSE():
    def __init__( self ):
        super().__init__()
        print( "Using BCE+MSE loss." )
    
    def forward( self, pred, target, device, weight1 = None ):
        if weight1 != None:
            weight1 = torch.tensor( weight1 )
        pred, conf = pred
        # Calculates the unweighted BCE loss.
        criterion1 = nn.BCELoss()
        criterion2 = nn.MSELoss()
        bce = criterion1( pred, target )
        
        with torch.no_grad():
            criterion = nn.BCELoss( reduction = "none" )
            confidence = criterion( conf, target )

        mse = criterion2( pred, confidence )
        loss = bce + weight1*mse

        return loss


# BCE Loss with logits Loss#########################
####################################################
class BCEwithLogits():
    def __init__( self ):
        super().__init__()
        print( "Using BCE with logits loss.." )
    
    def forward( self, pred, target, device, weight, target_mask ):
        # Pytorch version of a weighted BCE loss.
        # pred = pred.view( -1 )
        # target = target.view( -1 )
        apply_mask, target_mask = target_mask
        weight, pos_weight = weight
        weight = torch.tensor( weight )

        if pos_weight != None:
            # neg, pos = torch.bincount( target.int(), minlength = 2 ).float()
            # weight2 = int( ( neg + smooth )/( pos + smooth ) )
            pos_weight = torch.tensor( pos_weight )
            criterion = torch.nn.BCEWithLogitsLoss( weight = weight.to( device ), 
                                                    pos_weight = pos_weight.to( device ), 
                                                    reduction = "none" )
        
        else:
            criterion = torch.nn.BCEWithLogitsLoss( weight = weight.to( device ), reduction = "none" )
        if apply_mask:
            loss = criterion(pred, target)*target_mask
        else:
            loss = criterion(pred, target)
        loss = torch.mean( loss )
        return loss


class FocalLoss( nn.Module ):
    def __init__( self, weight = None, size_average = True ):
        super().__init__()
        print( "Using Focal loss." )
        """
        Focal loss down-weighs the well-classified examples and focuses on hard, misclassified examples.
        alpha --> weighting factor for balancing the contribution of positive and negative classes.
                range = 0-1
                    alpha = 0.5-0.75 - underrepreseted positive class.
                    alpha = near 0 - underrepreseted negative class.
        gamma --> focusing factor. Emphasises the importance of hard examples.
        """

    def forward( self, pred, target, device, weight ):
        alpha, gamma = weight
        # Flatten label and prediction tensors
        # pred = pred.view( -1 )
        # target = target.view( -1 )
        
        # First compute binary cross-entropy 
        bce = F.binary_cross_entropy( pred, target, reduction = 'none')
        p_t = torch.exp( -bce )
        focal_loss = alpha *( 1 - p_t )**gamma* bce
                       
        return torch.mean( focal_loss )


class SingularityEnhancedLoss( nn.Module ):
    def __init__( self ):
        super().__init__()
        print( "Using Singularity Enhanced loss." )

    def forward( self, pred, target, device, weight, target_mask ):
        # Default alpha = 0.83, beta = 3.
        alpha, beta = weight
        apply_mask, target_mask = target_mask
        
        # Flatten label and prediction tensors
        inputs = pred.view( -1 )
        targets = target.view( -1 )
        
        first_term = alpha*target*torch.log( pred )
        second_term = ( 1 - alpha )*( 1 - target )*torch.log( 1 - pred )*torch.pow( ( 1 + pred ), beta )

        if apply_mask:
            se_loss = torch.sum( -1*( first_term + second_term )*target_mask )
        else:
            se_loss = torch.sum( -1*( first_term + second_term ) )
                       
        return se_loss


# Interface Loss####################################
####################################################
class InterfaceLoss():
    def __init__( self ):
        super().__init__()
        print( "Using Interface loss." )
        self.se_loss = SingularityEnhancedLoss()
        # self.criterion = nn.BCEwithLogits()
        # self.criterion = torch.nn.BCEWithLogitsLoss( weight = torch.tensor( [1.0] ).to( "cuda" ), 
        #                                         pos_weight = torch.tensor( [10] ).to( "cuda" ) )
    
    def forward( self, pred, target, device, weight, target_mask ):
        # Calculates the unweighted BCE loss.
        pred1, pred2 = pred
        target1, target2 = target
        target_mask1 = [target_mask[0], target_mask[1][0]]
        target_mask2 = [target_mask[0], target_mask[1][1]]

        loss1 = self.se_loss( pred1, target1, device, weight, target_mask1 )
        loss2 = self.se_loss( pred2, target2, device, weight, target_mask2 )

        loss = loss1 + loss2

        return loss



class ConcentratorLoss( nn.Module ):
    def __init__( self, weight = None, size_average = True ):
        super().__init__()
        print( "Using Concentrator loss." )
        """
        Weighted loss function in which the weights are obtained based on the 
            difference from the target value.
        """

    def forward( self, pred, target, device, weight, target_mask ):
        alpha, gamma = weight
        apply_mask, target_mask = target_mask

        first_term = target*torch.log( pred )*torch.exp( target - pred )
        second_term = ( 1 - target )*torch.log( 1 - pred )*torch.exp( pred - target )

        if apply_mask:
            loss = -1*( first_term + second_term )*target_mask
        else:
            loss = -1*( first_term + second_term )
        loss = torch.mean( loss )

        return loss


## ------------------------------------------------------------------------------- ##
# Region based Loss functions#########################
######################################################

# Dice Loss#########################################
####################################################
class DiceLoss():
    def __init__( self, weight = None, size_average = True ):
        super().__init__()
        print( "Using Dice loss." )

    def forward( self, pred, target, device, smooth = 1 ):
        # Flatten label and prediction tensors.
        pred = pred.view( -1 )
        target = target.view( -1 )
        
        intersection = ( pred * target ).sum()
        dice = ( 2 * intersection + smooth )/( pred.sum() + target.sum() + smooth )
        
        return 1 - dice


# Tversky Loss#########################################
#######################################################
class TverskyLoss( nn.Module ):
    def __init__( self, weight = None, size_average = True ):
        super().__init__()
        print( "Using Tversky loss." )

    def forward( self, pred, target, device, weight, target_mask, smooth = 1 ):
        # Default: alpha = 0.5, beta = 0.5
        alpha, beta = weight
        # Flatten label and prediction tensors
        pred = pred.view( -1 )
        target = target.view( -1 )
        
        # True Positives, False Positives & False Negatives
        TP = ( pred * target ).sum()
        FP = ( ( 1 - target ) * pred ).sum()
        FN = ( targets * ( 1 - pred ) ).sum()

        alpha, beta = 0.7, 0.7
       
        Tversky = ( TP + smooth ) / ( TP + alpha*FP + beta*FN + smooth )
        
        return 1 - Tversky


## ------------------------------------------------------------------------------- ##
# Hybrid Loss functions#########################
################################################

# Dice-BCE Loss#####################################
####################################################
class DiceBCELoss( nn.Module ):
    def __init__( self, weight = None, size_average = True ):
        super().__init__()
        print( "Using Dice_BCE loss." )

    def forward( self, pred, target, device, weight, target_mask, smooth = 1 ):
        # Flatten label and prediction tensors.        
        pred = pred.view( -1 )
        target = target.view( -1 )
        
        intersection = ( pred * target ).sum()
        dice = 1 - ( 2 * intersection + smooth )/( pred.sum() + target.sum() + smooth )
        
        BCE = F.binary_cross_entropy( pred, target, reduction = "mean" )
        Dice_BCE = BCE + dice
        
        return Dice_BCE


# Representation Loss Functions (Dscript)#############################
######################################################################
class RepresentationLoss():
    def __init__( self ):
        super().__init__()
        print( "Using Representation loss." )
    
    def forward( self, pred, target, device, weight, target_mask ):
        lambda_, pos_weight = weight
        apply_mask, target_mask = target_mask
        # Implements a representation loss.
        if lambda_ >1 and lambda_ <0:
            raise ValueError( "Weight argument out of range (0-1)." )
        
        # criterion1 = nn.BCELoss()
        pos_weight = torch.tensor( pos_weight )
        criterion1 = torch.nn.BCEWithLogitsLoss( pos_weight = pos_weight.to( device ), reduction = "none" )
        if apply_mask:
            bce = torch.mean( criterion1( pred, target )*target_mask )
            representation_loss = torch.mean( torch.sigmoid( pred )*target_mask )
        else:
            bce = torch.mean( criterion1( pred, target ) )
            representation_loss = torch.mean( torch.sigmoid( pred ) )
        loss = lambda_*bce + (1 - lambda_)*representation_loss
        # loss = criterion1( pred, target ) + representation_loss
        return loss


# L1 regularized BCE Loss Functions#########################
############################################################
class L1RegularizedLoss():
    def __init__( self ):
        super().__init__()
        print( "Using L1 regularized hybrid BCE loss." )
    
    def forward( self, pred, target, device ):
        # Implements a regularized loss, combiing BCE and L1norm.
        criterion1 = nn.BCELoss()
        criterion2 = nn.L1Loss()
        loss = criterion1( pred, target ) + criterion2( pred, target )
        return loss


# L2 regularized BCE Loss Functions#########################
############################################################
class L2RegularizedLoss():
    def __init__( self ):
        super().__init__()
        print( "Using L2 regularized hybrid BCE loss." )
    
    def forward( self, pred, target, device ):
        # Implements a regularized loss, combiing BCE and L2norm.
        criterion1 = nn.BCELoss()
        criterion2 = nn.MSELoss()
        loss = criterion1( pred, target ) + criterion2( pred, target )
        return loss


# Confidence Loss Functions#################################
############################################################
class ConfidenceLoss():
    def __init__( self ):
        super().__init__()
        print( "Using Confidence Loss." )
    
    def forward( self, pred, target, device, lambda_ = None ):
        y_pred, confidence = pred
        criterion1 = nn.BCELoss( reduction = "none" )
        # criterion1 = nn.CrossEntropyLoss()
        target = target.type( torch.LongTensor ).to( device )

        conf_loss = -torch.log( confidence )
        loss = criterion1( y_pred, target ) + lambda_*conf_loss
        loss = torch.mean( loss )
        return loss


#Get Loss Functions###############################
##################################################
def get_loss_function( loss_ ):
    if loss_ == "ce":
        return CE()
    elif loss_ == "bce":
        return BCE()
    elif loss_ == "bce_with_logits":
        return BCEwithLogits()
    elif loss_ == "representation_loss":
        return RepresentationLoss()
    elif loss_ == "l1_regularized_loss":
        return L1RegularizedLoss()
    elif loss_ == "l2_regularized_loss":
        return L2RegularizedLoss()
    elif loss_ == "dice":
        return DiceLoss()
    elif loss_ == "dice_bce":
        return DiceBCELoss()
    elif loss_ == "focal":
        return FocalLoss()
    elif loss_ == "tversky":
        return TverskyLoss()
    elif loss_ == "se_loss":
        return SingularityEnhancedLoss()
    elif loss_ == "interface":
        return InterfaceLoss()
    elif loss_ == "concentrator":
        return ConcentratorLoss()
    elif loss_ == "conf_loss":
        return ConfidenceLoss()
    elif loss_ == "bce_mse":
        return BCE_MSE()
    else:
        return None

